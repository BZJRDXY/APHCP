"""Verify an incremental M3 Run and exercise a decision transition in isolated copies."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import sys

from holo_m3_decisions import resolve
from holo_m3_observations import review_records, key
from holo_m3_register import (
    backup, connect, database, file_ref, frozen_registry_root, index_review_results,
    read, recorded_assessments, require, row_counts, sha, utc_now, write,
)


def main(a):
    root=a.root.resolve(strict=True); runs=root/'runs'; pilot=a.pilot.resolve(strict=True); m2=a.m2.resolve(strict=True)
    seed=a.seed.resolve(strict=True); seed_db=database(frozen_registry_root(root,seed))
    staging=database(frozen_registry_root(root,pilot)); protected={str(p):sha(p) for p in (seed_db,staging)}
    code=Path(__file__).resolve().parent
    def invoke(*args):
        return json.loads(subprocess.check_output([sys.executable,'-B',str(code/'holo_m3_register.py'),*map(str,args)]))
    results=read(pilot/'review_results.json'); original=read(pilot/'registration_receipt.json'); replay=read(pilot/'rerun_receipt.json')
    require(replay['database_registration_reused'] and replay['prior_review_checks_reused']==len(results) and
            sha(pilot/'review_results.json')==original['results_sha256']==replay['results_sha256'], 'replay_mismatch')
    for row in results:
        for ref in row.get('files',{}).values():file_ref(ref['path'],(runs,),ref['sha256'])
    with connect(staging,readonly=True) as db:
        require(db.execute('PRAGMA integrity_check').fetchone()[0]=='ok' and not db.execute('PRAGMA foreign_key_check').fetchall()
                and db.execute('SELECT COUNT(*) FROM release').fetchone()[0]==0,'pilot_registry_check_failed')

    actual=next(row for row in results if row['disposition']=='REGISTERED'); sid=actual['system_id']; pb=actual['pdb_id']
    matched=invoke('lookup','--run',pilot,'--pdb',pb,'--verify-files')['systems']
    require(sid in {x.get('system_id') for x in matched},'PDB_lookup_missing_confirmed_system')
    system=invoke('lookup','--run',pilot,'--system',sid,'--verify-files')['systems'][0]
    require(system['disposition']=='REGISTERED','system_lookup_missing_confirmed_decision')

    previous=runs/(pilot.name+'-z-before-'+a.control_tag); later=runs/(pilot.name+'-a-after-'+a.control_tag); previous.mkdir(exist_ok=True); later.mkdir(exist_ok=True)
    prior_db=database(frozen_registry_root(root,previous,create=True))
    if not prior_db.exists():backup(seed_db,prior_db)
    prev_row=actual | {'disposition':'REVIEW','reasons':['SIMULATED_REVIEW_FOR_REGRESSION_ONLY']}
    outputs=previous/'review_results.json'; write(outputs,[prev_row]); seal=read(pilot/'frozen_staging.json')
    seal=seal | {'manifest':str(outputs),'manifest_sha256':sha(outputs),'test_only':True}
    index_review_results(prior_db.parent.parent,outputs,[prev_row],seal,outputs,utc_now(),0.)
    old=invoke('lookup','--run',previous,'--system',sid,'--verify-files')['systems'][0]
    require(old['disposition']=='REVIEW' and old['current_decision_id'] is not None, 'old_control_decision_missing')
    manifest=read(pilot/'registration_input.json'); n=next(x for x in manifest['review_candidates'] if key(x)==key(actual))
    n=n | {'supersedes':[old['current_decision_id']],
           'decision_basis':'Test only: remove the simulated hold using the normal fixed-source M3 assessment.'}
    new_manifest=later/'registration_input.json'; write(new_manifest,manifest | {'review_count':1,'review_candidates':[n]})
    cmd=('register','--root',root,'--run',later,'--manifest',new_manifest,'--git-head',seal['workspace_git_head'],
         '--seed-run',previous,'--jobs','1')
    answer=invoke(*cmd); replay_control=invoke(*cmd)
    require(answer['dispositions']=={'REGISTERED':1} and replay_control['database_registration_reused'], 'control_did_not_register')
    after=invoke('lookup','--run',later,'--system',sid,'--verify-files')['systems'][0]
    require(after['system_id']==sid and after['disposition']=='REGISTERED' and after['current_decision_id'].startswith(later.name)
            and after['files']==old['files'] and {r['disposition'] for r in after['provenance']}=={'REVIEW','REGISTERED'}, 'transition_lost_identity_evidence_or_files')
    with connect(prior_db,readonly=True) as pcon, connect(database(frozen_registry_root(root,later)),readonly=True) as db:
        prior_counts=row_counts(pcon); now_counts=row_counts(db)
        require(all(prior_counts[k]==now_counts[k] for k in ('entry_version','asset_file','work_item','system_member')),'transition_duplicated_assets')
        records=list(recorded_assessments(db,runs,system=sid)); n=len(records)
        require(n==2 and resolve(records)==resolve(records[::-1]) and latest_decision(records)==after['current_decision_id'], 'order_dependent_decision')
        catalog=review_records(db,runs)
        require(catalog[key(actual)][sid]['current_decision_id']==after['current_decision_id'],'catalog_effective_decision_differs')
    require(all(sha(p)==h for p,h in protected.items()),'protected_registry_changed')

    by_entry={}; routes=Counter(); sys_results=Counter(); observed_routes=Counter(); download_bytes=0; stored=Counter()
    for f in sorted((m2/'entry_results').glob('*.json')):
        e=read(f); p=e['pdb_id']; routes[e['stratum']]+=1
        proof=read(m2/'source_proofs'/f'{p}.json'); download_bytes+=sum(x['body_bytes'] for x in proof['attempts'])
        observed=Counter(s['status'] for s in read(m2/'system_results'/f'{p}.json')); sys_results.update(observed)
        if observed['EXPORTED_READBACK']:observed_routes[e['stratum']]+=1
        decisions=Counter(x['disposition'] for x in results if x['pdb_id']==p)
        by_entry[p]={'revision':e['revision'],'index_revision':e['index_source_revision'],'route':e['route'],'stratum':e['stratum'],
                     'compressed_bytes':e['source_bytes'],'export_statuses':dict(observed),'m3_outcomes':dict(decisions)}
    for name in ('raw','systems','assemblies'):
        stored[name]=sum(f.stat().st_size for f in (m2/name).rglob('*') if f.is_file())
    require(set(observed_routes)==set(routes) and sum(routes.values())==read(m2/'frozen.json')['plan_count'], 'pilot_routes_did_not_all_export')
    summary={'status':'PASS','verified_utc':utc_now(),'first_seconds':read(pilot/'increment_first.json')['elapsed_s'],
             'replay_seconds':read(pilot/'increment_receipt.json')['elapsed_s'],'entries':len(by_entry),'by_entry':by_entry,
             'strata':dict(routes),'entry_routes_with_readback':dict(observed_routes),'system_statuses':dict(sys_results),
             'm3_outcomes':original['dispositions'],'candidate_count':read(pilot/'candidate_evidence.json')['expected_candidates'],
             'registration_count':original['input_candidate_count'],'type_counts':read(pilot/'chemical_labels.receipt.json')['primary_classes'],
             'downloaded_bytes':download_bytes,'stored_bytes':dict(stored),'protected_registries_unchanged':True,'formal_release_count':0,
             'replay_preserved_files_and_results':True,'lookup_pdb':pb,'system_lookup':sid,'transition':{
                 'test_only':True,'before':old['current_decision_id'],'after':after['current_decision_id'],
                 'same_system':True,'same_files':True,'history_count':n,'reverse_order_pass':True,'catalog_pass':True}}
    write(pilot/'independent_verification.json',summary);print(json.dumps(summary,sort_keys=True))


def latest_decision(records):
    return resolve(records)['current_decision_id']


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','pilot','m2','seed'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--control-tag',default='v1')
    main(p.parse_args())
