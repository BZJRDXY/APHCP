"""Seal the final M2 first pass, independently sample it, and prepare a bounded M3 review list."""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime
import json
from pathlib import Path

from holo_candidate_batch import atomic_json, now, sha256
from holo_m2_blind_audit import audit


def get(p):
    return json.loads(p.read_text())


def table(p):
    with p.open(newline='',encoding='utf8') as f: return list(csv.DictReader(f))


def seal(run, entries, systems):
    dst=run/'blind_first_pass.json'
    if dst.exists():
        return get(dst)
    hashes={p:sha256(run/p) for p in ('frozen.json','build_inputs.json','entry_results.csv','system_results.csv','summary.json')}
    payload={'sealed_utc':now(),'inputs_and_tables_sha256':hashes,'entry_count':len(entries), 'system_count':len(systems),
             'entry_discovery_sha256':{},'entry_systems_sha256':{}}
    if not (run/'first_pass').is_dir(): (run/'first_pass').mkdir()
    for name in hashes:
        original=run/name; copy=run/'first_pass'/name
        if not copy.exists(): copy.write_bytes(original.read_bytes())
        assert sha256(copy)==hashes[name]
    for e in entries:
        p=e['pdb_id']
        for folder,out in [('discovery','entry_discovery_sha256'),('system_results','entry_systems_sha256')]:
            q=run/folder/f'{p}.json'
            if q.is_file(): payload[out][p]=sha256(q)
    atomic_json(dst,payload); return payload


def pick(evidence, entries, systems):
    targets={'route_index','route_native','route_mismatch','no_output','no_candidate','failure','pair','copy','crosslink','multiple','multi_model','altloc','glycan','peptide','negative','nonlocal','third_endpoint'}
    groups=defaultdict(list)
    for s in systems: groups[s['pdb_id']].append(s)
    out=[]; stats=Counter(); by={r['pdb_id']:r for r in evidence}
    available={item['pdb_id']:item for item in entries}
    views={}
    for p,e in available.items():
        ev=by[p]; f=set(ev['features']); g=groups[p]; tag=set()
        tag.add({'INDEX_MATCH':'route_index','NO_LIGAND_ANNOTATION':'route_native','INDEX_VERSION_MISMATCH':'route_mismatch'}[ev['group']]);
        if e['status']=='TECHNICAL_FAILURE': tag.add('failure')
        if int(e['structure_readback_pass'] or 0)==0: tag.add('no_output')
        if not g: tag.add('no_candidate')
        for src,dst in [('assembly_copy','copy'),('covalent_link','crosslink'),('repeated_ccd_instance','multiple'),('multi_model','multi_model'),('altloc_title_hint','altloc'),('branched_entity','glycan'),('polymer_peptide','peptide'),('additive_only_hint','negative')]:
            if src in f:tag.add(dst)
        if any('base_pairs' in s.get('member_grouping','') for s in g):tag.add('pair')
        if any('symmetry' in s.get('reason','') or 'symmetry' in s.get('member_grouping','') for s in g):tag.add('nonlocal')
        if any('endpoint' in s.get('reason','') and 'multi' in s.get('reason','') for s in g): tag.add('third_endpoint')
        views[p]=tag
    for _ in range(min(20,len(available))):
        def score(p):
            f=views[p]; useful=f&targets; return (sum((4 if t in ('failure','third_endpoint','nonlocal','no_output','no_candidate') else 2)/(.5+stats[t]) for t in useful)
                                                  + .25/(1+stats[by[p]['method']]), -by[p]['atom_count'], p)
        p=max(available,key=score); available.pop(p); tag=views[p]; stats.update(tag); stats[by[p]['method']]+=1
        out.append({'pdb_id':p,'reason':', '.join(sorted(tag))})
    return {'n':len(out),'selection_basis':'post-seal stratified diagnostic sample','coverage':dict(stats),'entries':out}


def m3_rows(run, entries, systems, audit_result):
    by={e['pdb_id']:e for e in entries}; good=[]; hold=Counter(); detailed=[]
    tested={x['disk']['candidate_key'] for x in audit_result['entries'] if x.get('disk') and x['audit_status']=='INDEPENDENT_READBACK_PASS'}
    for s in systems:
        if s['candidate_intake']!='PASS':
            hold[s.get('reason') or ('role_required' if s.get('role_qc')!='PASS' else 'partner_scope_required')]+=1; continue
        e=by[s['pdb_id']]; p=s['pdb_id']; detailed.append(s)
        complex_model=len(json.loads(e['models'] or '[]'))!=1
        branch=bool(s.get('member_grouping')) or len(json.loads(s['partner_instances'] or '[]'))!=1
        if complex_model or branch:
            hold['positive_candidate_but_complex_scope_requires_review']+=1; continue
        good.append((s['candidate_key'] not in tested,p,s))
    picked=[]; used=set(); per_entry={}; discoveries={}
    for _,p,s in sorted(good,key=lambda q:(q[0],q[1],q[2]['candidate_key'])):
        if p in used:continue
        entry=by[p]
        if p not in per_entry:
            per_entry[p]={z['candidate_key']:z for z in get(run/'system_results'/f'{p}.json')}
            discoveries[p]={z['candidate_key']:z['members'] for z in get(run/'discovery'/f'{p}.json')['native_observations']}
        row=per_entry[p][s['candidate_key']]; art=row['artifact']; parts=discoveries[p].get(s['candidate_key'],[])
        if not parts or len(parts)!=1 or parts[0].get('kind') not in ('component','polymer'):continue
        if parts[0].get('kind')=='polymer' and (parts[0].get('polymer_type','').lower().find('nucleotide')!=-1 or not art['connections']):continue
        if art['omitted_touching_connections'] or art.get('base_pairs'):continue
        used.add(p)
        picked.append({'pdb_id':p,'revision':s['revision'],'assembly':s['assembly'],'source_model':s['representative_model'],
                       'candidate':s['candidate_key'],'source_sha256':entry['source_sha256'],'frozen_source_proof':f'source_proofs/{p}.json',
                       'system_cif':s['system_path'],'system_sha256':s['system_sha256'],'mapping':s['mapping'],'role_evidence':row['role_basis'],
                       'independent_disk_sample':s['candidate_key'] in tested,'status':'REVIEW_BEFORE_PRODUCTION',
                       'remaining_checks':'Confirm biological role and full partner scope; then re-read pinned CIF/lineage and QC before registration.'})
        if len(picked)>=30:break
    return picked,detailed,hold


def run_manual_proof(run, rescued):
    extra=0
    for folder in run.glob('supplement_*'):
        note=folder/'imported.txt'
        if (folder/'transport.json').is_file() or not note.is_file():
            continue
        for p in note.read_text().split():
            proof=run/'source_proofs'/f'{p}.json'
            if proof.is_file() and get(proof).get('status')=='PASS' and p not in rescued:
                extra+=get(proof)['bytes']; rescued[p]+=1
    return extra


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('run',type=Path); args=ap.parse_args(); r=args.run.resolve()
    entries=table(r/'entry_results.csv'); systems=table(r/'system_results.csv'); frozen=get(r/'frozen.json'); metadata=get(r/'selection_evidence.json'); summary=get(r/'summary.json'); source_ids={x['pdb_id'] for x in frozen['files']}
    assert len(entries)==len(source_ids)==2000 and {x['pdb_id'] for x in entries}==source_ids and len(metadata)==2000 and len(systems)==len({(s['pdb_id'],s['candidate_key']) for s in systems})
    seal(r,entries,systems)
    if not (r/'audit_sample.json').exists(): atomic_json(r/'audit_sample.json',pick(metadata,entries,systems))
    previous=r/'independent_audit.json'
    if previous.is_file():
        result=get(previous)
        assert result['first_pass_sha256']==sha256(r/'blind_first_pass.json') and result['n']==20
    else:
        audit(r,frozen['files'])
        result=get(previous)
    good,det,hold=m3_rows(r,entries,systems,result)
    proofs=[get(x) for x in (r/'source_proofs').glob('*.json')]; rescued=Counter()
    transfers={}
    for f in sorted(r.glob('supplement_*/transport.json')):
        data=get(f); n=f.with_name('imported.txt'); names=n.read_text().split() if n.is_file() else []
        for v in data['results']:
            if v['status']=='CACHED': transfers.setdefault(v['pdb_id'],v['bytes'])
        rescued.update(names)
    manual=run_manual_proof(r, rescued)
    times=[(a['started_utc'],a['finished_utc']) for proof in proofs for a in proof.get('attempts',[]) if a.get('started_utc') and a.get('finished_utc')]
    intervals={}
    if times: intervals['elapsed_fetch_proof_clock_seconds']=(datetime.fromisoformat(max(x[1] for x in times))-datetime.fromisoformat(min(x[0] for x in times))).total_seconds()
    success=[e for e in entries if str(e['status'])!='TECHNICAL_FAILURE']; shown=str(len(good));
    atomic_json(r/'m3_first_review.json',{'created_utc':now(),'scope':'Preliminary candidates for individual review; no production approval',
        'run':str(r),'input_sha256':sha256(r/'frozen.json'),'build_versions':get(r/'build_inputs.json'),'verification_code_sha256':sha256(Path(__file__)),
        'review_count':len(good),'review_candidates':good,'all_initial_intake_pass_count':len(det),'held_candidate_reasons':dict(hold),
        'handoff_run_steps':['Use only pinned originals and frozen build versions.', 'Confirm intended ligand/partner role with independent source evidence for each candidate.', 'Read back atoms, experimental values, source instance and link identities before registration.', 'Keep unsupported and role-ambiguous candidates in REVIEW.'],
        'capacity_basis':{'entries':len(entries),'processed_without_technical_failure':len(success),'actual_system_and_mapping_bytes':summary['sizes']['systems']['bytes'],
            'actual_assembly_bytes':summary['sizes']['assemblies']['bytes'],'actual_compressed_original_bytes':sum(q['bytes'] for q in proofs if q['status']=='PASS'),
            'measured_bytes_per_selected_entry':round((summary['sizes']['systems']['bytes']+summary['sizes']['assemblies']['bytes'])/len(entries),1),
            'elapsed_build_invocation':summary['invocation'],**intervals}})
    atomic_json(r/'final_delivery_metrics.json',{'created_utc':now(),'run':str(r),'first_pass_sha256':sha256(r/'blind_first_pass.json'),
        'entry_status':dict(Counter(e['status'] for e in entries)),'no_system_output':sum(not int(e['structure_readback_pass'] or 0) for e in entries),
        'missing_seed_ids':sorted(source_ids-{e['pdb_id'] for e in entries}),'raw_proof':dict(Counter(q['status'] for q in proofs)),
        'remote_recorded_http_body_bytes':summary['network'].get('observed_http_body_bytes',0),
        'supplement_recorded_body_bytes_lower_bound':sum(transfers.values())+manual,'supplement_files_imported':len(rescued),
        'audit':result['status_counts'],'m3_review_count':len(good),'calendar':intervals})
    print('Final seal exists; independently audited',result['n'],'M3 review list',shown,flush=True)

if __name__=='__main__': main()
