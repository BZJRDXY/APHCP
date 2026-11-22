"""Read back the staged registry and its real files, plus a few scope/role controls."""
import argparse
from collections import Counter
from pathlib import Path
import subprocess
import sys

import holo_m3_register as m


def verify(run, subset=None):
    root=run.parent.parent; registration=m.read(run/'registration_receipt.json'); rerun=m.read(run/'rerun_receipt.json')
    path=Path(registration['registry_path']); rows=m.read(run/'review_results.json')
    assert registration['registry_counts']==rerun['registry_counts'] and registration['results_sha256']==rerun['results_sha256']
    assert registration['readback_files_reused_on_disk']==rerun['prior_review_checks_reused'] and rerun['database_registration_reused']
    keys={(r['batch'],r['pdb_id'],r['candidate_key']) for r in rows}; assert len(keys)==len(rows)
    by={x['pdb_id']:x for x in rows}; db_result={}
    with m.connect(path,readonly=True) as db:
        assert db.execute('PRAGMA integrity_check').fetchone()[0]=='ok' and not db.execute('PRAGMA foreign_key_check').fetchall()
        actual=Counter(r['status'] for r in db.execute("SELECT status FROM run_item WHERE stage='M3_REVIEW'"))
        assert actual==Counter(x['disposition'] for x in rows)
        assert db.execute('SELECT COUNT(*) FROM release').fetchone()[0]==0
        assert db.execute('SELECT COUNT(*) FROM batch').fetchone()[0]==1 and db.execute('SELECT COUNT(*) FROM run').fetchone()[0]==1
        db_result['work_item_kinds']={kind:n for kind,n in db.execute('SELECT kind,COUNT(*) FROM work_item GROUP BY kind')}
        db_result['registry_review_status']=dict(actual)
        own=set(db.execute('SELECT object_id FROM run_item WHERE stage=\'M3_REVIEW\' AND status=\'REGISTERED\'').fetchall())
        assert {v[0] for v in own}=={v['system_id'] for v in rows if v['disposition']=='REGISTERED'}
    exe=[sys.executable,'-B', str(run/'code/holo_m3_register.py'),'lookup','--run',str(run),'--verify-files']
    confirmed=by['1AVD']; held=by['1JPL']; checked={}
    for selector,arg,label in [('--pdb','1AVD','confirmed_by_pdb'),('--system',confirmed['system_id'],'confirmed_by_system'),
                               ('--system',held['system_id'],'held_by_system')]:
        viewed=__import__('json').loads(subprocess.check_output(exe+[selector,arg],text=True)); found=viewed['systems']
        expected=held if label=='held_by_system' else confirmed; assert found==[expected] and viewed['read_only_files_verified']
        assert viewed['registry']['entry']['id']==expected['pdb_id']
        files=found[0]['files']; block=m.gemmi.cif.read_file(files['system_cif']['path']).sole_block()
        assert len(block.find_values('_atom_site.id'))==expected['qc']['atoms']>0
        raw,rev=m.revision_from_cif(Path(files['source_mmcif']['path']),expected['pdb_id']); assert rev==expected['revision']
        mapping=m.read(files['source_mapping']['path']); assert mapping['artifact']['sha256']==files['system_cif']['sha256']
        if 'complete_assembly' in files:
            full=m.gemmi.cif.read_file(files['complete_assembly']['path']).sole_block()
            assert len(full.find_values('_atom_site.id'))==files['complete_assembly']['atoms']>expected['qc']['atoms']
        objects=viewed['registry']['objects']; system=next(x for x in objects['system'] if x['id']==expected['system_id'])
        assert system['protein_instances'] and system['ligand_instances']; assert objects['assembly'] and objects['chain'] and objects['ligand']
        checked[label]={'id':expected['system_id'],'pdb_id':expected['pdb_id'],'revision':rev,'atoms':expected['qc']['atoms'],
                        'outcome':expected['disposition'],'files':files,'reasons':expected['reasons']}
    asset=subprocess.check_output(exe+['--system',confirmed['system_id'],'--asset','system_cif'],text=True).strip()
    assert asset==confirmed['files']['system_cif']['path']
    env=__import__('os').environ.copy(); env['PYTHONPATH']=str(run/'code/src')
    official=__import__('json').loads(subprocess.check_output([sys.executable,'-B','-m','aphcp_control','--data-root',str(path.parents[1]),
                                            'status','pdb','1AVD'],env=env,text=True))
    assert '1AVD' in m.raw_json(official) and confirmed['system_id'] in m.raw_json(official)
    # Negative control from the same deposited experimental entry, not a different CCD in a different experiment.
    est=by['1AQU']; block,_=m.revision_from_cif(Path(est['files']['source_mmcif']['path']),'1AQU')
    st=m.gemmi.make_structure_from_block(block); model=next(x for x in st if x.num==est['source_model'])
    obs=[o for o in m.read(est['files']['discovery']['path'])['native_observations'] if o['assembly']==est['assembly'] and o['model']==model.num]
    target=lambda o:{r['ccd'] for v in o['members'] for r in v['residues']}
    roles=[]
    for o in obs:
        if target(o) in ({'EST'},{'A3P'}):
            status,e,conflict,contacts,title=m.source_role(block,o,st,model,est['assembly'])
            roles.append({'ccd':next(iter(target(o))),'role':status,'contact_count':len(contacts),'candidate_key':o['candidate_key'],
                          'comparison':'estradiol target / PAP (CCD A3P) cofactor named in same depositor title'})
    assert any(x['ccd']=='A3P' and x['role']!='PASS' and x['contact_count'] for x in roles)
    assert any(x['ccd']=='EST' and x['role']=='PASS' and x['contact_count'] for x in roles)
    # D-0018 exercises real source, map and coordinates, with no type annotation passed to registration preflight.
    prev=m.read(run/'frozen_staging.json'); manifest=m.read_manifest(Path(prev['manifest'])); key=next((a,b,c) for a,b,c,_ in manifest if b=='1AVD')
    inp=next(n for a,b,c,n in manifest if (a,b,c)==key); p=run.parent; old=p/'holo-m2-final2000-20261122a'; ref=p/'holo-m2-full-reconciliation-20261122a'
    base,_=m.pool_rows(ref,p/'holo-m2-chemtype-20261122a',{key})
    unknown=m.assess_candidate(key,inp,base[key],None,(old,ref,p/'holo-m2-blind200-20261122a'),old,None)
    assert unknown['type']['primary_class']=='UNKNOWN' and unknown['disposition']=='REGISTERED' and unknown['system_id']==confirmed['system_id']
    subset_receipt=m.read(subset/'registration_receipt.json') if subset else None
    if subset:
        assert subset_receipt['input_candidate_count']==3 and subset_receipt['dispositions']=={'REGISTERED':2,'REVIEW':1}
        small=m.read(subset/'review_results.json'); assert all(x==by[x['pdb_id']] for x in small)
    out={'registry':db_result,'lookup':checked,'d0018_unknown_type_control':{'id':unknown['system_id'],'disposition':unknown['disposition'],
         'primary_class':unknown['type']['primary_class']},'role_positive_negative_control':roles,'replicated_manifest':True,
         'existing_status_cli':'PASS','submanifest':subset_receipt}
    m.write(run/'verification.json',out);print(m.raw_json({k:v for k,v in out.items() if k!='lookup'}))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('run',type=Path);a.add_argument('--subset',type=Path);args=a.parse_args();verify(args.run,args.subset)
