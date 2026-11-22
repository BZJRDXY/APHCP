"""Freeze, run and resume the authorized APHCP M2 fixed-revision pilot."""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
from importlib.metadata import version
import json
from pathlib import Path
import re
import resource
import sqlite3
import time
import traceback

import gemmi
from plinder.core.release import PlinderRelease
from holo_candidate_batch import atomic_json, download_source, now, process_entry, revision_from_cif, sha256
from holo_native_candidates import inspect_source
from holo_source_mapping import assembly_instances, altloc, connection_edges, observed_endpoint, subchains


def load(path):
    return json.loads(path.read_text(encoding='utf8'))


def freeze(run, inp, args):
    nomination_path=inp/'m2_nominations.json'; nominations=load(nomination_path); plan=[]
    if not isinstance(nominations,list) or not nominations:
        raise ValueError('M2 nominations must be a nonempty list')
    for n in nominations:
        pdb=n['pdb_id'].upper(); f=n['m2_candidate_file']; rev=n['latest_revision']
        if not re.fullmatch(r'[0-9A-Z]{4}',pdb) or not re.fullmatch(r'\d+\.\d+',rev):
            raise ValueError(f'Invalid PDB ID or exact structure revision: {pdb} {rev}')
        low=pdb.lower(); base='pdb_0000'+low
        url=f'https://files-versioned.rcsb.org/pdb_versioned/data/entries/{low[1:3]}/{base}/{base}_xyz_v{rev.replace(".","-")}.cif.gz'
        if (f['source_url']!=url or f['revision']!=rev or f['file_type']!='pdbx_mmcif_gzip'
                or f['head_status']!=200 or f['compressed_bytes_p50']<=0):
            raise ValueError(f'Unverified archive identity or size: {pdb} {rev}')
        plan.append({'pdb_id':pdb,'revision':rev,'source_url':url,'expected_bytes':f['compressed_bytes_p50'],
            'm1_evidence':{k:n.get(k) for k in ('title','test_stratum','why_and_test','plinder_source_revision','annotation_present','assembly_ids','actual_index_examples')},
            'historical_comparison_not_authorized':n.get('optional_comparator_file')})
    if len(plan)!=len({r['pdb_id'] for r in plan}):
        raise ValueError('M2 nominations must have unique PDB IDs')
    code=Path(__file__).resolve().parent
    files=[nomination_path,inp/'source_role_rules.json',*[code/name for name in
        ('holo_m2_pilot.py','holo_candidate_batch.py','holo_native_candidates.py','holo_pilot.py','holo_source_mapping.py')]]
    signature={p.name:sha256(p) for p in files if p.is_file()}
    base={'batch_id':args.batch_id,'run_id':run.name,'git_head':args.git_head,
        'scope':'Exact revisions in the authorized nominations only','plinder_release':str(args.plinder_release),
        'environment':{'gemmi':gemmi.__version__,'python':__import__('platform').python_version(),'plinder':version('plinder')},
        'input_sha256':signature,'source_files':plan}
    frozen=run/'frozen_manifest.json'
    if frozen.is_file():
        old=load(frozen)
        if any(old.get(k)!=v for k,v in base.items()):raise SystemExit('Frozen inputs or code changed; use a new run directory')
        return old
    base['frozen_utc']=now(); atomic_json(frozen,base); return base


def strict_source(r, run):
    raw=run/'raw'; proof=run/'source_proofs'/f"{r['pdb_id']}.json"; ident={k:r[k] for k in ('pdb_id','revision','source_url','expected_bytes')}
    old=load(proof) if proof.is_file() else {}; source_meta={**ident,'source_sha256':old.get('sha256') if old.get('identity')==ident else None}
    source,status,h,n=download_source(source_meta,raw)
    if not old or status!='reused_verified':
        with gzip.open(source,'rb') as stream:
            while stream.read(1<<20):pass
        atomic_json(proof,{'identity':ident,'relative_path':source.relative_to(run).as_posix(),'sha256':h,'bytes':n,'first_download_utc':now(),'gzip_crc':'pass'})
    block,ver=revision_from_cif(source,r['pdb_id']); st=gemmi.make_structure_from_block(block)
    if ver!=r['revision'] or not len(st):raise ValueError(f'Expected {r["revision"]}; observed {ver}, models={len(st)}')
    return source,h,block,st


def write_native_projection(path, source, assembly_ids, selected):
    """Lossless source-ASU atom selection (all models); assembly IDs stay in the sidecar."""
    atom=source.find_mmcif_category('_atom_site.'); tags=list(atom.tags); col={s.split('.',1)[1]:i for i,s in enumerate(tags)}
    doc=gemmi.cif.Document(); b=doc.add_new_block('source_asu_observation'); b.set_pair('_entry.id',source.find_value('_entry.id'))
    loop=b.init_loop('_atom_site.',list(col)); count=0; models=Counter()
    for original in atom:
        row=list(original)
        if row[col['label_asym_id']] not in selected:continue
        loop.add_row(row); count+=1; models[row[col['pdbx_PDB_model_num']]]+=1
    if not count:raise ValueError('No selected source atoms')
    path.parent.mkdir(parents=True,exist_ok=True); doc.write_file(str(path)); disk=gemmi.cif.read(str(path)).sole_block(); structure=gemmi.make_structure_from_block(disk)
    reread=disk.find_mmcif_category('_atom_site.')
    matches=(list(reread.tags)==tags and len(reread)==count and all(list(x)==list(y) for x,y in zip((x for x in atom if x[col['label_asym_id']] in selected),reread)))
    if not matches or sum(x.count_atom_sites() for x in structure)!=count:raise ValueError('Exact native ASU source atom readback failed')
    return {'path':str(path),'sha256':sha256(path),'bytes':path.stat().st_size,'evidence':'Exact atom_site row equality (coordinates, IDs, altloc, occupation, model). Source ASU positions; use original mmCIF and explicit assembly operators for expansion. Assembly membership and contacts in sidecar. Projection omits struct_conn and unselected ASU labels; never a complete deposited assembly.', 'assembly_ids_represented_by_sidecar':assembly_ids,'source_label_asym_ids':sorted(selected),'models_and_atom_counts':dict(models),'atoms_readback':count}


def atom_evidence(block, st, rules):
    alt={}; assembly_list=[]; source_edges=list(connection_edges(block)); metrics=Counter(); audit=[]
    for model in st:
        for aid in [a.name for a in st.assemblies] or [None]:
            instances=assembly_instances(st,aid,model) if aid else {f'asu.{k}':(None,k,v) for k,v in subchains(model).items()}
            assembly_list.append({'model':model.num,'assembly':aid,'instances':[{'instance':key,'label_asym':label,'label_entity':str(rs[0][1].entity_id),'auth_asym':sorted({a for a,r in rs}),'operator':op.name if op else None,'residues':len(rs)} for key,(op,label,rs) in sorted(instances.items())]})
        for label, rr in subchains(model).items():
            for auth,res in rr:
                at=[{'atom':a.name,'altloc':altloc(a),'occ':float(a.occ)} for a in res if not a.element.is_hydrogen and (altloc(a) or 0<a.occ<.9999)]
                if at:alt[f'{model.num}/{label}/{auth}/{res.seqid}/{res.name}']=at
        labels=subchains(model)
        for edge in source_edges:
            matches=[[{'auth_asym':a,'auth_seq':str(res.seqid),'entity':str(res.entity_id),'ccd':res.name} for a,res in labels.get(end['label_asym'],[]) if observed_endpoint(end,res,a)] for end in edge['ends']]
            mode='BOTH_UNIQUE' if all(len(x)==1 for x in matches) else 'NONUNIQUE_OR_UNOBSERVED'
            metrics[mode]+=1; audit.append({'model':model.num,'edge':edge,'endpoint_matches':matches,'status':mode})
    return {'assemblies_and_members_nonwater':assembly_list,'altloc_or_partial_positive_occupancy_atoms':alt,'source_connection_checks':audit,'connection_states':dict(metrics)}


def verify_files(run,result):
    for x in [result['source'],result['native_projection']]+[{'path':s['derivative_relative_path'],'sha256':s['derivative_sha256']} for s in result['index_systems'] if s.get('derivative_sha256')]:
        path=run/x['path'];
        if not path.is_file() or sha256(path)!=x['sha256']: return False
    return True


def one(r, run, release, commit, rules, number, retry=False):
    t=now(); p=r['pdb_id']; qc=run/'qc'/f'{p}.json'; old=load(qc) if qc.is_file() else None
    if old and old['status']!='FAIL' and not retry and verify_files(run,old):return old,True
    if old and (retry or not verify_files(run,old) if old['status']!='FAIL' else retry):
        prior=run/'entries'/f'{p.lower()}.json'; prior.unlink(missing_ok=True)
    obs=run/'native'/f'{p}.cif'; src,h,b,st=strict_source(r,run)
    candidate={'pdb_id':p,'candidate_no':str(number),'primary_stratum':r['m1_evidence']['test_stratum'] or 'unassigned','selection_rationale':r['m1_evidence']['why_and_test'] or 'Authorized fixed-revision nomination'}
    old_run=process_entry(candidate,{**r,'source_sha256':h},run,release,commit,role_rules=rules)
    if old_run.get('source_parse_status')!='pass':raise ValueError(f'Legacy entry processing did not pass source parsing: {old_run.get("candidate_outcome")}')
    systems=old_run['system_results']
    native=inspect_source(b,st,p,r['revision'],rules)
    ae=atom_evidence(b,st,rules)
    target=[x for x in native if x['role']=='target_ligand']
    # Save all nonwater source label coordinates. Contacts, model/assembly relationships, and roles are separate observations.
    source_labels={label for model in st for label in subchains(model)}
    proj=write_native_projection(obs,b,[a.name for a in st.assemblies],source_labels); proj['path']=obs.relative_to(run).as_posix()
    for s in systems:
        if s.get('derivative_relative_path') and s['readback_status']=='pass':
            q=gemmi.cif.read(str(run/s['derivative_relative_path'])).sole_block()
            s['representative_connection_rows']=len(q.find_mmcif_category('_struct_conn.'))
            s['representative_model_ids']=[m.num for m in gemmi.make_structure_from_block(q)]
    tech=Counter(s['coordinate_validation_status'] for s in systems)
    if old_run.get('plinder_annotation_status')=='technical_failure':route='INDEX_LOOKUP_TECHNICAL_FALLBACK'
    elif old_run.get('plinder_annotation_status')=='absent':route='SOURCE_NO_ANNOTATION'
    elif not old_run.get('plinder_revision_match'):route='SOURCE_REVISION_DIFFERS'
    elif any(k!='coordinate_identity_pass' for k in tech):route='INDEX_WITH_SOURCE_FALLBACK'
    else:route='INDEX_MATCHED_PLUS_SOURCE_EVIDENCE'
    unresolved_targets=sum(bool(x['unresolved_connection_ids']) or not x['contact_instances_4A'] for x in target)
    relevant=[x for x in native if x['role']=='pending' and any(m['kind']!='polymer_modification' for m in x['members'])]
    boundary=not target and not relevant and bool(native) and all(x['role']=='non_target_component' for x in native)
    why=[]
    if route in ('INDEX_LOOKUP_TECHNICAL_FALLBACK','SOURCE_REVISION_DIFFERS'):why.append(route)
    if any(k!='coordinate_identity_pass' for k in tech):why.append('INDEX_GEOMETRY_OR_COMPONENT_REVIEW')
    if len(st)>1:why.append('REPRESENTATIVE_MODEL_IS_NOT_COMPLETE_ENSEMBLE')
    if unresolved_targets:why.append('TARGET_CONTACT_OR_CONNECTION_UNRESOLVED')
    if not target and relevant:why.append('TARGET_ROLE_PENDING')
    if any(x['role']=='target_ligand' and any(m['polymer_type'] and 'deoxyribonucleotide' in m['polymer_type'] for m in x['members']) for x in native):why.append('DNA_PARTNER_STRANDS_NOT_GROUPED_AS_DUPLEX')
    if boundary: status='NO_TARGET'
    elif why or not target:status='REVIEW'
    else:status='PASS'
    result={'pdb_id':p,'revision':r['revision'],'status':status,'review_reasons':why,'source':{'path':src.relative_to(run).as_posix(),'bytes':src.stat().st_size,'sha256':h,'gzip_crc':'pass','entry_revision':'pass'},'source_model_ids':[x.num for x in st],'assembly_ids':[x.name for x in st.assemblies],'native_projection':proj,'route':route,'index_systems':systems,'index_reported':old_run.get('plinder_system_count'), 'plinder_source_revision':old_run.get('plinder_source_revision'),'native_roles':dict(Counter(x['role'] for x in native)),'source_target_models':sorted({x['model'] for x in target}),'source_target_observations':len(target),'candidate_alias_ledger':old_run.get('candidate_ledger',[]),'source_observations':native,'source_qc':ae,'coupled_role_statement':'Local M2 role evidence; PASS is pilot observation/geometry only, not intake QC or a Holo release.', 'started_utc':t,'finished_utc':now()}
    atomic_json(qc,result); return result,False


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',type=Path)
    parser.add_argument('--plinder-release',type=Path,required=True)
    parser.add_argument('--git-head',required=True)
    parser.add_argument('--batch-id',required=True)
    parser.add_argument('--inputs',type=Path,help='Directory with m2_nominations.json and optional source_role_rules.json')
    parser.add_argument('--entries',nargs='+',help='Limit this invocation to these IDs in the frozen manifest')
    parser.add_argument('--smoke',nargs='+',help='Run these selected IDs first; pause the rest after a technical failure')
    parser.add_argument('--retry',action='store_true'); args=parser.parse_args()
    run=args.run.resolve(); inp=args.inputs.resolve() if args.inputs else Path(__file__).resolve().parent
    for x in ('raw','entries','qc','source_proofs'):(run/x).mkdir(parents=True,exist_ok=True)
    f=freeze(run,inp,args)
    sources={r['pdb_id']:(i,r) for i,r in enumerate(f['source_files'],1)}
    wanted={p.upper() for p in args.entries} if args.entries else set(sources)
    smoke=list(dict.fromkeys(p.upper() for p in (args.smoke or [])))
    if not wanted<=sources.keys() or not set(smoke)<=wanted:
        parser.error('--entries and --smoke must select IDs in this invocation of the frozen list')
    role_file=inp/'source_role_rules.json'; rules=load(role_file) if role_file.is_file() else {}
    if not isinstance(rules,dict):
        parser.error('source_role_rules.json must map PDB IDs to role evidence')
    release=PlinderRelease(data_dir=args.plinder_release)
    selected=smoke+[p for p in sources if p in wanted and p not in smoke]
    smoke_failed=[]; invoked=[]
    for p in selected:
        if smoke_failed and p not in smoke:
            break
        i,r=sources[p]; invoked.append(p)
        print(f'START {p} {r["revision"]}',flush=True)
        try:
            item,reused=one(r,run,release,args.git_head,rules.get(p,{}),i,args.retry)
        except Exception as ex:
            item={'pdb_id':p,'revision':r['revision'],'status':'FAIL','failure':f'{type(ex).__name__}: {ex}','traceback':traceback.format_exc().splitlines()[-12:],'finished_utc':now()}; atomic_json(run/'qc'/f'{p}.json',item); reused=False
        if p in smoke and item['status']=='FAIL':smoke_failed.append(p)
        print(json.dumps({'pdb_id':p,'state':item['status'],'reused':reused,'route':item.get('route'),'native':item.get('native_roles'),'index':Counter(s['coordinate_validation_status'] for s in item.get('index_systems',[])),'failure':item.get('failure')}),flush=True)
    a=[load(run/'qc'/f"{x['pdb_id']}.json") for x in f['source_files'] if (run/'qc'/f"{x['pdb_id']}.json").is_file()]; display=[]; idx=[]
    for x in a:
        display.append({'pdb_id':x['pdb_id'],'revision':x['revision'],'status':x['status'],'route':x.get('route'),'models':len(x.get('source_model_ids',[])),'assemblies':len(x.get('assembly_ids',[])),'plinder_source_revision':x.get('plinder_source_revision'),'index_reported':x.get('index_reported'),'index_attempted':len(x.get('index_systems',[])),'index_readback':sum(s['readback_status']=='pass' for s in x.get('index_systems',[])),'target_observations':x.get('source_target_observations'),'target_models':';'.join(map(str,x.get('source_target_models',[]))),'review':'|'.join(x.get('review_reasons',[])),'qc_json':f"qc/{x['pdb_id']}.json",'raw':x.get('source',{}).get('path'),'native_structure':x.get('native_projection',{}).get('path'),'failure':x.get('failure')})
        for s in x.get('index_systems',[]):
            idx.append({'pdb_id':x['pdb_id'],'system_id':s['system_id'],'write':s['reconstruction_status'],'readback':s['readback_status'],'mapping':s['coordinate_validation_status'],'difference':';'.join(sorted({d['kind'] for d in s.get('coordinate_validation',{}).get('differences',[])})),'mapping_error':s.get('mapping_error',{}).get('kind'),'derivative':s.get('derivative_relative_path'),'failure':s.get('failure')})
    for name,rows in [('entry_results.csv',display),('system_results.csv',idx)]:
        if rows:
            with (run/name).open('w',encoding='utf8',newline='') as out:
                w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    result={'batch_id':f['batch_id'],'run_id':f['run_id'],'planned':len(f['source_files']),'processed':len(a),'source_verified':sum('source' in x for x in a),'raw_bytes':sum(x.get('source',{}).get('bytes',0) for x in a),'states':dict(Counter(x['status'] for x in a)),'routes':dict(Counter(x.get('route') for x in a)),'index_reported':sum(x.get('index_reported') or 0 for x in a),'index_attempted':len(idx),'index_readback':sum(x['readback']=='pass' for x in idx),'index_mapping':dict(Counter(x['mapping'] for x in idx)), 'index_representative_bytes':sum(s.get('derivative_bytes',0) for r in a for s in r.get('index_systems',[])),'native_observation_roles':dict(sum((Counter(x.get('native_roles',{})) for x in a),Counter())), 'native_projection_bytes':sum(x.get('native_projection',{}).get('bytes',0) for x in a),'start_utc':f['frozen_utc'],'finished_utc':now(),'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'total_run_bytes_pre_summary':sum(p.stat().st_size for p in run.rglob('*') if p.is_file()),'entry_results':display}
    result['invocation']={'selected_ids':selected,'processed_ids':invoked,'smoke_ids':smoke,'smoke_failed_ids':smoke_failed}
    atomic_json(run/'summary.json',result);print(json.dumps({k:v for k,v in result.items() if k!='entry_results'}),flush=True)
    if smoke_failed or any(x['pdb_id'] in wanted and x['status']=='FAIL' for x in a):raise SystemExit(2)


if __name__=='__main__':
    import sys

    if '--v1-manifest' in sys.argv:
        from holo_m2_v1 import deliver
        raise SystemExit(deliver(sys.argv[1:]))
    main()
