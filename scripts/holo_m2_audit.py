"""M2 disk audits, readable native Assembly witnesses, and stage receipt."""
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json, math, os, sqlite3, sys
from pathlib import Path

import gemmi

run=Path(sys.argv[1]).resolve(); sys.path.insert(0,str(run/'bootstrap'))
from holo_candidate_batch import atomic_json, sha256, revision_from_cif
from holo_source_mapping import assembly_instances
from holo_m2_pilot import verify_files

read=lambda p:json.loads(p.read_text(encoding='utf8'))
state=read(run/'summary.json'); freeze=read(run/'frozen_manifest.json'); records=[read(run/'qc'/f"{p['pdb_id']}.json") for p in freeze['source_files']]
assert state['planned']==state['processed']==len(records)


def assembly_witness(q, block, model, instances, aid):
    atom=block.find_mmcif_category('_atom_site.'); names=[tag.split('.',1)[1] for tag in atom.tags]; pos={k:i for i,k in enumerate(names)}
    by_asym=defaultdict(list)
    for name,(op,label,_) in instances.items():by_asym[label].append((name,op))
    used={}; expected=[]; coords=('Cartn_x','Cartn_y','Cartn_z'); doc=gemmi.cif.Document(); b=doc.add_new_block('nonwater_assembly_witness'); b.set_pair('_entry.id',q['pdb_id']); loop=b.init_loop('_atom_site.',names)
    for atomrow in atom:
        row=list(atomrow); label=row[pos['label_asym_id']]
        if row[pos['pdbx_PDB_model_num']]!=str(model.num):continue
        if row[pos['label_comp_id']] in ('HOH','DOD'):continue
        for instance,operator in by_asym[label]:
            out=row.copy(); name='i'+instance.replace('.','_',1); used[name]=[instance,label,operator.name]; vals=operator.transform.apply(gemmi.Position(*(float(row[pos[k]]) for k in coords)))
            out[pos['label_asym_id']]=out[pos['auth_asym_id']]=name; out[pos['id']]=str(len(expected)+1)
            for k,v in zip(coords,vals):out[pos[k]]=f'{v:.7f}'
            loop.add_row(out); expected.append((out,(vals.x,vals.y,vals.z),row[pos['id']]))
    path=run/'native_assemblies'/f"{q['pdb_id']}-assembly{aid}-model{model.num}.cif"; path.parent.mkdir(exist_ok=True); doc.write_file(str(path))
    disk=gemmi.cif.read(str(path)).sole_block(); reread=disk.find_mmcif_category('_atom_site.'); gs=gemmi.make_structure_from_block(disk)
    if len(reread)!=len(expected) or len(gs)!=1 or gs[0].count_atom_sites()!=len(expected):raise ValueError('Assembly readback atom count mismatch')
    residual=0; lineage=[]
    for orig,new in zip(expected,reread):
        if orig[0]!=list(new):raise ValueError('Assembly exact row readback mismatch')
        delta=math.dist(orig[1], [float(new[pos[k]]) for k in coords]);residual=max(residual,delta); lineage.append(orig[2])
    if residual>1e-6 or sorted(used)!=sorted('i'+x.replace('.','_',1) for x in instances):raise ValueError('Assembly identity/precision mismatch')
    mapping={'pdb_id':q['pdb_id'],'source_model':model.num,'source_assembly':aid,'source_mmcif':q['source'],'cif':str(path.relative_to(run)),'sha256':sha256(path),'bytes':path.stat().st_size,'output_asym_to_source_instance_label_operator':used,'source_atom_site_id_by_output_atom_id':lineage,'atom_count':len(expected),'max_export_distance_A':residual,'scope':'All observed non-water instances for the first listed Assembly and model; new chain IDs; original atom occupancy / altloc retained; connection and contact graph separately in entry QC; not formal intake.'}
    loc=path.with_suffix('.json'); atomic_json(loc,mapping)
    return {k:v for k,v in mapping.items() if k not in ('source_atom_site_id_by_output_atom_id','output_asym_to_source_instance_label_operator')} | {'instance_count':len(used),'lineage':str(loc.relative_to(run))}


audit=[]; failure=[]; issues=Counter()
for q in records:
    pdb=q['pdb_id']; a={'pdb_id':pdb,'status':q['status']}; audit.append(a)
    if q['status']=='FAIL':failure.append(pdb); a['failure']=q['failure']; continue
    assert verify_files(run,q),pdb
    category=q['source_qc']; index=q['index_systems']; obs=q['source_observations']
    a.update(sha256_readback=True,revision=q['revision'],route=q['route'],native_role_counts=q['native_roles'],models=q['source_model_ids'],source_target_models=q['source_target_models'],assembly_ids=q['assembly_ids'],source_connections=category['connection_states'],review_reasons=q['review_reasons'])
    a['index_mapping_details']=[{'id':s['system_id'],'state':s['coordinate_validation_status'],'error':s.get('mapping_error'), 'failure':s.get('failure'), 'differences':[d['kind'] for d in s.get('coordinate_validation',{}).get('differences',[])], 'omitted_positive_altloc':s.get('coordinate_validation',{}).get('omitted_source_altloc_atoms'), 'occupancy_changes':s.get('coordinate_validation',{}).get('changed_occupancies'),'connections_in_derivative':s.get('representative_connection_rows'),'model_evidence':s.get('coordinate_validation',{}).get('model_evidence')} for s in index]
    issues.update(s['error']['kind'] for s in a['index_mapping_details'] if s['error'])
    issues.update(k for s in a['index_mapping_details'] for k in s['differences'])
    if q['route'].startswith('SOURCE_') and q['assembly_ids']:
        src=run/q['source']['path']; block,rev=revision_from_cif(src,pdb); st=gemmi.make_structure_from_block(block); aid=q['assembly_ids'][0]
        try:a['assembly_readback']=assembly_witness(q,block,st[0],assembly_instances(st,aid,st[0]),aid)
        except (ValueError,RuntimeError,KeyError) as err:a['assembly_review']=f'{type(err).__name__}: {err}'
    if pdb in ('5R42','5SBY','7OUN','7NM2','7NQC','7NMH','7K3O','7KX0','7N4V','7O1V','7OY8'):
        a['diagnostic_observations']={
          'model_1_samples':[{'role':r['role'],'identity':r['identity_status'],'members':[{'instance':m['instance'],'entity':m['entity'],'description':m['description'],'components':dict(Counter(x['ccd'] for x in m['residues']))} for m in r['members']], 'contacts':r['contact_instances_4A'],'connection_ids':[x['source']['id'] for x in r['connections']],'unresolved':r['unresolved_connection_ids']} for r in obs if r['model']==q['source_model_ids'][0]][:35],
          'positive_altloc_or_partial_occupancy_residue_count':len(category['altloc_or_partial_positive_occupancy_atoms']),
          'model_1_covalent_ends':[x for x in category['source_connection_checks'] if x['model']==q['source_model_ids'][0] and x['edge']['kind']=='covale'][:16]}
        if pdb=='5SBY':
            a['diagnostic_observations']['8BI_atoms']= {k:v for k,v in category['altloc_or_partial_positive_occupancy_atoms'].items() if k.endswith('/8BI')}
            a['diagnostic_observations']['8BI_mapped_atoms']=[x for s in index for group in s.get('coordinate_validation',{}).get('ligand_atom_map',{}).values() for x in group if x['ccd']=='8BI']
        if pdb in ('7NM2','7NQC'):
            a['diagnostic_observations']['all_source_model_observation_counts']={str(i):dict(Counter(r['role'] for r in obs if r['model']==i)) for i in q['source_model_ids']}

atomic_json(run/'diagnostic_audit.json',{'entries':audit,'issue_counts_by_index_system_or_difference':dict(issues),'witnesses_written':sum('assembly_readback' in x for x in audit)})
metadata=run/'run_receipt.json'; clock=datetime.now(timezone.utc); elapsed=(clock-datetime.fromisoformat(state['start_utc'])).total_seconds()
output={p.name:sha256(p) for p in [run/'summary.json',run/'entry_results.csv',run/'system_results.csv',run/'diagnostic_audit.json']}; out_hash=__import__('hashlib').sha256(json.dumps(output,sort_keys=True).encode()).hexdigest()
db_path=run.parent.parent/'registry'/'aphcp_registry.sqlite3'; batch=freeze['batch_id']; task=batch+'-pilot'; config=sha256(run/'frozen_manifest.json'); env=__import__('hashlib').sha256(json.dumps(freeze['environment'],sort_keys=True).encode()).hexdigest(); rid=freeze['run_id']
with sqlite3.connect(db_path) as db:
    db.execute('PRAGMA foreign_keys=ON'); db.execute('BEGIN IMMEDIATE')
    protected=['pdb_entry','entry_version','asset_file','work_item','system_member','run_item','release']; before={t:db.execute(f'SELECT COUNT(*) FROM {t}').fetchone()[0] for t in protected}
    result='SUCCEEDED' if not failure else 'PARTIAL'; stat='PILOT_COMPLETE' if not failure else 'PILOT_PARTIAL'
    db.execute('INSERT INTO batch(id,module_id,status,source_release,window_start,window_end,input_manifest,config_version,git_commit,candidate_count,next_step) VALUES (?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(id) DO UPDATE SET status=excluded.status,next_step=excluded.next_step',(batch,'holo',stat,'RCSB frozen coordinate versions; PLINDER 2026-09','2021-09-01','2021-10-01',str(run/'frozen_manifest.json'),'m2-202109-v1',freeze['git_head'],state['planned'],'Wait for project review'))
    db.execute('INSERT INTO task (id,batch_id,name,task_key,status) VALUES (?,?,?,?,?) ON CONFLICT(id) DO UPDATE SET status=excluded.status',(task,batch,'M2 fixed revision structure pilot',task,result))
    db.execute('''INSERT INTO run (id,batch_id,task_id,status,git_commit,code_version,config_version,config_hash,source_release,input_manifest_hash,output_manifest_hash,environment_lock_hash,host,start_time_utc,end_time_utc,elapsed_seconds,input_records,success_records,failed_records,uncertain_records,bytes_written,downloaded_bytes,peak_memory,error_summary,report_path) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(id) DO UPDATE SET status=excluded.status,end_time_utc=excluded.end_time_utc,output_manifest_hash=excluded.output_manifest_hash,failed_records=excluded.failed_records,uncertain_records=excluded.uncertain_records''',(rid,batch,task,result,freeze['git_head'],freeze['input_sha256']['holo_m2_pilot.py']+':'+sha256(Path(__file__)),'m2-202109-v1',config,'RCSB fixed, PLINDER 2026-09',config,out_hash,env,os.uname().nodename,state['start_utc'],clock.isoformat(),elapsed,state['planned'],state['states'].get('PASS',0)+state['states'].get('NO_TARGET',0),len(failure),state['states'].get('REVIEW',0),sum(p.stat().st_size for p in run.rglob('*') if p.is_file()),state['raw_bytes'],state['peak_rss_kib']*1024,','.join(failure) or None,str(metadata)))
    assert before=={t:db.execute(f'SELECT COUNT(*) FROM {t}').fetchone()[0] for t in protected}
    assert db.execute('SELECT status,output_manifest_hash FROM run WHERE id=?',(rid,)).fetchone()==(result,out_hash)
receipt={'batch_id':batch,'run_id':rid,'run_path':str(run),'status':stat,'start_utc':state['start_utc'],'finished_utc':clock.isoformat(),'wall_seconds_including_audit':round(elapsed,2),'raw_bytes':state['raw_bytes'],'run_disk_bytes':sum(p.stat().st_size for p in run.rglob('*') if p.is_file()),'native_witness_bytes':sum(x.get('assembly_readback',{}).get('bytes',0) for x in audit),'output_sha256':output,'registered':{'batch_id':batch,'task_id':task,'run_id':rid,'structural_registry_tables_unchanged':before},'technical_failure_ids':failure,'historical_files_downloaded':0}
atomic_json(metadata,receipt); print(json.dumps(receipt),flush=True)
