"""Register an existing V1 candidate manifest in an isolated APHCP registry; read-only on M2 assets."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import re
import socket
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))

import gemmi
from aphcp_control.registry import (RegistryDiskReserve, backup, connect, database, migrate,
                                    publish_registry, recover_empty_registry, require_backup_space, utc_now)
from aphcp_control.queries import pdb as pdb_query
from holo_candidate_batch import revision_from_cif
from holo_m3_decisions import resolve as resolve_decisions
from holo_m2_blind_structure import ContactContext, explicit_role_names, instances_for, normalize, role
from holo_role_evidence import assess, build_role_context, role_source
from holo_v1_export import ATOM, LABEL, MODEL, describe_instances, expanded_atoms, expanded_conn_context, export, output_label, rows, text

RULE='holo-m3-v1-staging-2'
STATE='M3_STAGING_UNPUBLISHED'
PDB=re.compile(r'^[A-Z0-9]{4}$')
KEY=re.compile(r'^[a-f0-9]{24}$')
BATCH=re.compile(r'^[a-zA-Z0-9][a-zA-Z0-9_-]{0,79}$')
RUN_ID=re.compile(r'^holo-m3-[a-z0-9-]+$')
QC_FIELDS=('coordinate_identity','experimental_fidelity','topology','partner_scope_qc')


def require(flag, reason):
    if not flag:raise ValueError(reason)


def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def read(path):return json.loads(Path(path).read_bytes())


def jsonl(path):
    with Path(path).open(encoding='utf8') as f:
        for line in f:yield json.loads(line)


def raw_json(x):return json.dumps(x, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def digest(x):return hashlib.sha256(raw_json(x).encode()).hexdigest()


def write(path,x):
    from holo_archive_runtime import durable_json
    durable_json(Path(path),x)


def read_only_path(path, allowed):
    p=Path(path).resolve(strict=True); require(any(p.is_relative_to(x) for x in allowed) and p.is_file(), 'asset_outside_frozen_runs');return p


def file_ref(path,roots,expected=None):
    p=read_only_path(path,roots); value=sha(p)
    require(expected is None or value==expected,'pinned_asset_sha256_differs')
    return {'path':str(p),'sha256':value,'bytes':p.stat().st_size}


def pool_rows(ref,chem, selected):
    def capture(path, k):
        found={}
        for x in jsonl(path):
            kk=k(x)
            if kk in selected:
                require(kk not in found,'duplicate_candidate_in_source_ledger');found[kk]=x
        return found
    key=lambda x:(str(x['batch']),x['pdb_id'],x['candidate_key'])
    base=capture(ref/'system_attribution.jsonl',key); a=capture(chem/'results/baseline_system_labels.jsonl',key)
    b=capture(chem/'results/supplemental_system_labels.jsonl',key)
    # Only a sealed historical increment may supersede baseline structure evidence.
    deltas={(d['pdb_id'],d['candidate_key']):d for d in read(ref/'historical_deltas.json')['deltas']
            if d.get('delta_intake',0)>0}
    return base, a|b, deltas


def historical_overlay(b, c, delta, runs):
    prior=(b or {}).get('historical_retest') or {}
    require(prior and all(delta.get(k)==v for k,v in prior.items()) and c, 'history_not_linked_to_baseline_and_chemistry')
    require(delta['revision']==b['revision'] and delta['source_sha256']==b['source_sha256'], 'history_changed_fixed_source')
    origin=Path(delta['result_root']).resolve(strict=True)
    require(origin.parent==runs and origin.name.startswith('holo-m2-') and c.get('evidence_layer')==delta['historical_layer'], 'history_outside_m2_evidence_root')
    require(all(c.get(k)==delta[k] for k in ('mapping_path','system_path','system_sha256','source_sha256','revision')), 'historical_chemistry_uses_different_artifact')
    discovery=read_only_path(c['discovery_path'],(origin,))
    evidence=origin/'system_results'/f"{b['pdb_id']}.json"
    proof=file_ref(evidence,(origin,)); found=[x for x in read(evidence) if x['candidate_key']==b['candidate_key']]
    require(len(found)==1,'history_candidate_missing_or_ambiguous'); h=found[0]
    require(h['revision']==b['revision'] and h['assembly']==b['assembly'] and h['representative_model']==b['model'] and
            sorted(h['partner_instances'])==sorted(b['member_instances']) and h['status']==delta['readback_after'] and
            h['role_qc']==delta['role_after'] and h['candidate_intake']==delta['candidate_intake_after'] and
            h['artifact']['sha256']==delta['system_sha256'] and origin/h['artifact']['path']==Path(delta['system_path']) and
            origin/h['mapping']==Path(delta['mapping_path']), 'history_result_and_seal_disagree')
    build=file_ref(origin/'build_inputs.json',(origin,)); pinned=read(build['path'])['code_files']['holo_v1_export.py']
    verifier=file_ref(origin/'bootstrap/holo_v1_export.py',(origin,),pinned)
    return b | {k:h[k] for k in QC_FIELDS} | {k:delta[k] for k in ('system_path','system_sha256','mapping_path','role_after')} | {
        'role_evidence_after':h['role_basis'],'intake_after':h['candidate_intake'],
        'discovery_path':str(discovery),'evidence_layer':delta['historical_layer'],'historical_result':proof,'historical_source_root':str(origin),'historical_build':build,'historical_verifier':verifier}


def pinned_exporter(ref):
    spec=importlib.util.spec_from_file_location('holo_frozen_verifier_'+ref['sha256'][:16],ref['path'])
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module.export


def memo(store, key, get):
    if key not in store:store[key]=get()
    return store[key]


def read_manifest(path, *, incremental=False):
    x=read(path); roster=x if isinstance(x,list) else x['review_candidates']
    require(roster and (isinstance(x,list) or x.get('review_count',len(roster))==len(roster)), 'empty_or_inconsistent_input_manifest')
    result=[]; ids=set()
    for n in roster:
        pb=n['pdb_id']; ck=n.get('candidate', n.get('candidate_key')); bid=str(n.get('batch','2000'))
        valid_batch = ('batch' in n and BATCH.fullmatch(bid)) if incremental else bid in {'200','2000'}
        require(PDB.fullmatch(pb) and KEY.fullmatch(ck) and valid_batch,'invalid_input_identity')
        if 'supersedes' in n:
            refs=n['supersedes']
            require(isinstance(refs,list) and refs and all(isinstance(x,str) and '/' in x for x in refs)
                    and len(refs)==len(set(refs)) and isinstance(n.get('decision_basis'),str)
                    and n['decision_basis'].strip(), 'invalid_decision_supersession')
        require((bid,pb,ck) not in ids,'duplicate_input_candidate'); ids.add((bid,pb,ck));result.append((bid,pb,ck,n))
    return result


def role_context(block, st, model, aid):
    title=text(block.find_value('_struct.title')); entity=role_source(block); names=explicit_role_names(block)
    polytype={text(z['entity_id']):text(z['type']) for z in rows(block,'_entity_poly.')}
    polylen=Counter(text(x['entity_id']) for x in rows(block,'_entity_poly_seq.'))
    ctx=ContactContext(instances_for(st,aid,model),{x for x,v in polytype.items() if v.startswith('polypeptide')},polylen)
    return title,names,polylen,ctx,entity,build_role_context(entity,names)


def source_role(observation, context):
    title,names,polylen,ctx,entity,aliases=context
    dist=ctx.contacts(observation['members']); n={**observation,'contact_instances_4A':sorted(dist)}
    outcome, evidence=role(n,title,names,polylen)
    addition=assess(n,outcome,title,names,polylen,entity,aliases)
    if addition:outcome,evidence=addition
    s=normalize(title); conflicts=[]
    for h in evidence.get('matches',[]):
        name=normalize(h['source_name']); name_re=r'\b'+re.escape(name)+r'\b'
        # Context belongs to the named partner, not to a different co-reported component.
        for m in re.finditer(name_re,s):
            left=s[max(0,m.start()-80):m.start()];right=s[m.end():m.end()+55]
            if re.search(r'\b(?:cryo protectant|cryoprotectant|detergent|(?:inactive )?cofactor|precipitant|buffer)(?: [a-z0-9]+){0,1} $',left) or re.match(
                    r'\s+(?:as (?:an? |the )?)?(?:cryo protectant|cryoprotectant|detergent|(?:inactive )?cofactor|precipitant|buffer)\b',right):
                conflicts.append({'kind':'named_component_has_competing_context','title':title,'matched_name':h['source_name']})
    return outcome,evidence,conflicts,dist,title


def context_gate(n):
    s=' '.join(normalize(h.get('text','')) for h in (n.get('role_evidence') or {}).get('matches',[]))
    # A catalog input may refer to unsuitable material. Such entries are counted but no asset, partner or role evidence is emitted.
    return bool(re.search(r'\b(?:exotoxin|botulinum neurotoxin|anthrax lethal factor|ricin toxin|saxitoxin biosynthesis)\b',s))


def check_nomination(n,b):
    for field in ('revision','assembly','source_sha256','system_sha256'):
        if n.get(field) is not None:require(str(n[field])==str(b[field]),f'nomination_{field}_mismatch')
    if n.get('source_model') is not None:require(int(n['source_model'])==b['model'],'nomination_model_mismatch')


def assess_candidate(key,n,b,c,roots,m2,cache,shared):
    batch,p,ck=key; work={'batch':batch,'pdb_id':p,'candidate_key':ck,'disposition':'REVIEW','reasons':[]}
    if context_gate(n):
        return work | {'disposition':'OUT_OF_SCOPE','reasons':['restricted_source_context_no_assets_staged']}
    if b is None:
        return work | {'disposition':'TECHNICAL_FAILURE','reasons':['not_found_in_frozen_candidate_evidence']}
    # Reconciliation supplies all primary facts; chemistry is joined by the same source identity.
    work.update(revision=b['revision'],assembly=b['assembly'],source_model=b['model'],m2_run=str(m2),evidence_layer=b['evidence_layer'],
                m2_route=b['route'], prior_role=b['role_after'],prior_program_intake=b['intake_after'],
                previous_qc={f:b.get(f) for f in QC_FIELDS})
    check_nomination(n,b)
    if context_gate({'role_evidence':b.get('role_evidence_after')}):
        return {'batch':batch,'pdb_id':p,'candidate_key':ck,'disposition':'OUT_OF_SCOPE','reasons':['restricted_source_context_no_assets_staged']}
    if not all(b.get(k) for k in ('source_path','system_path','mapping_path','system_sha256','discovery_path')):
        return work | {'reasons':['accepted_disk_structure_not_available']}
    get_ref=lambda path,expected=None:memo(shared,('file',str(path),expected),lambda:file_ref(path,roots,expected))
    raw=get_ref(b['source_path'],b['source_sha256']); cif=get_ref(b['system_path'],b['system_sha256']); mapping=get_ref(b['mapping_path'],b.get('mapping_sha256'))
    discovery=get_ref(b['discovery_path'],b.get('discovery_sha256'))
    refs={'source_mmcif':raw,'system_cif':cif,'source_mapping':mapping,'discovery':discovery}
    verifier=export
    if b.get('historical_result'):
        refs.update(historical_system_qc=b['historical_result'],historical_build=b['historical_build'],historical_verifier=b['historical_verifier'])
        m2=Path(b['historical_source_root'])
        verifier=memo(shared,('verifier',b['historical_verifier']['sha256']),lambda:pinned_exporter(b['historical_verifier']))
    if c:
        require(c['revision']==b['revision'] and c['source_sha256']==b['source_sha256'], 'chemistry_wrong_source_identity')
    ms=read(mapping['path'])
    proof_ref=n.get('frozen_source_proof') or ms.get('source_proof')
    if proof_ref:
        fp=get_ref(m2/proof_ref if not Path(proof_ref).is_absolute() else proof_ref)
        spr=memo(shared,('source_proof',fp['path']),lambda:read(fp['path'])); require(spr['status']==spr['gzip_crc']=='PASS' and spr['observed_revision']==b['revision'] and spr['sha256']==raw['sha256'] and spr['identity']['pdb_id']==p, 'source_proof_differs')
        refs['source_proof']=fp
    if cache and cache.get('refs_core')==refs:
        full=cache.get('files',{}).get('complete_assembly')
        if not full or get_ref(full['path'],full['sha256'])=={k:full[k] for k in ('path','sha256','bytes')}:
            return cache
    def read_raw():
        block,rev=revision_from_cif(Path(raw['path']),p)
        title=text(block.find_value('_struct.title'))
        gated=context_gate({'role_evidence':{'matches':[{'text':title}]}})
        return block,rev,title,None if gated else gemmi.make_structure_from_block(block)
    bsrc,rev,title,st=memo(shared,('struct',raw['path']),read_raw); require(rev==b['revision'],'source_model_revision_differs')
    if st is None:
        return {'batch':batch,'pdb_id':p,'candidate_key':ck,'disposition':'OUT_OF_SCOPE','reasons':['restricted_source_context_no_assets_staged']}
    model_num=b['model']; aid=b['assembly']; model=next((v for v in st if v.num==model_num),None)
    require(model is not None and aid is not None,'pinned_model_or_assembly_missing')
    ld=memo(shared,('discovery',discovery['path']),lambda:read(discovery['path'])); art=ms['artifact']
    for k,v in {'pdb_id':p,'revision':rev,'candidate_key':ck,'assembly':aid,'source_model':model_num,'source_sha256':raw['sha256']}.items():
        require(ms.get(k)==v,'mapping_identity_differs:'+k)
    require(art['sha256']==cif['sha256'] and art['bytes']==cif['bytes'] and art['model_ids']==[model_num], 'mapping_file_model_or_hash_differs')
    require(ld['pdb_id']==p and ld['revision']==rev,'discovery_entry_revision_differs')
    matches=[v for v in ld['native_observations'] if v['candidate_key']==ck]; require(len(matches)==1,'candidate_not_unique_in_discovery')
    obs=matches[0]; require(obs['assembly']==aid and obs['model']==model_num and obs['source_member_identity']==ms['partner_member_identities'],'discovery_and_mapping_partner_differ')
    require(set(b['member_instances'])=={m['instance'] for m in obs['members']}, 'reconciliation_and_mapping_instances_differ')
    ins, expanded=memo(shared,('instances',raw['path'],aid,model_num),lambda:describe_instances(st,aid,model))
    selected=ms['instances']; require(bool(art.get('output_instances'))==expanded and set(selected)<=set(ins), 'assembly_instance_expansion_differs')
    for k,v in selected.items():
        a=ins[k]; require(all(v.get(i)==a.get(i) for i in ('entity','label_asym','operator','auth_asym','matrix','vector','output_label_asym')), 'assembly_copy_lineage_differs')
    labels={output_label(v) for v in selected.values()}
    readback_kwargs={'expanded':ins if expanded else None,'verify':True}
    # Only the current verifier consumes this preparation; pinned historical exporters
    # keep their original preparation and readback behavior.
    if verifier is export:
        projected_count = art['atom_count'] + art['omitted_assembly_atom_rows'] if expanded else None

        def prepare_readback():
            original=memo(shared,('atom_rows',raw['path'],raw['sha256']),lambda:rows(bsrc,ATOM))
            if not expanded:
                return original,None
            # Older exports numbered atoms across all models. The sealed mapping's
            # projection size distinguishes that format from the model-scoped export.
            copies = Counter(v['label_asym'] for v in ins.values())
            selected_rows = [row for row in original if text(row[MODEL]) == str(model_num)]
            selected_count = sum(copies[text(row[LABEL])] for row in selected_rows)
            source_rows = selected_rows if projected_count == selected_count else original
            atoms=expanded_atoms(bsrc,ins,source_atoms=source_rows)
            require(len(atoms[0]) == projected_count, 'recorded_assembly_projection_size_differs')
            return atoms,expanded_conn_context(bsrc,ins,atoms,{str(model_num)})
        prepared, connections=memo(shared,('readback',RULE,raw['path'],raw['sha256'],model_num,aid,projected_count),prepare_readback)
        readback_kwargs.update(prepared=prepared,prepared_connections=connections)
    qc,atoms=verifier(Path(cif['path']),bsrc,p,aid,[model_num],labels,**readback_kwargs)
    require(qc['sha256']==cif['sha256'] and raw_json(qc['model_members'])==raw_json(art['model_members']) and qc['atom_count']==art['atom_count'],'derived_real_readback_differs')
    names={k:output_label(v) for k,v in selected.items()}
    protein=ms['binding_protein_instances']; members=obs['members']; req_instances=set(protein)|{m['instance'] for m in members}
    require(req_instances<=selected.keys() and not set(protein)&{m['instance'] for m in members},'missing_bound_or_partner_instance')
    ctx=memo(shared,('role',raw['path'],aid,model_num),lambda:role_context(bsrc,st,model,aid))
    result_role, role_evidence, conflicts, contacts,title=source_role(obs,ctx)
    require(set(contacts)==set(protein),'binding_chain_set_differs_on_fixed_source')
    residues={(text(x['label_asym_id']),text(x['auth_seq_id'])+text(x.get('pdbx_PDB_ins_code')),text(x['label_comp_id'])) for x in atoms}
    for m in members:
        require(all((names[m['instance']],r['auth_seq'],r['ccd']) in residues for r in m['residues']), 'output_target_residue_identity_differs')
    cls=c['classification'] if c else {'primary_class':'UNKNOWN','chemical_subclass':[],'annotation_status':'UNKNOWN',
                                          'status_reasons':['no_matching_frozen_annotation']}
    issues=(c or {}).get('partner_identity_evidence',{}).get('issue',[])
    partial_pos=cls.get('partial_observation',[])
    # Independent partner extent check works even without a type annotation.
    seq=defaultdict(set)
    for z in rows(bsrc,'_entity_poly_seq.'):
        seq[text(z['entity_id'])].add(int(text(z['num'])))
    declared_scope=[]
    for i,m in enumerate(members):
        if m['kind']=='polymer':
            positions={r['label_seq'] for r in m['residues'] if isinstance(r.get('label_seq'),int)}
            declared_scope.append({'instance':m['instance'],'declared_positions':len(seq[m['entity']]),'observed_positions':len(positions),
                                   'unobserved_positions':sorted(seq[m['entity']]-positions)})
    reasons=[]
    if any(b.get(f)!='PASS' for f in QC_FIELDS):reasons.append('prior_structure_or_scope_qc_not_passed')
    if result_role!='PASS' or b['role_after']!='PASS' or conflicts:reasons.append('target_role_evidence_insufficient_or_conflicting')
    partial_scope=bool(partial_pos or any(z['unobserved_positions'] for z in declared_scope))
    if partial_scope:reasons.append('deposited_partner_sequence_not_fully_observed')
    if issues:reasons.append('complete_partner_membership_requires_confirmation')
    dep_model_description=[text(bsrc.find_value(f)) for f in ('_struct.pdbx_model_details','_pdbx_nmr_representative.selection_criteria') if text(bsrc.find_value(f))]
    if re.search(r'\bminimized average structure\b',normalize(title)):
        reasons.append('only_deposited_nmr_average_needs_scope_decision')
    release_date=next((text(x['revision_date']) for x in rows(bsrc,'_pdbx_audit_revision_history.') if text(x.get('major_revision'))=='1' and text(x.get('minor_revision'))=='0'),None)
    methods=[text(z['method']) for z in rows(bsrc,'_exptl.')]
    assembly_path=m2/'assemblies'/p/f'a{aid}-m{model_num}.cif'
    def check_assembly():
        if not assembly_path.is_file():return None
        q=get_ref(assembly_path); d,aa=verifier(assembly_path,bsrc,p,aid,[model_num],{output_label(v) for v in ins.values()},
                                                  complete=True,**readback_kwargs)
        require(q['sha256']==d['sha256'],'complete_assembly_hash_differs')
        return q | {'model_scope':d['model_ids'],'observed_chains':len(d['model_members'][model_num]),'atoms':d['atom_count'],
                  'readback':'PASS','verifier_sha256':b.get('historical_verifier',{}).get('sha256')}
    whole=memo(shared,('complete',str(assembly_path),raw['path'],raw['sha256'],model_num,aid,
                       b.get('historical_verifier',{}).get('sha256',RULE)),check_assembly)
    files=refs|({'complete_assembly':whole} if whole else {})
    member_id=ms['partner_member_identities']
    identity={'pdb_id':p,'revision':rev,'source_sha256':raw['sha256'],'assembly':aid,'model_scope':[model_num],
              'binding_instances':sorted((k,ins[k]['label_asym'],ins[k]['operator']) for k in protein),
              'target_residue_instances':sorted(member_id),'altloc_scope':qc['altloc_policy']}
    system_id='holo:'+p+':'+digest(identity)[:28]
    result=work | {'system_id':system_id,'disposition':'REGISTERED' if not reasons else 'REVIEW','reasons':reasons,'identity':identity,
                  'title':title,'first_release_on':release_date,'methods':methods,'files':files,'refs_core':refs,
                  'output_instances':selected,'assembly_instances':ins,'direct_binding_protein_instances':sorted(protein),
                  'partner_instances':[{k:m[k] for k in ('instance','label_asym','entity','operator','kind','description')} for m in members],
                  'partner_scope_evidence':{'deposited_vs_observed':declared_scope,'concerns':issues},
                  'type':cls,'type_evidence':{k:c[k] for k in ('method_version','chosen_stage','chemical_identity','sources_of_labels') if k in c} if c else {},
                  'role':{'status': 'PASS' if result_role=='PASS' and not conflicts else 'REVIEW', 'source_evidence':role_evidence,'conflicts':conflicts,
                          'contact_distances_A':contacts,'reconciliation_role_evidence':b.get('role_evidence_after')},
                  'qc':{'coordinate_identity':'PASS','experimental_fidelity':'PASS','topology':'PASS',
                        'partner_scope':'REVIEW' if partial_scope or issues else 'PASS','atoms':qc['atom_count'],'source_model_scope':[model_num],
                        'model_representation_note':dep_model_description,'all_deposited_models':[z.num for z in st],
                        'named_altloc_atoms':qc['named_altloc_atoms'],'altloc_policy':qc['altloc_policy'],
                        'omitted_models':sorted(m.num for m in st if m.num != model_num),'explicit_connections':len(qc['connections']),
                        'omitted_touching_connections':len(qc['omitted_touching_connections']),'base_pairs':len(qc['base_pairs']),
                        'complete_assembly_available':whole is not None,'historical_verifier_sha256':b.get('historical_verifier',{}).get('sha256'),'transform_readback':qc.get('coordinate_transform_readback')}}
    return result



def review_group(group):
    shared={}; results=[]
    for position,key,n,b,c,delta,roots,m2,cached in group:
        p,ck=key[1:]
        try:
            if delta:
                b=historical_overlay(b,c,delta,m2.parent)
                roots=(*roots,Path(b['historical_source_root']))
            ans=assess_candidate(key,n,b,c,roots,m2,cached,shared)
        except (ValueError, OSError, KeyError, RuntimeError) as e:
            ans={'batch':key[0],'pdb_id':p,'candidate_key':ck,'disposition':'TECHNICAL_FAILURE','reasons':[str(e)],
                 'revision':(b or {}).get('revision')}
        reused=ans is cached
        if n.get('supersedes') and ans.get('system_id'):
            ans=ans | {'supersedes': sorted(n['supersedes']), 'decision_basis': n['decision_basis']}
        results.append((position,ans,reused))
    return results

def canonical_id(prefix,data):return f'{prefix}:{digest(data)[:28]}'


def index_review_results(registry_root, stamp, reports, seal, path_to_results, begin, elapsed):
    dbpath=database(registry_root); migrate(dbpath)
    label=registry_root.name; task=label+'-register'; run=label+'-run'; cv=RULE+':'+digest(seal)[:16]
    at=utc_now(); tally=Counter(r['disposition'] for r in reports)
    out_hash=sha(path_to_results)
    with connect(dbpath) as db:
        prev=db.execute('SELECT config_hash,output_manifest_hash FROM run WHERE id=?',(run,)).fetchone()
        if prev:
            require(prev['config_hash']==digest(seal) and prev['output_manifest_hash']==out_hash,'staged_run_inputs_or_results_changed')
            return row_counts(db),True
        for result in reports:
            if result.get('supersedes'):
                sid=result['system_id']
                old=[(rid,path,r) for rid,path,r in recorded_assessments(db, registry_root.parent.parent/'runs', system=sid)]
                resolve_decisions([*old,(run,str(path_to_results),result)])
        db.execute('INSERT INTO config_version VALUES (?,?,?)',(cv,digest(seal),str(stamp)))
        db.execute('''INSERT INTO batch(id,module_id,status,source_release,input_manifest,config_version,candidate_count,next_step)
                      VALUES (?,'holo',?,'wwPDB-fixed-mmCIF',?,?,?,?)''',(label,STATE,seal['manifest'],cv,len(reports),'暂存，等待发布授权'))
        db.execute('INSERT INTO task VALUES (?,?,?,?,?)',(task,label,'Holo V1 批次复核与登记',task,STATE))
        db.execute('''INSERT INTO run(id,batch_id,task_id,status,git_commit,code_version,config_version,config_hash,
                   source_release,input_manifest_hash,output_manifest_hash,environment_lock_hash,host,start_time_utc,end_time_utc,
                   elapsed_seconds,input_records,success_records,failed_records,uncertain_records,bytes_read,bytes_written,
                   downloaded_bytes,error_summary,report_path) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',
                   (run,label,task,STATE,seal['workspace_git_head'],seal['code_sha256'],cv,digest(seal),'wwPDB-fixed-mmCIF',
                    seal['manifest_sha256'],out_hash,digest(seal['runtime']),socket.gethostname(),begin,at,elapsed,len(reports),
                    tally['REGISTERED'],tally['TECHNICAL_FAILURE'],tally['REVIEW'],sum(f['bytes'] for x in reports for f in x.get('files',{}).values()),
                    path_to_results.stat().st_size,0,raw_json(tally),str(path_to_results)))
        def item(id,v,parent,kind,name,cc=None):
            db.execute('INSERT OR IGNORE INTO work_item(id,entry_version_id,parent_id,kind,label,component_id) VALUES (?,?,?,?,?,?)',
                       (id,v,parent,kind,name,cc))
        def stage(kind,id,p,key,status,codes=None,out=None,details=None):
            db.execute('INSERT INTO run_item(run_id,object_type,object_id,pdb_id,stage,status,error_code,error_message,output_uri,updated_at) VALUES (?,?,?,?,?,?,?,?,?,?)',
                       (run,kind,id,p,key,status,','.join(codes) if codes else None,raw_json(details) if details else None,out,at))
        seen_entry=set(); seen_assembly=set(); seen_system={}
        for row in reports:
            p=row['pdb_id']; db.execute('INSERT OR IGNORE INTO pdb_entry VALUES (?,?,?)',(p,row.get('title'),row.get('first_release_on')))
            outcome=row['disposition']; sid=row.get('system_id'); files=row.get('files',{})
            if not sid:
                stage('entry','candidate:'+row['candidate_key'],p,'M3_REVIEW',outcome,row['reasons'])
                continue
            ver=row['revision']; version=f'wwpdb:{p}:{ver}'; assembly=canonical_id('assembly',[version,row['assembly'],[row['source_model']]])
            db.execute('INSERT OR IGNORE INTO entry_version VALUES (?,?,?,?,?,?)',(version,p,'wwPDB','fixed-archive',p,ver))
            if version not in seen_entry:stage('entry',version,p,'FIXED_SOURCE','PASS',out=files['source_mmcif']['path']);seen_entry.add(version)
            item(assembly,version,None,'assembly',str(row['assembly']))
            if assembly not in seen_assembly:
                stage('assembly',assembly,p,'ASSEMBLY_IDENTITY','PASS',out=files.get('complete_assembly',{}).get('path'))
                seen_assembly.add(assembly)
            chain_ids={}
            for instance,data in sorted(row['assembly_instances'].items()):
                cid=canonical_id('chain',[assembly,instance,data['label_asym'],data['operator']]);chain_ids[instance]=cid
                item(cid,version,assembly,'chain',instance)
            for k,a in files.items():
                # Type and content identify unique assets; the per-observation receipt preserves each original locator.
                ak=canonical_id(k,[version,a['sha256']])
                db.execute('INSERT OR IGNORE INTO asset_file(id,entry_version_id,uri,sha256,bytes,status) VALUES (?,?,?,?,?,?)',
                           (ak,version,a['path'],a['sha256'],a['bytes'],'REUSED_VERIFIED'))
            ligand=canonical_id('ligand',[assembly,row['identity']['target_residue_instances']])
            ccs={x[-1] for x in row['identity']['target_residue_instances']}
            item(ligand,version,assembly,'ligand',';'.join(sorted(x['instance'] for x in row['partner_instances'])),next(iter(ccs)) if len(ccs)==1 else None)
            if sid in seen_system:
                # Aliases are retained in the detailed manifest; a distinct observation never shares the same key.
                require(seen_system[sid]['identity']==row['identity'] and seen_system[sid]['disposition']==outcome and
                        seen_system[sid]['files']['system_cif']['sha256']==files['system_cif']['sha256'],
                        'canonical_observation_has_conflicting_intake_or_artifact')
                stage('system_alias','candidate:'+row['candidate_key'],p,'CANONICAL_ALIAS',outcome,out=sid)
                continue
            seen_system[sid]=row
            item(sid,version,assembly,'system',row['candidate_key'])
            expected_members={(ligand,'ligand')}|{(chain_ids[c],'protein') for c in row['direct_binding_protein_instances']}
            stored={tuple(x) for x in db.execute('SELECT member_id,role FROM system_member WHERE system_id=?',(sid,))}
            require(not stored or stored==expected_members,'existing_canonical_system_member_identity_differs')
            db.executemany('INSERT OR IGNORE INTO system_member VALUES (?,?,?)',
                           [(sid,m,role) for m,role in sorted(expected_members)])
            stage('system',sid,p,'QC','PASS',out=files['system_cif']['path'],details=row['qc'])
            stage('system',sid,p,'ROLE',row['role']['status'],details=row['role'])
            stage('system',sid,p,'CHEMICAL_TYPE',row['type']['annotation_status'],details=row['type'])
            stage('system',sid,p,'M3_REVIEW',outcome,row['reasons'],out=files['system_cif']['path'])
        require(not db.execute('PRAGMA foreign_key_check').fetchall(),'staging_foreign_keys_failed')
        require(db.execute('PRAGMA integrity_check').fetchone()[0]=='ok','staging_registry_integrity_failed')
        require(db.execute('SELECT COUNT(*) FROM release').fetchone()[0]==0,'staging_registry_has_release')
        return row_counts(db),False


def row_counts(db):
    return {tbl:db.execute(f'SELECT COUNT(*) FROM {tbl}').fetchone()[0] for tbl in (
        'batch','task','run','pdb_entry','entry_version','asset_file','work_item','system_member','run_item','release')}


def frozen_registry_root(root, run, *, create=False):
    require(RUN_ID.fullmatch(run.name) and run.parent==root/'runs' and run.is_dir(),'not_an_isolated_m3_run')
    batch=root/'batches'/run.name
    require(batch.resolve()==batch and database(batch)!=database(root),'cannot_write_production_registry')
    target=database(batch)
    require(target.resolve()==target,'registry_path_is_a_symlink')
    if create:batch.mkdir(parents=True,exist_ok=True)
    return batch


def seed_registry(root,current,args):
    if not args.seed_run:return None
    source=args.seed_run.resolve(strict=True)
    require(source!=current,'cannot_seed_from_current_run')
    db=database(frozen_registry_root(root,source))
    with connect(db,readonly=True) as con:
        require(con.execute('SELECT COUNT(*) FROM release').fetchone()[0]==0,'seed_registry_not_unpublished')
        row=con.execute('SELECT output_manifest_hash FROM run WHERE id=?',(source.name+'-run',)).fetchone()
        if row:
            require(sha(source/'review_results.json') == row[0], 'seed_run_manifest_mismatch')
        else:
            for label in (source.name+'-observations', source.name):
                record = con.execute('SELECT config_hash,error_summary FROM run WHERE id=? AND status=?',
                                     (label, 'OBSERVATIONS_STAGED')).fetchone()
                if record:
                    require(digest(read(source/'frozen_catalog.json')) == record[0] and
                            read(source/'catalog_receipt.json') == json.loads(record[1]), 'seed_catalog_manifest_mismatch')
                    break
            require(record, 'seed_run_manifest_missing')
    return db.resolve(),sha(db)


def clone_sealed_registry(source, target):
    """Seed an isolated M3 DB from an immutable Run; never hard-link a writable DB."""
    source, target = Path(source), Path(target)
    if sys.platform != 'linux' or any(Path(f'{source}{suffix}').exists() for suffix in ('-wal', '-shm', '-journal')):
        return backup(source, target)
    target.parent.mkdir(parents=True, exist_ok=True)
    require_backup_space(source,target)
    with tempfile.TemporaryDirectory(prefix='.registry-seed-', dir=target.parent) as temp:
        staged=Path(temp)/target.name
        # Exclusive reflink is an optimization. SQLite's consistent backup remains the
        # portable fallback, including filesystems that do not support copy-on-write.
        copied=subprocess.run(['cp', '--reflink=always', '--', str(source), str(staged)],
                              capture_output=True, check=False)
        if copied.returncode:
            return backup(source, target)
        require(not any(Path(f'{source}{suffix}').exists() for suffix in ('-wal','-shm','-journal')),
                'sealed_catalog_became_writable_during_clone')
        with connect(staged,readonly=True) as saved:
            # Reflinks have a fresh page cache. mmap enables read-ahead during
            # the same full quick_check instead of reading one SQLite page at a time.
            saved.execute('PRAGMA mmap_size=536870912')
            require(saved.execute('PRAGMA quick_check').fetchone()[0]=='ok','cloned_catalog_integrity_failed')
        # Atomic, no replacement of an already staged destination.
        publish_registry(staged,target)
    return target


def incremental_inputs(root, inp, keys):
    spec=read(inp); require(spec.get('schema')=='aphcp-m3-increment-v1','unsupported_increment_manifest')
    evidence_ref=spec['candidate_evidence']; chem_ref=spec['chemistry']
    evidence_path=Path(evidence_ref['path']).resolve(); chem_path=Path(chem_ref['path']).resolve()
    file_ref(evidence_path,(root/'runs',),evidence_ref['sha256']); file_ref(chem_path,(root/'runs',),chem_ref['sha256'])
    bundle=read(evidence_path)
    source=Path(bundle['source_run']).resolve(strict=True)
    require(source.parent==root/'runs' and source.name.startswith('holo-m2-') and source.is_dir(),'invalid_increment_source_run')
    def keyed(rows):
        items={(x['batch'],x['pdb_id'],x['candidate_key']):x for x in rows}
        require(len(items)==len(rows),'duplicate_increment_evidence')
        return items
    base=keyed(bundle['candidates']); annotations=keyed(list(jsonl(chem_path)))
    require(BATCH.fullmatch(bundle['batch']) and all(k[0]==bundle['batch'] for k in base)
            and len(bundle['entries'])==bundle['expected_entries']==len(set(bundle['entries']))
            and len(base)==bundle['expected_candidates'] and keys<=base.keys()
            and all(k[1] in bundle['entries'] for k in base) and base.keys()==annotations.keys(),
            'increment_evidence_does_not_cover_frozen_input')
    return (source,),{bundle['batch']:source},base,annotations,{},(evidence_path,chem_path)


def timed_review_group(group):
    started=time.monotonic(); answers=review_group(group)
    return answers,round(time.monotonic()-started,3)


def review_groups(groups,jobs):
    work=[]
    for (batch,pdb), group in groups.items():
        slices=4 if jobs>=4 and len(group)>=256 else 1
        work.extend(((batch,pdb),group[i::slices]) for i in range(slices))
    started=time.monotonic(); pairs=[]; per_entry={}
    def add(results):
        for (key,_),(answers,elapsed) in zip(work,results,strict=True):
            pairs.extend(answers)
            info=per_entry.setdefault(key,{'batch':key[0],'pdb_id':key[1], 'candidates':0,
                                           'reused':0,'slice_elapsed_s':[]})
            info['candidates']+=len(answers); info['reused']+=sum(reused for _,_,reused in answers)
            info['slice_elapsed_s'].append(elapsed)
    if len(work)==1 or jobs==1:
        add(map(timed_review_group,(group for _,group in work)))
    else:
        with ProcessPoolExecutor(max_workers=min(jobs,len(work))) as pool:
            add(pool.map(timed_review_group,(group for _,group in work)))
    return sorted(pairs), {'review_wall_s':round(time.monotonic()-started,3),
                          'entries':[row for _,row in sorted(per_entry.items())]}


def execute(args):
    start=time.monotonic(); begin=utc_now(); root=args.root.resolve(); r=args.run.resolve(); ig=frozen_registry_root(root,r,create=not args.review_only)
    inp=args.manifest.resolve(); config=read(inp)
    incremental=isinstance(config,dict) and 'schema' in config
    manifest=read_manifest(inp,incremental=incremental); keys={(a,b,c) for a,b,c,_ in manifest}
    if incremental:
        roots,source_by_batch,base,annotations,deltas,ledger_paths=incremental_inputs(root,inp,keys)
    else:
        old=r.parent/'holo-m2-final2000-20261122a'; ref=r.parent/'holo-m2-full-reconciliation-20261122a'; chem=r.parent/'holo-m2-chemtype-20261122a'
        roots=tuple(z.resolve() for z in (old,ref,r.parent/'holo-m2-blind200-20261122a'))
        source_by_batch={'2000':roots[0],'200':roots[2]}
        ledger_paths=(ref/'system_attribution.jsonl',ref/'historical_deltas.json',chem/'results/baseline_system_labels.jsonl',chem/'results/supplemental_system_labels.jsonl')
        base,annotations,deltas=pool_rows(ref,chem,keys)
    deps=Path(__file__).resolve().parent
    runtime={k:importlib.metadata.version(k) for k in ['gemmi','numpy','scipy']} | {'python':sys.version.split()[0]}
    seed=seed_registry(root,r,args)
    seal={'stage':STATE,'rule':RULE,'manifest':str(inp),'manifest_sha256':sha(inp),'code_sha256':sha(__file__),
          'workspace_git_head':args.git_head,'runtime':runtime,'dependencies':{
          str(x.relative_to(deps)):sha(x) for x in sorted([*deps.glob('holo*.py'),*deps.glob('src/aphcp_control/*.py'),*deps.glob('migrations/*.sql')])},
          'inputs':{str(q):sha(q) for q in ledger_paths}}
    if incremental:
        build=read(next(iter(source_by_batch.values()))/'build_inputs.json')
        seal['incremental_source_code']={k:build[k] for k in ('algorithm','fingerprint','code_files')}
    if seed:seal['seed_registry']={'path':str(seed[0]),'sha256':seed[1]}
    seal_path=r/'frozen_staging.json'; data=r/'review_results.json'; current=read(data) if data.exists() else []
    if seal_path.exists():require(read(seal_path)==seal,'staging_input_or_code_changed_use_a_new_run')
    by_old={(x['batch'],x['pdb_id'],x['candidate_key']):x for x in current}
    groups=defaultdict(list)
    for position,(batch,p,k,n) in enumerate(manifest):
        key=(batch,p,k); b=base.get(key)
        delta=deltas.get((p,k)) if b and b.get('historical_retest') else None
        groups[(batch,p)].append((position,key,n,b,annotations.get(key),delta,roots,
                                  source_by_batch[batch],by_old.get(key)))
    pairs,timing=review_groups(groups,args.jobs)
    results=[x for _,x,_ in pairs]; reused=sum(y for _,_,y in pairs)
    counts=Counter(x['disposition'] for x in results); assert sum(counts.values())==len(manifest)
    if args.review_only:
        write(r/'diagnostic_review.json',results);print(raw_json({'dispositions':counts,'reasons':{x['pdb_id']:x['reasons'] for x in results if x['reasons']}}));return
    # After registration the run and its serialized evidence are stable; never overwrite a different registered result.
    if current and results!=current:raise ValueError('Immutable staging review differs; keep prior result and use a new run')
    if not current:write(data,results)
    if not seal_path.exists():write(seal_path,seal)
    clock=time.monotonic()
    if seed:
        recover_empty_registry(database(ig),(r/'registration_receipt.json', r/'phase_receipt.json',r/'catalog_receipt.json'))
        if not database(ig).exists(): clone_sealed_registry(seed[0],database(ig))
    counts_db,db_reuse=index_review_results(ig,seal_path,results,seal,data,begin,round(time.monotonic()-start,3))
    timing['registry_wall_s']=round(time.monotonic()-clock,3)
    write(r/('rerun_review_timing.json' if db_reuse else 'per_pdb_m3_timing.json'),timing)
    states={'input_candidate_count':len(manifest),'input_entry_count':len({x[1] for x in manifest}),'dispositions':counts,
            'retained_review_systems':sum(x['disposition']=='REVIEW' and 'system_id' in x for x in results),
            'unique_registered_systems':len({x['system_id'] for x in results if x['disposition']=='REGISTERED'}),
            'readback_files_reused_on_disk':sum(bool(x.get('system_id')) for x in results),
            'whole_assembly_verified':sum(bool(x.get('files',{}).get('complete_assembly')) for x in results),
            'prior_review_checks_reused':reused,'database_registration_reused':db_reuse,'coordinate_downloads':0,'coordinate_writes':0,
            'analysis_seconds':round(time.monotonic()-start,3),'registry_counts':counts_db,'registry_path':str(database(ig)),
            'results_path':str(data),'results_sha256':sha(data),'manifest_sha256':seal['manifest_sha256'],'status':STATE}
    write(r/('rerun_receipt.json' if db_reuse else 'registration_receipt.json'),states);print(raw_json(states))


def recorded_assessments(db, runs, *, system=None, pdb=None):
    criteria, value=('object_id',system) if system else ('pdb_id',pdb) if pdb else (None,None)
    suffix=f' AND run_item.{criteria}=?' if criteria else ''
    matched=db.execute("""SELECT DISTINCT run.id,run.report_path,run.output_manifest_hash FROM run
        JOIN run_item ON run.id=run_item.run_id WHERE run_item.stage='M3_REVIEW'"""+suffix,
        (value,) if criteria else ()).fetchall()
    for rec in matched:
        run_id=rec['id']; f=read_only_path(rec['report_path'],(runs,))
        require(f.name=='review_results.json' and f.parent.name+'-run'==run_id and sha(f)==rec['output_manifest_hash'], 'staging_result_file_changed')
        for x in read(f):
            if (system and x.get('system_id')!=system) or (pdb and x['pdb_id']!=pdb):continue
            sid=x.get('system_id')
            if sid:
                found=db.execute("SELECT status FROM run_item WHERE run_id=? AND object_id=? AND stage='M3_REVIEW'",(run_id,sid)).fetchone()
                require(found and found['status']==x['disposition'],'registry_and_science_receipt_disagree')
            yield run_id,str(f),x


def lookup(args):
    r=args.run.resolve(); root=r.parent.parent; batch=frozen_registry_root(root,r)
    with connect(database(batch),readonly=True) as db:
        value=args.system or args.pdb.upper()
        by_identity=defaultdict(list); verified=set()
        for run_id,f,x in recorded_assessments(db, root/'runs', system=args.system, pdb=args.pdb.upper() if args.pdb else None):
            sid=x.get('system_id')
            if sid and args.verify_files:
                for asset in x['files'].values():
                    proof=(asset['path'],asset['sha256'])
                    if proof not in verified:
                        file_ref(proof[0],(root/'runs',),proof[1]); verified.add(proof)
            ident=sid or ('candidate',x['batch'],x['pdb_id'],x['candidate_key'])
            by_identity[ident].append((run_id,f,x))
        require(by_identity,'staging_record_not_found')
        rows_=[resolve_decisions(group) for _,group in sorted(by_identity.items(),key=lambda x:str(x[0]))]
        if args.asset:
            require(len(rows_)==1 and args.asset in rows_[0].get('files',{}),'unique_system_asset_required')
            print(rows_[0]['files'][args.asset]['path']);return
        print(json.dumps({'status':STATE,'query':value,'registry':pdb_query(db,rows_[0]['pdb_id']),
                          'systems':rows_,'read_only_files_verified':args.verify_files},ensure_ascii=False))


def main():
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='command',required=True); run=sub.add_parser('register'); get=sub.add_parser('lookup')
    run.add_argument('--root',type=Path,required=True);run.add_argument('--run',type=Path,required=True);run.add_argument('--manifest',type=Path,required=True)
    run.add_argument('--seed-run',type=Path,help='Copy the registry of an earlier unpublished M3 run, once')
    run.add_argument('--git-head',required=True);run.add_argument('--jobs',type=int,default=8);run.add_argument('--review-only',action='store_true')
    get.add_argument('--run',type=Path,required=True);group=get.add_mutually_exclusive_group(required=True)
    group.add_argument('--pdb');group.add_argument('--system');get.add_argument('--verify-files',action='store_true')
    get.add_argument('--asset',choices=['system_cif','source_mmcif','source_mapping','complete_assembly'])
    args=p.parse_args()
    try:
        execute(args) if args.command=='register' else lookup(args)
    except RegistryDiskReserve as exc:
        raise SystemExit(73) from exc


if __name__=='__main__':main()
