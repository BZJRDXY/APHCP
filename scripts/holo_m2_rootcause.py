"""Read-only full M2 reconciliation and role-only replay over sealed batches."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
import json
import re
from pathlib import Path
import time

from holo_candidate_batch import atomic_json, now, revision_from_cif, sha256
from holo_v1_export import rows, text
from holo_m2_blind_structure import explicit_role_names, review_admission_evidence
from holo_role_evidence import EXPERIMENT, VERSION, assess, build_role_context, normalize, token

RUNS = {'200': 'holo-m2-blind200-20261122a', '2000': 'holo-m2-final2000-20261122a'}
EXPECTED = {'200': (200, 5425, 4212, 91), '2000': (2000, 41430, 38786, 419)}
SCAN_VERSION = 'v4'
HISTORICAL = (('admission','holo-m2-admission-retest-qqcYW9VT', {'7NDX','7C5Y','7S4W'}),
              ('covalent_copy','holo-m2-covalent-copy-20261122-Z2Fod80a',{'6TO3','6X8H'}))


def load(path):
    with path.open(encoding='utf8') as f:
        return json.load(f)


def csvrows(path):
    with path.open(newline='', encoding='utf8') as f:
        return list(csv.DictReader(f))


def dump(value):
    return json.dumps(value, ensure_ascii=False, separators=(',', ':'))


def source_metadata(block):
    entity = {text(x['id']): {k: text(x.get(k)) for k in ('type','pdbx_description','details','pdbx_mutation')}
              for x in rows(block, '_entity.')}
    exp = {k: [text(x) for x in block.find_values(k) if text(x)] for k in EXPERIMENT}
    extra = {k: text(block.find_value(k)) for k in ('_struct.title','_struct.pdbx_descriptor','_struct_keywords.text')}
    citations = [{k: text(x.get(k)) for k in ('id','title','pdbx_database_id_DOI','pdbx_database_id_PubMed')}
                 for x in rows(block, '_citation.') if text(x.get('id')).lower() == 'primary']
    conn = {text(x.get('id')): {k: text(v) for k,v in x.items()} for x in rows(block,'_struct_conn.')}
    scope = {text(a['id']):text(a['entity_id']) for a in rows(block,'_struct_asym.')}
    assembly_labels=defaultdict(set)
    for a in rows(block,'_pdbx_struct_assembly_gen.'):
        assembly_labels[text(a['assembly_id'])].update(text(a['asym_id_list']).replace(' ','').split(','))
    in_deposited_assembly=set().union(*assembly_labels.values()) if assembly_labels else set(scope)
    observed=defaultdict(set)
    atom=block.find('_atom_site.', ['label_asym_id','label_entity_id','occupancy','type_symbol'])
    for col in atom:
        try:
            occ=float(col[2])
        except ValueError:
            continue
        if occ>0 and col[3] not in ('H','D','T'):
            observed[text(col[1])].add(text(col[0]))
    return {'assembly_observed_entities':{a:sorted(k for k,v in observed.items() if v & labels) for a,labels in assembly_labels.items()},
            'scope':sorted({scope[l] for l in in_deposited_assembly if l in scope}),
            'observed_entities':{k:sorted(v) for k,v in observed.items()},
            'entity':entity, 'experiment':exp,'extra':extra,'citations':citations,'names':dict(explicit_role_names(block)),
            'polylen':dict(Counter(text(x['entity_id']) for x in rows(block,'_entity_poly_seq.'))),
            'methods':[text(x) for x in block.find_values('_exptl.method')], 'connections':conn}


def member_kind(member):
    kind = member['kind']
    if kind != 'polymer':
        return kind
    poly = str(member.get('polymer_type') or '').lower()
    return 'peptide_or_protein' if poly.startswith('polypeptide') else ('nucleic_acid' if 'ribonucleotide' in poly else 'other_polymer')


def join_json(row, name):
    return json.loads(row[name]) if row.get(name) else None


def selection_cause(obs, decision, meta):
    if decision['eligible_for_structure_review']:
        return 'admitted'
    members = obs['members']
    if not decision['models_with_contacts']:
        return 'no_compatible_direct_protein_contact'
    if all(x['kind']=='polymer_modification' for x in members):
        return 'modification_in_polymer_not_independent_partner'
    if any(member_kind(x)=='peptide_or_protein' and meta['polylen'].get(x['entity'], 0)>35 for x in members):
        return 'long_polypeptide_without_admission_evidence'
    if any(member_kind(x)=='peptide_or_protein' and meta['polylen'].get(x['entity'], 0)<2 for x in members):
        return 'single_residue_polymer_without_admission_evidence'
    if any(m['kind']=='component' for m in members):
        if all(r['positive_heavy_atoms'] == 1 for m in members for r in m['residues']):
            return 'single_heavy_atom_without_admission_evidence'
        if max((r['positive_heavy_atoms'] for m in members for r in m['residues']), default=0)<5:
            return 'small_component_below_five_heavy_atom_gate'
    return 'other_admission_gate_not_met'


def obs_short(obs):
    return [{k: m.get(k) for k in ('instance','label_asym','operator','entity','kind','description','polymer_type')} |
            {'residues':[(r['auth_asym'],r['auth_seq'],r['ccd'],r['positive_heavy_atoms']) for r in m['residues']]}
            for m in obs['members']]


def scan(inp):
    batch, root, entry, system, out, index_annotation = inp
    p=entry['pdb_id']; source=root/entry['raw']; revision=entry['revision']; rroot=Path(entry.get('result_root') or root)
    sha=sha256(source)
    if sha != entry['source_sha256']:
        raise ValueError(f'Raw SHA differs: {p}')
    block, reported=revision_from_cif(source,p)
    if reported!=revision:
        raise ValueError(f'Raw revision differs: {p} {reported} != {revision}')
    meta=source_metadata(block)
    ledger_path=rroot/'discovery'/f'{p}.json'; ledger=load(ledger_path)
    if ledger['revision']!=revision or ledger['pdb_id']!=p:
        raise ValueError(f'Ledger identity differs: {p}')
    keyed={r['candidate_key']:r for r in ledger['native_observations']}
    chosen=ledger['candidate_selection']; decisions={r['candidate_key']:r for r in chosen}
    eligible={r['candidate_key'] for r in chosen if r['eligible_for_structure_review']}
    keys={s['candidate_key'] for s in system}
    if keys != eligible:
        raise ValueError(f'System vs selected keys differ: {p} {len(keys^eligible)}')
    if len(keys)!=len(system) or len(system)!=int(entry['chosen_systems']):
        raise ValueError(f'Bad system counts: {p}')
    base={'batch':batch,'pdb_id':p,'revision':revision,'route':entry['route'],'stratum':entry['stratum'],
          'evidence_layer':entry.get('evidence_layer') or 'corrected_sealed_200','source_path':str(source),
          'source_sha256':sha,'discovery_path':str(ledger_path), 'build_fingerprint':ledger.get('signature'),
          'index_source_revision':entry.get('index_source_revision') or None}
    partition=Counter(); omitted=[]
    for d in chosen:
        o=keyed[d['candidate_key']]; cause=selection_cause(o,d,meta); partition[cause]+=1
        if cause!='admitted':
            omitted.append({'candidate_key':d['candidate_key'],'assembly':d['assembly'],'model':d['model'],
                            'observed_models':d['observed_models'],'models_with_contacts':d['models_with_contacts'],
                            'direct_contacts':o['contact_instances_4A'],'reason':cause,'baseline_reason':d['reason'],
                            'members':obs_short(o),'role_input':{k:o[k] for k in ('members','contact_instances_4A','role_qc','role_basis')},
                            'current_admission_evidence':review_admission_evidence(o,meta['extra']['_struct.title'],meta['names'],meta['polylen'])})
    summary=base | {'candidate_count':len(system),'observations':len(keyed),'loci':len(chosen),
                    'selection_reasons':dict(partition),'omitted':omitted,'methods':meta['methods'],
                    'entity_counts':dict(Counter(v['type'] for v in meta['entity'].values())),
                    'source_connection_kinds':dict(Counter(x['conn_type_id'] for x in meta['connections'].values())),
                    'baseline_status':entry['status'], 'baseline_reason':entry.get('review_reason'),
                    'readbacks':int(entry['structure_readback_pass']), 'intake':int(entry['locally_supported_candidates'])}
    result=[]
    flags=('ligand_is_artifact','ligand_is_cofactor_like','ligand_is_ion','ligand_is_covalent','ligand_molecule_type')
    source_flags={(r['system_id'],r['ligand_id']):{k:r.get(k) for k in flags} for r in index_annotation}
    for s in system:
        k=s['candidate_key']; o=keyed[k]; folder=Path(s.get('result_root') or rroot)
        eid={m['entity'] for m in o['members']}
        codes={r['ccd'] for m in o['members'] if m['kind']=='component' for r in m['residues']}
        result.append({**base,'candidate_key':k,'assembly':o['assembly'],'model':o['model'],
            'model_scope':join_json(s,'model_scope'),'kinds':sorted({member_kind(m) for m in o['members']}),
            'entity_ids':sorted(eid),'ccds':sorted(codes),'direct_contacts':o['contact_instances_4A'],
            'role_input':{k:o[k] for k in ('members','contact_instances_4A')},
            'member_instances':[m['instance'] for m in o['members']], 'member_grouping':s.get('member_grouping') or None,
            'partner_scope_qc':s.get('partner_scope_qc') or None, 'prior_role':s['role_qc'],'prior_role_basis':o['role_basis'],
            'prior_status':s['status'],'prior_intake':s['candidate_intake'],'coordinate_identity':s['coordinate_identity'],
            'experimental_fidelity':s['experimental_fidelity'],'topology':s['topology'],'export_reason':s.get('reason') or None,
            'export_detail':s.get('detail') or None,'mapping_error_evidence':join_json(s,'mapping_error_evidence'),
            'system_path':str(folder/s['system_path']) if s.get('system_path') else None,
            'mapping_path':str(folder/s['mapping']) if s.get('mapping') else None,
            'system_sha256':s.get('system_sha256') or None,'index_matches':o['index_aliases'],
            'same_version_index_flags':[source_flags[(r['system_id'],r['ligand_id'])] for r in o['index_aliases']
                if r['status']=='SOURCE_INSTANCE_MATCHED' and (r['system_id'],r['ligand_id']) in source_flags],
            'source_connection_kinds':sorted({x['source']['kind'] for x in o['connections']}),
            'unresolved_connections':o['unresolved_connection_ids']})
    summary['observed_entities']=meta['observed_entities']; summary['deposited_assembly_entities']=meta['scope']
    summary['assembly_observed_entities']=meta['assembly_observed_entities']
    summary['declared_entity_types']={k:v['type'] for k,v in meta['entity'].items()}
    compact = {'entry':summary,'systems':result,'meta':meta}
    target=out/f'scan_{SCAN_VERSION}'/f'{p}.json'; temp=target.with_suffix('.tmp')
    temp.write_text(dump(compact),encoding='utf8');temp.replace(target)
    return p


def tasks(root, out):
    overall=set(); pending=[]; baseline={}
    for batch, name in RUNS.items():
        run=root/name; table=run/('handoff' if batch=='2000' else '')
        ep=table/'entry_results.csv'; sp=table/'system_results.csv'
        entries=csvrows(ep); system=csvrows(sp); by=defaultdict(list)
        for s in system:
            by[s['pdb_id']].append(s)
        ids={e['pdb_id'] for e in entries}
        totals=(len(entries),len(system),sum(s['status']=='EXPORTED_READBACK' for s in system),
                sum(s['candidate_intake']=='PASS' for s in system))
        assert not ids&overall and len(ids)==len(entries) and set(by)<=ids and totals==EXPECTED[batch],(batch,totals)
        overall.update(ids)
        baseline[batch]={'entry_table':str(ep),'entry_sha256':sha256(ep),'system_table':str(sp),'system_sha256':sha256(sp),'counts':totals}
        idx=load(run/'index_snapshot.json')['pdb']
        pending.extend((batch,run,e,by[e['pdb_id']],out,idx[e['pdb_id']]['annotation']) for e in entries)
    return pending,baseline


def diagnosis(x, meta):
    tags=[]; notes=x['prior_role_basis']; fam=set(x['kinds']); title=meta['extra']['_struct.title']; t=normalize(title)
    companions=[]
    for m in x['role_input']['members']:
        ent=meta['entity'].get(m['entity'],{}); choices=[]
        if m['kind']=='component':
            for r in m['residues']:
                cc=r['ccd']; choices+=list(meta['names'].get(cc, []))+[('literal CCD',cc)]
        else:
            choices=[('_entity.pdbx_description',m.get('description'))]
        hit=[(f,n) for f,n in choices if len(normalize(n))>=3 and re.search(token(normalize(n)),t)]
        variants=[(f,n) for f,n in choices if len(normalize(n))>=7 and normalize(n).endswith('ic acid')
                  and re.search(token(normalize(n)[:-7]+'ate'),t)]
        section=[]
        for f,n in choices:
            a=normalize(n)
            if len(a)<4:
                continue
            for exp, strings in meta['experiment'].items():
                if any(re.search(token(a),normalize(z)) for z in strings):
                    section.append(exp)
            if any(re.search(token(a),normalize(c['title'])) for c in meta['citations']):
                section.append('_citation.title')
        if ent.get('details'):
            section.append('_entity.details')
        companions.append({'entity':m['entity'],'title_hits':sorted(set(tuple(z) for z in hit)),
                           'possible_spelling_variants':variants,'other_fields':sorted(set(section))})
    if x['prior_role']=='PASS':
        tags.append('prior_title_role_pass')
    else:
        if len(x['entity_ids'])>1 or len(x['role_input']['members'])>1:
            tags.append('multi_member_same_role_not_established')
        if 'nucleic_acid' in fam:
            tags.append('nucleic_role_or_complete_partner_not_established')
        if 'polymer_modification' in fam:
            tags.append('polymer_modification_needs_context')
        if any(member_kind(m)=='peptide_or_protein' and meta['polylen'].get(m['entity'],0)>35 for m in x['role_input']['members']):
            tags.append('peptide_vs_protein_length_gate')
        if 'competing environment' in notes['reason']:
            tags.append('role_conflict_explicitly_detected_by_legacy')
        if any(m['title_hits'] for m in companions):
            tags.append('title_mentions_candidate_without_sufficient_role_syntax')
        if any(m['possible_spelling_variants'] for m in companions):
            tags.append('possible_title_alias_not_resolved')
        if any(m['other_fields'] for m in companions):
            tags.append('non_title_field_has_clue_but_role_not_established')
        for key in ('ligand_is_artifact','ligand_is_cofactor_like','ligand_is_ion'):
            if any(flag[key] is True for flag in x['same_version_index_flags']):
                tags.append('index_'+key.removeprefix('ligand_is_')+'_hint_role_not_established')
        if not tags:
            tags.append('no_specific_contextual_role_evidence_found')
    return tags,companions


def topology_detail(x, meta):
    q=x['mapping_error_evidence'] or {}; ids=q.get('ids',[]) + ([q['connection']] if q.get('connection') else [])
    conns={i:meta['connections'].get(i) for i in ids if (q.get('category') or '_struct_conn.')=='_struct_conn.'}
    sym={tuple(r.get(f'ptnr{i}_symmetry','') for i in (1,2)) for r in conns.values() if r}
    sym.update([tuple(q['symmetry'])] if q.get('symmetry') else [])
    categories=sorted({r['conn_type_id'] for r in conns.values() if r})
    reason=x['export_reason']; subtype=None
    if reason=='connection_endpoint_unresolved':
        subtype='endpoint_not_observed_in_export_model' if q.get('matches')==0 else 'endpoint_ambiguous_in_export_model'
    elif reason=='expanded_connection_instance_unresolved':
        subtype='mapped_endpoint_has_zero_assembly_copies' if q.get('matches') is not None and any(not a for a in q['matches']) else 'copy_assignment_not_unique_or_not_supported'
    elif reason=='unresolved_partner_connection':
        subtype='nonlocal_partner_endpoints_not_uniquely_resolved' if any(a!=b for a,b in sym) else 'partner_endpoint_unresolved_local_or_unreported_symmetry'
    elif reason=='expanded_topology_pending':
        subtype='legacy_global_copy_gate_explicit_base_pair' if q.get('category')=='_ndb_struct_na_base_pair.' else 'legacy_global_copy_gate_explicit_connection'
        if q.get('category')=='_struct_conn.':
            categories=sorted({r['conn_type_id'] for r in meta['connections'].values()})
    return subtype,conns,sorted(sym),categories


def zero_cause(e):
    if e['observations']:
        return 'observations_exist_none_admitted'
    ts=e['declared_entity_types']; scope=set(e['deposited_assembly_entities']); observed=set(e['observed_entities']);
    non_water={k for k,t in ts.items() if t in ('non-polymer','branched')}
    if non_water:
        if not non_water & observed:
            return 'declared_nonpolymer_or_branch_not_observed_as_positive_heavy_atoms'
        if not non_water & observed & scope:
            return 'observed_nonpolymer_or_branch_outside_deposited_assembly'
        return 'independent_entity_observed_but_not_emitted_needs_mapping_review'
    if sum(ts.get(i)=='polymer' for i in scope&observed)>1:
        assembly=e['assembly_observed_entities']
        if assembly and not any(sum(ts.get(i)=='polymer' for i in entities)>1 for entities in assembly.values()):
            return 'observed_polymer_entities_deposited_in_separate_assemblies'
        return 'multiple_observed_polymer_entities_not_emitted_needs_context_review'
    return 'only_one_observed_polymer_entity_or_no_independent_partner_in_scope'


def load_historical(root):
    overlay={}; controls={}; checks=[]
    for layer,dirname,ids in HISTORICAL:
        run=root/dirname; er={x['pdb_id']:x for x in csvrows(run/'entry_results.csv')}; sr=csvrows(run/'system_results.csv')
        assert set(er)==ids and all(s['pdb_id'] in ids for s in sr)
        for p,e in er.items():
            control={'scan_only_prior_review': e.get('review_reason'), 'targeted_candidates':int(e['chosen_systems']),
                     'targeted_readbacks':int(e['structure_readback_pass']), 'root':str(run), 'historical_layer':layer}
            controls[p]=control
        for s in sr:
            p=s['pdb_id']; k=s['candidate_key']; e=er[p]
            assert s['revision']==e['revision'] and (p,k) not in overlay and sha256(run/s['system_path'])==s['system_sha256']
            overlay[p,k]={'pdb_id':p,'candidate_key':k,'revision':s['revision'],'historical_layer':layer,'result_root':str(run),
                          'source_sha256':e['source_sha256'],'readback_after':s['status'],'role_after':s['role_qc'],
                          'candidate_intake_after':s['candidate_intake'],'system_sha256':s['system_sha256'],
                          'system_path':str(run/s['system_path']),'mapping_path':str(run/s['mapping'])}
        checks.append({'layer':layer,'entry_csv_sha256':sha256(run/'entry_results.csv'),
                       'system_csv_sha256':sha256(run/'system_results.csv'),'systems':len(sr)})
    return overlay,controls,checks


def m3_controls(out, m3):
    prior=load(out.parent/RUNS['2000']/'m3_first_review.json')['review_candidates']
    current={(x['pdb_id'],x['candidate_key']):x for x in m3}
    checked=[]
    for row in prior:
        key=row['pdb_id'],row['candidate']; x=current[key]
        assert row['source_sha256']==x['source_sha256'] and row['system_sha256']==x['system_sha256'] and sha256(Path(x['system_path']))==x['system_sha256']
        checked.append({'pdb_id':key[0],'candidate_key':key[1],'sha256':x['system_sha256'],'status':'UNCHANGED_PROGRAM_PASS_NO_PRODUCTION_APPROVAL'})
    assert len(checked)==30
    return checked


def aggreg(out):
    start=time.monotonic(); before=load(out/'baseline.json'); historical,controls,hverify=load_historical(out.parent)
    by=defaultdict(Counter); slice_counts=defaultdict(Counter); transition=defaultdict(Counter); examples=defaultdict(list)
    omitted_by=defaultdict(Counter); candidates_by_pdb={}; eligible_changes=[]; changed_groups=defaultdict(Counter)
    n=Counter(); m3=[]; drops=[]; historical_out=[]; hits=set()
    targets=[out/x for x in ('entry_attribution.jsonl','system_attribution.jsonl','role_changes.jsonl','zero_candidate_entries.jsonl')]
    with targets[0].with_suffix('.tmp').open('w',encoding='utf8') as ef, targets[1].with_suffix('.tmp').open('w',encoding='utf8') as sf, targets[2].with_suffix('.tmp').open('w',encoding='utf8') as cf, targets[3].with_suffix('.tmp').open('w',encoding='utf8') as zf:
        for f in sorted((out/f'scan_{SCAN_VERSION}').glob('*.json')):
            frame=load(f); e=frame['entry']; meta=frame['meta']; p=e['pdb_id']; batch=e['batch']; b=by[batch]
            b['entries']+=1; n['raw_bytes_reused']+=Path(e['source_path']).stat().st_size
            b['candidates']+=len(frame['systems']); stage_counts=Counter(); reasons=Counter(); e_roles=Counter(); by_ent=defaultdict(Counter)
            title=meta['extra']['_struct.title']; context=build_role_context(meta,meta['names'])
            for o in e['omitted']:
                v=o.pop('role_input'); new=assess(v,v['role_qc'],title,meta['names'],meta['polylen'],meta,context)
                if new:
                    o['replay_role_assertion']={'status':new[0],'evidence':new[1]}
                if o['current_admission_evidence'] or (new and new[0]=='PASS'):
                    new_origin='role_evidence' if new and new[0]=='PASS' else 'explicit_source_review_discovery'
                    eligible_changes.append({'batch':batch,'pdb_id':p,'revision':e['revision'],'candidate_key':o['candidate_key'],
                         'assembly':o['assembly'],'model':o['model'],'reason':new_origin,'historically_verified':(p,o['candidate_key']) in historical,
                         'evidence':(new or (None,o['current_admission_evidence']))[1]})
                if (p,o['candidate_key']) in historical:
                    x=historical[p,o['candidate_key']]; assert x['revision']==e['revision'] and x['source_sha256']==e['source_sha256']
                    historical_out.append(x | {'baseline_stage':'not_admitted','delta_candidates':1,'delta_readbacks':1})
                    hits.add((p,o['candidate_key']));o['historical_retest']=x
            for x in frame['systems']:
                b['readbacks_before']+=x['prior_status']=='EXPORTED_READBACK'; b['intake_before']+=x['prior_intake']=='PASS'
                role, role_basis=x['prior_role'],x['prior_role_basis']; tags, clue=diagnosis(x,meta)
                assertion=assess(x['role_input'],role,title,meta['names'],meta['polylen'],meta,context)
                subtype,conns,syms,conn_kinds=topology_detail(x,meta)
                readback=x['prior_status']=='EXPORTED_READBACK'; scope=x['partner_scope_qc']=='PASS'
                if assertion:
                    role,role_basis=assertion
                    if readback:
                        group=(batch,p,role_basis['rule'],role)
                        changed_groups[group][str((tuple(x['entity_ids']),tuple(x['ccds'])))]+=1
                        assert sha256(Path(x['system_path']))==x['system_sha256'],x['system_path']
                    cf.write(dump({k:v for k,v in x.items() if k!='role_input'} | {'role_after':role,'evidence_after':role_basis})+'\n')
                intake=readback and scope and role=='PASS' and all(x[z]=='PASS' for z in ('coordinate_identity','experimental_fidelity','topology'))
                if x['prior_intake']=='PASS' and not intake:
                    drops.append((batch,p,x['candidate_key']))
                scope_cause=(x['member_grouping'] if str(x['member_grouping']).startswith('REVIEW') else
                    'nucleic_partner_scope_requires_independent_review_even_if_pairing_mapped' if set(x['kinds'])&{'nucleic_acid','other_polymer'} else 'partner_scope_review')
                stage,primary=('intake','program_pass') if intake else (
                    ('structure_export',subtype or x['export_reason'] or 'unclassified_export_hold') if not readback else
                    ('role','confirmed_preparation_role_not_target') if role=='NON_TARGET_EVIDENCED' else
                    ('partner_scope',scope_cause) if not scope else ('role',tags[0]))
                main_before='intake' if x['prior_intake']=='PASS' else 'structure_export' if not readback else 'role' if scope else 'role_and_scope'
                stage_counts[stage]+=1; reasons[primary]+=1; e_roles[role]+=1; transition[batch][x['prior_role']+' -> '+role]+=1
                b['role_pending_before_readback']+=readback and x['prior_role']!='PASS'
                b['role_pending_after_readback']+=readback and role=='REVIEW'
                b['non_target_after_readback']+=readback and role=='NON_TARGET_EVIDENCED'
                b['role_pass_after_readback']+=readback and role=='PASS'; b['intake_after']+=intake
                b['scope_and_role_review_after_readback']+=readback and not scope and role=='REVIEW'
                b['structure_holds']+=not readback
                for kind in x['kinds']:
                    slice_counts[(batch,'kind',kind)][primary]+=1
                for kind in sorted(set(x['source_connection_kinds'])|set(conn_kinds)) or ['unknown_or_no_record']:
                    slice_counts[(batch,'conn',kind)][primary]+=1
                slice_counts[(batch,'route',x['route'])][primary]+=1
                for tag in tags:
                    slice_counts[(batch,'role_clue_all_readback' if readback else 'role_clue_no_readback',tag)][primary]+=1
                for entity in x['entity_ids']:
                    by_ent[entity][primary]+=1
                if len(examples[(batch,primary)])<3:
                    examples[(batch,primary)].append([p,x['candidate_key']])
                key=p,x['candidate_key']; linked=historical.get(key)
                if linked:
                    assert linked['revision']==x['revision'] and linked['source_sha256']==x['source_sha256']; hits.add(key)
                    historical_out.append(linked | {'baseline_stage':main_before, 'delta_candidates':0,
                          'delta_readbacks':int(not readback),'delta_intake':int(linked['candidate_intake_after']=='PASS' and x['prior_intake']!='PASS')})
                x.pop('role_input')
                x |= {'stage_before':main_before,'stage_after':stage,'primary_reason':primary,'role_clues':clue,'role_subcauses':tags,
                      'role_after':role,'role_evidence_after':role_basis,'intake_after':'PASS' if intake else 'EXCLUDED' if role=='NON_TARGET_EVIDENCED' else 'REVIEW',
                      'source_connection_rows':conns, 'connection_symmetry':syms,'error_connection_kinds':conn_kinds,
                      'historical_retest':linked}
                sf.write(dump(x)+'\n')
                if intake:
                    m3.append({k:x[k] for k in ('batch','pdb_id','revision','candidate_key','assembly','model','model_scope',
                              'source_path','source_sha256','system_path','system_sha256','mapping_path','role_evidence_after',
                              'coordinate_identity','experimental_fidelity','topology','partner_scope_qc','evidence_layer')})
            if not e['candidate_count']:
                e['entry_stage']='discovery'; e['primary_reason']=zero_cause(e)
                omitted_by[(batch,e['stratum'])].update(e['selection_reasons']); b['zero_candidates']+=1
                slice_counts[(batch,'zero_entry',e['stratum'])][e['primary_reason']]+=1
                zf.write(dump(e | {'source_title':title})+'\n')
            elif not e['readbacks']:
                e['entry_stage']='structure_export'; e['primary_reason']=reasons.most_common(1)[0][0]
            elif stage_counts['intake']:
                e['entry_stage']='intake';e['primary_reason']='at_least_one_program_pass'
            else:
                e['entry_stage']='readback_not_intake'; e['primary_reason']=reasons.most_common(1)[0][0]
            b['entry_'+e['entry_stage']]+=1
            for k,v in reasons.items(): slice_counts[(batch,'primary','all')][k]+=v
            e['system_stages_after']=dict(stage_counts); e['system_reasons_after']=dict(reasons);e['role_counts_after']=dict(e_roles)
            e['historical_control']=controls.get(p); e['source_title']=title; e['by_deposited_entity']=dict(by_ent)
            ef.write(dump(e)+'\n')
            candidates_by_pdb[p] = e['candidate_count']
    assert not drops and hits==set(historical) and len(candidates_by_pdb)==2200,(drops,set(historical)-hits)
    for batch,(en,sy,rb,ip) in EXPECTED.items():
        b=by[batch]; assert (b['entries'],b['candidates'],b['readbacks_before'],b['intake_before'])==(en,sy,rb,ip)
    for p in targets: p.with_suffix('.tmp').replace(p)
    atomic_json(out/'historical_deltas.json',{'evidence':hverify,'deltas':historical_out,
                'already_in_2000_baseline':'9Q6G point_symmetry_retest via handoff; do not add again'})
    atomic_json(out/'structure_retest_selection.json',eligible_changes)
    atomic_json(out/'m3_program_shortlist.json',{'scope':'Program QC handoff only; no production approval',
         'original_30_controls':m3_controls(out,m3),'candidates':m3})
    atomic_json(out/'attribution_summary.json',{'written_utc':now(),'rules':VERSION,'baseline':before,'by_batch':dict(by),
        'role_transition':dict(transition),'slice_counts':[{'batch':b,'dimension':d,'value':v,'reasons':c} for (b,d,v),c in sorted(slice_counts.items())],
        'zero_candidate_omission_locus_counts':[{'batch':b,'stratum':s,'counts':c} for (b,s),c in sorted(omitted_by.items())],
        'role_changed_groups':[{'batch':b,'pdb_id':p,'rule':r,'status':s,'by_entity_ccd':c} for (b,p,r,s),c in sorted(changed_groups.items())],
        'example_keys':[{'batch':b,'reason':r,'keys':v} for (b,r),v in sorted(examples.items())],
        'new_admission_loci':len(eligible_changes),'already_targeted_admission_loci':sum(l['historically_verified'] for l in eligible_changes),
        'reused_raw_bytes':n['raw_bytes_reused'],'new_coordinate_download_bytes':0,'replayed_role_pass_regressions':len(drops),
        'elapsed_s':round(time.monotonic()-start,3),'table_sha256':{p.name:sha256(p) for p in targets},
        'code_sha256':{p.name:sha256(p) for p in Path(__file__).parent.glob('holo*.py')}})
    print(dump(dict(by)),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out',type=Path); parser.add_argument('--jobs',type=int,default=24)
    parser.add_argument('--mode',choices=('all','scan','replay'),default='all')
    args=parser.parse_args(); out=args.out.resolve(); (out/f'scan_{SCAN_VERSION}').mkdir(parents=True,exist_ok=True)
    if args.mode in ('all','scan'):
        start=time.monotonic(); all_tasks,baseline=tasks(out.parent,out); baseline_path=out/'baseline.json'
        baseline_index={k:sha256(out.parent/RUNS[k]/'index_snapshot.json') for k in RUNS}
        for k in baseline: baseline[k]['index_snapshot_sha256']=baseline_index[k]
        if baseline_path.is_file() and load(baseline_path)!=baseline:
            current=load(baseline_path)
            if not (all(baseline[k]['entry_sha256']==current[k]['entry_sha256'] and baseline[k]['system_sha256']==current[k]['system_sha256'] for k in baseline) and
                    ('index_snapshot_sha256' not in current['200'] or current['200']['index_snapshot_sha256']==baseline_index['200']) and
                    ('index_snapshot_sha256' not in current['2000'] or current['2000']['index_snapshot_sha256']==baseline_index['2000'])):
                raise ValueError('Baseline identity differs')
        atomic_json(baseline_path,baseline)
        todo=[t for t in all_tasks if not (out/f'scan_{SCAN_VERSION}'/f"{t[2]['pdb_id']}.json").is_file()]
        print('scan pending',len(todo),'reuse',len(all_tasks)-len(todo),flush=True)
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            for n,p in enumerate(pool.map(scan,todo,chunksize=1),1):
                if n%100==0 or n==len(todo): print(n,p,round(time.monotonic()-start,2),flush=True)
        atomic_json(out/'scan_receipt.json',{'new':len(todo),'reused':len(all_tasks)-len(todo),'elapsed_s':round(time.monotonic()-start,3)})
    if args.mode in ('all','replay'):
        aggreg(out)


if __name__=='__main__':
    main()
