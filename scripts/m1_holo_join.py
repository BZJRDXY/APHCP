"""Join an APHCP M1 monthly metadata snapshot to the fixed PLINDER release."""
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import resource
import sys
import time
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

from m1_holo_config import read_config, sha256, window

run = Path(sys.argv[1])
config = read_config(run)
start, end, cutoff = window(config)
start_day, end_day = start.isoformat(), end.isoformat()
cutoff_day = cutoff.isoformat() if cutoff else None
release = config['plinder_release']
default_root = Path('/srv/aphcp/upstream/PLINDER')
root = Path(config.get('plinder_index', default_root / release / 'index'))
out=run/'inventory'; out.mkdir(exist_ok=True)
started=datetime.now(timezone.utc); t0=time.monotonic()
load=lambda p:json.loads(p.read_text(encoding='utf-8'))
def dump(name,data):
    path=out/name; tmp=path.with_suffix('.tmp'); tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf8'); tmp.replace(path)
def csvout(name,rows,cols=None):
    path=out/name; tmp=path.with_suffix('.tmp')
    with tmp.open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols or (list(rows[0]) if rows else [])); w.writeheader()
        w.writerows({k:json.dumps(v,ensure_ascii=False,separators=(',',':')) if isinstance(v,(list,dict)) else v for k,v in row.items()} for row in rows)
    tmp.replace(path)

def index_rows(name, cols, filt):
    file=root/(name+'.parquet'); receipt=load(file.with_suffix('.parquet.source.json')); meta=pq.ParquetFile(file)
    assert file.stat().st_size==receipt['bytes']==receipt['expected_bytes'] and f'/PLINDER-{release}/' in receipt['url'],name
    rows=ds.dataset(str(file),format='parquet').to_table(columns=cols,filter=filt).to_pylist()
    index[name]={'rows':meta.metadata.num_rows,'columns':len(meta.schema_arrow),'rows_selected':len(rows),'bytes':file.stat().st_size,
        'source_receipt':receipt,'selected_schema':{c:str(meta.schema_arrow.field(c).type) for c in cols}}
    return rows

def by(rows,k='entry_pdb_id'):
    out=defaultdict(list)
    for r in rows:out[r[k]].append(r)
    return out

def rev(a):
    a=a or {}; m=a.get('major_revision'); n=a.get('minor_revision')
    return f'{m}.{n}' if m is not None and n is not None else None

def when(date):
    if not date or not cutoff_day:
        return 'UNKNOWN'
    day = date[:10]
    return 'BEFORE' if day < cutoff_day else ('ON' if day == cutoff_day else 'AFTER')

def labels_from_entity(row, kind):
    return row.get(f'rcsb_{kind}_entity_container_identifiers') or {}

search=load(run/'search'/'summary.json'); entries=load(run/'data'/'entries.json'); ids=search['entry_ids']; lower=[i.lower() for i in ids]
assert search['concordant'] and [e['rcsb_id'] for e in entries]==ids and len(ids)==len(set(ids))
seen=set(lower); fields=ds.field('entry_pdb_id'); id_filter=fields.isin(lower); date=ds.field('entry_release_date'); index={}
meta=index_rows('entry_metadata','entry_pdb_id entry_release_date entry_determination_method entry_failed_assembly_ids'.split(),id_filter | ((date >= start_day) & (date < end_day)) | date.is_null() | (date==''))
M=by(meta)
period={r['entry_pdb_id'] for r in meta if start_day <= (r['entry_release_date'] or '')[:10] < end_day}
undated={r['entry_pdb_id'] for r in meta if not r['entry_release_date']}
union=seen|period|undated
F=fields.isin(sorted(union)); src=by(index_rows('entry_sources','entry_pdb_id source_mmcif_major_revision source_mmcif_minor_revision'.split(),F))
anno_cols='entry_pdb_id system_biounit_id system_id ligand_id system_protein_chains_asym_id system_ligand_chains_asym_id ligand_asym_id ligand_instance ligand_instance_chain ligand_instance_chains ligand_residue_numbers ligand_ccd_code ligand_unique_ccd_code ligand_molecule_type ligand_bird_id ligand_is_proper ligand_is_covalent ligand_covalent_linkages ligand_is_ion ligand_is_cofactor ligand_is_cofactor_like ligand_is_artifact ligand_is_monosaccharide ligand_is_oligosaccharide ligand_is_mononucleotide ligand_is_oligonucleotide ligand_is_monopeptide ligand_is_oligopeptide ligand__is_multi_residue ligand_num_interactions system_pocket_Pfam system_pocket_UniProt'.split()
ann=by(index_rows('annotation_table',anno_cols,F)); annotated_period={pid for pid, rows in ann.items() if rows} & period; extra=sorted(period-seen)
extra_rows=[{**M[p][0],'annotation_rows':len(ann[p]),'candidate_systems':len({r['system_id'] for r in ann[p]}),'source_revisions':src[p], 'rcsb_protein_eligibility':'NOT_CHECKED'} for p in extra]
B=by(index_rows('entry_biounit_chains','entry_pdb_id biounit_id chain_instance chain_asym_id chain_role'.split(),id_filter)); C=by(index_rows('entry_chains','entry_pdb_id chain_asym_id chain_auth_id chain_entity_id chain_type chain_uniprot_ids'.split(),id_filter))
pockets=index_rows('ligand_pocket_residues','entry_pdb_id system_id ligand_id chain_instance is_pli'.split(),id_filter & (ds.field('is_pli')==True))
contacts=defaultdict(set)
for r in pockets:
    if r['chain_instance']: contacts[r['system_id'],r['ligand_id']].add(r['chain_instance'])
del pockets

out_entries=[]; out_pairs=[]; out_chains=[]; assemblies=[]; entity_rows=[]; problems=[]
Mcount=Counter(); methods=Counter(); method_resolution=Counter(); family_types=Counter(); family_entities=Counter(); modes=Counter(); classes=Counter(); role_status=Counter(); trail=Counter(); joint=Counter(); rev_count=Counter(); dates=Counter(); revisions_seen=set(); source_rev_seen=set(); pdb_ligands=set(); pdb_proteins=set(); pl_assemblies=set(); systems=set(); ww_assemblies=set(); family_values=defaultdict(Counter)
atom_rows=[]
for e in entries:
    pdb=e['rcsb_id']; pid=pdb.lower(); a=e.get('rcsb_accession_info') or {}; ci=e.get('rcsb_entry_container_identifiers') or {}; info=e.get('rcsb_entry_info') or {}; title=(e.get('struct') or {}).get('title'); history=e.get('pdbx_audit_revision_history') or []; pub=e.get('rcsb_primary_citation') or {}; dates[(a.get('initial_release_date') or 'UNKNOWN')[:10]]+=1
    expected={}; observed={}; asym_entities=defaultdict(list); tags=set(); path_evidence=[]; comp=[]; family_labels=set(); c_reasons=set(); method_list=sorted({m['method'] for m in e.get('exptl') or [] if m.get('method')}); methods.update(method_list or ['UNKNOWN'])
    if info.get('resolution_combined'):method_resolution.update(method_list)
    if ci.get('non_polymer_entity_ids'):tags.add('A_NONPOLYMER')
    if ci.get('branched_entity_ids'):tags.add('B_BRANCHED')
    entry_prd=[r for r in e.get('pdbx_molecule_features') or [] if r.get('prd_id')]; model_history=[r for r in history if r.get('data_content_type')=='Structure model' and rev(r)]
    if entry_prd:c_reasons.add('ENTRY_PRD'); path_evidence+=entry_prd
    if len(ci.get('polymer_entity_ids') or [])>1:trail['multiple_polymer_entities']+=1
    for kind,listed in [('polymer','polymer_entity_ids'),('nonpolymer','non_polymer_entity_ids'),('branched','branched_entity_ids')]:
        expected[kind]=set(ci.get(listed) or []); observed[kind]=set()
        for entity in (e.get(kind+'_entities') or []):
            if not entity:continue
            id_block=labels_from_entity(entity,kind); eid=id_block.get('entity_id'); ids_asym=id_block.get('asym_ids') or []; desc=None; code=None; pol=(entity.get('entity_poly') or {}); f={}
            if eid:observed[kind].add(eid)
            if id_block.get('prd_id'):
                c_reasons.add(f'{kind.upper()}_PRD'); path_evidence.append({'entity_id':eid,'entity_kind':kind,'prd_id':id_block['prd_id']})
            if kind=='polymer':
                pe=entity.get('rcsb_polymer_entity') or {}; desc=pe.get('pdbx_description'); polymer_type=pol.get('rcsb_entity_polymer_type'); length=pol.get('rcsb_sample_sequence_length')
                if polymer_type and polymer_type!='Protein':
                    c_reasons.add('NONPROTEIN_POLYMER'); path_evidence.append({'entity_id':eid,'polymer_type':polymer_type,'description':desc})
                elif polymer_type=='Protein' and re.search(r'\b(?:cyclic peptide|peptide (?:inhibitor|substrate)|insulin|amyloid[- ]beta peptide|peptide ligand)\b',desc or '',re.I):
                    c_reasons.add('NAMED_PEPTIDE'); path_evidence.append({'entity_id':eid,'description':desc,'length':length})
                elif polymer_type=='Protein' and length and length<=30 and len(expected[kind])>1:
                    c_reasons.add('SHORT_POLYMER_REVIEW'); path_evidence.append({'entity_id':eid,'description':desc,'length':length})
                if polymer_type=='Protein':
                    f['family_annotations']=[x for x in entity.get('rcsb_polymer_entity_annotation') or [] if x.get('annotation_id') and str(x.get('type','')).upper() in {'PFAM','CATH','SCOP','SCOP2','SCOP2B','ECOD','INTERPRO','PANTHER','SUPERFAMILY'}]
                    f['uniprot_ids']=id_block.get('uniprot_ids') or []; f['ec']=pe.get('pdbx_ec') or []
                    family_entities['protein_entities']+=1
                    if f['family_annotations']:family_entities['has_family_classification']+=1
                    for tp in set(f['family_annotations'][i]['type'].upper() for i in range(len(f['family_annotations']))):family_entities['class_'+tp]+=1
                    if f['uniprot_ids']:family_entities['has_uniprot']+=1
                    if f['ec']:family_entities['has_ec']+=1
                    for x in f['family_annotations']:
                        tp=x['type'].upper(); family_labels.add(tp); family_values[tp][x['annotation_id']]+=1
                    if f['uniprot_ids']:family_labels.add('UNIPROT')
                    if f['ec']:family_labels.add('EC')
            elif kind=='nonpolymer':
                v=entity.get('pdbx_entity_nonpoly') or {}; code=v.get('comp_id'); desc=v.get('name')
                if code:comp.append(code)
            else:
                v=entity.get('pdbx_entity_branch') or {}; desc=v.get('type'); f['component_count']=v.get('rcsb_branched_component_count'); f['monomers']=id_block.get('chem_comp_monomers')
            er={'pdb_id':pdb,'entity_kind':kind,'entity_id':eid,'label_asym_ids':ids_asym,'auth_asym_ids':id_block.get('auth_asym_ids'),'ccd':code,'prd':id_block.get('prd_id'),'description':desc,'polymer_type':pol.get('rcsb_entity_polymer_type'),'polymer_length':pol.get('rcsb_sample_sequence_length'),'family_evidence':f}
            entity_rows.append(er)
            for label in ids_asym:asym_entities[label].append(er)
    family_types.update(family_labels or ['NO_FAMILY_FUNCTION_OR_ID_EVIDENCE'])
    if not any(x not in ('EC','UNIPROT') for x in family_labels):trail['entry_family_class_unknown']+=1
    if c_reasons:tags.add('C_PARTNER_CLUE')
    if c_reasons=={'SHORT_POLYMER_REVIEW'}:trail['C_short_only']+=1
    trail.update(c_reasons)
    if not tags:tags.add('INSUFFICIENT')
    joint['+'.join(sorted(tags))]+=1
    source=src[pid]; rrev=rev(a); prev=None
    if len(source)==1:
        rr=source[0]; n=rr['source_mmcif_major_revision']; m=rr['source_mmcif_minor_revision']; prev=f'{n}.{m}' if n is not None and m is not None else None
    if rrev:revisions_seen.add((pdb,rrev))
    if prev:source_rev_seen.add((pdb,prev))
    pl=ann[pid]; status='NO_ANNOTATION' if not pl else ('REVISION_UNKNOWN' if not rrev or not prev else ('MATCH' if rrev==prev else 'DIFFERENT'))
    rev_count[status]+=1
    lookup=by(B[pid], 'biounit_id'); chainref=by(C[pid],'chain_asym_id'); assembly_objs={}; matrix_members={}
    for obj in e.get('assemblies') or []:
        if not obj:continue
        bid=(obj.get('rcsb_assembly_container_identifiers') or {}).get('assembly_id')
        if not bid:continue
        assert bid not in assembly_objs,(pdb,bid)
        assembly_objs[bid]=obj; gen=obj.get('pdbx_struct_assembly_gen') or []; matrix_members[bid]={label for g in gen for label in (g.get('asym_id_list') or [])}
        ww_assemblies.add((pdb,bid))
        assemblies.append({'pdb_id':pdb,'assembly_id':bid,'models':ci.get('model_ids'),'generators':gen,'size':obj.get('rcsb_assembly_info'),'oligomeric':obj.get('pdbx_struct_assembly')})
    related=set(); count_pl_systems=set(); row_status=Counter()
    for row in pl:
        sid, lid=row['system_id'],row['ligand_id']; assert sid and lid and (sid,lid) not in related,(sid,lid); related.add((sid,lid)); systems.add(sid); count_pl_systems.add(sid)
        aid=row['system_biounit_id']; pl_assemblies.add((pdb,aid)); by_inst=by(lookup[aid],'chain_instance')
        system_p=row['system_protein_chains_asym_id'] or []; system_l=row['system_ligand_chains_asym_id'] or []
        instances=sorted(set(row['ligand_instance_chains'] or []) | ({row['ligand_instance_chain']} if row['ligand_instance_chain'] else set()))
        instances_by_role={'PROTEIN':system_p,'LIGAND':instances}
        reasons=set(); matched_ligand_members=[]
        if aid not in (ci.get('assembly_ids') or []): reasons.add('ASSEMBLY_ID_ABSENT_RCSB')
        if not lookup[aid]:reasons.add('ASSEMBLY_CHAINS_ABSENT_INDEX')
        if not instances:reasons.add('NO_LIGAND_INSTANCE_REFERENCE')
        if not system_p:reasons.add('NO_PROTEIN_INSTANCE_REFERENCE')
        for role, refs in instances_by_role.items():
            for label in refs:
                collection=by_inst[label]; ch=collection[0] if len(collection)==1 else {}; asym=ch.get('chain_asym_id'); source_chains=chainref[asym] if asym else []; w=asym_entities[asym] if asym else []
                if len(collection)!=1:reasons.add('INDEX_CHAIN_MISSING_OR_AMBIGUOUS')
                if len(w)!=1:reasons.add('RCSB_ENTITY_FOR_ASYM_MISSING_OR_AMBIGUOUS')
                if asym and matrix_members.get(aid) and asym not in matrix_members[aid]:reasons.add('RCSB_ASSEMBLY_GENERATOR_LACKS_ASYM')
                if role=='LIGAND' and not source_chains:role_status['ligand_ref_without_entry_chains']+=1
                if len(collection)==1:
                    (pdb_proteins if role=='PROTEIN' else pdb_ligands).add((pdb,prev,aid,label))
                    if role=='LIGAND':matched_ligand_members.append(label)
                out_chains.append({'pdb_id':pdb,'source_revision':prev,'assembly_id':aid,'system_id':sid,'ligand_id':lid,'role_index':role,'instance_id':label,'plinder_asym_id':asym,'plinder_assembly_chain_matches':len(collection),'plinder_auth_entity_chains':source_chains,'rcsb_asym_entity_matches':[{k:x[k] for k in ('entity_id','entity_kind','auth_asym_ids','ccd','prd')} for x in w],'rcsb_assembly_generator_contains_asym':asym in matrix_members[aid] if asym and matrix_members.get(aid) else None})
        for rr in reasons:role_status[rr]+=1
        row_status.update(reasons)
        ct=sorted(contacts[sid,lid]); mode='CROSS_CHAIN_CONTACT' if len(ct)>1 else ('SINGLE_CHAIN_CONTACT' if ct else 'NO_INDEXED_PLI_CHAIN_EVIDENCE'); modes[mode]+=1
        if row['ligand_is_covalent'] is True:modes['COVALENT_FLAG']+=1
        if row['ligand_covalent_linkages']:modes['COVALENT_ENDPOINT_EVIDENCE']+=1
        category=row['ligand_molecule_type'] or 'UNKNOWN'; classes[category]+=1
        row.update(pdb_id=pdb,source_revision=prev,rcsb_latest_revision=rrev,rcsb_model_ids=ci.get('model_ids'),index_binding_mode=mode,pli_chain_instances=ct,index_metadata_review=sorted(reasons),structure_qc='NOT_RUN')
        out_pairs.append(row)
    gaps=[]
    for kind in expected:
        if expected[kind]!=observed[kind]:gaps.append(f'{kind}_entity_ids_metadata_gap')
    if set(ci.get('assembly_ids') or [])!=set(assembly_objs):gaps.append('assembly_metadata_gap')
    if e.get('_confirmed_unavailable'):gaps.append('entry_confirmed_unavailable')
    elif e.get('_missing'):gaps.append('entry_metadata_missing')
    else:
        if not info.get('polymer_entity_count_protein') or not method_list or not a.get('initial_release_date') or info.get('structure_determination_methodology')!='experimental' or not start_day <= a.get('initial_release_date','')[:10] < end_day:gaps.append('search_predicate_review')
    if not rrev:gaps.append('latest_revision_missing')
    if not model_history:gaps.append('structure_model_history_missing')
    if model_history and rrev and (latest:=max(model_history,key=lambda h:(h['major_revision'],h['minor_revision'],h.get('ordinal') or 0))) and (rev(latest)!=rrev or (latest.get('revision_date') or '')[:10]!=(a.get('revision_date') or '')[:10]):gaps.append('current_revision_history_difference')
    if row_status:problems.append({'pdb_id':pdb,'index_rows_review':dict(sorted(row_status.items()))})
    m=M[pid]; rel_day=(a.get('initial_release_date') or '')[:10]
    out_entries.append({'pdb_id':pdb,'title':title,'first_release_date':a.get('initial_release_date'),'deposition_date':a.get('deposit_date'),'publication_year':pub.get('year'),'publication_date':None,'pubmed_id':pub.get('pdbx_database_id_PubMed'),'methods':method_list,'resolution_angstrom':info.get('resolution_combined'),'atom_count':info.get('deposited_atom_count'),'model_ids':ci.get('model_ids'),'model_count':info.get('deposited_model_count'),'rcsb_latest_revision':rrev,'rcsb_latest_revision_date':a.get('revision_date'),'latest_revision_vs_cutoff':when(a.get('revision_date')),'revision_history':history,'selected_revision':None,'selected_revision_vs_cutoff':'UNKNOWN','rcsb_assembly_ids':ci.get('assembly_ids'),'protein_entity_count':info.get('polymer_entity_count_protein'),'source_entity_ids':{k:sorted(x) for k,x in expected.items()},'discoveries':sorted(tags),'nonpolymer_ccds':sorted(set(comp)),'partner_evidence':path_evidence,'family_evidence_types':sorted(family_labels),'plinder_index_presence':bool(pl),'plinder_annotated_assemblies':sorted({r['system_biounit_id'] for r in pl}),'plinder_annotation_rows':len(pl),'plinder_candidate_system_count':len(count_pl_systems) if pl else None,'plinder_source_revision_if_available':prev,'plinder_revision_status':status,'plinder_metadata_release_day':m[0]['entry_release_date'][:10] if len(m)==1 and m[0]['entry_release_date'] else None,'plinder_release_day_state':('MATCH' if len(m)==1 and m[0]['entry_release_date'] and rel_day==m[0]['entry_release_date'][:10] else ('DIFFERENT' if len(m)==1 and m[0]['entry_release_date'] and rel_day else 'UNKNOWN')),'metadata_review':gaps})

assert sum(len(ann[i]) for i in seen)==len(out_pairs) and len(out_entries)==len(ids)
valid_order=ids==[r['pdb_id'] for r in sorted(out_entries,key=lambda r:((r['first_release_date'] or '')[:10],r['pdb_id']))] if not any(r['first_release_date'] is None for r in out_entries) else None
boundary = [r for r in out_entries if (r['first_release_date'] or '')[:10] == cutoff_day]
assert valid_order is not False, 'Search and Data ordering disagree'
boundary_check = search['boundary_counts'].get(cutoff_day, {})
if valid_order is not None and boundary_check.get('usable'):
    assert len(boundary) == (boundary_check['count'] or 0), 'Search and Data cutoff counts disagree'
for name,rows in [('candidate_entries',out_entries),('plinder_ligand_rows',out_pairs),('indexed_chain_links',out_chains),('rcsb_assemblies',assemblies),('rcsb_entities',entity_rows)]:
    csvout(name+'.csv',rows)
    if name in ('candidate_entries','plinder_ligand_rows'):
        (out/(name+'.jsonl')).write_text(''.join(json.dumps(row,ensure_ascii=False,separators=(',',':'))+'\n' for row in rows),encoding='utf8')
if cutoff_day:
    csvout(f'boundary_{cutoff_day}.csv',[{'pdb_id':r['pdb_id'],'first_release_date':r['first_release_date']} for r in boundary],['pdb_id','first_release_date'])
csvout('plinder_period_only.csv',extra_rows,['entry_pdb_id','entry_release_date','entry_determination_method','entry_failed_assembly_ids','annotation_rows','candidate_systems','source_revisions','rcsb_protein_eligibility'])
dump('index_sources.json',index); dump('index_review.json',problems)
counts={'rcsb_entries':len(entries),'rcsb_entity_rows':len(entity_rows),'release_dates':dates,'cutoff_entry_count':len(boundary),'data_assembly_order_verified':valid_order,'methods':methods,'methods_with_resolution':method_resolution,'has_resolution':sum(bool(r['resolution_angstrom']) for r in out_entries),'has_atom_count':sum(r['atom_count'] is not None for r in out_entries),'has_model_count':sum(r['model_count'] is not None for r in out_entries),'has_publication_year':sum(r['publication_year'] is not None for r in out_entries),'has_deposition_date':sum(bool(r['deposition_date']) for r in out_entries),'has_revision_history':sum(bool(r['revision_history']) for r in out_entries),'revision_vs_boundary':Counter(r['latest_revision_vs_cutoff'] for r in out_entries),'metadata_review':Counter(g for r in out_entries for g in r['metadata_review']),'discovery_combinations':joint,'discovery_margins':Counter(x for r in out_entries for x in r['discoveries']),'discovery_details':trail,'family_entry_evidence':family_types,'protein_entity_evidence':family_entities,'family_top_values':{k:dict(v.most_common(12)) for k,v in family_values.items()},'rcsb_known_entry_revisions':len(revisions_seen),'plinder_known_entry_revisions_for_rcsb':len(source_rev_seen),'rcsb_assembly_pairs':len(ww_assemblies),'rcsb_assembly_entries':len({p for p,_ in ww_assemblies}),'plinder_annotation_rows':len(out_pairs),'plinder_revision_status':rev_count,'plinder_annotation_entries':sum(bool(r['plinder_index_presence']) for r in out_entries),'plinder_no_annotation_entries':rev_count['NO_ANNOTATION'],'plinder_same_period_metadata_entries':len(period),'plinder_same_period_annotated_entries':len(annotated_period),'plinder_only_calendar_subset':len(annotated_period-seen),'plinder_only_calendar_metadata':len(extra),'all_source_revision_comparisons':Counter('UNKNOWN' if not r['rcsb_latest_revision'] or not r['plinder_source_revision_if_available'] else ('MATCH' if r['rcsb_latest_revision']==r['plinder_source_revision_if_available'] else 'DIFFERENT') for r in out_entries),'unannotated_discovery_combinations':Counter('+'.join(r['discoveries']) for r in out_entries if not r['plinder_index_presence']),'plinder_undated_annotated_entries_global':len(set(ann)&undated),'plinder_release_day_state':Counter(r['plinder_release_day_state'] for r in out_entries),'plinder_assembly_pairs':len(pl_assemblies),'index_metadata_review_rows':sum(bool(r['index_metadata_review']) for r in out_pairs),'index_metadata_review_entries':len(problems),'index_reference_counts':role_status,'plinder_distinct_protein_assembly_instances':len(pdb_proteins),'plinder_distinct_ligand_assembly_chain_instances':len(pdb_ligands),'plinder_distinct_ligand_ids':len({r['ligand_id'] for r in out_pairs}),'plinder_candidate_systems':len(systems),'indexed_ligand_row_classes':classes,'pli_binding_evidence':modes,'plinder_proper_flag':Counter(str(r['ligand_is_proper']) for r in out_pairs)}
dump('inventory_summary.json',{'started_utc':started.isoformat(),'finished_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':round(time.monotonic()-t0,3),'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'python':sys.version.split()[0],'pyarrow':pa.__version__,'counts':counts,'entry_metadata_sha256':sha256(run/'data'/'entries.json')})
print(json.dumps(counts,ensure_ascii=False),flush=True)
