"""Evidence and structure stages for the isolated blind M2 batch."""
from collections import Counter, defaultdict
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import csv
import glob
import json
from pathlib import Path
import re
import time
from sys import getsizeof

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from holo_candidate_batch import atomic_json, now, revision_from_cif, sha256
from holo_native_candidates import columns, inspect_source, key_for
from holo_role_evidence import assess, build_role_context, role_source
from holo_source_mapping import altloc, assembly_instances, MappingError, xyz
from holo_v1_export import CONN, PAIR, COVALENT, choose_connected, describe_instances, deposited_copy_matches, export, output_label, rows, text

VERSION = 'blind-m2-assembly-evidence-2'
INDEXCOLS = ('entry_pdb_id system_biounit_id system_id ligand_id ligand_instance_chains ligand_residue_numbers '
             'ligand_unique_ccd_code ligand_ccd_code ligand_member_asym_ids ligand_protein_chains_asym_id '
             'system_protein_chains_asym_id ligand_is_proper ligand_is_covalent ligand_molecule_type '
             'ligand_is_cofactor_like ligand_is_artifact ligand_is_ion ligand_covalent_linkages').split()
ROLE_CUE = re.compile(r'\b(in complex with|complexed with|complex with|in the presence of|bound to|bound with|bound by|binding of|with bound|soaked with)\b')
COMPONENT_CONTEXT = re.compile(r'\b(substrate|inhibitor|ligand|agonist|antagonist|activator|drug)\b')
GENERIC = {'alpha','beta','protein','peptide','dna','rna','ligand','sugar','ion','water','a','b','c','h','atom','complex','amino acid'}
CONTEXT_NON_TARGET = {'HOH','DOD'}


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def normalize(s):
    return ' '.join(re.findall('[a-z0-9]+',str(s or '').lower()))


def explicit_role_names(block):
    out = defaultdict(list)
    fields = ('name','pdbx_synonyms')
    for r in rows(block,'_chem_comp.'):
        for field in fields:
            for name in text(r.get(field)).split(';'):
                out[text(r['id'])].append((field,name.strip()))
    for r in rows(block,'_pdbx_entity_nonpoly.'):
        out[text(r.get('comp_id'))].append(('entity_nonpoly.name',text(r.get('name'))))
    return out


def find_named_role(title, names, match_ccd=None):
    t=normalize(title)
    if match_ccd and len(match_ccd)>=3 and any(ch.isdigit() for ch in match_ccd):
        names=[*names,('literal CCD',match_ccd)]
    for field, raw in names:
        key=normalize(raw)
        if len(key)<5 or key in GENERIC or not re.search(r'\b'+re.escape(key)+r'\b',t):
            continue
        preceded=any(re.search(r'\b'+re.escape(key)+r'\b', t[m.end():m.end()+215]) for m in ROLE_CUE.finditer(t))
        # Prefix before "bound", e.g., "glycine bound". No role inferred from contact or index flags.
        postfix=bool(re.search(r'\b'+re.escape(key)+r'\s+(?:bound|as (?:a |the )?(?:substrate|inhibitor|ligand))\b',t))
        if preceded or postfix:
            return {'field':'_struct.title','source_name_field':field,'source_name':raw,'text':title}
    return None


def review_admission_evidence(row, title, names, polylen):
    """Keep explicitly reported partners for review when a size filter excludes them."""
    if len(row['members']) != 1 or not row['contact_instances_4A']:
        return None
    member = row['members'][0]
    if (member['kind'] == 'polymer'
            and str(member.get('polymer_type') or '').startswith('polypeptide')
            and polylen.get(member['entity'], 0) > 35):
        description = normalize(member.get('description')).split()
        title_text = normalize(title)
        for cue in ROLE_CUE.finditer(title_text):
            peptide = re.match(r'\s*(?:an? |the )?peptide (?:of|from|derived from) ([a-z0-9]+)\b',
                               title_text[cue.end():])
            if peptide and len(peptide[1]) >= 4 and peptide[1] in description:
                return {'kind': 'source_named_peptide', 'field': '_struct.title', 'text': title,
                        'matched_entity_description': member.get('description'),
                        'declared_length': polylen[member['entity']]}
            # SEQRES may describe an entire construct while only a named tail is observed.
            fragment = re.match(r'\s*(?:an? |the )?([a-z][a-z0-9]{3,}) (tail|fragment|peptide)\b', title_text[cue.end():])
            observed=[r['label_seq'] for r in member['residues'] if r['positive_heavy_atoms']]
            if (fragment and fragment[1] in description and 2 <= len(observed) <= 35
                    and all(isinstance(x,int) for x in observed)
                    and sorted(set(observed))==list(range(min(observed),max(observed)+1))):
                return {'kind':'source_named_observed_fragment', 'field':'_struct.title + _atom_site.label_seq_id',
                        'text':title,'matched_entity_description':member.get('description'),
                        'declared_length':polylen[member['entity']],'positive_observed_label_seq':sorted(set(observed))}
    if member['kind'] == 'component' and len(member['residues']) == 1:
        residue = member['residues'][0]
        if residue['positive_heavy_atoms'] == 1:
            for field, source_name in names.get(residue['ccd'], ()):
                name = normalize(source_name)
                if not name.endswith(' ion'):
                    continue
                match = find_named_role(title, [(field, name.removesuffix(' ion'))])
                if match:
                    return match | {'kind': 'source_named_single_atom_ion', 'ccd': residue['ccd'],
                                    'source_name': source_name}
    return None


def role(row, title, names, polylen):
    members=row['members']
    if not members or any(m['kind']=='polymer_modification' for m in members):
        return 'REVIEW',{'reason':'Modification needs contextual target-role evidence.'}
    if any(m['kind']=='polymer' and not str(m.get('polymer_type','')).startswith('polypeptide') for m in members):
        return 'REVIEW',{'reason':'Nucleic-acid role and complete partner membership require independent evidence.'}
    hits=[]
    for m in members:
        if m['kind']=='polymer':
            if polylen.get(m['entity'],999)>35:
                continue
            match=find_named_role(title,[('_entity.pdbx_description',m.get('description'))])
        elif m['kind']=='branched':
            match=find_named_role(title,[('_entity.pdbx_description',m.get('description'))])
        elif m['kind']=='component' and len(m['residues'])==1:
            cc=m['residues'][0]['ccd']; match=find_named_role(title,names[cc],cc)
        else:
            match=None
        if match:
            hits.append(match)
    if len(hits)==len(members) and row['contact_instances_4A']:
        # A co-reported component could be a substrate, inhibitor, or an environmental additive.
        # Do not promote named detergent/solvent/cofactor by mere co-reporting in a title.
        if any(re.search(r'detergent|digitonin|glycerol|ethylene glycol|maltoside|sulfate|chloride|nanodisc|phosphate buffer|polyethylene glycol|hydrogen sulfate|octyl glucoside|he(m|me)|iron sulfur|cholesterol|phosphatid',normalize(x['source_name'])) for x in hits):
            return 'REVIEW',{'reason':'Named component has competing environment/cofactor and target interpretations.','matches':hits}
        return 'PASS',{'reason':'Source depositor title explicitly names matching observed partner; positive compatible protein contact.','matches':hits}
    return 'REVIEW',{'reason':'No unambiguous depositor target-role statement for the entire observed object.'}


def points(rr, op, subset=None):
    out=defaultdict(list)
    for author,res in rr:
        if subset is not None and (author,str(res.seqid),res.name) not in subset:
            continue
        for atom in res:
            if atom.occ>0 and not atom.element.is_hydrogen:
                out[altloc(atom) or ''].append(xyz(atom.pos))
    m=np.asarray(op.transform.mat.tolist()) if op else None
    v=np.asarray(op.transform.vec.tolist()) if op else None
    return {a:(np.asarray(x) @ m.T + v if op else np.asarray(x)) for a,x in out.items()}


class ContactContext:
    """The geometry gate is per deposited model and compatible altloc (same or blank)."""
    def __init__(self, instances, protein_entities, polylen):
        self.instances=instances
        self.proteins={k for k,(_,_,rr) in instances.items() if str(rr[0][1].entity_id) in protein_entities}
        self.receptors={k for k in self.proteins if polylen.get(str(instances[k][2][0][1].entity_id),999)>35} or self.proteins
        self.grouped={k:points(instances[k][2],instances[k][0]) for k in self.receptors}
        self.trees={}
        for k,g in self.grouped.items():
            if g:
                lo=np.min(np.concatenate(list(g.values())),axis=0)
                hi=np.max(np.concatenate(list(g.values())),axis=0)
                self.trees[k]=(lo,hi,cKDTree(np.concatenate(list(g.values()))),{a:cKDTree(p) for a,p in g.items()})

    def contacts(self,members):
        coords=defaultdict(list); own=set()
        for m in members:
            k=m['instance']; own.add(k)
            op,label,rr=self.instances[k]
            select={(r['auth_asym'],r['auth_seq'],r['ccd']) for r in m['residues']}
            source=self.grouped[k] if k in self.grouped and len(select)==len(rr) else points(rr,op,select)
            for a,p in source.items():
                coords[a].append(p)
        qry={a:np.concatenate(p) for a,p in coords.items()}
        hits={}
        if not qry:
            return hits
        allcoords=np.concatenate(list(qry.values())); lo=allcoords.min(axis=0); hi=allcoords.max(axis=0)
        for k,(plo,phi,alltree,trees) in self.trees.items():
            if k in own or np.any(plo-hi>4.) or np.any(lo-phi>4.):
                continue
            values=[]
            for a, arr in qry.items():
                search=[alltree] if not a else [trees[x] for x in ('',a) if x in trees]
                for tree in search:
                    x=tree.query(arr,k=1,distance_upper_bound=4.00000001)[0].min()
                    if x<=4.:
                        values.append(float(x))
            if values:
                hits[k]=round(min(values),6)
        return hits


def instances_for(st,a,model):
    if a:
        return assembly_instances(st,a,model)
    from holo_source_mapping import subchains
    return {f'asu.{k}':(None,k,v) for k,v in subchains(model).items()}


def same_source_aliases(native, index, biounit, revision_match, contexts):
    if not revision_match:
        return {}, []
    lookup=defaultdict(list)
    for r in native:
        # A source-native group may subsume more than one index member; record containment distinctly.
        for m in r['members']:
            lookup[(r['assembly'],r['model'],m['instance'])].append(r)
    bindings=defaultdict(list); result=[]
    for i in index:
        aid=str(i['system_biounit_id']); ids=set(i['ligand_instance_chains'] or []); partners=set(i['ligand_protein_chains_asym_id'] or [])
        states=[]
        for (a,model), at in contexts.items():
            if a!=aid or not ids or not ids<=at.instances.keys() or not partners<=at.proteins:
                continue
            if any(biounit.get((aid,k))!=at.instances[k][1] for k in ids|partners):
                continue
            ccd=Counter((i['ligand_ccd_code'] or '').split('-'))
            seen=set()
            for one in ids:
                for r in lookup[(a,model,one)]:
                    key=r['candidate_key']
                    if key in seen:
                        continue
                    seen.add(key)
                    by=[m for m in r['members'] if m['instance'] in ids]
                    if {m['instance'] for m in by}!=ids:
                        continue
                    identity=Counter(v['ccd'] for m in by for v in m['residues'])
                    labelseq=[v['label_seq'] for m in by for v in m['residues']]
                    same=(identity==ccd)
                    rawnums=i['ligand_residue_numbers'] or []
                    # For partial polymer/branched maps, or inconsistent cardinality, do not promote.
                    if not same or (len(rawnums) != sum(identity.values())):
                        continue
                    # Where CIF label sequence is defined for every member, demand its exact set.
                    if all(x is not None for x in labelseq) and Counter(labelseq)!=Counter(rawnums):
                        continue
                    match={'candidate_key':key,'model':model,'same_complete_partner':len(by)==len(r['members']),
                           'index_contact_alignment':sorted(partners & set(r['contact_instances_4A']))}
                    states.append(match)
        unique={x['candidate_key']:x for x in states}
        chosen=next(iter(unique.values())) if len(unique)==1 else None
        r={'assembly':aid,'system_id':i['system_id'],'ligand_id':i['ligand_id'],'instance_ids':sorted(ids),
           'proper':i['ligand_is_proper'],'ccd':i['ligand_ccd_code'],'index_protein_instances':sorted(partners),
           'status':'SOURCE_INSTANCE_MATCHED' if chosen and chosen['same_complete_partner'] else 'SOURCE_MEMBER_CONTAINED' if chosen else 'AMBIGUOUS_OR_UNMATCHED',
           'matches':list(unique.values())}
        if chosen:
            bindings[chosen['candidate_key']].append(r)
        result.append(r)
    return bindings,result


def nucleic_evidence(native, source, contexts, paired):
    """Group observed DNA/RNA members only where deposited pair instances resolve."""
    def is_na(r):
        return any(m['kind']=='polymer' and any(z in str(m.get('polymer_type','')).lower() for z in ('deoxyribonucleotide','ribonucleotide')) for m in r['members'])
    grouped=defaultdict(list)
    for r in native:
        grouped[(r['assembly'],r['model'])].append(r)
    out=[]
    for ctx, rr in grouped.items():
        labels=defaultdict(list)
        for r in rr:
            if is_na(r):
                for m in r['members']:
                    labels[m['label_asym']].append((r,m['operator'],m))
        link=defaultdict(set); evidence=defaultdict(list); any_pair=defaultdict(set)
        source_instances = {key:{'label_asym':label,'operator':op.name if op else None}
                            for key,(op,label,_) in contexts[ctx].instances.items()}
        row_grouping = {}
        for a in paired:
            if text(a.get('model_number')) not in ('',str(ctx[1])):
                continue
            ends = [text(a.get(f'{side}_label_asym_id')) for side in ('i','j')]
            left,right = (labels[e] for e in ends)
            for side in (left,right):
                for r,_,m in side:
                    any_pair[r['candidate_key']].add(text(a.get('pair_number')) or 'source base pair')
            syms = [text(a.get(f'{side}_symmetry')) for side in ('i','j')]
            is_nonlocal = any(sym not in ('','1_555') for sym in syms)
            requested = deposited_copy_matches(source,source_instances,syms if is_nonlocal else ['1_555']*2,ends)
            for choices in requested:
                if any(len(opts)!=1 for opts in choices):
                    continue
                key_l,key_r = (side[0][0] for side in choices)
                if key_l==key_r:
                    continue
                lmatches=[(r,m) for r,_,m in left if m['instance']==key_l]
                rmatches=[(r,m) for r,_,m in right if m['instance']==key_r]
                if len(lmatches)!=1 or len(rmatches)!=1:
                    continue
                (r,m),=lmatches; (other,n),=rmatches
                if other is r or not all(any(text(a.get(f'{side}_label_seq_id')) == str(z['label_seq'])
                                and text(a.get(f'{side}_label_comp_id')) == z['ccd'] for z in x['residues'])
                                            for side,x in [('i',m),('j',n)]):
                    continue
                l,k=r['candidate_key'],other['candidate_key']; link[l].add(k); link[k].add(l)
                key=frozenset((l,k)); evidence[key].append(a)
                row_grouping[key]=row_grouping.get(key,False) or is_nonlocal
        by={r['candidate_key']:r for r in rr}; visited=set()
        for k,r in by.items():
            if k in visited:
                continue
            q=list(link[k]); connected={k}
            while q:
                n=q.pop()
                if n not in connected:
                    connected.add(n); q.extend(link[n]-connected)
            visited |= connected
            if len(connected)>1:
                group=[by[x] for x in sorted(connected)]
                members=[m for x in group for m in x['members']]; identities=[m for x in group for m in x['source_member_identity']]
                r=r | {'members':members,'source_member_identity':identities,
                       'candidate_key':key_for(r['entry'],r['revision'],r['model'],r['assembly'],identities),
                       'member_grouping':('explicit_deposited_symmetry_base_pairs' if any(v and pair<=connected for pair,v in row_grouping.items())
                                          else 'explicit_same_operation_local_base_pairs'),
                       'member_grouping_rows':sum(len(v) for pair,v in evidence.items() if pair<=connected),
                       'connections':[m for x in group for m in x['connections']],
                       'unresolved_connection_ids':sorted({m for x in group for m in x['unresolved_connection_ids']})}
            elif is_na(r):
                r['member_grouping']='REVIEW_UNRESOLVED_SOURCE_PAIRS' if any_pair[k] else 'REVIEW_NO_EXPLICIT_SOURCE_PAIR_EVIDENCE'
            out.append(r)
    return out


def prepare_index(run, plan, index):
    import pyarrow.parquet as pq
    names=('annotation_table','entry_sources','entry_biounit_chains')
    files={n:index/f'{n}.parquet' for n in names}
    stamp={n:{'bytes':f.stat().st_size,'mtime_ns':f.stat().st_mtime_ns} for n,f in files.items()}
    dst=run/'index_snapshot.json'
    if dst.exists():
        old=read(dst)
        if old['source_tables']==stamp and len(old['pdb'])==len(plan):
            return old['pdb'], 'reused'
    ids=[r['pdb_id'].lower() for r in plan]; selected={x.upper():{'annotation':[],'source_revision':None,'biounit_chains':[]} for x in ids}
    fields={'annotation_table':INDEXCOLS,'entry_sources':['entry_pdb_id','source_mmcif_major_revision','source_mmcif_minor_revision'],
            'entry_biounit_chains':['entry_pdb_id','biounit_id','chain_instance','chain_asym_id','chain_role']}
    for n in names:
        for r in pq.read_table(files[n],columns=fields[n],filters=[('entry_pdb_id','in',ids)]).to_pylist():
            one=selected[r['entry_pdb_id'].upper()]
            if n=='annotation_table':
                one['annotation'].append(r)
            elif n=='entry_sources':
                major,minor=r['source_mmcif_major_revision'],r['source_mmcif_minor_revision']
                if major is not None and minor is not None:
                    rev=f'{major}.{minor}'
                    if one['source_revision'] not in (None,rev):
                        raise ValueError(f'Inconsistent PLINDER entry_sources for {r["entry_pdb_id"]}')
                    one['source_revision']=rev
            else:
                one['biounit_chains'].append(r)
    atomic_json(dst,{'created_utc':now(),'source_tables':stamp,'pdb':selected})
    return selected,'built'


def member_locus(r):
    return (r['assembly'],tuple(sorted((m['instance'], m['entity'],m['kind'],
                             tuple((v['auth_asym'],v['auth_seq'],v['ccd']) for v in m['residues'])) for m in r['members'])))


def one(run, spec, index, fingerprint, shard=None, split_above=None):
    start=time.monotonic(); p=spec['pdb_id']; ver=spec['revision']; proof=read(run/'source_proofs'/f'{p}.json')
    rr=run/'entry_results'/f'{p}.json'; ledger=run/'discovery'/f'{p}.json'; sysfile=run/'system_results'/f'{p}.json'
    result={'pdb_id':p,'revision':ver,'index_source_revision':index['source_revision'],'index_rows':len(index['annotation']),
            'indexed_systems':len({r['system_id'] for r in index['annotation']}),'signature':fingerprint}
    same=bool(index['annotation']) and index['source_revision']==ver
    stratum='INDEX_MATCH' if same else 'NO_LIGAND_ANNOTATION' if not index['annotation'] else 'INDEX_VERSION_MISMATCH'
    result['stratum']=stratum
    if proof['status']!='PASS':
        result.update(status='TECHNICAL_FAILURE',failure_stage='fetch',error=proof['error'])
        atomic_json(rr,result);return p,'failed_fetch'
    raw=run/proof['relative_path']; result.update(raw=proof['relative_path'],source_sha256=proof['sha256'],source_bytes=proof['bytes'])
    if rr.is_file() and prior_reusable(run, rr, ledger, sysfile, result):
        return p,'reused_verified'
    started=now()
    try:
        source, got=revision_from_cif(raw,p)
        if got!=ver:
            raise ValueError(f'Fixed revision {ver}, got {got}')
        st=gemmi.make_structure_from_block(source)
        if not len(st):
            raise ValueError('No readable model')
        models=[m.num for m in st]; ais=[a.name for a in st.assemblies] or [None]
        title=text(source.find_value('_struct.title')); names=explicit_role_names(source)
        types=columns(source,'_entity_poly.','entity_id','type'); protein={k for k,v in types.items() if str(v).startswith('polypeptide')}
        polylen=Counter(str(text(x.get('entity_id'))) for x in rows(source,'_entity_poly_seq.'))
        role_meta=role_source(source); role_context=build_role_context(role_meta,names)
        declared=text(source.find_value('_pdbx_nmr_representative.conformer_id'))
        rep=int(declared) if declared.isdigit() and int(declared) in models else models[0]
        native=inspect_source(source,st,p,ver,{},calculate_distances=False,preserve_protein_context=True)
        contexts={(aid,m.num):ContactContext(instances_for(st,aid,m),protein,polylen) for m in st for aid in ais}
        native=nucleic_evidence(native,source,contexts,rows(source,PAIR))
        for n in native:
            ctx=contexts[n['assembly'],n['model']]
            n['contact_distances_A']=ctx.contacts(n['members']);n['contact_instances_4A']=sorted(n['contact_distances_A'])
            n['contact_definition']='heavy atoms with positive occupancy, same model and identical or blank altloc, <= 4 A; excludes optional short receptor heuristics'
            n['protein_assembly_instances']=sorted(ctx.proteins)
            n['role_qc'], n['role_basis']=role(n,title,names,polylen)
            addition=assess(n,n['role_qc'],title,names,polylen,role_meta,role_context)
            if addition:
                n['role_qc'], n['role_basis']=addition
            n['role']={'PASS':'target_candidate_with_depositor_evidence',
                       'NON_TARGET_EVIDENCED':'explicit_experimental_preparation_role'}.get(n['role_qc'],'pending')
        bit={(str(x['biounit_id']),x['chain_instance']):x['chain_asym_id'] for x in index['biounit_chains']}
        aliases, idx=same_source_aliases(native,index['annotation'],bit,same,contexts)
        idx_summary=dict(Counter(x['status'] for x in idx))
        route='INDEX_SOURCE_MATCH_WITH_NATIVE_COVERAGE' if aliases else 'NATIVE_ONLY_INDEX_UNRESOLVED' if same else 'NATIVE_ONLY_NO_ANNOTATION' if not index['annotation'] else 'NATIVE_ONLY_DIFFERENT_REVISION'
        source_by_key={n['candidate_key']:n for n in native}
        by_locus=defaultdict(list)
        for n in native:
            by_locus[member_locus(n)].append(n)
            n['index_aliases']=[{'system_id':r['system_id'],'ligand_id':r['ligand_id'],'status':r['status']} for r in aliases.get(n['candidate_key'],[])]
        decisions=[]
        for scope,nn in by_locus.items():
            contacted=[x for x in nn if x['contact_instances_4A']]
            preferred=next((x for x in contacted if x['model']==rep),None) or (contacted[0] if contacted else nn[0]);
            important=any(x['role_qc']=='PASS' for x in nn)
            key=preferred['candidate_key']; m=preferred['members']; contexts_to_check=[x for x in nn if x['model']==preferred['model']]
            index_proper=any(i['proper'] for x in nn for i in aliases.get(x['candidate_key'],[]))
            short=all(x['kind']!='polymer' or not str(x.get('polymer_type','')).startswith('polypeptide') or 2<=polylen.get(x['entity'],999)<=35 for x in m)
            chemistry=any(x['kind'] in ('component','branched') and any(r['positive_heavy_atoms']>=5 for r in x['residues']) for x in m)
            partner=any(x['kind']=='branched' or x['kind']=='polymer' and short for x in m)
            role_only_modification=all(x['kind']=='polymer_modification' for x in m)
            admissible=bool(contacted and not role_only_modification and (important or index_proper or chemistry or partner or preferred['connections']))
            why='named depositor role claim' if important else 'matched proper index clue plus source contact; role pending' if index_proper else 'native chemical/polymer/connection evidence with source contact; role pending'
            review_evidence = None
            if not admissible and not role_only_modification:
                review_evidence = review_admission_evidence(preferred, title, names, polylen)
                if review_evidence:
                    admissible = True
                    why = 'Explicit source partner clue and source contact; admitted for structure review only.'
            if not contacted:
                why='No compatible positive-occupancy direct protein contact in observed models.'
            elif not admissible:
                why='Observed but insufficient standalone partner evidence; retained for role and omission review.'
            decision = {'candidate_key':key,'assembly':preferred['assembly'],'model':preferred['model'],
                        'observed_models':sorted({n['model'] for n in nn}),
                        'models_with_contacts':sorted({n['model'] for n in contacted}),
                        'eligible_for_structure_review':admissible,'reason':why}
            if review_evidence:
                decision['review_admission_evidence'] = review_evidence
            decisions.append(decision)
        atomic_json(ledger,{'pdb_id':p,'revision':ver,'signature':fingerprint,'finished_utc':now(),'route':route,'stratum':stratum,
                            'models':models,'representative_model':rep,'source_title':title,'index_mapping':idx,
                            'unmapped_index_policy':'Native discovery was performed on all source contexts. Different/unknown index source revisions were not used for assembly/instance assignment.',
                            'native_observations':native,'candidate_selection':decisions})
        if split_above is not None and sum(d['eligible_for_structure_review'] for d in decisions)>=split_above:
            return p,'needs_slices'
        exported=[]; assembly_cache={}; mapping={}; complete={}; pending_export=[]
        for choice in decisions:
            n=source_by_key[choice['candidate_key']]; key=n['candidate_key']
            if not choice['eligible_for_structure_review'] or (shard and int(key[:16],16)%shard[1]!=shard[0]):
                continue
            aid=n['assembly']; model=n['model']; ctx=contexts[aid,model]; used={m['instance'] for m in n['members']}; binding=set(n['contact_instances_4A']); labels=used | binding
            for edge in n['connections']:
                if edge['source']['kind'] in COVALENT:
                    labels.update(edge['instances'])
            out={'pdb_id':p,'candidate_key':key,'revision':ver,'assembly':aid,'representative_model':model,
                 'model_scope':choice['observed_models'],'partner_instances':sorted(used),'direct_binding_chains':sorted(binding),
                 'assembly_protein_instances':n['protein_assembly_instances'],'role_qc':n['role_qc'],'role_basis':n['role_basis'],
                 'coordinate_identity':'NOT_RUN','experimental_fidelity':'NOT_RUN','topology':'NOT_RUN','candidate_intake':'REVIEW',
                 'index_aliases':n['index_aliases'],'member_grouping':n.get('member_grouping'),
                 'connections_from_source':n['connections'],'unresolved_connections':n['unresolved_connection_ids']}
            if 'review_admission_evidence' in choice:
                out['review_admission_evidence'] = choice['review_admission_evidence']
            pending_export.append((choice,n,ctx,binding,labels,out))
        # One source-row materialization per used context, shared by every projected system.
        from holo_v1_export import ATOM, expanded_atoms, expanded_conn_context
        used_ctx={(n['assembly'],n['model']) for _,n,_,_,_,_ in pending_export if n['assembly']}
        source_atoms=rows(source,ATOM) if used_ctx else []
        base_by_model=defaultdict(list)
        for arow in source_atoms:
            base_by_model[int(text(arow['pdbx_PDB_model_num']))].append(arow)
        for choice,n,ctx,binding,labels,out in pending_export:
            aid=n['assembly']; model=n['model']; key=n['candidate_key']; m=st[models.index(model)]
            def ensure_context():
                ck=(aid,model)
                if ck in assembly_cache:
                    v=assembly_cache[ck]
                    if isinstance(v,Exception):raise v
                    return v
                try:
                    instance_map, expansion=describe_instances(st,aid,m)
                    exp=instance_map if expansion else None
                    atoms=expanded_atoms(source,exp,source_atoms=base_by_model[model]) if exp else base_by_model[model]
                    links=expanded_conn_context(source,exp,atoms,{str(model)}) if exp else None
                    assembly_cache[ck]=(instance_map,exp,atoms,links)
                    return assembly_cache[ck]
                except MappingError as exc:
                    assembly_cache[ck]=exc; raise
            try:
                if not aid:
                    raise MappingError('assembly_missing','No explicit Assembly; ASU observation withheld from V1.')
                instance_map,exp,prepared,prepared_connections=ensure_context()
                if n['unresolved_connection_ids'] and not exp:
                    raise MappingError('unresolved_partner_connection','Source partner connection does not map unambiguously.',ids=n['unresolved_connection_ids'])
                scope={output_label(instance_map[k]) for k in labels}
                if not exp:
                    scope=choose_connected(source,scope)
                    if not scope<={output_label(t) for t in instance_map.values()}:
                        raise MappingError('topology_crosses_assembly_scope','Covalent closure would add a non-member of this Assembly.')
                system_path=run/'systems'/p/key/'system.cif'
                artifact,_=export(system_path,source,p,aid,[model],scope,expanded=exp,prepared=prepared,prepared_connections=prepared_connections)
                artifact['omitted_models']=sorted(set(models)-{model})
                if artifact['path'].startswith(str(run)):
                    artifact['path']=str(system_path.relative_to(run))
                mapping_data=artifact['output_instances'] if exp else {k:v for k,v in instance_map.items() if output_label(v) in scope}
                atomic_json(system_path.parent/'source_mapping.json',{'pdb_id':p,'revision':ver,'candidate_key':key,'assembly':aid,
                                 'source_model':model,'source_proof':f'source_proofs/{p}.json','source_sha256':proof['sha256'],
                                 'instances':mapping_data,'partner_member_identities':n['source_member_identity'],
                                 'binding_protein_instances':sorted(binding),'artifact':artifact})
                out.update(coordinate_identity='PASS',experimental_fidelity='PASS',topology='PASS',artifact=artifact,
                           mapping=str((system_path.parent/'source_mapping.json').relative_to(run)),
                           status='EXPORTED_READBACK',partner_scope_qc='REVIEW' if str(n.get('member_grouping','')).startswith('REVIEW') or
                              any(x['kind']=='polymer' and not str(x.get('polymer_type','')).startswith('polypeptide') for x in n['members']) else 'PASS')
                if out['role_qc']=='PASS' and out['partner_scope_qc']=='PASS':
                    out['candidate_intake']='PASS'
                elif out['role_qc']=='NON_TARGET_EVIDENCED':
                    out['candidate_intake']='EXCLUDED'
                if exp or len(models)>1:
                    ckey=(aid,model)
                    if ckey not in complete:
                        try:
                            whole={output_label(x) for x in instance_map.values()}; dst=run/'assemblies'/p/f'a{aid}-m{model}.cif'
                            c,_=export(dst,source,p,aid,[model],whole,complete=True,expanded=exp,prepared=prepared,prepared_connections=prepared_connections)
                            c['omitted_models']=sorted(set(models)-{model}); c['path']=str(dst.relative_to(run))
                            evidence_file=dst.with_suffix('.evidence.json')
                            atomic_json(evidence_file,c)
                            complete[ckey]={'status':'PASS','artifact':{k:c[k] for k in ('path','bytes','sha256')},
                                           'shared_evidence':str(evidence_file.relative_to(run)),
                                           'shared_evidence_sha256':sha256(evidence_file)}
                        except (MappingError,ValueError,RuntimeError) as exc:
                            complete[ckey]={'status':'REVIEW','reason':getattr(exc,'kind',type(exc).__name__),'detail':str(exc)}
                    out['complete_assembly']=complete[ckey]
                if len(choice['observed_models'])>1:
                    out['experimental_models_retrieval']={'source_model_ids':models,'source_raw':result['raw'],
                        'export_models':[model],'other_model_contact_ids':choice['models_with_contacts'],
                        'note':'Other deposited models retained in fixed original; representative CIF covers the one declared source model.'}
            except MappingError as exc:
                out.update(status='REVIEW',topology='REVIEW',reason=exc.kind,detail=str(exc),mapping_error_evidence=exc.evidence)
            except (ValueError,RuntimeError,KeyError) as exc:
                out.update(status='EXPORT_TECHNICAL_FAILURE',reason=type(exc).__name__,detail=str(exc))
            exported.append(out)
        atomic_json(sysfile,exported)
        usable=sum(r['candidate_intake']=='PASS' for r in exported); built=sum(r['status']=='EXPORTED_READBACK' for r in exported)
        result.update(status='HAS_QC_PASS_CANDIDATE' if usable else 'REVIEW',route=route,
                      models=models,assemblies=ais,observations=len(native),distinct_candidate_loci=len(decisions),
                      chosen_systems=len(exported),structure_readback_pass=built,locally_supported_candidates=usable,
                      system_status_counts=dict(Counter(x['status'] for x in exported)),index_mapping=idx_summary,
                      models_with_observations=sorted({x['model'] for x in native}),started_utc=started,finished_utc=now(),
                      elapsed_seconds=round(time.monotonic()-start,3),
                      review_reason=None if usable else ('Candidate roles or topology require review.' if exported else 'No confident standalone structure selection; check for missed targets.'))
        result['discovery_sha256']=sha256(ledger);result['systems_sha256']=sha256(sysfile)
        atomic_json(rr,result)
        return p,'built'
    except Exception as exc:
        import traceback
        result.update(status='TECHNICAL_FAILURE',failure_stage='discovery_or_export',error=f'{type(exc).__name__}: {exc}',traceback=traceback.format_exc(),elapsed_seconds=round(time.monotonic()-start,3))
        atomic_json(rr,result)
        return p,'technical_failure'


def prior_reusable(run, file, ledger, sysfile, signature):
    old=read(file)
    if (old['status']=='TECHNICAL_FAILURE' or old['signature']!=signature['signature'] or old.get('source_sha256')!=signature['source_sha256']
            or not ledger.is_file() or not sysfile.is_file() or sha256(ledger)!=old['discovery_sha256'] or sha256(sysfile)!=old['systems_sha256']):
        return False
    checked_assembly=set()
    for a in read(sysfile):
        if a['status']=='EXPORT_TECHNICAL_FAILURE':
            return False
        if 'artifact' in a and (not (run/a['mapping']).is_file() or not artifact_ok(run,a['artifact'])):
            return False
        assembly=a.get('complete_assembly',{})
        if assembly.get('status')=='PASS':
            entry=assembly['artifact']; key=entry['path']
            if key not in checked_assembly:
                if not artifact_ok(run,entry): return False
                if assembly.get('shared_evidence'):
                    evidence=run/assembly['shared_evidence']
                    if not evidence.is_file() or sha256(evidence)!=assembly['shared_evidence_sha256']: return False
                checked_assembly.add(key)
    return True


def artifact_ok(run,a):
    p=run/a['path']
    return p.is_file() and p.stat().st_size==a['bytes'] and sha256(p)==a['sha256']


def summarize(run,plan,invocation):
    entries=[]; flat=[]; raw=[]; missing=[]
    for r in plan:
        p=r['pdb_id']; path=run/'entry_results'/f'{p}.json'; sec=run/'system_results'/f'{p}.json'; proof=run/'source_proofs'/f'{p}.json'
        if path.exists():
            ent=read(path); entries.append(ent)
            if sec.exists() and ent['status']!='TECHNICAL_FAILURE':
                flat+=read(sec)
        else:
            missing.append(p)
        if proof.is_file():
            raw.append(read(proof))
    def csv_out(name, source, preferred):
        tmp=run/f'{name}.tmp'; fields=list(dict.fromkeys([*preferred,*(k for x in source for k in x if k!='traceback')]))
        with tmp.open('w',newline='',encoding='utf8') as f:
            writer=csv.DictWriter(f,fields,extrasaction='ignore'); writer.writeheader()
            for r in source:
                writer.writerow({k:json.dumps(v,ensure_ascii=False,separators=(',',':')) if isinstance(v,(list,dict)) else v for k,v in r.items()})
        tmp.replace(run/name)
    csv_out('entry_results.csv',entries,['pdb_id','revision','stratum','status','route','index_source_revision','observations','chosen_systems','structure_readback_pass','locally_supported_candidates'])
    # Keep compact aggregate exports; per-Entry JSON supplies atom paths/connection evidence.
    compact=[{k:v for k,v in x.items() if k not in ('artifact','complete_assembly','connections_from_source','role_basis')} |
             {'system_path':x.get('artifact',{}).get('path'),'system_sha256':x.get('artifact',{}).get('sha256'),
              'complete_assembly_status':x.get('complete_assembly',{}).get('status'),
              'complete_assembly_path':x.get('complete_assembly',{}).get('artifact',{}).get('path')} for x in flat]
    csv_out('system_results.csv',compact,['pdb_id','candidate_key','status','coordinate_identity','experimental_fidelity','topology','role_qc','candidate_intake','system_path'])
    sums=defaultdict(Counter)
    for x in entries:
        c=sums[x['stratum']]; c['entries']+=1;c[x['status']]+=1
        for k in ('observations','distinct_candidate_loci','chosen_systems','structure_readback_pass','locally_supported_candidates','index_rows','indexed_systems'):
            c[k]+=x.get(k,0)
    counters=Counter(); sizes={}
    for folder in ('raw','discovery','systems','assemblies','entry_results','system_results','source_proofs'):
        fs=[p for p in (run/folder).rglob('*') if p.is_file()];sizes[folder]={'files':len(fs),'bytes':sum(p.stat().st_size for p in fs)}
    for s in raw:
        for a in s.get('attempts',[]):counters['observed_http_body_bytes']+=a['body_bytes']
    atomic_json(run/'summary.json',{'updated_utc':now(),'run_id':run.name,'frozen_entries':len(plan),'finished_entries':len(entries),'missing_entries':missing,
        'counts_by_stratum':sums,'entry_status':dict(Counter(x['status'] for x in entries)), 'entry_routes':dict(Counter(x.get('route') for x in entries)),
        'system_status':dict(Counter(x['status'] for x in flat)),'system_roles':dict(Counter(x['role_qc'] for x in flat)),
        'system_candidate_intake':dict(Counter(x['candidate_intake'] for x in flat)),
        'system_reasons':dict(Counter(x['reason'] for x in flat if x.get('reason'))),
        'assembly_status':dict(Counter(x['complete_assembly']['status'] for x in flat if x.get('complete_assembly'))),
        'raw_verified':sum(s['status']=='PASS' for s in raw),'network':counters,'sizes':sizes,'invocation':invocation})


def build(run,plan,selected,index,jobs,fetch_counts,start):
    rows,cache=prepare_index(run,plan,index)
    bootstrap=Path(__file__).resolve().parent
    signatures={f.name:sha256(f) for f in bootstrap.glob('holo_*.py')}
    fingerprint=__import__('hashlib').sha256(json.dumps(signatures,sort_keys=True).encode()).hexdigest()
    atomic_json(run/'build_inputs.json',{'code_files':signatures,'fingerprint':fingerprint,'algorithm':VERSION,'checked_utc':now()})
    counts=Counter()
    # Use the same worker budget for ordinary entries and proven large-entry slices.
    from holo_m2_parallel_tail import assemble, dispatch, promote
    slots=jobs; tail=run/'parallel_tail'; slices=4
    split_above=256 if slots>=slices else None
    targets={r['pdb_id']:r for r in plan if r['pdb_id'] in selected}
    complete={}; started={}; total=len(targets)
    with ProcessPoolExecutor(max_workers=slots) as pool:
        pending={pool.submit(one,run,r,rows[p],fingerprint,split_above=split_above):(p,None)
                 for p,r in targets.items()}
        while pending:
            finished,_=wait(pending,return_when=FIRST_COMPLETED)
            for task in finished:
                pdb,part=pending.pop(task)
                try:
                    p,status=task.result()
                    if p!=pdb:
                        raise ValueError('Structure worker returned a different PDB ID')
                except Exception as e:
                    status='worker_failure'; print(pdb,repr(e),flush=True)
                if part is None and status=='needs_slices':
                    started[pdb]=time.monotonic(); complete[pdb]=[]
                    for i in range(slices):
                        pending[pool.submit(dispatch,run,tail,targets[pdb],rows[pdb],fingerprint,i,slices)]=(pdb,i)
                    continue
                if part is not None:
                    complete[pdb].append(status in ('built','reused_verified'))
                    if len(complete[pdb])!=slices:
                        continue
                    try:
                        if not all(complete[pdb]):
                            raise ValueError('One or more candidate slices did not finish')
                        assemble(run,tail,pdb,slices,fingerprint)
                        promote(run,tail,pdb,round(time.monotonic()-started[pdb],3))
                        status='built_sliced'
                    except (AssertionError,OSError,ValueError) as exc:
                        atomic_json(run/'entry_results'/f'{pdb}.json',{
                            'pdb_id':pdb,'revision':targets[pdb]['revision'],'signature':fingerprint,
                            'status':'TECHNICAL_FAILURE','failure_stage':'candidate_slice_merge',
                            'error':str(exc) or type(exc).__name__})
                        status='slice_failure'
                counts[status]+=1; n=sum(counts.values())
                if n%10==0 or n==total or status not in ('built','reused_verified'):
                    print('structure',n,'/',total,pdb,dict(counts),round(time.monotonic()-start,2),flush=True)
    summarize(run,plan,{'finished_utc':now(),'fetch':fetch_counts,'index_snapshot':cache,'structure':counts,'elapsed_s':round(time.monotonic()-start,3)})

