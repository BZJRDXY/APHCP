"""Read-only APHCP chemical annotations: validate identity first, then classify once per key."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import time

import gemmi

from holo_candidate_batch import revision_from_cif, sha256
from holo_v1_export import rows, text
from holo_source_mapping import connection_edges

VERSION = 'aphcp-chemtype-v1'
REF = 'holo-m2-full-reconciliation-20261122a'
TYPE_ROOTS = {
    'LIPID': {'lipid'}, 'STEROID': {'steroid'}, 'STEROL': {'sterol'}, 'FATTY_ACID': {'fatty acid'},
    'NUCLEOTIDE': {'nucleotide', 'nucleoside phosphate'}, 'NUCLEOSIDE': {'nucleoside'},
    'CARBOHYDRATE': {'carbohydrate'}, 'POLYOL': {'polyol', 'alditol'},
    'MONOSACCHARIDE': {'monosaccharide'}, 'OLIGOSACCHARIDE': {'oligosaccharide'},
    'AMINO_ACID': {'amino acid'}, 'PEPTIDE': {'peptide'},
    '_KETONE_BODY_BRANCH': {'ketone body'},
}
METALS = frozenset('LI NA K RB CS FR BE MG CA SR BA RA AL GA IN TL SN PB BI SC TI V CR MN FE CO NI CU ZN Y ZR NB MO TC RU RH PD AG CD HF TA W RE OS IR PT AU HG LA CE PR ND PM SM EU GD TB DY HO ER TM YB LU AC TH PA U NP PU'.split())
COV = {'covale','covale_base','covale_phosphate','covale_sugar','disulf','modres'}


def load(path):
    return json.loads(Path(path).read_bytes())


def jsonl(path):
    with Path(path).open(encoding='utf8') as f:
        yield from map(json.loads, f)


def pack(x):
    return json.dumps(x, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def digest(x):
    return hashlib.sha256(pack(x).encode()).hexdigest()


def write(path, x):
    path=Path(path); tmp=path.with_suffix('.tmp')
    tmp.write_text(pack(x)+'\n',encoding='utf8'); tmp.replace(path)


def cif_blocks(path):
    with gzip.open(path,'rt',encoding='utf8') as f:
        lines=[]; key=None
        for line in f:
            if line.startswith('data_'):
                if key is not None:
                    yield key,lines
                lines=[]; key=line[5:].strip()
            if key is not None:
                lines.append(line)
        if key is not None:
            yield key,lines


def ccd_table(path,codes):
    result={}
    for k, lines in cif_blocks(path):
        if k not in codes:
            continue
        b=gemmi.cif.read_string(''.join(lines)).sole_block(); d=rows(b,'_chem_comp.')[0]
        aa={text(r['atom_id']):text(r['type_symbol']).upper() for r in rows(b,'_chem_comp_atom.')}
        bonds=rows(b,'_chem_comp_bond.'); keys=sorted({text(r['descriptor']) for r in rows(b,'_pdbx_chem_comp_descriptor.')
               if text(r['type']).lower()=='inchikey' and re.fullmatch(r'[A-Z]{14}-[A-Z]{10}-[A-Z]',text(r['descriptor']))})
        result[k]={'name':text(d.get('name')),'type':text(d.get('type')).lower(),'formula':text(d.get('formula')),
                   'release':text(d.get('pdbx_release_status')), 'modified':text(d.get('pdbx_modified_date')),
                   'formal_charge':text(d.get('pdbx_formal_charge')),'parent':text(d.get('mon_nstd_parent_comp_id')),
                   'elements':dict(Counter(aa.values())),'inchikeys':keys,
                   'carbon_bond':any(tuple(sorted([aa.get(text(x['atom_id_1']),''),aa.get(text(x['atom_id_2']),'')])) in {('C','C'),('C','H'),('C','D')} for x in bonds),
                   'bond_graph':{'atoms':aa,'bonds':[[text(v['atom_id_1']),text(v['atom_id_2']),text(v['value_order'])] for v in bonds]},
                   'source':'wwPDB CCD _chem_comp / _chem_comp_atom / _chem_comp_bond / _pdbx_chem_comp_descriptor'}
    return result


def bird_table(path):
    result={}
    for k, lines in cif_blocks(path):
        b=gemmi.cif.read_string(''.join(lines)).sole_block(); d=rows(b,'_pdbx_reference_molecule.')[0]
        if text(d.get('release_status'))!='REL':
            continue
        s=defaultdict(list)
        for row in rows(b,'_pdbx_reference_entity_poly_seq.'):
            s[text(row['ref_entity_id'])].append((int(text(row['num'])),text(row['mon_id'])))
        seqs=[ [q for _,q in sorted(v)] for v in s.values() ]; ccd=text(d.get('chem_comp_id'))
        links=rows(b,'_pdbx_reference_entity_poly_link.')
        result[k]={'name':text(d.get('name')),'type':text(d.get('type')),'representation':text(d.get('represent_as')),
                    'ccd':ccd,'seqs':seqs,'link_rows':len(links),'entity_parts':len(rows(b,'_pdbx_reference_entity_list.')),
                    'sequence_links':sorted(sorted([[int(text(x[f'entity_seq_num_{i}'])),text(x[f'comp_id_{i}']),text(x[f'atom_id_{i}'])]
                     for i in (1,2)])+[text(x.get('value_order')).lower()] for x in links if all(text(x.get(f'entity_seq_num_{i}')).isdigit() for i in (1,2)))}
    return result


def entity_context(block):
    def rt(cat, k):
        return {text(r[k]):{i:text(v) for i,v in r.items()} for r in rows(block,cat)}
    def sequences(category,number,code):
        out=defaultdict(lambda:defaultdict(set))
        for r in rows(block,category):
            out[text(r['entity_id'])][int(text(r[number]))].add(text(r[code]))
        return {k:{n:sorted(v) for n,v in sorted(a.items())} for k,a in out.items()}
    branchlinks=defaultdict(list)
    for r in rows(block,'_pdbx_entity_branch_link.'):
        def end(i):
            return [int(text(r[f'entity_branch_list_num_{i}'])),text(r[f'comp_id_{i}']),text(r[f'atom_id_{i}'])]
        branchlinks[text(r['entity_id'])].append(sorted([end(1),end(2)])+[text(r.get('value_order'))])
    bird=defaultdict(set)
    for r in rows(block,'_pdbx_molecule.'):
        if text(r.get('prd_id')).startswith('PRD_'):
            bird[text(r.get('asym_id'))].add(text(r['prd_id']))
    labelmap=rt('_struct_asym.','id')
    return {'entities':rt('_entity.','id'),'poly':rt('_entity_poly.','entity_id'), 'branch':rt('_pdbx_entity_branch.','entity_id'),
            'polyseq':sequences('_entity_poly_seq.','num','mon_id'),'branchseq':sequences('_pdbx_entity_branch_list.','num','comp_id'),
            'branchlinks':{k:sorted(v) for k,v in branchlinks.items()},'nonpoly':rt('_pdbx_entity_nonpoly.','entity_id'),
            'label_entity':{k:v['entity_id'] for k,v in labelmap.items()},'prd':dict(bird),
            'local_ccd':rt('_chem_comp.','id')}


def labels_for_polymer(poly):
    p=(poly or '').lower()
    if p=='polydeoxyribonucleotide':return ['NUCLEIC_ACID','DNA']
    if p=='polyribonucleotide':return ['NUCLEIC_ACID','RNA']
    if 'deoxyribonucleotide/polyribonucleotide hybrid' in p:return ['NUCLEIC_ACID','DNA_RNA_HYBRID']
    if p.startswith('polypeptide'):return ['POLYPEPTIDE']
    if p=='peptide nucleic acid':return ['PEPTIDE_NUCLEIC_ACID']
    return []


def connected(vertices,edges):
    if not vertices:return False
    visited={min(vertices)}
    while True:
        before=len(visited)
        for a,b in edges:
            if a in visited or b in visited: visited.update((a,b))
        if len(visited)==before:return visited==set(vertices)


def member_descriptor(m,c,model):
    eid=m['entity']; lab=m['label_asym']; kind=m['kind']; residues=m['residues']
    if c['label_entity'].get(lab)!=eid:
        return {'status':'AMBIGUOUS','reason':'source_label_entity_disagrees','kind':kind}
    base={'kind':kind,'prd':sorted(c['prd'].get(lab,[]))}
    if kind=='component' and len(residues)==1:
        code=residues[0]['ccd']
        declared=c['nonpoly'].get(eid,{}).get('comp_id')
        if declared and declared!=code:
            return base | {'status':'AMBIGUOUS','reason':'declared_component_conflicts_with_observation','ccd':code}
        embedded=c['local_ccd'].get(code,{})
        ids=c['nonpoly_atom_signatures'].get((lab,str(model),code),{})
        return base | {'status':'SINGLE_CCD','ccd':code,'local_type':embedded.get('type'),
                       'embedded_formula':embedded.get('formula'),'entity_comp_id':declared,'source_heavy_atoms':ids}
    if kind not in {'polymer','branched'}:
        return base | {'status':'AMBIGUOUS','reason':'standalone_identity_not_resolved'}
    seq=c['polyseq' if kind=='polymer' else 'branchseq'].get(eid,{})
    branch_links=c['branchlinks'].get(eid,[])
    ptype=c['poly'].get(eid,{}).get('type','') if kind=='polymer' else c['branch'].get(eid,{}).get('type','')
    obs=defaultdict(set)
    for residue in residues:
        pos=residue['label_seq']
        if kind=='branched' and pos is None:
            # In branching nomenclature atom_site label_seq_id can be null; scheme links auth to branch num.
            k=(lab,residue['auth_asym'],residue['auth_seq'],residue['ccd']); pos=c['branch_positions'].get(k)
        if pos is not None:
            obs[pos].add(residue['ccd'])
    ref_len=len(seq); missing=sorted(set(seq)-set(obs)); conflict=sorted(n for n,o in obs.items() if n not in seq or not o<=set(seq[n]))
    verified=bool(seq and not conflict and sum(map(len,obs.values()))==len(residues))
    return base | {'status':'ENTITY_SEQUENCE' if verified else 'ENTITY_SEQUENCE_REVIEW','polymer_type':ptype,'source_length':ref_len,
                  'sequence':list(seq.values()),'sequence_positions':list(seq),'observed_positions':sorted(obs),
                  'missing_declared_positions':missing,'discordant_positions':conflict,'branch_links':branch_links,
                  'branch_connected':kind!='branched' or connected(set(seq),[(r[0][0],r[1][0]) for r in branch_links]),
                  'observed_count':len(residues)}


def branch_scheme(c,block):
    schemes=defaultdict(lambda:[defaultdict(set),defaultdict(set)])
    at=defaultdict(set); links=defaultdict(list)
    for r in rows(block,'_pdbx_branch_scheme.'):
        lab=text(r.get('asym_id')); pdb=text(r.get('pdb_asym_id')); auth=text(r.get('auth_asym_id')) or pdb
        if text(r.get('num')).isdigit():
            n=int(text(r['num'])); mon=text(r.get('mon_id'))
            for scheme,prefix,chain in zip(schemes[lab],('pdb','auth'),(pdb,auth)):
                seq=text(r.get(f'{prefix}_seq_num')); comp=text(r.get(f'{prefix}_mon_id')) or mon
                if seq and chain:scheme[(lab,chain,seq,comp)].add(n)
    if schemes:
        for r in rows(block,'_atom_site.'):
            lab=text(r.get('label_asym_id'))
            if lab in schemes:
                k=(lab,text(r.get('auth_asym_id')),text(r.get('auth_seq_id'))+text(r.get('pdbx_PDB_ins_code')),text(r.get('label_comp_id')))
                at[lab].add(k)
        for edge in connection_edges(block):
            if edge['kind'] not in COV:continue
            a,b=edge['ends']; lab=a['label_asym']
            if lab in schemes and b['label_asym']==lab:
                links[lab].append([((lab,v['auth_asym'],(v['auth_seq'] or '')+(v['insertion'] or ''),v['ccd']),v['atom']) for v in [a,b]])
    positions={}
    for lab,options in schemes.items():
        raw=c['branchlinks'].get(c['label_entity'][lab],[])
        reference={tuple(sorted([tuple(a),tuple(b)])) for a,b,_ in raw}
        ranked=[]
        for option in options:
            unique={k:next(iter(v)) for k,v in option.items() if len(v)==1}
            matches=wrong=0
            for ends in links[lab]:
                if all(k in unique for k,_ in ends):
                    observed=tuple(sorted((unique[k],k[3],atom) for k,atom in ends))
                    matches+=observed in reference;wrong+=observed not in reference
            if reference and wrong:continue
            ranked.append(((matches,len(set(unique)&at[lab])),unique))
        if ranked:
            top=[ix for rank,ix in ranked if rank==max(rank for rank,_ in ranked)]
            for k in set.intersection(*[set(x) for x in top]):
                vals={ix[k] for ix in top}
                if len(vals)==1:positions[k]=next(iter(vals))
    c['branch_positions']=positions


def ranked_edge(a,b):
    return sorted([a,b],key=pack)


def partner_record(o,c):
    members=o['members']; desc=[member_descriptor(m,c,o['model']) for m in members]
    mi=defaultdict(list)
    for i,m in enumerate(members):
        for r in m['residues']:
            if r['label_seq'] is None and m['kind']=='branched':
                pos=c['branch_positions'].get((m['label_asym'],r['auth_asym'],r['auth_seq'],r['ccd']))
            else:pos=r['label_seq']
            mi[m['instance']].append((i,str(pos) if pos is not None else '',r['ccd'],r['auth_asym'],r['auth_seq']))
    edges=[]; extr=[]; cannot_map=[]
    for x in o.get('connections',[]):
        edge=x['source']; ends=edge['ends']; labs=x['instances']; kind=edge['kind']
        if kind not in COV:
            continue
        if len(ends)!=2 or len(labs)!=2:
            cannot_map.append(edge['id']);continue
        options=[]
        for inst,en in zip(labs,ends):
            ins=en.get('insertion') or ''; auth=(en.get('auth_seq') or '')+ins
            if any(m['instance']==inst and m['kind']=='component' for m in members):
                matches=[(i,1,en['atom']) for i,pos,ccd,ch,a in mi.get(inst,()) if ccd==en['ccd'] and
                         ((en.get('label_seq') and pos==str(en['label_seq'])) or (ch==en.get('auth_asym') and a==auth))]
            else:
                matches=[(i,int(pos),en['atom']) for i,pos,ccd,ch,a in mi.get(inst,()) if pos.lstrip('-').isdigit() and ccd==en['ccd'] and
                         ((en.get('label_seq') and pos==str(en['label_seq'])) or (ch==en.get('auth_asym') and a==auth))]
            options.append(matches)
        if all(len(v)==1 for v in options):
            a,b=options[0][0],options[1][0]
            edges.append({'id':edge['id'],'kind':kind,'ends':[list(a),list(b)]})
        elif any(len(v)>1 for v in options) or all(l in mi for l in labs):
            cannot_map.append(edge['id'])
        elif any(options):
            extra_i=0 if not options[0] else 1
            foreign=ends[extra_i]['label_asym']; foreign_eid=c['label_entity'].get(foreign)
            protein=c['poly'].get(foreign_eid,{}).get('type','').lower().startswith('polypeptide')
            is_receptor=protein and labs[extra_i] in o.get('protein_assembly_instances',[]) and (
                         labs[extra_i] in o.get('contact_instances_4A',[]) or labs[extra_i] in o.get('v1_binding_protein_instances',[]))
            # A short protein participant may itself be an ungrouped partner.
            is_receptor=is_receptor and (len(c['polyseq'].get(foreign_eid,{}))>35 or labs[extra_i] in o.get('v1_binding_protein_instances',[]))
            extr.append({'id':edge['id'],'foreign_instance':labs[extra_i],'foreign_label_asym':foreign,'foreign_entity':foreign_eid,
                         'outside_member_reason':'attached_to_protein_context' if is_receptor else 'additional_covalent_component_not_in_candidate_members'})
        elif any(l in mi for l in labs):
            cannot_map.append(edge['id'])
    grouping=o.get('member_grouping') or ''
    issues=[]
    if any(d['status']=='AMBIGUOUS' for d in desc):issues.append('member_identity_ambiguous')
    if any(d['status']=='ENTITY_SEQUENCE_REVIEW' for d in desc):issues.append('member_sequence_not_verified')
    if any(d.get('branch_connected') is False for d in desc):issues.append('branched_graph_not_verified')
    if any(x['outside_member_reason'].startswith('additional') for x in extr):issues.append('external_covalent_component_not_classified_as_partner')
    if cannot_map:issues.append('explicit_covalent_endpoints_not_mapped')
    if grouping.lower().startswith('review') or (grouping.startswith('explicit_') and not o.get('member_grouping_rows')):issues.append('complete_partner_group_not_verified')
    if len(members)>1 and not grouping.startswith('explicit_') and not connected(set(range(len(members))),
              [(x['ends'][0][0],x['ends'][1][0]) for x in edges]):issues.append('multi_member_connectivity_not_verified')
    if o.get('unresolved_connection_ids'):
        known_covalent=set(c['covalent_ids'])
        if known_covalent.intersection(o['unresolved_connection_ids']):issues.append('unresolved_partner_covalent_connection')
    # Identity chemistry excludes Assembly/model (instance details are carried separately).
    node=[{k:d[k] for k in ('kind','ccd','polymer_type','sequence','sequence_positions','branch_links') if k in d} for d in desc]
    internal=[(a,b,x['kind']) for x in edges for a,b in [x['ends']]]
    nontrivial=[]
    for a,b,k in internal:
        if a[0]==b[0] and desc[a[0]]['kind'] in {'polymer','branched'} and abs(a[1]-b[1])==1 and (({a[2],b[2]}=={'C','N'} and desc[a[0]]['polymer_type'].lower().startswith('polypeptide')) or
             ({a[2],b[2]} in ({"O3'",'P'},{'O3*','P'}) and 'nucleotide' in desc[a[0]]['polymer_type'].lower())):
            continue  # ordinary adjacent backbone link does not define a different compound key
        if a[0]==b[0] and desc[a[0]]['kind']=='branched' and desc[a[0]].get('branch_connected'):
            continue  # branch topology already authoritative in its branch table
        nontrivial.append((a,b,k))
    if len(desc)==1 and desc[0]['status']=='SINGLE_CCD' and not issues:
        chemical_id='ccd:'+desc[0]['ccd']
    else:
        order=sorted(range(len(node)),key=lambda i:(pack(node[i]),i)); map_idx={i:j for j,i in enumerate(order)}
        graph=sorted([[k,ranked_edge([map_idx[a[0]],*a[1:]],[map_idx[b[0]],*b[1:]])] for a,b,k in nontrivial],key=pack)
        # Avoid merging two different ambiguous groups or permutations of indistinguishable nodes with arbitrary positional assignment.
        twins=bool(graph and len({pack(v) for v in node})!=len(node))
        identity={'nodes':[node[i] for i in order],'covalent':graph,'association':'paired_nucleic_group' if grouping.startswith('explicit_') else None}
        if twins:issues.append('indistinguishable_graph_nodes_source_scoped')
        if issues:
            identity['source_scope']=[o['entry'],o['revision'],o['assembly'],o['model'],
                sorted((m['entity'],m['label_asym'],m['instance']) for m in members),sorted(set(issues))]
        chemical_id=('review:' if issues else 'partner:')+digest(identity)[:24]
    return {'identity':chemical_id,'issue':sorted(set(issues)), 'members':desc,'covalent_edges':edges,
            'external_covalent':extr,'unmapped_covalent':cannot_map,'grouping':grouping or None,
            'grouping_source_rows':o.get('member_grouping_rows',[]),
            'source_partner_scope_qc':o.get('partner_scope_qc'),
            'member_instances':[{'instance':m['instance'],'label_asym':m['label_asym'],'source_entity':m['entity'],'operator':m.get('operator'),
              'source_residues':[[r['auth_asym'],r['auth_seq'],r['label_seq'],r['ccd']] for r in m['residues']] if len(m['residues'])<25 else None}
              for m in members]}


def ccd_basic(d):
    if not d or d['release']!='REL' or not d['elements']:
        return [], 'ccd_definition_unavailable_or_not_released'
    t=d['type']; e=d['elements']; heavy=sum(n for a,n in e.items() if a not in {'H','D','T'})
    if heavy==1 and sum(e.values())==1:
        element=next(iter(e)); charge=d['formal_charge']
        if charge and re.fullmatch(r'[+-]?\d+',charge) and int(charge):
            return ['ION','METAL_ION' if element in METALS else 'MONATOMIC_NONMETAL_ION'],None
        return ['MONATOMIC_COMPONENT'],None
    if 'saccharide' in t:
        return ['CARBOHYDRATE','SACCHARIDE_COMPONENT'],None
    if 'peptide' in t:
        if 'linking' in t or 'terminus' in t:
            return ['SMALL_MOLECULE','AMINO_ACID_RELATED_COMPONENT'],None
        return ['SMALL_MOLECULE','PEPTIDE_LIKE_COMPONENT'],None
    if 'rna' in t or 'dna' in t:
        return ['SMALL_MOLECULE','NUCLEIC_ACID_RELATED_COMPONENT'],None
    if d['carbon_bond']:
        return ['SMALL_MOLECULE','ORGANIC_SMALL_MOLECULE'],None
    if 'C' not in e:
        if re.fullmatch(r'[+-]?\d+',d['formal_charge'] or '') and int(d['formal_charge']):
            return ['ION','POLYATOMIC_INORGANIC_ION'],None
        return ['SMALL_MOLECULE','INORGANIC_COMPONENT'],None
    return ['SMALL_MOLECULE'],None


def bird_matches(d,bird):
    found=[]
    for pid in d.get('prd', ()):
        b=bird.get(pid)
        if not b:continue
        if d.get('ccd'):
            if b['ccd']==d['ccd'] and b['representation'].lower().startswith('single'):
                found.append((pid,b,'exact deposited CCD'))
        elif d['status']=='ENTITY_SEQUENCE' and b['seqs'] and b['entity_parts']==1 and len(b['seqs'])==1 and d.get('sequence') and len(b['seqs'][0])==len(d['sequence']):
            if all(x in choices for x,choices in zip(b['seqs'][0],d['sequence'])):
                if d['kind']=='branched' and d['branch_connected'] and b['link_rows']==len(b['sequence_links']) and b['sequence_links']==d['branch_links']:
                    found.append((pid,b,'deposited PRD, exact ordered CCD sequence and branch atom graph'))
                elif d['kind']=='polymer':
                    found.append((pid,b,'deposited PRD and full ordered CCD sequence; PRD topology not independently verified'))
    return found


def base_class(p,ccd,bird):
    tags=set(); basis=[]; body=[]; found=[]; cycle_rows=[]
    for i,d in enumerate(p['members']):
        t=set(); branch_reason=None
        if d['status']=='SINGLE_CCD' and d.get('ccd'):
            cc=d['ccd']; cd=ccd.get(cc); labels,branch_reason=ccd_basic(cd)
            t.update(labels)
            if branch_reason is None:
                local=(d.get('local_type') or '').lower()
                check={'source':'CCD','component':cc,'reference_type':cd['type'],'embedded_type':local,'CCD_modified':cd['modified']}
                if local and local!=cd['type']:
                    observed=d['source_heavy_atoms']; expected=heavy_atom_signature(cd['bond_graph']['atoms'].items())
                    formula_alias_match=bool(d['embedded_formula']) and sorted(d['embedded_formula'].split())==sorted(cd['formula'].split())
                    consistent=(d['entity_comp_id']==cc and formula_alias_match and observed==expected and bool(expected['count']))
                    check['type_disagreement']={'resolved':consistent,'source_entity_comp_id':d['entity_comp_id'],
                        'source_formula':d['embedded_formula'],'CCD_formula':cd['formula'],'source_heavy_atom_signature':observed,'CCD_heavy_atom_signature':expected}
                    if not consistent:
                        branch_reason='embedded_vs_reference_chem_comp_type_conflict';t.clear()
                basis.append(check)
        elif d['status']=='ENTITY_SEQUENCE':
            if d['kind']=='polymer':
                t.update(labels_for_polymer(d['polymer_type']))
                if 'POLYPEPTIDE' in t:
                    short=2<=d['source_length']<=35
                    t.add('PEPTIDE' if short else 'LONG_POLYPEPTIDE')
                    if short:
                        start,end=d['sequence_positions'][0],d['sequence_positions'][-1]
                        for x in p['covalent_edges']:
                            a,b=x['ends']
                            if a[0]==b[0]==i and abs(a[1]-b[1])>1 and x['kind'] in {'covale','modres'}:
                                # A deposited non-adjacent covalent edge closes a path along the sequence backbone.
                                cycle_rows.append(x['id']); t.add('CYCLIC_PEPTIDE')
                                left,right=sorted([a,b],key=lambda v:v[1])
                                t.add('HEAD_TO_TAIL' if (left[1],right[1],left[2],right[2])==(start,end,'N','C')
                                      else 'SIDE_CHAIN_OR_MODIFIED_MACROCYCLE')
                            elif a[0]==b[0]==i and x['kind']=='disulf':
                                t.add('DISULFIDE_CROSSLINKED')
                basis.append({'source':'fixed_mmCIF','fields':['_entity_poly.type','_entity_poly_seq','_struct_conn'],'source_polymer_type':d['polymer_type']})
            elif d['kind']=='branched' and d['branch_connected'] and 'saccharide' in d['polymer_type'].lower():
                t.update(['CARBOHYDRATE','BRANCHED_ENTITY'])
                if d['source_length']>1:t.add('OLIGOSACCHARIDE' if d['source_length']<=10 else 'GLYCAN')
                basis.append({'source':'fixed_mmCIF','fields':['_pdbx_entity_branch','_pdbx_entity_branch_list','_pdbx_entity_branch_link'],
                              'source_branch_type':d['polymer_type'],'links':d['branch_links']})
        for pid,b,match in bird_matches(d,bird):
            found.append({'prd':pid,'match':match,'type':b['type'],'link_rows':b['link_rows']})
            kind=b['type'].lower(); mapped={'polypeptide':'POLYPEPTIDE','oligopeptide':'PEPTIDE','peptide-like':'PEPTIDE_LIKE_COMPONENT',
                                          'peptide':'PEPTIDE','oligosaccharide':'OLIGOSACCHARIDE'}.get(kind)
            # Use PRD type only when the claimed label applies to verified identity,
            # not as evidence of cyclization or of an unverified chemical graph.
            if mapped and branch_reason is None and 'not independently verified' not in match:
                t.add(mapped)
        if t and not branch_reason: body.append(sorted(t))
        else: body.append(['UNKNOWN'])
        if branch_reason:basis.append({'reason':branch_reason})
    tags.update(t for member in body for t in member if t!='UNKNOWN')
    if any('UNKNOWN' in m for m in body):p_issues=['one_or_more_members_not_classified']
    else:p_issues=[]
    if len(body)==2 and all('DNA' in m for m in body) and p['grouping_source_rows'] and str(p['grouping']).startswith('explicit_'):
        tags.add('DNA_DUPLEX')
    if len(body)==2 and all('RNA' in m for m in body) and p['grouping_source_rows'] and str(p['grouping']).startswith('explicit_'):
        tags.add('RNA_PAIRED_GROUP')
    if len(body)>1 and len({tuple(t) for t in body})>1:tags.add('MIXED_MEMBER_CHEMISTRY')
    reasons=sorted(set(p_issues+p['issue']+[v['reason'] for v in basis if 'reason' in v]))
    return {'labels':sorted(tags), 'member_labels':body,
            'reason':reasons or ([] if tags else ['insufficient_source_type_evidence']), 'evidence':basis, 'bird_verified':found,
            'cycle_connection_ids':sorted(set(cycle_rows))}


def chebi_ontology(path,keys):
    index=defaultdict(list); terms={}; header={}
    def save(t):
        if not t.get('id'):return
        terms[t['id']]={k:t[k] for k in ('name','parents','obsolete') if k in t}
        for key in t['keys'] & keys:
            index[key].append(t['id'])
    state={}; chapter='header'
    with gzip.open(path,'rt',encoding='utf8') as f:
        for l in f:
            if l.startswith('['):
                save(state); state={'parents':[],'keys':set()};chapter=l.strip(); continue
            if chapter=='header':
                if l.startswith(('data-version:', 'date:')): k,v=l.split(':',1);header[k]=v.strip()
                continue
            if chapter!='[Term]':continue
            if l.startswith(('id: ','name: ')):
                k,v=l.split(': ',1);state[k]=v.strip()
            elif l.startswith('is_a: '):state['parents'].append(l.split()[1])
            elif l.startswith('is_obsolete: true'):state['obsolete']=True
            elif 'inchi_key_string ' in l:
                m=re.search(r'"([A-Z]{14}-[A-Z]{10}-[A-Z])"',l)
                if m:state['keys'].add(m[1])
    save(state)
    named={name:label for label,names in TYPE_ROOTS.items() for name in names}
    roots={k:named[v['name']] for k,v in terms.items() if v.get('name') in named and not v.get('obsolete')}
    if not {'LIPID','NUCLEOTIDE','STEROID','CARBOHYDRATE'}<=set(roots.values()):
        raise ValueError(f'Missing ontology roots: {roots}')
    @lru_cache(maxsize=None)
    def tags(term):
        item=terms.get(term,{})
        if item.get('obsolete'):return frozenset()
        return frozenset(([roots[term]] if term in roots else [])+[x for par in item.get('parents',[]) for x in tags(par)])
    idx={}
    for key,ids in sorted(index.items()):
        active=sorted(set(x for x in ids if not terms[x].get('obsolete')))
        if active:
            sets=[tags(t) for t in active]; common=set.intersection(*(set(x) for x in sets))
            domain_exclusions={'LIPID':'ChEBI ketone-body ancestry alone is not evidence of structural lipid class'} if (
                'LIPID' in common and any('_KETONE_BODY_BRANCH' in a for a in sets)) else {}
            domain_labels=sorted(v for v in common if not v.startswith('_') and v not in domain_exclusions)
            idx[key]={'terms':[{'id':t,'name':terms[t].get('name'),'labels':sorted(v)} for t,v in zip(active,sets)],
                      'common_labels':domain_labels,'excluded_in_scope':domain_exclusions,'conflicting_labels':len(set(sets))>1}
    return idx, {'header':header,'roots':{k:{'label':v,'name':terms[k]['name']} for k,v in sorted(roots.items())},'terms':len(terms)}


def ch_add(base, cd, ontology):
    keys=cd.get('inchikeys',[]) if cd else []
    if len(keys)!=1:
        return {'status':'CCD_FULL_INCHIKEY_NOT_UNIQUE_OR_MISSING','added':[]}
    key=keys[0]; data=ontology.get(key)
    if not data:return {'status':'NO_EXACT_CHEBI_IDENTITY','full_inchikey':key,'added':[]}
    allow=set(data['common_labels'])-set(base['labels'])
    if not base['labels']:allow=set() # unresolved reference type stays review, not cross-source override
    return {'status':'EXACT_FULL_INCHIKEY','full_inchikey':key, **data,'added':sorted(allow)}


def heavy_atom_signature(items):
    aa=sorted(set((a,e) for a,e in items if e not in {'D','T','H'}))
    return {'count':len(aa),'fingerprint':digest(aa)[:32]}


def read_context(source,p,revision,expected):
    block, actual=revision_from_cif(Path(source),p)
    if actual!=revision or sha256(Path(source))!=expected:
        raise ValueError(f'Frozen source changed for {p} {revision}')
    c=entity_context(block);branch_scheme(c,block)
    positions=defaultdict(set); ligand_labels={k for k,e in c['label_entity'].items() if e in c['nonpoly']}
    if ligand_labels:
        for fields in block.find('_atom_site.',['label_asym_id','pdbx_PDB_model_num','label_comp_id','label_atom_id','type_symbol']):
            if fields[0] in ligand_labels and fields[4] not in {'H','D','T'}:
                positions[tuple(fields[i] for i in range(3))].add((text(fields[3]),text(fields[4]).upper()))
    c['nonpoly_atom_signatures']={key:heavy_atom_signature(vals) for key,vals in positions.items()}
    c['covalent_ids']={edge['id'] for edge in connection_edges(block) if edge['kind'] in COV}
    return block,c


def entry_worker(args):
    p, records=args; first=records[0]
    _,c=read_context(first['source_path'],p,first['revision'],first['source_sha256'])
    if any((v['source_path'],v['revision'],v['source_sha256'])!=(first['source_path'],first['revision'],first['source_sha256']) for v in records):
        raise ValueError(f'Conflicting frozen sources for {p}')
    discoveries={}
    for s in records:
        path=s['discovery_path']; k=s['candidate_key']
        if path not in discoveries:
            ledger=load(path)
            if ledger['pdb_id']!=p or ledger['revision']!=first['revision']:
                raise ValueError(f'Wrong discovery ledger: {path}')
            discoveries[path]={o['candidate_key']:o for o in ledger['native_observations']}
        o=discoveries[path][k]; actual={m['entity'] for m in o['members']}
        if s.get('entity_ids') and set(s['entity_ids'])!=actual:
            raise ValueError(f'Observation does not match baseline attribution: {p} {k}')
        s['partner']=partner_record(o,c)
        for key in ('assembly','model'):
            if key in s and s[key]!=o[key]:raise ValueError(f'{p} {k} wrong {key}')
            s[key]=o[key]
    return p,records


def layer_rows(r,base):
    hist=load(r/'historical_deltas.json')['deltas']; validate=load(r/'validation.json')['structure_retest']['systems']
    existing={(x['pdb_id'],x['candidate_key']):x for x in base}
    prov={x['pdb_id']:x for x in jsonl(r/'entry_attribution.jsonl')}
    made=[]
    for h in hist+validate:
        p=h['pdb_id']; k=h['candidate_key']; p0=prov[p]; old=existing.get((p,k))
        assert h['revision']==p0['revision'] and h['source_sha256']==p0['source_sha256']
        rr=Path(h['result_root']) if 'result_root' in h else r/'new_structure'/p
        layer=h.get('historical_layer','new_structure_reconciliation')
        made.append({**(old or {z:p0[z] for z in ('batch','pdb_id','stratum','revision','source_path','source_sha256')}),
                     'evidence_layer':layer,'candidate_key':k, 'candidate_is_increment':old is None,
                     'prior_status':'EXPORTED_READBACK','role_after':h.get('role_after',h.get('role')),
                     'intake_after':h.get('candidate_intake_after',h.get('intake')),
                     'system_path':h['system_path'],'system_sha256':h['system_sha256'],'mapping_path':h['mapping_path'],
                     'discovery_path':str(rr/'discovery'/f'{p}.json')})
    assert len(hist)==14 and len(validate)==1 and sum(x['candidate_is_increment'] for x in made)==9
    return made,prov


def from_sources(baseline,overlays,jobs):
    rows_by=defaultdict(list)
    for s in baseline+overlays:rows_by[s['pdb_id']].append(s)
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        for p,records in pool.map(entry_worker, sorted(rows_by.items()),chunksize=2):
            yield p,records
from contextlib import contextmanager

@contextmanager
def jsonl_output(path):
    path=Path(path); tmp=path.with_suffix('.tmp')
    with tmp.open('w',encoding='utf8',newline='\n') as f:
        yield lambda row:f.write(pack(row)+'\n')
    tmp.replace(path)


def collect_evidence(run,ref,jobs):
    start=time.monotonic(); baseline=list(jsonl(ref/'system_attribution.jsonl')); overlay,entries=layer_rows(ref,baseline)
    expected={(s['batch'],s['pdb_id'],s['candidate_key']) for s in baseline}
    assert len(baseline)==len(expected)==46855 and len(entries)==2200
    assert sum(1 for x in entries.values() if x['batch']=='200')==200
    incoming=Counter((s['batch'],s['pdb_id'],s['candidate_key']) for s in overlay)
    assert len(incoming)==len(overlay)==15 and len(expected|set(incoming))==46864
    assert all(s['pdb_id'] in entries for s in baseline)
    ccds={k for x in baseline for k in x['ccds']}; assert len(ccds)==1028,len(ccds)
    shards=run/'identity-v4'; shards.mkdir(exist_ok=True)
    src_names=('entry_attribution.jsonl','system_attribution.jsonl','historical_deltas.json','validation.json')
    fingerprints={name:sha256(ref/name) for name in src_names}
    dependency_names=('holo_m2_chemtype.py','holo_candidate_batch.py','holo_v1_export.py','holo_source_mapping.py')
    extraction_rule='source-identity-v4' # update if extraction changes; classifier changes do not invalidate this checkpoint
    seal={'inputs':fingerprints,'identity_rule':extraction_rule}
    seal_path=shards/'seal.json'
    if seal_path.exists() and load(seal_path)!=seal:raise ValueError('Identity checkpoint input changed; use new checkpoint directory')
    if not seal_path.exists():write(seal_path,seal)
    expected_by=defaultdict(Counter)
    for x in baseline+overlay:expected_by[x['pdb_id']][(x['candidate_key'],x.get('candidate_is_increment'),x['evidence_layer'])]+=1
    reused=set(); todo={p for p in expected_by if not (shards/f'{p}.jsonl').is_file()}
    for p in sorted(set(expected_by)-todo):
        found=Counter((x['candidate_key'],x.get('candidate_is_increment'),x['evidence_layer']) for x in jsonl(shards/f'{p}.jsonl'))
        if found!=expected_by[p]:raise ValueError(f'Incomplete shard for {p}')
        reused.add(p)
    fields=('batch','pdb_id','candidate_key','stratum','route','revision','source_path','source_sha256','discovery_path','evidence_layer',
            'candidate_is_increment','assembly','model','model_scope','ccds','role_after','prior_role','intake_after','prior_status',
            'coordinate_identity','experimental_fidelity','topology','partner_scope_qc','prior_intake','system_path','system_sha256','mapping_path','partner')
    done=0
    for p,records in from_sources([x for x in baseline if x['pdb_id'] in todo],[x for x in overlay if x['pdb_id'] in todo],jobs):
        with jsonl_output(shards/f'{p}.jsonl') as emit:
            for x in records:emit({k:x[k] for k in fields if k in x})
        done+=1
        if done%100==0:print(f'identity: {done}/{len(todo)} new entry shards',flush=True)
    assert done==len(todo)
    receipt={'identity':seal,'baseline_ccds':sorted(ccds),'baseline':len(baseline),'overlays':len(overlay),'unique_after_overlay':len(expected|set(incoming)),
             'entries':len(entries),'candidate_entries':len(expected_by),'source_shards_reused':len(reused),'source_shards_built':done,
             'elapsed_s':round(time.monotonic()-start,2), 'code':{n:sha256(Path(__file__).parent/n) for n in dependency_names}}
    write(run/'identity_receipt.json',receipt)
    return receipt,entries


# Reference additions cannot overwrite source identity or the primary type.
GENERIC={'SMALL_MOLECULE','ORGANIC_SMALL_MOLECULE','INORGANIC_COMPONENT','MIXED_MEMBER_CHEMISTRY'}

def primary_for(p, base):
    if p['issue'] or any(x==['UNKNOWN'] for x in base['member_labels']):return 'UNKNOWN'
    categories=set()
    for labs in base['member_labels']:
        for tags,kind in (({'ION','MONATOMIC_COMPONENT'},'ION_OR_ATOMIC_COMPONENT'),({'NUCLEIC_ACID','PEPTIDE_NUCLEIC_ACID'},'NUCLEIC_ACID'),
                          ({'PEPTIDE'},'PEPTIDE'),({'POLYPEPTIDE'},'PROTEIN_PARTNER'),({'CARBOHYDRATE'},'CARBOHYDRATE'),
                          ({'SMALL_MOLECULE'},'SMALL_MOLECULE')):
            if tags.intersection(labs):
                categories.add(kind);break
        else:return 'UNKNOWN'
    if len(categories)==1:return next(iter(categories))
    return 'MIXED_COVALENT_PARTNER' if p['covalent_edges'] and len(categories)>1 else 'UNKNOWN'


def graph_polyol(d):
    """Narrow RDKit trial: saturated acyclic hydrocarbon skeleton with >=2 free OH, no other heterochemistry."""
    if not d or not {'C','O'}<=d['elements'].keys() or set(d['elements'])-{'C','O','H'} or d['elements']['O']<2 or d['formal_charge'] not in {'0','+0','-0'}:
        return False
    from rdkit import Chem
    src=d['bond_graph']; g=Chem.RWMol(); idx={k:g.AddAtom(Chem.Atom(v)) for k,v in src['atoms'].items() if v!='H'}
    for a,b,order in src['bonds']:
        if a not in idx or b not in idx:continue
        if order.lower()!='sing':return False
        g.AddBond(idx[a],idx[b],Chem.BondType.SINGLE)
    mol=g.GetMol()
    try:Chem.SanitizeMol(mol)
    except ValueError:return False
    if len(Chem.GetMolFrags(mol))!=1 or mol.GetRingInfo().NumRings():return False
    return all(at.GetDegree()==1 and at.GetTotalNumHs()==1 and at.GetNeighbors()[0].GetAtomicNum()==6
               for at in mol.GetAtoms() if at.GetAtomicNum()==8)


def annotation(p,base, add, extras):
    primary=primary_for(p,base); labels=set(base['labels'])|set(add)|set(extras)
    sub=sorted(labels-{'SMALL_MOLECULE','POLYPEPTIDE','ION','NUCLEIC_ACID'}-({'CARBOHYDRATE'} if primary=='CARBOHYDRATE' else set()))
    reasons=sorted(set(base['reason']))
    missing=[{'member_index':i,'deposited_count':v['source_length'],'observed_count':v['observed_count'],
              'missing_positions':v['missing_declared_positions']} for i,v in enumerate(p['members']) if v.get('missing_declared_positions')]
    subtype='KNOWN'
    if primary=='UNKNOWN':subtype='WHOLE_PARTNER_UNRESOLVED'
    elif primary=='SMALL_MOLECULE' and not (labels-GENERIC):
        subtype='NOT_RESOLVED_AT_TRIAL_SUBTYPE_GRANULARITY'; reasons.append('no_validated_fine_subclass')
    if primary=='UNKNOWN' and not reasons:reasons.append('complete_partner_primary_not_justified')
    if missing:reasons.append('deposited_entity_sequence_partly_observed')
    stat='UNKNOWN' if not labels else 'PARTIAL' if primary=='UNKNOWN' or subtype!='KNOWN' or missing else 'LABELED'
    return {'primary_class':primary,'chemical_subclass':sub if primary!='UNKNOWN' else [],'annotation_status':stat,
            'known_member_labels':base['member_labels'],'whole_partner_unresolved':primary=='UNKNOWN',
            'subtype_status':subtype,'partial_observation':missing,'status_reasons':reasons}


def classifier(ccd,bird,ontology):
    @lru_cache(maxsize=None)
    def rd_trial(code):return graph_polyol(ccd.get(code))
    @lru_cache(maxsize=None)
    def exact_chebi(code):
        d=ccd.get(code);return ch_add({'labels':ccd_basic(d)[0]},d,ontology)
    @lru_cache(maxsize=32768)
    def basic_signature(sig):
        q=json.loads(sig); a=base_class(q,ccd,{})
        return a, base_class(q,ccd,bird) if any(d.get('prd') for d in q['members']) else a
    def classify(p):
        sig={k:p[k] for k in ('members','covalent_edges','grouping','grouping_source_rows','issue')}
        sig['grouping_source_rows']=bool(sig['grouping_source_rows'])
        sig['covalent_edges']=[{**v,'id':i} for i,v in enumerate(sig['covalent_edges'])]
        a,b=basic_signature(pack(sig)); ref={'status':'NOT_A_VERIFIED_SINGLE_CCD','added':[]}; rd=[]; rd_status='INELIGIBLE'
        if p['identity'].startswith('ccd:') and b['labels']:
            code=p['identity'][4:]; exact=exact_chebi(code); ref={**exact,'added':sorted(set(exact['added'])-set(b['labels']))}
            if 'POLYOL' not in (set(b['labels'])|set(ref['added'])):
                rd_status='NO_MATCH'
                if rd_trial(code):rd=['POLYOL'];rd_status='MATCH_SATURATED_ACYCLIC_ALKANEPOLYOL'
        bird_add=sorted(set(b['labels'])-set(a['labels']))
        stages={'mmcif_ccd':annotation(p,a,[],[]),'plus_bird':annotation(p,b,[],[]),
                'plus_chebi':annotation(p,b,ref['added'],[]),'plus_rdkit':annotation(p,b,ref['added'],rd)}
        return stages,{'ccd_mmcif_evidence':b['evidence'],'bird_verified':b['bird_verified'],'bird_added_labels':bird_add,
                      'chebi':ref,'rdkit':{'status':rd_status,'added':rd},
                      'cycle_connection_ids':[p['covalent_edges'][i]['id'] for i in b['cycle_connection_ids']]}
    return classify


def single_ccd(code):
    return {'identity':'ccd:'+code,'members':[{'status':'SINGLE_CCD','kind':'component','ccd':code,'prd':[]}],
            'covalent_edges':[],'issue':[],'grouping':None,'grouping_source_rows':[]}


def source_catalog(run):
    formats={
        'components.cif.gz':('CCD','https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz',121031032,'Sat, 14 Nov 2026 17:40:24 GMT'),
        'prd-all.cif.gz':('BIRD','https://files.wwpdb.org/pub/pdb/data/bird/prd/prd-all.cif.gz',288703,'Sat, 14 Nov 2026 17:40:02 GMT'),
        'chebi-257.obo.gz':('ChEBI 257','https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz',47744176,'Thu, 05 Nov 2026 18:11:13 GMT')}
    out={}
    sources=run/'sources'
    for name,(ref,url,n,modified) in formats.items():
        f=sources/name
        if f.stat().st_size!=n:raise ValueError(f'Reference snapshot length changed: {f}')
        with gzip.open(f,'rb') as inp:
            while inp.read(1024*1024):pass
        out[ref]={'path':str(f),'bytes':n,'sha256':sha256(f),'gzip_verified':True,'source':url,'header_last_modified_utc':modified}
    return out


def historical_control(run,pdb):
    if pdb=='1HEW':
        root=run.parent/'holo-pilot-20261120c'; raw=root/'raw/1hew-v2.1.cif.gz'; old=root/'reconstructed/1hew/system.cif'
        proof=sha256(raw); b,c=read_context(raw,pdb,'2.1',proof); scheme=rows(b,'_pdbx_branch_scheme.')
        assert len(scheme)==3 and all(text(x['pdb_mon_id'])=='NAG' for x in scheme)
        # Observation from the deposited source coordinate rows, not from expected type.
        lab=text(scheme[0]['asym_id']); author=text(scheme[0]['pdb_asym_id']); entity=text(scheme[0]['entity_id'])
        seen={(text(x['auth_seq_id']),text(x['label_comp_id'])) for x in rows(b,'_atom_site.') if text(x['label_asym_id'])==lab}
        member={'kind':'branched','label_asym':lab,'entity':entity,'instance':f'1.{lab}','operator':'1',
                'residues':[{'auth_asym':author,'auth_seq':i,'ccd':cc,'label_seq':None} for i,cc in sorted(seen)]}
        con=[{'source':x,'instances':[member['instance']]*2} for x in connection_edges(b) if all(v['label_asym']==lab for v in x['ends'])]
        ds=[x for x in connection_edges(gemmi.cif.read_file(str(old)).sole_block()) if all(v['ccd']=='NAG' for v in x['ends'])]
        assert len(ds)==len(con)==2 and {tuple(sorted(v['atom'] for v in x['ends'])) for x in ds}=={('C1','O4')}
        cand={'entry':pdb,'revision':'2.1','model':1,'assembly':'1','members':[member],'connections':con}
        return partner_record(cand,c),{'fixed_source':str(raw),'sha256':proof,'readback':str(old),'branch_positions':sorted(seen),'deposited_and_output_link_count':len(ds)}
    a=run.parent/'holo-v1-acceptance-20261121a'; q=load(a/'qc'/f'{pdb}.json'); src=q['source']
    b,c=read_context(src['path'],pdb,q['revision'],src['sha256']); art=q['artifacts']['system']
    assert sha256(Path(art['path']))==art['sha256']; out=gemmi.cif.read_file(art['path']).sole_block()
    members=[]; con={}; wanted=set(q['partner_instances'])
    for o in q['observed_roles']:
        if {m['instance'] for m in o['members']}<=wanted:
            members.extend(o['members']); con.update({x['source']['id']:x for x in o['connections']})
    assert len(members)==len(wanted) and {m['instance'] for m in members}==wanted
    base={'entry':pdb,'revision':q['revision'],'model':q['representative_source_model'],'assembly':q['assembly'],
          'members':members,'connections':list(con.values()),'protein_assembly_instances':q['binding_protein_instances'],
          'v1_binding_protein_instances':q['binding_protein_instances']}
    proof={'fixed_source':src['path'],'sha256':src['sha256'],'readback':art['path'],'readback_sha256':art['sha256']}
    if 'partner_group' in q['checks']:
        cat=q['checks']['partner_group']['source_and_disk_category']; original=rows(b,cat); actual=rows(out,cat)
        assert len(original)==len(actual)==len(q['checks']['partner_group']['pairs']) and len(wanted)==2
        base.update(member_grouping='explicit_same_operation_local_base_pairs',member_grouping_rows=len(actual))
        proof.update(pair_category=cat,source_and_readback_pair_count=len(actual))
    proof['output_connections']=[x['id'] for x in connection_edges(out)]
    assert {v['source']['id'] for v in base['connections'] if v['source']['kind'] in COV}<=set(proof['output_connections'])
    return partner_record(base,c),proof


def controls(run,ccd,classify):
    # Targets and rejected labels below were specified by manual chemistry/identity inspection.
    # Development samples have no vote in the independent score.
    specs=[
      ('development','ATP','nucleotide','SMALL_MOLECULE',['NUCLEOTIDE'],['NUCLEIC_ACID'],'The CCD graph has adenosine and a triphosphate.'),
      ('development','GOL','polyol','SMALL_MOLECULE',['POLYOL'],['SACCHARIDE_COMPONENT'],'Three aliphatic hydroxyls; no saccharide ring.'),
      ('development','1HEW','oligosaccharide','CARBOHYDRATE',['OLIGOSACCHARIDE'],[],'Three NAG on one connected branch, not three independent ligands.'),
      ('development','7OUN','macrocyclic_peptide','PEPTIDE',['CYCLIC_PEPTIDE','SIDE_CHAIN_OR_MODIFIED_MACROCYCLE'],['HEAD_TO_TAIL'],'14-residue entity; ALA1 CA to CYS13 SG plus two disulfides.'),
      ('development','6VTX','dna_duplex','NUCLEIC_ACID',['DNA','DNA_DUPLEX'],['RNA'],'Two separate 12-nt DNA chains with 12 source/writeback pair rows.'),
      ('holdout','GTP','nucleotide','SMALL_MOLECULE',['NUCLEOTIDE'],['NUCLEIC_ACID'],'Guanosine connected to triphosphate.'),
      ('holdout','FAD','nucleotide','SMALL_MOLECULE',['NUCLEOTIDE'],['NUCLEIC_ACID'],'Diphosphate-linked dinucleotide despite its possible cofactor role.'),
      ('holdout','CLR','steroid','SMALL_MOLECULE',['LIPID','STEROID','STEROL'],['FATTY_ACID'],'Cholesterol steroid skeleton and one hydroxyl.'),
      ('holdout','EDO','polyol','SMALL_MOLECULE',['POLYOL'],['NUCLEOTIDE'],'Ethylene glycol, a saturated diol.'),
      ('holdout','PGO','polyol','SMALL_MOLECULE',['POLYOL'],['NUCLEOTIDE'],'S-1,2-propanediol, a saturated diol.'),
      ('holdout','MPD','polyol','SMALL_MOLECULE',['POLYOL'],['SACCHARIDE_COMPONENT'],'2-methyl-2,4-pentanediol; two free hydroxyls.'),
      ('holdout','IMD','ordinary_small_molecule','SMALL_MOLECULE',[],['POLYOL','NUCLEOTIDE','LIPID'],'Imidazole heteroaromatic ring lacks sugar/phosphate.'),
      ('holdout','ACN','ordinary_small_molecule','SMALL_MOLECULE',[],['POLYOL','NUCLEOTIDE','LIPID'],'Acetone contains a ketone, no hydroxyls.'),
      ('holdout','ZN','metal_ion','ION_OR_ATOMIC_COMPONENT',['METAL_ION'],['POLYOL'],'Formal ion charge from CCD; role unspecified.'),
      ('holdout','SO4','polyatomic_ion','ION_OR_ATOMIC_COMPONENT',['POLYATOMIC_INORGANIC_ION'],['METAL_ION'],'Sulfate tetrahedral anion; not a metal or a role label.'),
      ('holdout','BGC','saccharide','CARBOHYDRATE',['SACCHARIDE_COMPONENT'],['OLIGOSACCHARIDE'],'One beta-D-glucopyranose ring from CCD.'),
      ('holdout','5R42','short_peptide','PEPTIDE',['PEPTIDE'],['CYCLIC_PEPTIDE'],'Four-residue peptide; covalent protein context and role uncertainty do not change its chemistry.'),
      ('holdout','1CKT/ff3dcb5ff66c752d35c8ea9b','dna_duplex','NUCLEIC_ACID',['DNA','DNA_DUPLEX'],['RNA'],'Two full 16-residue DNA sequences; source explicitly pairs them.'),
      ('holdout','11MH/b93ea66e7ba7014629a63049','oligosaccharide','CARBOHYDRATE',['OLIGOSACCHARIDE'],[],'Two linked NAG, not the same identity as the 1HEW tri-saccharide.'),
      ('holdout','6KZU/3ddc6b97612c2751ce0dc48d','macrocyclic_peptide','PEPTIDE',['CYCLIC_PEPTIDE'],['HEAD_TO_TAIL'],'Side-chain staple between modified CCD positions 1 and 5.'),
      ('holdout','1A9N/8e1739bee4c7d5582f9e085d','scope_ambiguity','UNKNOWN',[],['RNA_PAIRED_GROUP'],'An observed RNA chain, but complete pairing/partner boundary unresolved.'),
      ('holdout','1GAG/214c674064a88d2958a152c1','scope_ambiguity','UNKNOWN',[],['PEPTIDE'],'Peptide linked to CCD 112; this candidate excludes that component.'),
      ('holdout','2E4H/5ab8be0bdbe5c0a5eaaa2d63','observed_protein_fragment','PROTEIN_PARTNER',['LONG_POLYPEPTIDE'],['PEPTIDE'],'Deposited length 36 vs only 8 observed; retain both extents and distinguish observed fragment.'),
      ('fresh_holdout','PLM','fatty_acid','SMALL_MOLECULE',['LIPID','FATTY_ACID'],['POLYOL'],'Palmitic acid, straight-chain saturated carboxylic acid.'),
      ('fresh_holdout','OLC','lipid_ester','SMALL_MOLECULE',['LIPID'],['POLYOL'],'Glyceryl oleate; hydroxyls do not make this a simple unfunctionalized alkanepolyol.'),
      ('fresh_holdout','CIT','hydroxy_carboxylic_acid','SMALL_MOLECULE',[],['NUCLEOTIDE','POLYOL'],'Citric acid; only one free aliphatic OH, not a polyol.'),
      ('fresh_holdout','DMS','ordinary_small_molecule','SMALL_MOLECULE',[],['LIPID','POLYOL','NUCLEOTIDE'],'Dimethyl sulfoxide does not represent any of the three targeted moieties.'),
      ('fresh_holdout','EPE','ordinary_small_molecule','SMALL_MOLECULE',[],['LIPID','POLYOL','NUCLEOTIDE'],'HEPES sulfonate with just one alcohol; no polyol.'),
    ]
    evaluated=[]; ids=set()
    for subset,sample,family,primary,required,excluded,note in specs:
        if '/' in sample:
            p,k=sample.split('/'); x=next(x for x in jsonl(run/'identity-v4'/f'{p}.jsonl') if x['candidate_key']==k)
            partner=x['partner']; b,_=read_context(x['source_path'],p,x['revision'],x['source_sha256'])
            proof={'fixed_source':x['source_path'],'sha256':x['source_sha256'],'candidate':k,'discovery':x['discovery_path'],
                   'members':[(m['kind'],m.get('source_length'),m.get('observed_count'),m.get('polymer_type')) for m in partner['members']],
                   'source_covalent_ids':[v['id'] for v in connection_edges(b) if v['kind'] in COV]}
            if x['system_path']:
                assert sha256(Path(x['system_path']))==x['system_sha256']
                out=gemmi.cif.read_file(x['system_path']).sole_block()
                proof['readback']=x['system_path'];proof['readback_atom_rows']=len(out.find_values('_atom_site.id'))
        elif sample in {'1HEW','7OUN','6VTX','5R42'}:
            partner,proof=historical_control(run,sample)
        else:
            partner=single_ccd(sample); d=ccd[sample]; proof={k:d[k] for k in ('name','type','formula','formal_charge','inchikeys')}
        if partner['identity'] in ids:raise ValueError(f'Same chemical control counted twice: {sample}')
        ids.add(partner['identity']); stages,evidence=classify(partner); answers={}
        for stage,v in stages.items():
            missed=sorted(set(required)-set(v['chemical_subclass'])); wrong=sorted(set(excluded)&set(v['chemical_subclass']))
            if v['primary_class']!=primary:wrong.append(f'primary_class:{v["primary_class"]}')
            if sample=='2E4H/5ab8be0bdbe5c0a5eaaa2d63' and not v['partial_observation']:missed.append('observed sequence extent')
            answers[stage]={**v,'missing_expected':missed,'wrong_checked':wrong,
                            'evaluation':'MISCLASSIFIED' if wrong else 'MISSED' if missed else 'CORRECT'}
        evaluated.append({'split':subset,'sample':sample,'identity':partner['identity'],'family':family,'expected_primary':primary,
                          'checked_required':required,'checked_forbidden':excluded,'manual_basis':note,'proof':proof,'stages':answers,
                          'source_additions':evidence})
    def measures_for(split):
        return score_controls([x for x in evaluated if x['split']==split])
    return {'cases':evaluated,'holdout':measures_for('holdout'),'fresh_holdout':measures_for('fresh_holdout'),
            'manual_case_count':len(evaluated),'independent_count':sum(x['split']!='development' for x in evaluated)}


def score_controls(holdout):
    measures={}; prior=None
    for stag in holdout[0]['stages']:
        rows_=[v['stages'][stag] for v in holdout]; statuses=[v['evaluation'] for v in rows_]; by_family=defaultdict(Counter)
        for x,v in zip(holdout,rows_):by_family[x['family']][v['evaluation']]+=1
        by_stage={k:statuses.count(k) for k in ('CORRECT','MISCLASSIFIED','MISSED')}
        new=[]; regress=[]; added_tags=0
        if prior:
            for x,v in zip(holdout,rows_):
                q=x['stages'][prior]; added_tags+=len((set(v['chemical_subclass'])-set(q['chemical_subclass']))&set(x['checked_required']))
                if q['evaluation']!='CORRECT' and v['evaluation']=='CORRECT':new.append(x['sample'])
                if q['evaluation']=='CORRECT' and v['evaluation']!='CORRECT':regress.append(x['sample'])
        measures[stag]={**by_stage,'by_class':by_family,'newly_correct_cases':new,'lost_correct_cases':regress,'new_correct_checked_tags':added_tags,
                        'whole_partner_unknown':sum(v['primary_class']=='UNKNOWN' for v in rows_),
                        'fully_unknown':sum(v['annotation_status']=='UNKNOWN' for v in rows_)}
        prior=stag
    return measures


def annotation_count():
    return {'status':Counter(),'primary':Counter(),'subclass':Counter(),'reasons':Counter(),
            'whole_partner_unresolved':0,'unknown_structure_readback_retained':0}


def count_annotation(result, anno, row=None):
    result['status'][anno['annotation_status']]+=1;result['primary'][anno['primary_class']]+=1
    result['subclass'].update(anno['chemical_subclass']);result['reasons'].update(anno['status_reasons'])
    result['whole_partner_unresolved']+=int(anno['whole_partner_unresolved'])
    if row and anno['primary_class']=='UNKNOWN' and row.get('prior_status')=='EXPORTED_READBACK':
        result['unknown_structure_readback_retained']+=1


def execute(run,ref,jobs,controls_only=False):
    start=time.monotonic(); identity,entries=collect_evidence(run,ref,jobs)
    print(f"identity check: {identity['source_shards_built']} new / {identity['source_shards_reused']} reused",flush=True)
    sources=source_catalog(run); sr=run/'sources'
    controls_ccd={'ATP','GOL','GTP','FAD','CLR','EDO','PGO','MPD','IMD','ACN','ZN','SO4','BGC','PLM','OLC','CIT','DMS','EPE'}
    observed=set(identity['baseline_ccds']); supplemental_codes={m['ccd'] for f in sorted((run/'identity-v4').glob('*.jsonl')) for x in jsonl(f) if x.get('candidate_is_increment')
                        for m in x['partner']['members'] if m.get('ccd')}-observed
    ccd=ccd_table(sr/'components.cif.gz',observed|controls_ccd|supplemental_codes)
    bird=bird_table(sr/'prd-all.cif.gz'); ontology,ontology_scope=chebi_ontology(sr/'chebi-257.obo.gz',
                       {k for d in ccd.values() for k in d['inchikeys']})
    from rdkit import rdBase
    sources['methods']={'version':VERSION,'rdkit':rdBase.rdkitVersion,'gemmi':gemmi.__version__,'chebi':ontology_scope,
                        'peptide_length_working_boundary':35,'covalent_role_policy':'Chemical type never updates experimental role, QC or intake.'}
    write(run/'source_receipt.json',sources)
    cls=classifier(ccd,bird,ontology); audit=controls(run,ccd,cls); write(run/'controls.json',audit)
    measures=audit['holdout']; trial=measures['plus_rdkit']; cheb=measures['plus_chebi']; basic=measures['plus_bird']
    if any(v['MISCLASSIFIED'] for m in (measures,audit['fresh_holdout']) for v in m.values()) or cheb['lost_correct_cases'] or trial['lost_correct_cases']:
        raise ValueError(f'Holdout regression/misclassification: see {run}/controls.json')
    method='plus_rdkit' if trial['CORRECT']>cheb['CORRECT'] else 'plus_chebi' if cheb['CORRECT']>basic['CORRECT'] else 'plus_bird'
    print(f'controls: {audit["independent_count"]} independent; method {method}; {measures[method]["CORRECT"]} correct',flush=True)
    if controls_only:return
    results=run/'results';results.mkdir(exist_ok=True)
    coverage=defaultdict(annotation_count); method_counts=defaultdict(Counter); stages_count=defaultdict(Counter)
    identities={}; union_ids=set(); base_ids=set(); new_ids=set(); deltas=Counter(); part_status=defaultdict(Counter)
    entries_by_batch=defaultdict(Counter); entries_states=defaultdict(Counter); rollup=Counter(); resources=defaultdict(Counter)
    variant_sources={'bird':'bird_added_labels','chebi':'chebi','rdkit':'rdkit'}
    with jsonl_output(results/'ccd_labels.jsonl') as emit:
        for k in sorted(observed|supplemental_codes):
            ss,refs=cls(single_ccd(k)); d=ccd.get(k,{})
            emit({'ccd':k,'reference':{i:v for i,v in d.items() if i not in {'bond_graph'}},'annotation':ss[method],
                  'source_evidence':refs,'stages':ss,'version':VERSION,'baseline_observed':k in observed})
            if k not in observed:
                count_annotation(coverage['supplemental_only_ccds'],ss[method]);continue
            count_annotation(coverage['unique_baseline_ccd'],ss[method])
            for name,row in ss.items():stages_count[f'baseline_ccd/{name}'][row['annotation_status']]+=1
            for item,v in (('chebi',refs['chebi']['added']),('rdkit',refs['rdkit']['added'])):
                if v:resources[f'unique_baseline_ccd/{item}'].update(v)
            if refs['chebi']['status']=='EXACT_FULL_INCHIKEY':rollup['ccd_exact_chebi']+=1
            elif refs['chebi']['status']!='EXACT_FULL_INCHIKEY':rollup['ccd_'+refs['chebi']['status']]+=1
            if refs['chebi'].get('excluded_in_scope'):rollup['ccd_chebi_label_outside_project_scope']+=1
            if refs['chebi'].get('conflicting_labels'):rollup['ccd_chebi_matching_terms_disagree']+=1
    from contextlib import ExitStack
    with ExitStack() as stack:
        emits={k:stack.enter_context(jsonl_output(results/f'{k}.jsonl')) for k in ('baseline_system_labels','supplemental_system_labels')}
        for shard in sorted((run/'identity-v4').glob('*.jsonl')):
            for x in jsonl(shard):
                p=x.pop('partner'); ss,evidence=cls(p); ann=ss[method]; cid=p['identity']; scope='baseline' if 'candidate_is_increment' not in x else 'new_candidate' if x['candidate_is_increment'] else 'overlay_existing'
                key=(x['batch'],x['pdb_id'],x['candidate_key']); rollup[f'{scope}_rows']+=1
                if scope=='baseline':
                    assert key not in base_ids;base_ids.add(key);entries_states[x['pdb_id']][ann['annotation_status']]+=1
                else:
                    deltas[(scope,x['evidence_layer'])]+=1
                    if scope=='new_candidate':assert key not in base_ids and key not in new_ids;new_ids.add(key)
                if scope!='overlay_existing':union_ids.add(key)
                for level in {scope,f'{scope}/batch_{x["batch"]}'}:count_annotation(coverage[level],ann,x)
                for name,a in ss.items():stages_count[f'{scope}/{name}'][a['annotation_status']]+=1
                for q,field in variant_sources.items():
                    val=evidence[field] if q=='bird' else evidence[field]['added']
                    if val:method_counts[scope][q+'_systems_with_additions']+=1;resources[f'{scope}/{q}'].update(val)
                for record in evidence['ccd_mmcif_evidence']:
                    if 'type_disagreement' in record:
                        method_counts[scope]['CCD_type_metadata_resolved' if record['type_disagreement']['resolved'] else 'CCD_type_metadata_still_unresolved']+=1
                if evidence['bird_verified']:
                    method_counts[scope]['bird_any_match']+=1
                    method_counts[scope]['bird_exact_graph_or_ccd']+=any('not independently verified' not in i['match'] for i in evidence['bird_verified'])
                row={**x,'chemical_identity':cid,'chemical_identity_scope':'source_scoped_review' if cid.startswith('review:') else 'deduplicated_full_partner',
                     'partner_identity_evidence':p,'classification':ann,'sources_of_labels':evidence,'method_version':VERSION,'chosen_stage':method,
                     'classification_changes_existing_decision':False}
                emits['baseline_system_labels' if scope=='baseline' else 'supplemental_system_labels'](row)
                i=identities.setdefault(cid,{'chemical_identity':cid,'identity_scope':row['chemical_identity_scope'],'baseline_systems':0,
                       'evidence_layers':Counter(),'primary_profiles':Counter(),'status_profiles':Counter(),'example':None,'stages':None,
                       'source_addition_statuses':Counter(),'basis_prd':set(),'baseline_stages':None})
                i['evidence_layers'][x['evidence_layer']]+=1
                i['baseline_systems']+=scope=='baseline'; i['primary_profiles'][ann['primary_class']]+=1
                i['status_profiles'][ann['annotation_status']]+=1
                i['basis_prd'].update(v['prd'] for v in evidence['bird_verified'])
                i['source_addition_statuses'].update(q for q,field in variant_sources.items() if (evidence[field] if q=='bird' else evidence[field]['added']))
                part_status[cid][f'{scope}/{ann["annotation_status"]}']+=1
                rank={'UNKNOWN':0,'PARTIAL':1,'LABELED':2}
                if scope=='baseline' and (i['baseline_stages'] is None or rank[ann['annotation_status']]>rank[i['baseline_stages'][method]['annotation_status']]):
                    i['baseline_stages']=ss
                if i['example'] is None or rank[ann['annotation_status']]>rank[i['stages'][method]['annotation_status']]:
                    i['example']={k:x[k] for k in ('pdb_id','revision','candidate_key','evidence_layer','source_path','discovery_path')}
                    i['stages']=ss
    assert len(base_ids)==46855 and len(new_ids)==9 and len(union_ids)==46864 and rollup['overlay_existing_rows']==6
    assert {p for _,p,_ in union_ids}<=entries.keys()
    with jsonl_output(results/'partner_labels.jsonl') as emit:
        for k,value in sorted(identities.items()):
            value['basis_prd']=sorted(value['basis_prd']);value['layer_statuses']=part_status[k]; emit(value)
            level='unique_ambiguous_group' if k.startswith('review:') else 'unique_complete_partner'
            count_annotation(coverage[level],value['stages'][method]); rollup[level]+=1
            if value['baseline_systems']:
                v=value['baseline_stages'][method]
                count_annotation(coverage['baseline_'+level],v)
            if len(value['primary_profiles'])>1:rollup['identities_with_differing_primary_profile']+=1
    with jsonl_output(results/'entry_dispositions.jsonl') as emit:
        by_new=Counter(p for _,p,_ in new_ids)
        for p,x in sorted(entries.items()):
            status=entries_states[p]; disp='ANNOTATED_OR_EXPLICIT_UNKNOWN' if status or by_new[p] else 'NO_CANDIDATE_IN_BASELINE_OR_RECOVERY'
            entries_by_batch[x['batch']][disp]+=1
            emit({'pdb_id':p,'batch':x['batch'],'revision':x['revision'],'source_path':x['source_path'],'source_sha256':x['source_sha256'],
                  'candidate_evidence':x['discovery_path'],'dataset_route':x['route'],'status':disp,'baseline_annotations':status,
                  'unique_candidate_increment':by_new[p],'previous_candidate_count':x['candidate_count'],
                  'no_candidate_evidence':{k:x.get(k) for k in ('baseline_status','baseline_reason','selection_reasons')} if not status else None,
                  'establishes_no_target':False})
    assert sum(sum(v.values()) for v in entries_by_batch.values())==2200
    artifacts={f.name:{'bytes':f.stat().st_size,'sha256':sha256(f)} for f in sorted(results.glob('*.jsonl'))}
    summary={'version':VERSION,'method':method,'input':{k:v for k,v in identity.items() if k!='baseline_ccds'},'coverage':coverage,
             'rollup':rollup,'supplemental_only_ccds':sorted(supplemental_codes),'entries_by_batch':entries_by_batch,'supplemental_counts':{'/'.join(k):v for k,v in sorted(deltas.items())},
             'stage_results':stages_count,'sources_of_labels':method_counts,'sources_of_labels_by_tag':resources,
             'validation':{k:v for k,v in audit.items() if k!='cases'},'source_snapshots':sources,'artifacts':artifacts,
             'reference_resource_bytes':sum(sources[r]['bytes'] for r in ('CCD','BIRD','ChEBI 257')),
             'coordinate_downloads':0,'structures_written_or_changed':0,'role_qc_or_intake_changes':0,'elapsed_s':round(time.monotonic()-start,2)}
    write(run/'summary.json',summary)
    print(f'annotated: {len(base_ids)} baseline / {len(new_ids)} incremental candidates / {len(entries)} entries; {summary["elapsed_s"]} s',flush=True)


def annotate_batch(evidence_path, reference_run, output, jobs):
    """Classify a frozen candidate batch with the validated historical method and references."""
    start=time.monotonic(); src=load(evidence_path); records=src['candidates']; batch=src['batch']
    ids={(r['pdb_id'],r['candidate_key']) for r in records}
    if (len(ids)!=len(records) or len(records)!=src['expected_candidates'] or
            len(set(src['entries']))!=len(src['entries']) or len(src['entries'])!=src['expected_entries'] or
            any(x['batch']!=batch or x['pdb_id'] not in src['entries'] for x in records)):
        raise ValueError('Frozen chemistry evidence count or identity mismatch')
    ref=reference_run.resolve(strict=True); study=load(ref/'summary.json'); receipt=load(ref/'source_receipt.json')
    method=study['method']
    if study['version']!=VERSION or method not in {'plus_bird','plus_chebi','plus_rdkit'}:
        raise ValueError('Unsupported validated chemistry method')
    files={}
    for name in ('CCD','BIRD','ChEBI 257'):
        r=receipt[name]; path=Path(r['path']).resolve(strict=True)
        if not path.is_relative_to(ref/'sources') or path.stat().st_size!=r['bytes'] or sha256(path)!=r['sha256']:
            raise ValueError(f'Chemical reference changed: {name}')
        files[name]=path
    from rdkit import rdBase
    seal={'candidate_evidence_sha256':sha256(evidence_path),'code_sha256':sha256(Path(__file__)),
          'reference_receipt_sha256':sha256(ref/'source_receipt.json'),'method':method,
          'version':VERSION,'gemmi':gemmi.__version__,'rdkit':rdBase.rdkitVersion}
    output=Path(output); proof=output.with_suffix('.receipt.json')
    if proof.is_file():
        previous=load(proof)
        if previous['seal']!=seal or previous['results_sha256']!=sha256(output):
            raise ValueError('Sealed chemistry outputs changed; use a new Run')
        return previous | {'reused':True,'elapsed_s':round(time.monotonic()-start,3)}
    made=[x for _,rows_ in from_sources(records,[],jobs) for x in rows_] if records else []
    codes={m['ccd'] for x in made for m in x['partner']['members'] if m.get('ccd')}
    ccd=ccd_table(files['CCD'],codes) if codes else {}
    bird=bird_table(files['BIRD']) if made else {}
    ontology,_=chebi_ontology(files['ChEBI 257'],{k for d in ccd.values() for k in d['inchikeys']}) if ccd else ({},None)
    classify=classifier(ccd,bird,ontology); classes=Counter(); coverage=Counter()
    with jsonl_output(output) as emit:
        for row in sorted(made,key=lambda x:(x['pdb_id'],x['candidate_key'])):
            partner=row.pop('partner'); stages,labels=classify(partner); ann=stages[method]
            classes[ann['primary_class']]+=1; coverage[ann['annotation_status']]+=1
            emit({k:row[k] for k in ('batch','pdb_id','candidate_key','revision','source_sha256')} |
                 {'chemical_identity':partner['identity'],'partner_identity_evidence':partner,
                  'classification':ann,'sources_of_labels':labels,'method_version':VERSION,'chosen_stage':method,
                  'classification_changes_existing_decision':False})
    summary={'seal':seal,'batch':batch,'entries':len(src['entries']),'candidates':len(made),'primary_classes':dict(classes),
             'annotation_states':dict(coverage),'results_sha256':sha256(output),'elapsed_s':round(time.monotonic()-start,3)}
    write(proof,summary)
    return summary | {'reused':False}


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('run',type=Path);a.add_argument('--jobs',type=int,default=16)
    a.add_argument('--collect-only',action='store_true');a.add_argument('--controls-only',action='store_true')
    a.add_argument('--batch-evidence',type=Path);a.add_argument('--reference-run',type=Path)
    args=a.parse_args();r=args.run;ref=r.parent/REF
    if args.batch_evidence or args.reference_run:
        if not (args.batch_evidence and args.reference_run) or args.collect_only or args.controls_only:
            a.error('Incremental mode requires both --batch-evidence and --reference-run')
        print(pack(annotate_batch(args.batch_evidence,args.reference_run,r/'chemical_labels.jsonl',args.jobs)),flush=True)
    elif args.collect_only:
        z,_=collect_evidence(r,ref,args.jobs);print(pack({k:v for k,v in z.items() if k!='baseline_ccds'}),flush=True)
    else:execute(r,ref,args.jobs,args.controls_only)
