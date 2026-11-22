"""Project deposited three-atom coordination angles without inventing covalent bonds."""
import argparse
from collections import Counter, defaultdict
from math import acos, degrees, dist
from pathlib import Path
import re

import gemmi

from holo_candidate_batch import atomic_json, now, revision_from_cif, sha256
from holo_m2_blind_audit import rd
from holo_m2_topology_readback import rows, lookup
from holo_v1_export import deposited_copy_matches

ANGLE='_pdbx_struct_conn_angle.'
TRACE='_aphcp_coordination_angle_source.'
KEYS=(('label_asym_id','label_asym_id'),('label_comp_id','label_comp_id'),('label_seq_id','label_seq_id'),
      ('auth_seq_id','auth_seq_id'),('PDB_ins_code','pdbx_PDB_ins_code'),('label_atom_id','label_atom_id'),('label_alt_id','label_alt_id'))


def normalized(value):
    return None if value in (False,None,'','.','?') else value


def residue_atom_signature(source,prefix=''):
    return tuple(normalized(source.get(prefix+(f if prefix else a))) for f,a in KEYS)


def replace_category(block, category, records):
    if records:
        keys=list(records[0]); block.set_mmcif_category(category,{k:[x[k] for x in records] for k in keys})


def angle_xyz(a,b,c):
    a,b,c=([float(row[f'Cartn_{d}']) for d in 'xyz'] for row in (a,b,c))
    u=[i-j for i,j in zip(a,b)]; v=[i-j for i,j in zip(c,b)]
    den=dist(a,b)*dist(c,b)
    return degrees(acos(max(-1.0,min(1.0,sum(x*y for x,y in zip(u,v))/den)))) if den else float('inf')


def project(run, parent, p):
    proof=rd(run/'source_proofs'/f'{p}.json'); raw=run/proof['relative_path']; source,rev=revision_from_cif(raw,p)
    assert sha256(raw)==proof['sha256'] and rev==proof['observed_revision']
    originals=rows(source,ANGLE); source_angles=lookup(originals,'id'); source_atoms=rows(source,'_atom_site.')
    source_atom_index=defaultdict(list)
    for a in source_atoms:
        source_atom_index[(a['pdbx_PDB_model_num'],residue_atom_signature(a))].append(a)
    output_receipts=[]; misses=Counter()
    for rec in rd(run/'system_results'/f'{p}.json'):
        if not rec.get('artifact'): continue
        art=rec['artifact']; path=run/art['path']; mp=rd(run/rec['mapping']); assert sha256(path)==art['sha256']
        doc=gemmi.cif.read(str(path)); block=doc.sole_block(); model=str(rec['representative_model']); originals_by_output=rows(block,'_aphcp_atom_source.')
        output_atoms=lookup(rows(block,'_atom_site.'),'id'); remap={ (r['source_atom_site_id'],r['assembly_instance_id']):r['atom_site_id'] for r in originals_by_output }
        cols=[]; traces=[]; seen=set(); finite=[]; nonlocal_found=0
        for row in originals:
            syms=[row[f'ptnr{j}_symmetry'] for j in (1,2,3)]
            if not all(syms): misses['missing_deposited_symmetry']+=1; continue
            ends=[source_atom_index[(model,residue_atom_signature(row,f'ptnr{j}_'))] for j in (1,2,3)]
            if any(len(a)!=1 for a in ends): misses['unresolved_source_atom_or_altloc']+=1; continue
            ends=[a[0] for a in ends]
            labels=[a['label_asym_id'] for a in ends]
            if originals_by_output:
                choices=deposited_copy_matches(source,mp['instances'],syms,labels)
                proposals=[[(remap.get((a['id'],matches[0][0])),matches[0][0]) for a,matches in zip(ends,group)] for group in choices if all(len(s)==1 for s in group)]
            else:
                proposals=[[(a['id'],next((k for k,v in mp['instances'].items() if v['label_asym']==a['label_asym_id']),None)) for a in ends]] if set(syms)=={'1_555'} else []
            if not any(all(i in output_atoms and copy for i,copy in proposal) for proposal in proposals):
                misses['endpoint_outside_declared_output_or_ambiguous']+=1
                continue
            for proposal in proposals:
                ids=tuple(i for i,_ in proposal)
                if any(i not in output_atoms for i in ids) or any(not k for _,k in proposal) or len(set(ids))!=3: continue
                endpoints=[output_atoms[i] for i in ids]
                key=(row['id'],ids)
                if key in seen: continue
                seen.add(key)
                target=float(row['value']); error=abs(angle_xyz(*endpoints)-target)
                if error>0.12:
                    misses['deposited_geometry_not_reproduced']+=1; continue
                finite.append(error); new=dict(row); new['id']=str(len(cols)+1)
                trace={'output_angle_id':new['id'],'source_angle_id':row['id'],'source_model_id':model}
                for j,(a,(atom_id,copy),o,sym) in enumerate(zip(ends,proposal,endpoints,syms),1):
                    stem=f'ptnr{j}_'; new[stem+'label_asym_id']=o['label_asym_id']; new[stem+'auth_asym_id']=o['auth_asym_id']; new[stem+'symmetry']='1_555'
                    trace.update({stem+'source_atom_site_id':a['id'],stem+'output_atom_site_id':atom_id,
                                  stem+'source_symmetry':sym,stem+'assembly_instance_id':copy})
                cols.append(new); traces.append(trace); nonlocal_found+=len(set(syms))>1
        if not cols: continue
        original_links=rows(block,'_struct_conn.'); replace_category(block,ANGLE,cols); replace_category(block,TRACE,traces)
        out=parent/p/rec['candidate_key']/'system.cif'; out.parent.mkdir(parents=True,exist_ok=True); doc.write_file(str(out))
        disk=gemmi.cif.read(str(out)).sole_block(); read_atoms=lookup(rows(disk,'_atom_site.'),'id')
        assert read_atoms==output_atoms and original_links==rows(disk,'_struct_conn.')
        assert cols==rows(disk,ANGLE) and traces==rows(disk,TRACE)
        for row,t in zip(cols,traces):
            a=source_angles[t['source_angle_id']]; tr=[]
            for j in (1,2,3):
                stem=f'ptnr{j}_'; out_id=t[stem+'output_atom_site_id']; m=read_atoms[out_id]; tr.append(m)
                assert residue_atom_signature(row,stem)==residue_atom_signature(m), (j,row,m)
                assert a[stem+'symmetry']==t[stem+'source_symmetry']
                if originals_by_output: assert remap[t[stem+'source_atom_site_id'],t[stem+'assembly_instance_id']]==out_id
            assert abs(angle_xyz(*tr)-float(a['value']))<=0.12
        output_receipts.append({'base_file':art['path'],'file':str(out.relative_to(run)),'base_sha256':art['sha256'],'file_sha256':sha256(out),
                                'angles':len(cols),'nonlocal_angles':nonlocal_found,'max_deviation_degrees':max(finite),
                                'structure_atoms_unchanged':True,'covalent_topology_unchanged':True})
    return {'pdb_id':p,'exact_revision':rev,'source_sha256':proof['sha256'],'source_angles':len(originals),'outputs':output_receipts,'nonwritten_row_encounters':dict(misses)}


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('run',type=Path); ap.add_argument('entries',nargs='+'); args=ap.parse_args()
    run=args.run.resolve(); dst=run/'angle_review'; dst.mkdir(exist_ok=True); result=[]
    for p in args.entries:
        assert re.fullmatch('[0-9][A-Z0-9]{3}',p)
        try:
            record=project(run,dst,p)
        except (KeyError,ValueError,AssertionError,RuntimeError) as e:
            record={'pdb_id':p,'status':'REVIEW','reason':f'{type(e).__name__}: {str(e)[:500]}'}
        print(p,record.get('source_angles'),sum(a['angles'] for a in record.get('outputs',[])),record.get('reason',''),flush=True)
        result.append(record)
    atomic_json(dst/'receipt.json',{'utc':now(),'scope':'Explicit deposited coordination angles, not extra bonds','entries':result})


if __name__=='__main__': main()
