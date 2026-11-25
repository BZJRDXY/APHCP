"""Compare sampled source-mapped pocket residues with PLINDER's recorded neighbors."""
import argparse
import json
from pathlib import Path
import sqlite3

import gemmi
import numpy as np


def check(obs, annotations):
    identity = obs['identity']
    target = identity['target_residue_instances']
    if len(target) != 1:
        return {'state':'NOT_SINGLE_RESIDUE'}
    proteins = {x[0] for x in identity['binding_instances']}
    matches = [x for x in annotations if set(x.get('ligand_protein_chains_asym_id') or []) == proteins
               and set(x.get('ligand_instance_chains') or []) == {target[0][0]}]
    if not matches:
        return {'state':'NO_CHAIN_COMPARABLE_ANNOTATION'}
    mapped = json.loads(Path(obs['files']['source_mapping']['path']).read_bytes())
    instances = mapped['instances']
    protein_labels = {instances[k].get('output_label_asym', instances[k]['label_asym']): k for k in proteins}
    ligand = target[0]
    ligand_label = instances[ligand[0]].get('output_label_asym', instances[ligand[0]]['label_asym'])
    block = gemmi.cif.read_file(obs['files']['system_cif']['path']).sole_block()
    a = block.get_mmcif_category('_atom_site.')
    n = len(a['Cartn_x'])
    ligand_xyz, protein_xyz, protein_ids = [], [], []
    for j in range(n):
        elem = (a.get('type_symbol') or ['']*n)[j]
        if str(elem).upper() in ('H','D'):
            continue
        label = a['label_asym_id'][j]
        xyz = tuple(float(a[k][j]) for k in ('Cartn_x','Cartn_y','Cartn_z'))
        if label == ligand_label and a['label_comp_id'][j] == ligand[6]:
            ins = (a.get('pdbx_PDB_ins_code') or [None]*n)[j]
            auth = str(a['auth_seq_id'][j]) + (str(ins) if ins not in (None,'.','?','') else '')
            if auth == str(ligand[5]):
                ligand_xyz.append(xyz)
        elif label in protein_labels:
            seq = (a.get('label_seq_id') or [None]*n)[j]
            if seq not in (None,'.','?',''):
                protein_xyz.append(xyz)
                protein_ids.append((protein_labels[label],str(seq)))
    if not ligand_xyz or not protein_xyz:
        return {'state':'NO_MAPPED_ATOMS', 'ligand_atoms':len(ligand_xyz), 'protein_atoms':len(protein_xyz)}
    la, pa = np.asarray(ligand_xyz), np.asarray(protein_xyz)
    near = set()
    for start in range(0,len(pa),2048):
        distances = ((pa[start:start+2048,None,:]-la[None,:,:])**2).sum(axis=2)
        for index in np.nonzero((distances <= 36.0).any(axis=1))[0]:
            near.add(protein_ids[start+int(index)])
    plinder = set()
    for match in matches:
        for item in match.get('ligand_neighboring_residues') or []:
            bits = item.rsplit('_',4)
            if len(bits) == 5 and bits[0] in proteins:
                plinder.add((bits[0],bits[1]))
    return {'state':'COMPARED', 'cutoff_A':6.0,'ligand_heavy_atoms':len(la),
            'mapped_near_residues':len(near),'plinder_neighbor_residues':len(plinder),
            'intersection':len(near & plinder), 'only_mapped':sorted(near-plinder)[:15],
            'only_plinder':sorted(plinder-near)[:15]}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--database',type=Path,required=True)
    p.add_argument('--sample',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args = p.parse_args()
    db = sqlite3.connect(args.database.resolve().as_uri()+'?mode=ro', uri=True)
    sample = json.loads(args.sample.read_text(encoding='utf-8'))
    results = []
    for group, data in sample['strata'].items():
        for item in data['selected']:
            sid = item['system_id']
            row = db.execute('select observation.document,result.annotations,result.match_level from observation join result using(system_id) where system_id=?',(sid,)).fetchone()
            if row[2] != 'SAME_REVISION_INSTANCE_CHAIN_MATCH':
                results.append({'group':group,'pdb_id':item['pdb_id'],'system_id':sid,'state':'NO_SAME_REVISION_CHAIN_MATCH'})
                continue
            try:
                result = check(json.loads(row[0]),json.loads(row[1]))
            except (OSError, ValueError, KeyError, RuntimeError) as exc:
                result = {'state':'CHECK_FAILED','reason':str(exc)}
            results.append({'group':group,'pdb_id':item['pdb_id'],'system_id':sid,**result})
    with args.output.open('x',encoding='utf-8') as f:
        json.dump(results,f,indent=2,ensure_ascii=False)
    for x in results:
        print(json.dumps({k:v for k,v in x.items() if k not in ('system_id','only_mapped','only_plinder')},ensure_ascii=False))


if __name__ == '__main__':
    main()
