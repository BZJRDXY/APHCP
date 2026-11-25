"""Check whether fixed multi-chain assignments are equivalent under actual Assembly symmetry."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys
import time

import gemmi
import numpy as np
from scipy.spatial import cKDTree

import apo_c_pair
import apo_c_pilot
import apo_c_triplet_regression
import holo_source_mapping
from apo_c_pair import align_coordinates, ca, get_sequence, revision, sequence_map
from apo_c_pilot import now, read_json, sha, write_json
from apo_c_triplet_regression import (component_contacts, positive_heavy, residue_index, target_members,
                                      unobserved_records)
from holo_source_mapping import assembly_instances


def matrix(operator):
    return np.array(operator.transform.mat.tolist()), np.array(operator.transform.vec.tolist())


def validate_assignment(assignment, required):
    holo = [x[0] for x in assignment]
    candidate = [x[1] for x in assignment]
    if len(holo) != len(set(holo)) or set(holo) != set(required):
        raise ValueError('assignment does not cover each required Holo chain exactly once')
    if len(candidate) != len(set(candidate)):
        raise ValueError('one candidate chain copy cannot cover multiple required Holo chains')


def atom_map(rows, operator):
    result = {}
    for author, residue in rows:
        for atom in positive_heavy(residue, operator):
            key = (author, str(residue.seqid), residue.name, atom['name'], atom['element'], atom['altloc'],
                   round(atom['occupancy'], 4))
            if key in result:
                raise ValueError(f'ambiguous atom key in Assembly: {key}')
            result[key] = np.array(atom['xyz'])
    return result


def assignment_evidence(assignment, hs, cs, candidate_block, hi, ci, ligand):
    holo_points, candidate_points, pocket_holo, pocket_candidate = [], [], [], []
    site_rows, site_atoms, keyed_site = [], [], {}
    ligand_tree = cKDTree(ligand)
    for hkey, ckey in assignment:
        hop, hlabel, hrows = hi[hkey]
        cop, clabel, crows = ci[ckey]
        hres, cres = residue_index(hrows), residue_index(crows)
        unobserved = unobserved_records(candidate_block, '_pdbx_unobs_or_zero_occ_residues.', clabel, cs[0].num)
        mapping, seq = sequence_map(get_sequence(hs, hrows), get_sequence(cs, crows))
        for hn, cn in mapping.items():
            if hn not in hres or cn not in cres:
                continue
            hp, cp = ca(hres[hn][1], hop), ca(cres[cn][1], cop)
            if hp is not None and cp is not None:
                holo_points.append(hp)
                candidate_points.append(cp)
        for hn, (hauth, hr) in sorted(hres.items()):
            hatoms = positive_heavy(hr, hop)
            if not hatoms or np.min(ligand_tree.query(np.array([a['xyz'] for a in hatoms]))[0]) > 4.0:
                continue
            cn = mapping.get(hn)
            row = {'holo_instance': hkey, 'holo_label_seq': hn, 'holo_author_seq': str(hr.seqid),
                   'holo_component': hr.name, 'candidate_instance': ckey, 'candidate_label_seq': cn,
                   'candidate_deposited_unobserved_residue_records': unobserved.get(cn, [])}
            if cn in cres:
                cauth, cr = cres[cn]
                catoms = positive_heavy(cr, cop)
                names = {a['name'] for a in catoms}
                row.update({'candidate_author_chain': cauth, 'candidate_author_seq': str(cr.seqid),
                            'candidate_component': cr.name, 'candidate_has_CA': 'CA' in names,
                            'missing_relative_holo_heavy_atoms': sorted({a['name'] for a in hatoms} - names)})
                hp, cp = ca(hr, hop), ca(cr, cop)
                if hp is not None and cp is not None:
                    pocket_holo.append(hp)
                    pocket_candidate.append(cp)
                for atom in catoms:
                    record = atom | {'label_seq': cn, 'author_chain': cauth,
                                     'author_seq': str(cr.seqid), 'residue': cr.name}
                    site_atoms.append(record)
                    key = (hkey, hn, atom['name'], atom['altloc'])
                    if key in keyed_site:
                        raise ValueError(f'ambiguous mapped site atom {key}')
                    keyed_site[key] = np.array(atom['xyz'])
            else:
                row['candidate_has_CA'] = False
            site_rows.append(row)
    if len(holo_points) < 3 or len(pocket_holo) < 3:
        raise ValueError('insufficient chain or site CA atoms')
    R, t, rmsd = align_coordinates(np.array(holo_points), np.array(candidate_points))
    _, _, local_rmsd = align_coordinates(np.array(pocket_holo), np.array(pocket_candidate))
    projected = ligand @ R + t
    contacts = component_contacts(ci, {x[1] for x in assignment}, site_atoms, projected)
    return {'assignment': assignment, 'global_ca_count': len(holo_points), 'global_ca_rmsd_A': rmsd,
            'pocket_ca_count': len(pocket_holo), 'local_pocket_ca_rmsd_A': local_rmsd,
            'site_residues': site_rows, 'actual_components': contacts,
            '_site_atoms': keyed_site, '_projected_target': projected}


def symmetry_proof(ci, first, second, result_a, result_b):
    qkey = first[0][0]
    from_key = dict(first)[qkey]
    to_key = dict(second)[qkey]
    op1, label1, _ = ci[from_key]
    op2, label2, _ = ci[to_key]
    if label1 != label2:
        raise ValueError('the two chain copies do not arise from the same deposited subchain')
    R1, t1 = matrix(op1)
    R2, t2 = matrix(op2)
    G = R2 @ np.linalg.inv(R1)
    g = t2 - G @ t1
    instance_map = {}
    max_atom_error, atom_count = 0.0, 0
    for key, (operator, label, rows) in ci.items():
        rot, trans = matrix(operator)
        matches = [candidate for candidate, (op2, label2, _) in ci.items()
                   if label == label2 and np.allclose(matrix(op2)[0], G @ rot, atol=1e-6)
                   and np.allclose(matrix(op2)[1], G @ trans + g, atol=1e-6)]
        if len(matches) != 1:
            raise ValueError(f'Assembly symmetry image is not unique for {key}: {matches}')
        target = matches[0]
        instance_map[key] = target
        a, b = atom_map(rows, operator), atom_map(ci[target][2], ci[target][0])
        if set(a) != set(b):
            raise ValueError(f'atom inventory differs under symmetry for {key}')
        if a:
            errors = [np.linalg.norm(G @ a[x] + g - b[x]) for x in a]
            max_atom_error = max(max_atom_error, float(max(errors)))
            atom_count += len(errors)
    a, b = result_a['_site_atoms'], result_b['_site_atoms']
    if set(a) != set(b):
        raise ValueError('mapped site atom identities differ between assignments')
    site_error = max((float(np.linalg.norm(G @ a[k] + g - b[k])) for k in a), default=0.0)
    p, q = result_a['_projected_target'], result_b['_projected_target']
    target_error = float(np.max(np.linalg.norm(p @ G.T + g - q, axis=1)))
    def contact_index(result, transform):
        rows = result['actual_components']['nonwater'] + result['actual_components']['nearby_waters']
        return {(transform.get(x['instance'], x['instance']), x['author_chain'], x['author_seq'], x['component'], x['kind']):
                x['min_distance_to_mapped_site_atoms_A'] for x in rows}
    x, y = contact_index(result_a, instance_map), contact_index(result_b, {})
    same_contacts = set(x) == set(y) and all(abs(x[k] - y[k]) <= 0.001 for k in x)
    return {'rotation': G.tolist(), 'translation': g.tolist(), 'instance_map': instance_map,
            'symmetry_squared_rotation_error': float(np.linalg.norm(G @ G - np.eye(3))),
            'symmetry_squared_translation_error_A': float(np.linalg.norm(G @ g + g)),
            'assembly_atoms_compared': atom_count, 'assembly_atom_max_error_A': max_atom_error,
            'mapped_site_atoms_compared': len(a), 'mapped_site_atom_max_error_A': site_error,
            'observed_holo_target_probe_max_error_A': target_error,
            'nearby_component_identity_and_distance_equivalent': same_contacts,
            'nearby_component_count': len(x)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', required=True)
    p.add_argument('--case', required=True)
    args = p.parse_args()
    run = Path(args.run)
    entries = [x for x in read_json(run / 'inputs.json') if x['id'] == args.case]
    if len(entries) != 1:
        p.error('case must match one fixed input')
    entry, old = entries[0], entries[0]['baseline']
    if sha(Path(entry['baseline_path'])) != entry['baseline_sha256']:
        raise ValueError('baseline evidence changed')
    if sha(Path(entry['sample']['source_uri'])) != entry['sample']['identity']['source_sha256']:
        raise ValueError('Holo fixed original changed')
    if sha(Path(old['candidate_source']['uri'])) != old['candidate_source']['sha256']:
        raise ValueError('candidate fixed original changed')
    assignments = old['near_equivalent_assignments']
    if len(assignments) != 2:
        raise ValueError('this bounded check requires exactly two recorded candidate assignments')
    required = [x[0] for x in entry['sample']['identity']['binding_instances']]
    for assignment in assignments:
        validate_assignment(assignment, required)
    base = run / 'multichain' / args.case
    base.mkdir(parents=True, exist_ok=True)
    n = 1
    while (base / f'attempt{n:03d}').exists():
        n += 1
    out = base / f'attempt{n:03d}'
    out.mkdir()
    deps = {}
    for path in [Path(__file__), Path(apo_c_pair.__file__), Path(apo_c_pilot.__file__),
                 Path(apo_c_triplet_regression.__file__), Path(holo_source_mapping.__file__)]:
        digest = sha(path)
        dest = run / 'code' / f'{path.stem}-{digest[:12]}.py'
        if not dest.exists():
            shutil.copy2(path, dest)
        deps[path.name] = digest
    result = {'case': args.case, 'started_utc': now(), 'dependency_sha256': deps}
    start = time.monotonic()
    try:
        hb = gemmi.cif.read(entry['sample']['source_uri']).sole_block()
        cb = gemmi.cif.read(old['candidate_source']['uri']).sole_block()
        hv, _ = revision(hb, entry['sample']['pdb'])
        cv, _ = revision(cb, old['request']['candidate'])
        if hv != entry['sample']['identity']['revision'] or cv != old['candidate_source']['revision']:
            raise ValueError('fixed source revision differs')
        hs, cs = gemmi.make_structure_from_block(hb), gemmi.make_structure_from_block(cb)
        hi = assembly_instances(hs, str(entry['sample']['identity']['assembly']), hs[0], include_water=True)
        ci = assembly_instances(cs, old['request']['assembly'], cs[0], include_water=True)
        members, ligand = target_members(hb, hi, entry['sample']['identity'])
        a = assignment_evidence(assignments[0], hs, cs, cb, hi, ci, ligand)
        b = assignment_evidence(assignments[1], hs, cs, cb, hi, ci, ligand)
        proof = symmetry_proof(ci, assignments[0], assignments[1], a, b)
        site_a = [(x['holo_instance'], x['holo_label_seq'], x['candidate_label_seq'], x.get('candidate_has_CA'),
                   x.get('missing_relative_holo_heavy_atoms'), x.get('candidate_deposited_unobserved_residue_records'))
                  for x in a['site_residues']]
        site_b = [(x['holo_instance'], x['holo_label_seq'], x['candidate_label_seq'], x.get('candidate_has_CA'),
                   x.get('missing_relative_holo_heavy_atoms'), x.get('candidate_deposited_unobserved_residue_records'))
                  for x in b['site_residues']]
        proof['residue_mapping_and_missing_atom_evidence_equal'] = site_a == site_b
        passed = (proof['symmetry_squared_rotation_error'] < 1e-6 and
                  proof['symmetry_squared_translation_error_A'] < 1e-6 and
                  proof['assembly_atom_max_error_A'] < 0.001 and proof['mapped_site_atom_max_error_A'] < 0.001 and
                  proof['observed_holo_target_probe_max_error_A'] < 0.001 and
                  proof['nearby_component_identity_and_distance_equivalent'] and
                  proof['residue_mapping_and_missing_atom_evidence_equal'])
        for item in (a, b):
            item.pop('_site_atoms')
            item.pop('_projected_target')
        for name in ['holo_system.cif', 'candidate_assembly.cif']:
            src, dst = Path(entry['baseline_path']).parent / name, out / name
            shutil.copy2(src, dst)
            check = gemmi.read_structure(str(dst))
            if not check or not check[0].count_atom_sites():
                raise ValueError(f'derivative cannot be read: {name}')
        result.update({'execution_status': 'MEASURED', 'holo_system_id': entry['sample']['id'],
                       'holo_status': entry['sample']['holo_status'], 'holo_identity': entry['sample']['identity'],
                       'candidate_source': old['candidate_source'], 'baseline_evidence': entry['baseline_path'],
                       'holo_target_members': members, 'assignments': [a, b], 'symmetry_proof': proof,
                       'mapping_conclusion': 'OBSERVED_ATOMS_SYMMETRY_EQUIVALENT' if passed else 'EQUIVALENCE_NOT_ESTABLISHED',
                       'formal_pair_status': 'PENDING_INDEPENDENT_EVIDENCE_REVIEW',
                       'limits': ['Equivalence concerns this fixed Assembly, target and observed atoms only.',
                                  'Missing site atoms and experimental unbound evidence remain separate questions.',
                                  'Nonblank altloc correlation across residues is not assumed.']})
    except Exception as exc:
        result.update({'execution_status': 'UNRESOLVED', 'error': f'{type(exc).__name__}: {exc}'})
    result['elapsed_s'] = round(time.monotonic() - start, 3)
    result['finished_utc'] = now()
    write_json(out / 'evidence.json', result)
    print(json.dumps({'case': args.case, 'path': str(out), 'status': result['execution_status'],
                      'mapping': result.get('mapping_conclusion'), 'error': result.get('error')}), flush=True)
    if result['execution_status'] != 'MEASURED':
        sys.exit(1)


if __name__ == '__main__':
    main()
