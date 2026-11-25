"""Independent fixed-source checks for a small Apo/Holo evidence regression."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import gemmi
import numpy as np
from scipy.spatial import cKDTree


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for part in iter(lambda: f.read(1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()


def data(path):
    return json.loads(path.read_text(encoding='utf8'))


def last_evidence(run, case, multi=False):
    folder = run / ('multichain' if multi else 'pairs') / case
    paths = sorted(folder.glob('attempt*/evidence.json'))
    if not paths:
        raise ValueError(f'no evidence for {case}')
    d = data(paths[-1])
    if d.get('execution_status') != 'MEASURED':
        raise ValueError(f'latest evidence unresolved: {paths[-1]}')
    return paths[-1], d


def atom_rows(block):
    cat = block.get_mmcif_category('_atom_site.')
    if not cat:
        raise ValueError('no atom_site category')
    for i in range(len(cat['id'])):
        row = {k: v[i] for k, v in cat.items()}
        try:
            occ = float(row['occupancy'])
        except (TypeError, ValueError):
            continue
        if occ <= 0 or str(row['type_symbol']).upper() in {'H', 'D'}:
            continue
        if str(row.get('pdbx_PDB_model_num', '1')) != '1':
            continue
        row['xyz'] = np.array([float(row['Cartn_x']), float(row['Cartn_y']), float(row['Cartn_z'])])
        yield row


def min_dist(a, b):
    if not a or not b:
        raise ValueError('cannot compute distance without both atom sets')
    return float(np.min(cKDTree(np.array([x['xyz'] for x in a])).query(np.array([x['xyz'] for x in b]))[0]))


def candidate_site(block, evidence):
    label = evidence['candidate_source']['label_asym']
    positions = {str(x['candidate_label_seq']) for x in evidence['mapped_site']
                 if x.get('candidate_label_seq') is not None and x.get('candidate_observed')}
    return [x for x in atom_rows(block) if x['label_asym_id'] == label and str(x['label_seq_id']) in positions]


def near(actual, expected, tolerance=0.005):
    return abs(actual - expected) <= tolerance


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', required=True)
    p.add_argument('--output-name', default='independent_source_audit.json')
    args = p.parse_args()
    run = Path(args.run)
    if Path(args.output_name).name != args.output_name:
        raise ValueError('output name must be a simple filename')
    result = {'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'checks': {}}
    checks = result['checks']
    evidence = {}
    for case in ['3ABX_3A64', '2SRC_1FMK', '1AKE_4AKE', '1AKE_4AKE_chainB', '1AKE_1ANK', '5MO0_2XTT',
                 '1GAG_1I44', '1GAG_1I44_composite_target_research', '2SRC_1A07']:
        path, e = last_evidence(run, case)
        evidence[case] = e
        checks[case + '_fixed_source_hashes'] = (
            digest(Path(e['holo_source']['path'])) == e['holo_source']['sha256'] and
            digest(Path(e['candidate_source']['path'])) == e['candidate_source']['sha256'])
        checks[case + '_derivative_readback'] = all(
            gemmi.read_structure(v['path'])[0].count_atom_sites() == v['readback_atoms'] and
            digest(Path(v['path'])) == v['sha256'] for v in e['derivatives'].values())

    e = evidence['3ABX_3A64']
    h = gemmi.cif.read(e['holo_source']['path']).sole_block()
    c = gemmi.cif.read(e['candidate_source']['path']).sole_block()
    ccd = h.get_mmcif_category('_chem_comp_atom.')
    expected = {ccd['atom_id'][i] for i, v in enumerate(ccd['comp_id'])
                if v == 'RCB' and ccd['type_symbol'][i].upper() not in {'H', 'D'}}
    observed = {x['label_atom_id'] for x in atom_rows(h) if x['label_asym_id'] == 'C'
                and x['label_comp_id'] == 'RCB' and str(x['auth_seq_id']) == '1001'}
    checks['3ABX_RCB_32_of_43_from_atom_site_and_CCD'] = len(expected) == 43 and len(observed) == 32 and observed <= expected
    mg = [x for x in atom_rows(c) if x['label_comp_id'] == 'MG' and x['auth_asym_id'] == 'A'
          and str(x['auth_seq_id']) == '396']
    checks['3A64_Mg_actual_site_distance_9_6545A'] = near(min_dist(mg, candidate_site(c, e)), 9.6545)
    checks['3ABX_refinement_start_3A64'] = '3A64' in gemmi.cif.as_string(h.find_value('_refine.pdbx_starting_model'))

    e = evidence['2SRC_1FMK']
    h, c = (gemmi.cif.read(e[k]['path']).sole_block() for k in ['holo_source', 'candidate_source'])
    def ptr(block):
        return [x for x in atom_rows(block) if x['label_comp_id'] == 'PTR' and x['auth_asym_id'] == 'A'
                and str(x['auth_seq_id']) == '527']
    checks['2SRC_and_1FMK_PTR527_present'] = len(ptr(h)) == 16 and len(ptr(c)) == 16
    checks['2SRC_refinement_start_1FMK'] = '1FMK' in gemmi.cif.as_string(h.find_value('_refine.pdbx_starting_model'))
    checks['2SRC_target_ANP_31_heavy_atoms'] = len([x for x in atom_rows(h) if x['label_comp_id'] == 'ANP'
                                                   and x['label_asym_id'] == 'B']) == 31
    holo_positions = {str(x['holo_label_seq']) for x in e['mapped_site']}
    holo_site = [x for x in atom_rows(h) if x['label_asym_id'] == 'A' and str(x['label_seq_id']) in holo_positions]
    holo_ptr_distance = min_dist(ptr(h), holo_site)
    candidate_ptr_distance = min_dist(ptr(c), candidate_site(c, e))
    checks['2SRC_PTR527_distance_to_actual_site'] = (near(holo_ptr_distance, 24.5169) and
        near(e['holo_modifications']['selected_observed'][0]['min_distance_to_actual_site_atoms_A'], holo_ptr_distance))
    checks['1FMK_PTR527_distance_to_actual_site'] = (near(candidate_ptr_distance, 24.4625) and
        near(e['candidate_modifications']['selected_observed'][0]['min_distance_to_actual_site_atoms_A'], candidate_ptr_distance))

    e = evidence['1AKE_4AKE']
    c = gemmi.cif.read(e['candidate_source']['path']).sole_block()
    b36 = [x for x in atom_rows(c) if x['auth_asym_id'] == 'B' and str(x['auth_seq_id']) == '36'
           and x['label_comp_id'] == 'ARG']
    checks['4AKE_chain_B_real_contact_to_mapped_site'] = near(min_dist(b36, candidate_site(c, e)), 2.8763)
    checks['4AKE_large_motion_not_automatic_exclusion'] = e['global_ca_rmsd_A'] > 7 and e['formal_pair_status'].startswith('PENDING')
    e = evidence['1AKE_4AKE_chainB']
    c = gemmi.cif.read(e['candidate_source']['path']).sole_block()
    a53 = [x for x in atom_rows(c) if x['auth_asym_id'] == 'A' and str(x['auth_seq_id']) == '53'
           and x['label_comp_id'] == 'MET']
    checks['4AKE_chain_B_is_separate_observation'] = (e['candidate_source']['author_chain'] == 'B' and
                                                       near(min_dist(a53, candidate_site(c, e)), 2.8763))

    e = evidence['1AKE_1ANK']
    c = gemmi.cif.read(e['candidate_source']['path']).sole_block()
    for comp, expected_distance in [('ANP', 2.5138), ('AMP', 2.5576)]:
        atoms = [x for x in atom_rows(c) if x['label_comp_id'] == comp and x['auth_asym_id'] == 'A']
        checks['1ANK_' + comp + '_real_site_contact'] = near(min_dist(atoms, candidate_site(c, e)), expected_distance)

    e = evidence['5MO0_2XTT']
    c = gemmi.cif.read(e['candidate_source']['path']).sole_block()
    arg = [x for x in atom_rows(c) if x['label_comp_id'] == 'ARG' and x['auth_asym_id'] == 'A'
           and str(x['auth_seq_id']) == '29']
    checks['2XTT_inhibitor_chain_ARG29_real_site_contact'] = near(min_dist(arg, candidate_site(c, e)), 2.6357)

    partial, whole = evidence['1GAG_1I44'], evidence['1GAG_1I44_composite_target_research']
    h = gemmi.cif.read(partial['holo_source']['path']).sole_block()
    c = gemmi.cif.read(partial['candidate_source']['path']).sole_block()
    atoms = list(atom_rows(h))
    p107 = [x for x in atoms if x['label_asym_id'] == 'B' and x['label_comp_id'] == 'PHE'
            and str(x['auth_seq_id']) == '107' and x['label_atom_id'] == 'CZ']
    n112 = [x for x in atoms if x['label_asym_id'] == 'E' and x['label_comp_id'] == '112'
            and str(x['auth_seq_id']) == '101' and x['label_atom_id'] == 'NS']
    checks['1GAG_covalent_endpoints_1_4453A'] = near(min_dist(p107, n112), 1.4453)
    conn = h.get_mmcif_category('_struct_conn.')
    checks['1GAG_source_covalent_connection_present'] = any(
        conn['conn_type_id'][i] == 'covale' and conn['id'][i] == 'covale6'
        for i in range(len(conn.get('id', []))))
    checks['1GAG_scope_extension_adds_15_residues'] = (
        partial['holo_contact_residue_count'] == 16 and whole['holo_contact_residue_count'] == 31 and
        len(whole['baseline_site_difference']['added_holo_label_seq']) == 15)
    checks['1GAG_partial_flags_external_covalent'] = partial['holo_target_scope_status'] == 'EXTERNAL_COVALENT_CONNECTION_REQUIRES_REVIEW'
    unobs = c.get_mmcif_category('_pdbx_unobs_or_zero_occ_residues.')
    absent = {str(unobs['auth_seq_id'][i]) for i in range(len(unobs.get('id', [])))
              if unobs['auth_asym_id'][i] == 'A' and str(unobs['PDB_model_num'][i]) == '1'}
    checks['1I44_unobserved_1167_to_1170'] = {'1167', '1168', '1169', '1170'} <= absent
    missing_loop = [row for row in partial['mapped_site'] if row['holo_author_seq'] in {'1167', '1168', '1169', '1170'}]
    checks['1I44_missing_loop_flag_is_unobserved'] = len(missing_loop) == 4 and all(
        any(x.get('occupancy_flag') == '1' and x.get('deposited_observation_status') == 'UNOBSERVED'
            for x in row.get('candidate_deposited_unobserved_residue_records', [])) for row in missing_loop)
    acp = [x for x in atom_rows(c) if x['label_comp_id'] == 'ACP' and x['auth_asym_id'] == 'A']
    checks['1I44_ACP_contacts_composite_site'] = near(min_dist(acp, candidate_site(c, whole)), 2.6607)
    full_holo = gemmi.cif.read(whole['derivatives']['holo_system.cif']['path']).sole_block()
    checks['1GAG_composite_holo_file_contains_component_112'] = any(
        x['label_comp_id'] == '112' for x in atom_rows(full_holo))
    sub = data(sorted((run / 'subsites' / '112_ACP').glob('attempt*/evidence.json'))[-1])
    hlabel = 'A'
    clabel = whole['candidate_source']['label_asym']
    hrows = list(atom_rows(h))
    crows = list(atom_rows(c))
    X, Y = [], []
    for row in sub['holo_member_contact_residues']:
        if row.get('candidate_label_seq') is None:
            continue
        ha = [a for a in hrows if a['label_asym_id'] == hlabel and str(a['label_seq_id']) == str(row['holo_label_seq'])
              and a['label_atom_id'] == 'CA']
        ca = [a for a in crows if a['label_asym_id'] == clabel and str(a['label_seq_id']) == str(row['candidate_label_seq'])
              and a['label_atom_id'] == 'CA']
        if len(ha) == 1 and len(ca) == 1:
            X.append(ha[0]['xyz'])
            Y.append(ca[0]['xyz'])
    X, Y = np.array(X), np.array(Y)
    x0, y0 = X.mean(axis=0), Y.mean(axis=0)
    u, _, vt = np.linalg.svd((X - x0).T @ (Y - y0))
    flip = np.eye(3)
    flip[-1, -1] = np.linalg.det(u @ vt)
    rot, trans = u @ flip @ vt, y0 - x0 @ u @ flip @ vt
    independent_rmsd = float(np.sqrt(np.mean(np.sum((X @ rot + trans - Y) ** 2, axis=1))))
    member = [a for a in hrows if a['label_asym_id'] == 'E' and a['label_comp_id'] == '112'
              and str(a['auth_seq_id']) == '101']
    probe = np.array([a['xyz'] for a in member]) @ rot + trans
    acp_xyz = np.array([a['xyz'] for a in acp])
    distances = cKDTree(probe).query(acp_xyz)[0]
    checks['1GAG_112_local_frame_independently_recomputed'] = (
        len(X) == 17 and near(independent_rmsd, 1.1131) and near(float(np.min(distances)), 0.0784)
        and int(np.sum(distances <= 3.0)) == 31)
    site_positions = {str(row['candidate_label_seq']) for row in sub['holo_member_contact_residues']
                      if row.get('candidate_label_seq') is not None}
    actual_site = [a for a in crows if a['label_asym_id'] == clabel and str(a['label_seq_id']) in site_positions]
    checks['1I44_ACP_actual_nucleotide_site_contact'] = near(min_dist(acp, actual_site), 2.6607)

    e = evidence['2SRC_1A07']
    checks['1A07_site_not_mapped'] = (len(e['mapped_site']) == 16 and
                                      all(x['candidate_label_seq'] is None for x in e['mapped_site']))
    cb = gemmi.cif.read(e['candidate_source']['path']).sole_block()
    ep = cb.get_mmcif_category('_entity_poly.')
    seq = str(ep['pdbx_seq_one_letter_code_can'][0])
    checks['1A07_construct_107_residues'] = len(''.join(seq.split())) == 107

    path, multi = last_evidence(run, '5MYQ_1RAV', multi=True)
    source_entry = next(x for x in data(run / 'inputs.json') if x['id'] == '5MYQ_1RAV')
    checks['5MYQ_1RAV_fixed_source_hashes'] = (
        digest(Path(source_entry['sample']['source_uri'])) == multi['holo_identity']['source_sha256'] and
        digest(Path(multi['candidate_source']['uri'])) == multi['candidate_source']['sha256'])
    baseline_dir = Path(multi['baseline_evidence']).parent
    checks['5MYQ_1RAV_derivatives_match_fixed_baseline'] = all(
        digest(path.parent / name) == digest(baseline_dir / name) and
        gemmi.read_structure(str(path.parent / name))[0].count_atom_sites() > 0
        for name in ['holo_system.cif', 'candidate_assembly.cif'])
    cb = gemmi.cif.read(multi['candidate_source']['uri']).sole_block()
    ops = cb.get_mmcif_category('_pdbx_struct_oper_list.')
    i = ops['id'].index('2')
    mat = np.array([[float(ops[f'matrix[{x}][{y}]'][i]) for y in (1, 2, 3)] for x in (1, 2, 3)])
    vec = np.array([float(ops[f'vector[{x}]'][i]) for x in (1, 2, 3)])
    proof = multi['symmetry_proof']
    checks['1RAV_recorded_operator_2_matches_proof'] = (
        np.allclose(mat, proof['rotation'], atol=1e-7) and np.allclose(vec, proof['translation'], atol=1e-7))
    checks['1RAV_equivalent_mapping_has_observed_scope_and_missing_residues'] = (
        multi['mapping_conclusion'] == 'OBSERVED_ATOMS_SYMMETRY_EQUIVALENT' and
        proof['residue_mapping_and_missing_atom_evidence_equal'] and
        proof['nearby_component_identity_and_distance_equivalent'] and
        all(a['pocket_ca_count'] == 16 and len(a['site_residues']) == 20 for a in multi['assignments']))
    checks['1RAV_symmetry_is_involution'] = np.allclose(mat @ mat, np.eye(3), atol=1e-7) and np.allclose(mat @ vec + vec, 0, atol=1e-7)
    missing = [row for row in multi['assignments'][0]['site_residues']
               if row['holo_author_seq'] in {'37', '38', '39', '40'} and not row.get('candidate_has_CA')]
    checks['1RAV_four_site_residues_are_zero_occupancy'] = len(missing) == 4 and all(
        any(x.get('occupancy_flag') == '0' and x.get('deposited_observation_status') == 'ZERO_OCCUPANCY'
            for x in row.get('candidate_deposited_unobserved_residue_records', [])) for row in missing)
    for key in checks:
        checks[key] = bool(checks[key])
    result['finished_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    result['all_pass'] = all(bool(x) for x in checks.values())
    result['note'] = 'Fixed-source and coordinate checks only; not an independent human scientific adjudication.'
    out = run / 'receipts' / args.output_name
    with out.open('x', encoding='utf8') as f:
        json.dump(result, f, indent=2)
    print(json.dumps({'path': str(out), 'checks': len(checks), 'passed': sum(bool(x) for x in checks.values()),
                      'all_pass': result['all_pass']}), flush=True)
    if not result['all_pass']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
