"""Independently check fixed candidate-union triage against deposited mmCIF data."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import gemmi
import numpy as np
from scipy.spatial import cKDTree


def load(path):
    return json.loads(path.read_text(encoding='utf8'))


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def atom_rows(block):
    cat = block.get_mmcif_category('_atom_site.')
    for i in range(len(cat.get('id', []))):
        row = {k: v[i] for k, v in cat.items()}
        try:
            occupancy = float(row['occupancy'])
        except (TypeError, ValueError):
            continue
        if occupancy <= 0 or str(row['type_symbol']).upper() in {'H', 'D'} or str(row.get('pdbx_PDB_model_num', '1')) != '1':
            continue
        row['xyz'] = np.array([float(row['Cartn_x']), float(row['Cartn_y']), float(row['Cartn_z'])])
        alt = row.get('label_alt_id')
        row['alt'] = '' if alt in (None, False, '', '.', '?') else str(alt)
        yield row


def direct_contacts(d, component):
    obs = next(x for x in d['entities'][0]['assembly_observations'] if x['assembly'] == '1')
    block = gemmi.cif.read(d['candidate_source']['path']).sole_block()
    gen = block.get_mmcif_category('_pdbx_struct_assembly_gen.')
    matching = [i for i, a in enumerate(gen.get('assembly_id', [])) if str(a) == '1']
    if len(matching) != 1 or str(gen['oper_expression'][matching[0]]) != '1':
        raise ValueError('direct source-coordinate check requires the recorded identity Assembly operation')
    labels = set(str(gen['asym_id_list'][matching[0]]).split(','))
    op = block.get_mmcif_category('_pdbx_struct_oper_list.')
    idx = op['id'].index('1')
    mat = np.array([[float(op[f'matrix[{x}][{y}]'][idx]) for y in (1, 2, 3)] for x in (1, 2, 3)])
    vec = np.array([float(op[f'vector[{x}]'][idx]) for x in (1, 2, 3)])
    if not np.allclose(mat, np.eye(3), atol=1e-7) or not np.allclose(vec, 0, atol=1e-7):
        raise ValueError('operation 1 is not identity')
    component_labels = {x['label_asym'] for x in obs['actual_component_contacts']['nonwater']
                        if x['component'] == component and x['atoms_within_4A_of_mapped_site'] > 0}
    if obs['label_asym'] not in labels or not component_labels or not component_labels <= labels:
        raise ValueError('selected chain or component is not in the declared Assembly')
    seqs = {str(x['candidate_label_seq']) for x in obs['mapped_site'] if x.get('candidate_label_seq') is not None}
    atoms = list(atom_rows(block))
    site = [x for x in atoms if x['label_asym_id'] == obs['label_asym'] and str(x['label_seq_id']) in seqs]
    ligand = [x for x in atoms if x['label_comp_id'] == component and x['label_asym_id'] in component_labels]
    if not site or not ligand:
        raise ValueError('no actual site or component atoms')
    tree = cKDTree(np.array([x['xyz'] for x in site]))
    unambiguous, both_alt, minimum = 0, 0, None
    for a in ligand:
        for j in tree.query_ball_point(a['xyz'], 4.0):
            b = site[j]
            distance = float(np.linalg.norm(a['xyz'] - b['xyz']))
            if a['alt'] and b['alt']:
                both_alt += 1
            else:
                unambiguous += 1
                minimum = distance if minimum is None else min(minimum, distance)
    return {'component': component, 'unambiguous_atom_pairs_within_4A': unambiguous,
            'both_nonblank_altloc_pairs_requiring_compatibility_check': both_alt,
            'minimum_unambiguous_distance_A': round(minimum, 4) if minimum is not None else None,
            'component_labels': sorted(component_labels), 'site_label_asym': obs['label_asym']}


def direct_nonwater_near_site(d, observation):
    block = gemmi.cif.read(d['candidate_source']['path']).sole_block()
    gen = block.get_mmcif_category('_pdbx_struct_assembly_gen.')
    matching = [i for i, a in enumerate(gen.get('assembly_id', [])) if str(a) == str(observation['assembly'])]
    if len(matching) != 1 or str(gen['oper_expression'][matching[0]]) != '1':
        raise ValueError('direct nonwater check requires one identity Assembly operation')
    labels = set(str(gen['asym_id_list'][matching[0]]).split(','))
    ops = block.get_mmcif_category('_pdbx_struct_oper_list.')
    i = ops['id'].index('1')
    rot = np.array([[float(ops[f'matrix[{x}][{y}]'][i]) for y in (1, 2, 3)] for x in (1, 2, 3)])
    trans = np.array([float(ops[f'vector[{x}]'][i]) for x in (1, 2, 3)])
    if not np.allclose(rot, np.eye(3), atol=1e-7) or not np.allclose(trans, 0, atol=1e-7):
        raise ValueError('operation 1 is not identity')
    atoms = list(atom_rows(block))
    seqs = {str(x['candidate_label_seq']) for x in observation['mapped_site'] if x.get('candidate_label_seq') is not None}
    site = [x for x in atoms if x['label_asym_id'] == observation['label_asym'] and str(x['label_seq_id']) in seqs]
    components = [x for x in atoms if x['label_asym_id'] in labels and x['label_asym_id'] != observation['label_asym']
                  and x['group_PDB'] == 'HETATM' and x['label_comp_id'] not in {'HOH', 'DOD'}]
    if not site:
        raise ValueError('selected site has no observed atoms')
    if not components:
        return {'component_atom_count': 0, 'within_4A': 0, 'minimum_distance_A': None}
    distances = cKDTree(np.array([x['xyz'] for x in site])).query(np.array([x['xyz'] for x in components]))[0]
    return {'component_atom_count': len(components), 'within_4A': int(np.sum(distances <= 4.0)),
            'minimum_distance_A': round(float(np.min(distances)), 4)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', required=True)
    p.add_argument('--triage-name', default='material_triage_v3')
    p.add_argument('--output-name', default='independent_material_audit.json')
    args = p.parse_args()
    run = Path(args.run)
    if Path(args.triage_name).name != args.triage_name or Path(args.output_name).name != args.output_name:
        raise ValueError('triage and output names must be simple names')
    directory = run / 'rank_comparison' / args.triage_name
    files = sorted(directory.glob('*.json'))
    if len(files) != 21:
        raise ValueError(f'expected 21 triage records, found {len(files)}')
    rows = {p.stem: load(p) for p in files}
    checks = {}
    for name, row in rows.items():
        source = row.get('candidate_source') or {}
        checks[name + '_fixed_source_read_and_hash'] = (
            row['status'] == 'SOURCE_AND_ENTITY_READ' and bool(row.get('entities')) and
            bool(source.get('path')) and digest(Path(source['path'])) == source['sha256'])
    src = load(sorted((run / 'pairs/2SRC_1FMK').glob('attempt*/evidence.json'))[-1])
    hb = gemmi.cif.read(src['holo_source']['path']).sole_block()
    ref = hb.get_mmcif_category('_struct_ref_seq.')
    hi = next(i for i, acc in enumerate(ref['pdbx_db_accession']) if acc == 'P12931' and ref['pdbx_strand_id'][i] == 'A')
    qb, qe = int(ref['seq_align_beg'][hi]), int(ref['seq_align_end'][hi])
    ub, ue = int(ref['db_align_beg'][hi]), int(ref['db_align_end'][hi])
    if qe - qb != ue - ub:
        raise ValueError('Holo source reference alignment is not an equal-length segment')
    target_uniprot = {ub + x['holo_label_seq'] - qb for x in src['mapped_site']}
    if len(target_uniprot) != 16:
        raise ValueError('unexpected Holo ATP contact-residue count')
    domain_candidates = ['1A08', '1A09', '1A1A', '1A1B', '1A1C', '1A1E', '1HCS', '1HCT']
    domain_details = []
    for pdb in domain_candidates:
        d = rows['2SRC_' + pdb]
        b = gemmi.cif.read(d['candidate_source']['path']).sole_block()
        c = b.get_mmcif_category('_struct_ref_seq.')
        ranges = [(int(c['db_align_beg'][i]), int(c['db_align_end'][i]))
                  for i, acc in enumerate(c.get('pdbx_db_accession', [])) if acc == 'P12931']
        overlap = sorted(x for x in target_uniprot if any(a <= x <= z for a, z in ranges))
        title = gemmi.cif.as_string(b.find_value('_struct.title'))
        checks['2SRC_' + pdb + '_SH2_reference_range_excludes_ATP_site'] = (
            bool(ranges) and not overlap and 'SH2' in title.upper() and
            d['entities'][0]['site_residues_covered_by_construct_sequence'] == 0)
        domain_details.append({'pdb': pdb, 'source_reference_ranges': ranges, 'overlap_uniprot_positions': overlap})
    contacts = {}
    for pdb in ['1E4V', '1E4Y', '3HPQ', '3HPR', '4X8L']:
        d = rows['1AKE_' + pdb]
        value = direct_contacts(d, 'AP5')
        contacts['1AKE_' + pdb] = value
        checks['1AKE_' + pdb + '_same_target_has_actual_compatible_site_contact'] = (
            value['unambiguous_atom_pairs_within_4A'] > 0)
    for comp in ['ADP', 'AMP']:
        value = direct_contacts(rows['1AKE_2ECK'], comp)
        contacts['1AKE_2ECK_' + comp] = value
        checks['1AKE_2ECK_' + comp + '_actual_site_contact'] = value['unambiguous_atom_pairs_within_4A'] > 0
    value = direct_contacts(rows['3ABX_3A9B'], 'BGC')
    contacts['3ABX_3A9B_BGC'] = value
    checks['3ABX_3A9B_cellobiose_member_contacts_site'] = value['unambiguous_atom_pairs_within_4A'] > 0
    d = rows['3ABX_3A9B']
    conn = gemmi.cif.read(d['candidate_source']['path']).sole_block().get_mmcif_category('_struct_conn.')
    checks['3A9B_BGC_members_have_recorded_covalent_link'] = any(
        conn['conn_type_id'][i] == 'covale' and conn['ptnr1_label_comp_id'][i] == 'BGC'
        and conn['ptnr2_label_comp_id'][i] == 'BGC' for i in range(len(conn.get('id', []))))
    value = direct_contacts(rows['3ABX_3VOF'], 'BGC')
    contacts['3ABX_3VOF_BGC'] = value
    checks['3ABX_3VOF_BGC_contacts_mapped_site'] = value['unambiguous_atom_pairs_within_4A'] > 0
    for pdb, complete in [('4X8H', 34), ('4X8M', 33)]:
        d = rows['1AKE_' + pdb]
        e = d['entities'][0]
        o = e['assembly_observations'][0]
        direct = direct_nonwater_near_site(d, o)
        contacts['1AKE_' + pdb + '_direct_nonwater'] = direct
        checks['1AKE_' + pdb + '_site_observation_and_no_nonwater_contact'] = (
            o['site_observed_residues'] == 34 and o['site_complete_relative_holo_residues'] == complete and
            direct['within_4A'] == 0 and
            not any(x['kind'] != 'water' and x['atoms_within_4A_of_mapped_site'] > 0
                    for x in o['actual_component_contacts']['nonwater']))
    d = rows['1AKE_6RZE']
    growth = d['candidate_deposited_metadata']['crystal_growth_details'] or ''
    o = d['entities'][0]['assembly_observations'][0]
    b6 = gemmi.cif.read(d['candidate_source']['path']).sole_block()
    checks['6RZE_AP5A_experimental_context_and_site_mutation'] = (
        'AP5A' in growth.upper() and not d['entities'][0]['common_uniprot'] and
        any(x['holo_author_seq'] == '119' and x.get('same_residue_component') is False for x in o['mapped_site']))
    checks['6RZE_no_positive_occupancy_AP5_coordinates'] = not any(
        a['label_comp_id'] == 'AP5' for a in atom_rows(b6))
    yo = rows['2SRC_1YOJ']['entities'][0]['assembly_observations']
    missing = {x['instance']: [y['holo_author_seq'] for y in x['mapped_site'] if y['observed_positive_heavy_atoms'] == 0]
               for x in yo}
    checks['1YOJ_chain_observations_have_distinct_missing_sites'] = missing.get('1.A') == ['404'] and missing.get('1.B') == ['276']
    yi = rows['2SRC_1YI6']['entities'][0]['assembly_observations']
    yi_direct = [direct_nonwater_near_site(rows['2SRC_1YI6'], x) for x in yi]
    contacts['2SRC_1YI6_direct_nonwater'] = yi_direct
    checks['1YI6_both_chain_sites_observed_without_nonwater_contact'] = len(yi) == 2 and all(
        x['site_observed_residues'] == 16 and x['site_complete_relative_holo_residues'] == 16 and
        not any(y['kind'] != 'water' and y['atoms_within_4A_of_mapped_site'] > 0
                for y in x['actual_component_contacts']['nonwater']) for x in yi) and all(
                    x['within_4A'] == 0 for x in yi_direct)
    byo = gemmi.cif.read(rows['2SRC_1YOJ']['candidate_source']['path']).sole_block()
    unobs = byo.get_mmcif_category('_pdbx_unobs_or_zero_occ_residues.')
    source_unobserved = {(str(unobs['auth_asym_id'][i]), str(unobs['auth_seq_id'][i]), str(unobs['occupancy_flag'][i]))
                         for i in range(len(unobs.get('id', [])))}
    checks['1YOJ_chain_specific_missing_residues_recorded_unobserved'] = (
        ('A', '406', '1') in source_unobserved and ('B', '278', '1') in source_unobserved)
    for key in checks:
        checks[key] = bool(checks[key])
    result = {'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
              'checks': checks, 'all_pass': all(checks.values()), 'target_uniprot_positions_2SRC': sorted(target_uniprot),
              'SH2_domain_source_reference_checks': domain_details, 'independent_actual_contact_checks': contacts,
              'limits': ['Structure/contact checks are independent code, not independent human scientific adjudication.',
                         'Sequence and source-range overlap do not establish a fully accepted Pair.',
                         'Both-nonblank altloc contact pairs need separate compatibility evidence.']}
    out = run / 'rank_comparison' / args.output_name
    with out.open('x', encoding='utf8') as f:
        json.dump(result, f, indent=2)
    print(json.dumps({'path': str(out), 'checks': len(checks), 'passed': sum(checks.values()),
                      'all_pass': result['all_pass']}), flush=True)
    if not result['all_pass']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
