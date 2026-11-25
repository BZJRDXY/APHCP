"""Fixed-case evidence run: compare deposited observations, not just projected ligands."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
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
import holo_source_mapping
from apo_c_pair import align_coordinates, ca, get_sequence, revision, sequence_map
from apo_c_pilot import now, read_json, sha, write_json
from holo_source_mapping import assembly_instances, connection_edges


def positive_heavy(residue, operator):
    out = []
    for atom in residue:
        if atom.occ <= 0 or atom.element.is_hydrogen:
            continue
        point = operator.transform.apply(atom.pos)
        alt = '' if atom.altloc in ('\0', ' ', '.', '?') else atom.altloc
        out.append({'name': atom.name, 'element': atom.element.name, 'occupancy': float(atom.occ),
                    'altloc': alt, 'xyz': [point.x, point.y, point.z]})
    return out


def residue_index(rows):
    result = {}
    for author, residue in rows:
        if residue.label_seq is None:
            continue
        if residue.label_seq in result:
            raise ValueError(f'multiple residues at label_seq {residue.label_seq}')
        result[residue.label_seq] = (author, residue)
    return result


def ccd_heavy_atoms(block, component):
    cat = block.get_mmcif_category('_chem_comp_atom.')
    if not cat:
        return None
    return sorted({cat['atom_id'][i] for i, name in enumerate(cat.get('comp_id', []))
                   if name == component and cat['type_symbol'][i].upper() not in {'H', 'D'}})


def target_members(block, instances, identity):
    members, all_points = [], []
    for inst, label, _, _, author, number, component in identity['target_residue_instances']:
        if inst not in instances:
            raise ValueError(f'target instance missing: {inst}')
        operator, _, rows = instances[inst]
        found = [r for chain, r in rows if chain == author and str(r.seqid) == str(number) and r.name == component]
        if len(found) != 1:
            raise ValueError(f'target member is not unique: {inst}/{author}/{number}/{component}')
        residue = found[0]
        atoms = positive_heavy(residue, operator)
        if not atoms:
            raise ValueError(f'target member has no observed heavy atoms: {inst}/{number}')
        expected = ccd_heavy_atoms(block, component)
        observed = sorted({x['name'] for x in atoms})
        is_polymer = residue.entity_type == gemmi.EntityType.Polymer
        members.append({'instance': inst, 'label_asym': label, 'author_chain': author,
                        'author_seq': str(number), 'component': component,
                        'entity_type': str(residue.entity_type), 'observed_heavy_atom_names': observed,
                        'observed_heavy_atoms': atoms,
                        'ccd_heavy_atom_names': expected, 'missing_ccd_heavy_atom_names': (
                            sorted(set(expected) - set(observed)) if expected is not None and not is_polymer else None),
                        'ccd_comparison_note': ('Individual polymer-residue CCD templates do not define whole partner completeness.'
                                                if is_polymer else 'CCD atom names as embedded in the fixed source mmCIF.')})
        all_points.extend(x['xyz'] for x in atoms)
    return members, np.array(all_points)


def connections_for_target(block, identity):
    target = {(str(x[1]), str(x[4]), str(x[5]), str(x[6])) for x in identity['target_residue_instances']}
    def key(end):
        insertion = end.get('insertion')
        suffix = '' if insertion in (None, '', '.', '?') else str(insertion)
        return (str(end.get('label_asym')), str(end.get('auth_asym')),
                str(end.get('auth_seq')) + suffix, str(end.get('ccd')))
    edges = []
    for edge in connection_edges(block):
        membership = [key(end) in target for end in edge['ends']]
        if any(membership):
            edges.append(edge | {'target_endpoint_membership': membership})
    return edges


def modifications(block, structure, label, operator, site_atoms):
    cat = block.get_mmcif_category('_pdbx_struct_mod_residue.')
    records = []
    if cat:
        for i in range(len(next(iter(cat.values())))):
            row = {k: v[i] for k, v in cat.items()}
            if str(row.get('label_asym_id') or '') == label:
                records.append(row)
    observed = []
    tree = cKDTree(np.array([a['xyz'] for a in site_atoms])) if site_atoms else None
    for chain in structure[0]:
        for res in chain:
            if res.subchain == label and res.name in {'PTR', 'TPO', 'SEP', 'MSE', 'CSO'}:
                atoms = positive_heavy(res, operator)
                item = {'author_chain': chain.name, 'author_seq': str(res.seqid),
                        'label_seq': res.label_seq, 'component': res.name,
                        'positive_heavy_atoms': len(atoms)}
                if tree is not None and atoms:
                    distances, indices = tree.query(np.array([a['xyz'] for a in atoms]))
                    i = int(np.argmin(distances))
                    item.update({'min_distance_to_actual_site_atoms_A': round(float(distances[i]), 4),
                                 'closest_modified_atom': atoms[i],
                                 'closest_site_atom': site_atoms[int(indices[i])]})
                observed.append(item)
    return {'declared': records, 'selected_observed': observed}


def unobserved_records(block, category, label, model):
    cat = block.get_mmcif_category(category)
    out = {}
    if not cat:
        return out
    for i in range(len(next(iter(cat.values())))):
        row = {k: v[i] for k, v in cat.items()}
        if str(row.get('label_asym_id') or '') != label:
            continue
        if row.get('PDB_model_num') not in (None, '.', '?', str(model)):
            continue
        try:
            seq = int(row['label_seq_id'])
        except (KeyError, TypeError, ValueError):
            continue
        flag = str(row.get('occupancy_flag'))
        row['deposited_observation_status'] = {'0': 'ZERO_OCCUPANCY', '1': 'UNOBSERVED'}.get(flag, 'UNSPECIFIED')
        out.setdefault(seq, []).append(row)
    return out


def atom_alt_compatibility(a, b):
    if a['altloc'] and b['altloc']:
        return 'cross_residue_altloc_compatibility_unverified'
    return 'one_or_both_blank_altloc'


def deposited_metadata(block):
    def value(key):
        raw = block.find_value(key)
        if raw is None or raw in ('.', '?'):
            return None
        return gemmi.cif.as_string(raw)
    return {'title': value('_struct.title'), 'experimental_method': value('_exptl.method'),
            'refinement_starting_model': value('_refine.pdbx_starting_model'),
            'crystal_growth_details': value('_exptl_crystal_grow.pdbx_details')}


def component_contacts(instances, selected, site_atoms, projection_points):
    if not site_atoms:
        return {'status': 'NO_MAPPED_OBSERVED_SITE_ATOMS', 'nonwater': [], 'nearby_waters': []}
    coords = np.array([x['xyz'] for x in site_atoms])
    tree = cKDTree(coords)
    target_tree = cKDTree(projection_points) if projection_points is not None else None
    nonwater, water = [], []
    for instance, (operator, label, residues) in instances.items():
        if instance in selected:
            continue
        for author, residue in residues:
            atoms = positive_heavy(residue, operator)
            if not atoms:
                continue
            query = np.array([a['xyz'] for a in atoms])
            distances, indices = tree.query(query)
            idx = int(np.argmin(distances))
            nearest = site_atoms[int(indices[idx])]
            minimum = float(distances[idx])
            kind = ('water' if residue.is_water() else 'polymer' if residue.entity_type == gemmi.EntityType.Polymer
                    else 'single_atom' if len(atoms) == 1 else 'nonpolymer')
            row = {'instance': instance, 'label_asym': label, 'author_chain': author,
                   'author_seq': str(residue.seqid), 'component': residue.name, 'kind': kind,
                   'min_distance_to_mapped_site_atoms_A': round(minimum, 4),
                   'atoms_within_4A_of_mapped_site': int(np.sum(distances <= 4.0)),
                   'closest_component_atom': atoms[idx],
                   'closest_site_atom': nearest,
                   'altloc_compatibility': atom_alt_compatibility(atoms[idx], nearest)}
            if target_tree is not None:
                td = target_tree.query(query)[0]
                row['min_distance_to_globally_projected_observed_holo_target_A'] = round(float(np.min(td)), 4)
            if kind == 'water':
                if minimum <= 4.0:
                    water.append(row)
            elif kind != 'polymer' or minimum <= 6.0:
                nonwater.append(row)
    nonwater.sort(key=lambda x: x['min_distance_to_mapped_site_atoms_A'])
    water.sort(key=lambda x: x['min_distance_to_mapped_site_atoms_A'])
    return {'status': 'MEASURED', 'nonwater': nonwater, 'nearby_waters': water,
            'meaning': 'Distances are to actual mapped candidate-site atoms; proximity alone does not prove competitive occupancy.'}


def initialize(args):
    config = read_json(Path(args.config))
    baseline = Path(config['baseline_run'])
    samples = read_json(baseline / 'samples.json')
    run = Path(args.run)
    if run.exists():
        raise FileExistsError(f'refusing to overwrite Run {run}')
    entries = []
    for spec in config['cases']:
        cases = [x for x in samples if x['id'] == spec['system_id']]
        if len(cases) != 1:
            raise ValueError(f'baseline sample is not unique: {spec["system_id"]}')
        bp = baseline / spec['baseline']
        data = read_json(bp)
        if data.get('execution_status') != 'INSPECTED':
            raise ValueError(f'baseline evidence not inspected: {bp}')
        holo_path = Path(cases[0]['source_uri'])
        candidate_path = Path(data['candidate_source']['uri'])
        if sha(holo_path) != cases[0]['identity']['source_sha256']:
            raise ValueError(f'Holo checksum differs: {holo_path}')
        if sha(candidate_path) != data['candidate_source']['sha256']:
            raise ValueError(f'candidate checksum differs: {candidate_path}')
        entries.append({'id': spec['id'], 'sample': cases[0], 'baseline_path': str(bp),
                        'baseline_sha256': sha(bp), 'baseline': data,
                        'source_evidence': spec['source_evidence']})
    run.mkdir(parents=True)
    (run / 'code').mkdir()
    for path in [Path(__file__), Path(apo_c_pair.__file__), Path(apo_c_pilot.__file__), Path(holo_source_mapping.__file__)]:
        shutil.copy2(path, run / 'code' / f'{path.stem}-{sha(path)[:12]}.py')
    write_json(run / 'config.json', config)
    write_json(run / 'inputs.json', entries)
    write_json(run / 'run.json', {'started_utc': now(), 'purpose': 'Fixed observations and evidence-range regression',
                                 'baseline_run': str(baseline), 'config_sha256': sha(Path(args.config)),
                                 'code_sha256': sha(Path(__file__)), 'cases': len(entries),
                                 'formal_release': False})
    print(json.dumps({'run': str(run), 'cases': len(entries)}), flush=True)


def inspect_single(entry, out):
    sample, old = entry['sample'], entry['baseline']
    if len(sample['identity']['binding_instances']) != 1:
        raise ValueError('multi-chain case needs its own explicit copy-equivalence analysis')
    hpath = Path(sample['source_uri'])
    cpath = Path(old['candidate_source']['uri'])
    if sha(hpath) != sample['identity']['source_sha256'] or sha(cpath) != old['candidate_source']['sha256']:
        raise ValueError('fixed original changed since initialization')
    hb, cb = gemmi.cif.read(str(hpath)).sole_block(), gemmi.cif.read(str(cpath)).sole_block()
    hpdb, cpdb = sample['pdb'], old['candidate_pdb']
    hrev, hdate = revision(hb, hpdb)
    crev, cdate = revision(cb, cpdb)
    if hrev != sample['identity']['revision'] or crev != old['candidate_source']['revision']:
        raise ValueError('source revision differs from fixed input')
    hs, cs = gemmi.make_structure_from_block(hb), gemmi.make_structure_from_block(cb)
    hi = assembly_instances(hs, str(sample['identity']['assembly']), hs[0], include_water=True)
    ci = assembly_instances(cs, str(old['candidate_assembly']), cs[0], include_water=True)
    hkey, ckey = sample['identity']['binding_instances'][0][0], old['candidate_binding_instance']
    hop, hlabel, hrows = hi[hkey]
    cop, clabel, crows = ci[ckey]
    hres, cres = residue_index(hrows), residue_index(crows)
    members, ligand_points = target_members(hb, hi, sample['identity'])
    ligand_tree = cKDTree(ligand_points)
    actual_holo_site, holo_site_atoms = set(), []
    for n, (author, res) in hres.items():
        atoms = positive_heavy(res, hop)
        if atoms and np.min(ligand_tree.query(np.array([a['xyz'] for a in atoms]))[0]) <= 4.0:
            actual_holo_site.add(n)
            holo_site_atoms.extend(a | {'label_seq': n, 'author_chain': author,
                                        'author_seq': str(res.seqid), 'residue': res.name} for a in atoms)
    baseline_site = {x['holo_label_seq'] for x in old['site_residues']}
    site_difference = {'added_holo_label_seq': sorted(actual_holo_site - baseline_site),
                       'removed_holo_label_seq': sorted(baseline_site - actual_holo_site)}
    if actual_holo_site != baseline_site and not entry.get('research_scope_extension'):
        raise ValueError(f'baseline Holo site differs from direct source: {sorted(actual_holo_site ^ baseline_site)}')

    mapping, sequence_evidence = sequence_map(get_sequence(hs, hrows), get_sequence(cs, crows))
    candidate_unobserved = unobserved_records(cb, '_pdbx_unobs_or_zero_occ_residues.', clabel, cs[0].num)
    candidate_unobserved_atoms = unobserved_records(cb, '_pdbx_unobs_or_zero_occ_atoms.', clabel, cs[0].num)
    hx, cy = [], []
    for hn, cn in mapping.items():
        if hn in hres and cn in cres:
            a, b = ca(hres[hn][1], hop), ca(cres[cn][1], cop)
            if a is not None and b is not None:
                hx.append(a)
                cy.append(b)
    projection, global_rmsd = None, None
    if len(hx) >= 3:
        rotation, translation, global_rmsd = align_coordinates(np.array(hx), np.array(cy))
        projection = ligand_points @ rotation + translation
    rows, site_atoms = [], []
    for hn in sorted(actual_holo_site):
        hauth, hr = hres[hn]
        cn = mapping.get(hn)
        row = {'holo_label_seq': hn, 'holo_author_chain': hauth, 'holo_author_seq': str(hr.seqid),
               'holo_component': hr.name, 'candidate_label_seq': cn,
               'candidate_deposited_unobserved_residue_records': candidate_unobserved.get(cn, []),
               'candidate_deposited_unobserved_atom_records': candidate_unobserved_atoms.get(cn, [])}
        if cn in cres:
            cauth, cr = cres[cn]
            hatoms, catoms = positive_heavy(hr, hop), positive_heavy(cr, cop)
            hnms, cnms = {a['name'] for a in hatoms}, {a['name'] for a in catoms}
            row.update({'candidate_author_chain': cauth, 'candidate_author_seq': str(cr.seqid),
                        'candidate_component': cr.name, 'same_residue_component': cr.name == hr.name,
                        'candidate_observed_heavy_atom_names': sorted(cnms),
                        'candidate_missing_relative_holo_heavy_atom_names': sorted(hnms - cnms),
                        'candidate_backbone_atoms': sorted(cnms & {'N', 'CA', 'C', 'O'}),
                        'candidate_observed': bool(catoms)})
            for atom in catoms:
                site_atoms.append(atom | {'label_seq': cn, 'author_chain': cauth,
                                          'author_seq': str(cr.seqid), 'residue': cr.name})
        else:
            row['candidate_observed'] = False
        rows.append(row)
    component_evidence = component_contacts(ci, {ckey}, site_atoms, projection)
    target_instances = {x[0] for x in sample['identity']['target_residue_instances']}
    holo_context = component_contacts(hi, {hkey} | target_instances, holo_site_atoms, ligand_points)
    baseline_dir = Path(entry['baseline_path']).parent
    for name in ['holo_system.cif', 'candidate_assembly.cif']:
        src, dst = baseline_dir / name, out / name
        if name == 'holo_system.cif' and entry.get('research_scope_extension'):
            apo_c_pair.write_assembly(hs, str(sample['identity']['assembly']), dst)
        else:
            if not src.is_file():
                raise ValueError(f'baseline derivative missing: {src}')
            shutil.copy2(src, dst)
        st = gemmi.read_structure(str(dst))
        if not st or not st[0].count_atom_sites():
            raise ValueError(f'derivative cannot be read: {dst}')
    derivatives = {name: {'path': str(out / name), 'sha256': sha(out / name),
                          'readback_atoms': gemmi.read_structure(str(out / name))[0].count_atom_sites()}
                   for name in ['holo_system.cif', 'candidate_assembly.cif']}
    edges = connections_for_target(hb, sample['identity'])
    external_covalent = [edge for edge in edges if edge['kind'] == 'covale'
                         and sum(edge['target_endpoint_membership']) == 1]
    return {'case': entry['id'], 'system_id': sample['id'], 'holo_status': sample['holo_status'],
            'research_scope_extension': entry.get('research_scope_extension'),
            'holo_identity': sample['identity'], 'holo_qc': sample['qc'],
            'holo_source': {'path': str(hpath), 'sha256': sha(hpath), 'revision': hrev,
                            'first_public_model_date': hdate, 'model': hs[0].num},
            'candidate_source': {'path': str(cpath), 'sha256': sha(cpath), 'revision': crev,
                                 'first_public_model_date': cdate, 'model': cs[0].num,
                                 'assembly': old['candidate_assembly'], 'instance': ckey,
                                 'author_chain': old['candidate_author_chain'], 'label_asym': clabel},
            'deposited_metadata': {'holo': deposited_metadata(hb), 'candidate': deposited_metadata(cb)},
            'holo_target_members': members, 'holo_target_connections': edges,
            'holo_target_scope_status': ('EXTERNAL_COVALENT_CONNECTION_REQUIRES_REVIEW' if external_covalent
                                         else 'NO_EXTERNAL_COVALENT_CONNECTION_RECORDED'),
            'external_covalent_connections': external_covalent,
            'holo_contact_definition': 'Observed positive-occupancy target and receptor heavy atoms within 4 A.',
            'holo_contact_residue_count': len(actual_holo_site),
            'source_site_matches_baseline': actual_holo_site == baseline_site,
            'baseline_site_difference': site_difference,
            'mapped_site': rows, 'mapped_site_observed_residues': sum(x['candidate_observed'] for x in rows),
            'mapped_site_complete_relative_holo': all(x['candidate_observed'] and x.get('same_residue_component') and
                 not x.get('candidate_missing_relative_holo_heavy_atom_names') for x in rows),
            'sequence_alignment': sequence_evidence, 'sifts_common_uniprot': old.get('sifts_common_uniprot'),
            'holo_modifications': modifications(hb, hs, hlabel, hop, holo_site_atoms),
            'candidate_modifications': modifications(cb, cs, clabel, cop, site_atoms),
            'holo_other_components': holo_context,
            'actual_candidate_components': component_evidence,
            'global_alignment_ca_count': len(hx), 'global_ca_rmsd_A': round(global_rmsd, 4) if global_rmsd is not None else None,
            'projected_target_is_diagnostic_only': True,
            'baseline_geometric_screen': old.get('geometric_screen'), 'baseline_site_observation': old.get('site_observation_status'),
            'baseline_evidence': entry['baseline_path'], 'source_evidence': entry['source_evidence'],
            'formal_pair_status': 'PENDING_INDEPENDENT_EVIDENCE_REVIEW', 'derivatives': derivatives,
            'limitations': ['Nearby component/site-atom distances are not by themselves a competition or occupancy decision.',
                            'A superposed Holo target is a probe only; large motion requires alternate evidence.',
                            'Sequence alignment is not a residue-level SIFTS file.',
                            'Nonblank altlocs on separate residues are not assumed mutually compatible.',
                            'Experimental unbound evidence and Holo partner completeness require separate adjudication.']}


def run_case(args):
    run = Path(args.run)
    input_path = Path(args.input_file) if args.input_file else run / 'inputs.json'
    if not input_path.resolve().is_relative_to(run.resolve()):
        raise ValueError('alternate input file must be inside the research Run')
    matches = [x for x in read_json(input_path) if x['id'] == args.case]
    if len(matches) != 1:
        raise ValueError('case must match one frozen input')
    entry = matches[0]
    if sha(Path(entry['baseline_path'])) != entry['baseline_sha256']:
        raise ValueError('baseline evidence changed since initialization')
    base = run / 'pairs' / entry['id']
    base.mkdir(parents=True, exist_ok=True)
    n = 1
    while (base / f'attempt{n:03d}').exists():
        n += 1
    out = base / f'attempt{n:03d}'
    out.mkdir()
    code_hash = sha(Path(__file__))
    snapshot = run / 'code' / f'{Path(__file__).stem}-{code_hash[:12]}.py'
    if not snapshot.exists():
        shutil.copy2(Path(__file__), snapshot)
    start, clock = now(), time.monotonic()
    result = {'started_utc': start, 'case': entry['id'], 'code_sha256': code_hash}
    try:
        result.update(inspect_single(entry, out))
        result['execution_status'] = 'MEASURED'
    except Exception as exc:
        result.update({'execution_status': 'UNRESOLVED', 'error': f'{type(exc).__name__}: {exc}'})
    result['elapsed_s'] = round(time.monotonic() - clock, 3)
    result['finished_utc'] = now()
    write_json(out / 'evidence.json', result)
    print(json.dumps({'case': entry['id'], 'status': result['execution_status'], 'path': str(out),
                      'error': result.get('error'), 'elapsed_s': result['elapsed_s']}), flush=True)
    if result['execution_status'] != 'MEASURED':
        sys.exit(1)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('stage', choices=['init', 'case'])
    p.add_argument('--run', required=True)
    p.add_argument('--config')
    p.add_argument('--case')
    p.add_argument('--input-file', help='explicit, immutable research scope-extension input inside this Run')
    args = p.parse_args()
    if args.stage == 'init':
        if not args.config:
            p.error('init requires --config')
        initialize(args)
    else:
        if not args.case:
            p.error('case requires --case')
        run_case(args)


if __name__ == '__main__':
    main()
