"""Revalidate the pinned Holo pilot without obtaining new source data."""

import argparse
import csv
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import gemmi
from plinder.core import PlinderSystem
from plinder.core.release import PlinderRelease

from holo_candidate_batch import atomic_json, parse_system_id, revision_from_cif, sha256, verify_system_source
from holo_native_candidates import assemble_candidate_ledger
from holo_source_mapping import MappingError


def table(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def json_from(path):
    return json.loads(path.read_text(encoding='utf-8'))


def export(path, rows):
    if not rows:
        return
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def brief(value):
    return json.dumps(value, ensure_ascii=False, separators=(',', ':'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('roles', type=Path)
    parser.add_argument('plinder_release', type=Path)
    parser.add_argument('--entries', nargs='+')
    args = parser.parse_args()
    old, out = args.baseline.resolve(), args.output.resolve()
    out.mkdir(exist_ok=False, parents=True)
    started = datetime.now(timezone.utc).isoformat()
    rules = json_from(args.roles)['entries']
    release = PlinderRelease(data_dir=args.plinder_release)
    candidates = table(old / 'inputs/candidates_display.csv')
    source_manifest = {r['pdb_id']: r for r in table(old / 'inputs/sources.csv')}
    baseline_entries = {r['pdb_id']: r for r in table(old / 'receipts_final/entry_results.csv')}
    baseline_systems = defaultdict(list)
    for row in table(old / 'receipts_final/system_results.csv'):
        baseline_systems[row['pdb_id']].append(row)
    prior_rechecks = {r['system_id']: r for r in json_from(old / 'recheck_final.json')['systems']}
    if (len(candidates) != len(source_manifest) or len(candidates) != len(baseline_entries)
        or sum(map(len, baseline_systems.values())) != 291):
        raise SystemExit('Unexpected prior denominators')
    entry_rows, system_rows, unresolved29, native_rows, observations = [], [], [], [], []
    for candidate in candidates:
        pdb = candidate['pdb_id']
        if args.entries and pdb not in args.entries:
            continue
        entry = baseline_entries[pdb]
        raw_root = (old / 'results' / entry['entry_receipt']).resolve().parent.parent
        old_receipt = json_from((old / 'results' / entry['entry_receipt']).resolve())
        source = raw_root / old_receipt['source_relative_path']
        assert sha256(source) == entry['original_sha256']
        block, revision = revision_from_cif(source, pdb)
        assert revision == source_manifest[pdb]['revision']
        structure = gemmi.make_structure_from_block(block)
        receipt_by_system = {r['system_id']: r for r in old_receipt['system_results']}
        systems, ledger_inputs = [], []

        for old_system in baseline_systems[pdb]:
            system_id = old_system['system_id']
            path = (old / 'results' / old_system['derivative_path']).resolve()
            previous_pass = (old_system['coordinate_comparison'] == 'coordinate_identity_pass' or
                             prior_rechecks.get(system_id, {}).get('status') == 'coordinate_identity_pass')
            system = {'system_id': system_id, 'source_path': str(source), 'derivative': str(path),
                      'old_coordinate_status': old_system['coordinate_comparison'], 'old_failure': old_system['problem'],
                      'old_recheck': prior_rechecks.get(system_id), 'previous_pass': previous_pass,
                      'derivative_checksum_ok': False, 'readback': 'not_done', 'formal_qc': 'not_run'}
            assembly, _, _ = parse_system_id(system_id)
            annotations = receipt_by_system[system_id]['annotation_rows']
            result = None
            try:
                if sha256(path) != old_system['derivative_sha256']:
                    raise MappingError('derivative_checksum_mismatch', 'Baseline derivative hash changed')
                system['derivative_checksum_ok'] = True
                derivative = gemmi.cif.read(str(path)).sole_block()
                dstructure = gemmi.make_structure_from_block(derivative)
                system['readback'] = 'pass' if len(dstructure) and dstructure[0].count_atom_sites() else 'failed'
                result = verify_system_source(structure, block, derivative, system_id, annotations)
                system.update(result)
            except MappingError as exc:
                category = ('unsupported_branch' if 'unsupported' in exc.kind else 'source_observation_unresolved')
                system.update(status=category, reason=exc.kind, message=str(exc), diagnostic=exc.evidence)
            except Exception as exc:
                system.update(status='technical_failure', reason=type(exc).__name__, message=str(exc))

            ledger_inputs.append((system_id, assembly, annotations, result))
            if system['status'] == 'confirmed_derivative_difference':
                try:
                    rebuild_dir = out / 'reconstructed' / pdb / system_id
                    plinder = PlinderSystem(system_id=system_id, release=release, source_mmcif=source,
                                            reconstruction_dir=rebuild_dir)
                    rebuilt = Path(plinder.system_cif)
                    replay = verify_system_source(structure, block, gemmi.cif.read(str(rebuilt)).sole_block(),
                                                  system_id, annotations)
                    system['rebuild'] = {'path': str(rebuilt.relative_to(out)), 'sha256': sha256(rebuilt),
                                         'status': replay['status'],
                                         'same_differences': replay['differences'] == system['differences']}
                except Exception as exc:
                    system['rebuild'] = {'status': 'technical_failure', 'message': str(exc)}
            systems.append(system)

        absent = entry['plinder_annotation'] == 'absent'
        ledger = assemble_candidate_ledger(block, structure, pdb, revision, ledger_inputs, rules.get(pdb, {}), direct=absent)
        needs_native = any(s['status'] != 'coordinate_identity_pass' for s in systems)
        if absent or needs_native:
            for observation in ledger:
                key = observation['candidate_key']
                if absent:
                    native_rows.append({'pdb_id': pdb, 'revision': revision, 'model': observation['model'],
                        'assembly': observation['assembly'], 'candidate_key': key,
                        'members': '; '.join(f"{m['instance']}@{m['operator']}:{m['entity']}:" +
                            (m['description'] if m['kind'] == 'polymer' else '+'.join(
                                f"{r['ccd']}({r['auth_asym']}:{r['auth_seq']})" for r in m['residues']))
                            for m in observation['members']),
                        'role': observation['role'], 'role_detail': observation['role_detail'],
                        'role_evidence': observation['role_evidence'], 'identity_status': observation['identity_status'],
                        'protein_candidates_A': brief(observation['protein_instance_candidates']),
                        'source_connection_ids': ';'.join(str(link['source']['id']) for link in observation['connections']),
                        'unresolved_connections': ';'.join(observation['unresolved_connection_ids']),
                        'formal_qc': 'not_run', 'entry_receipt': f'entries/{pdb.lower()}.json'})
        role_counts = Counter(row['role'] for row in ledger if 'source_mmcif' in row['routes'])
        system_counts = Counter(row['status'] for row in systems)
        root_receipt = f'entries/{pdb.lower()}.json'
        atomic_json(out / root_receipt, {'pdb_id': pdb, 'revision': revision, 'source': str(source),
                    'source_sha256': entry['original_sha256'], 'source_models': [m.num for m in structure],
                    'source_assemblies': [a.name for a in structure.assemblies],
                    'index_annotation': entry['plinder_annotation'], 'systems': systems,
                    'candidate_ledger': ledger, 'candidate_count': len(ledger)})
        observations.extend(ledger)
        entry_rows.append({'pdb_id': pdb, 'revision': revision, 'source_parse': 'pass', 'source_models': len(structure),
            'index_annotation': entry['plinder_annotation'], 'systems': len(systems),
            'readback_pass': sum(s['readback'] == 'pass' for s in systems),
            'coordinate_pass': system_counts['coordinate_identity_pass'],
            'confirmed_difference': system_counts['confirmed_derivative_difference'],
            'other_systems': len(systems) - system_counts['coordinate_identity_pass'] - system_counts['confirmed_derivative_difference'],
            'source_branch_status': ('observations_recorded' if role_counts else 'no_observed_ligand_candidate') if absent or needs_native else 'not_needed',
            'target_ligand': role_counts['target_ligand'], 'non_target_component': role_counts['non_target_component'],
            'pending_role': role_counts['pending'], 'formal_qc': 'not_run', 'receipt': root_receipt})
        for s in systems:
            row = {'pdb_id': pdb, 'system_id': s['system_id'], 'old_coordinate_status': s['old_coordinate_status'],
                   'old_combined_pass': s['previous_pass'], 'old_failure': s['old_failure'], 'readback': s['readback'],
                   'source_model': s.get('source_model'), 'models_in_source': len(structure),
                   'new_status': s['status'], 'conclusion': ('mapping_fixed' if not s['previous_pass'] and s['status'] == 'coordinate_identity_pass'
                        else ('retained_pass' if s['previous_pass'] and s['status'] == 'coordinate_identity_pass' else s['status'])),
                   'reason': s.get('reason') or ';'.join(d['kind'] for d in s.get('differences', [])),
                   'evidence': brief(s.get('differences') or s.get('diagnostic') or s.get('message') or {
                       'source_operators': s.get('source_operators'), 'max_distance_A': s.get('max_distance_A'),
                       'tolerance_A': s.get('precision', {}).get('distance_tolerance_A'),
                       'model_evidence': s.get('model_evidence')}),
                   'rebuilt': s.get('rebuild', {}).get('status'), 'formal_qc': 'not_run', 'receipt': root_receipt}
            system_rows.append(row)
            if not s['previous_pass']:
                unresolved29.append(row)
        print(brief({'entry': pdb, 'systems': len(systems), 'statuses': system_counts, 'source_branch': role_counts}), flush=True)

    if not args.entries and (len(entry_rows), len(system_rows), len(unresolved29)) != (100, 291, 29):
        raise SystemExit('Full run did not cover the exact baseline denominators')
    export(out / 'entry_results.csv', entry_rows)
    export(out / 'system_results.csv', system_rows)
    export(out / 'previous_29_results.csv', unresolved29)
    export(out / 'source_no_annotation_observations.csv', native_rows)
    atomic_json(out / 'candidate_ledger.json', observations)
    summary = {'started_utc': started, 'finished_utc': datetime.now(timezone.utc).isoformat(),
        'entries': len(entry_rows), 'systems': len(system_rows), 'readback_pass': sum(r['readback'] == 'pass' for r in system_rows),
        'old_passes_retained': sum(r['conclusion'] == 'retained_pass' for r in system_rows),
        'old_passes_changed': [r['system_id'] for r in system_rows if r['old_combined_pass'] and r['conclusion'] != 'retained_pass'],
        'systems_by_status': dict(Counter(r['new_status'] for r in system_rows)),
        'previous_29_by_conclusion': dict(Counter(r['conclusion'] for r in unresolved29)),
        'no_annotation_entries': sum(r['index_annotation'] == 'absent' for r in entry_rows),
        'native_observations_by_role': dict(Counter(r['role'] for r in native_rows)),
        'roles_file_sha256': sha256(args.roles), 'baseline': str(old), 'plinder_release': str(args.plinder_release),
        'versions': {name: version(name) for name in ('gemmi', 'plinder', 'biotite', 'numpy', 'scipy')},
        'input_sha256': {str(path.name): sha256(path) for path in (
            old / 'inputs/candidates_display.csv', old / 'inputs/sources.csv',
            old / 'receipts_final/entry_results.csv', old / 'receipts_final/system_results.csv',
            old / 'recheck_final.json', Path(__file__), Path(__file__).with_name('holo_candidate_batch.py'),
            Path(__file__).with_name('holo_source_mapping.py'),
            Path(__file__).with_name('holo_native_candidates.py'))},
        'new_downloads': 0, 'formal_qc': 'not_run'}
    atomic_json(out / 'summary.json', summary)
    print(brief(summary), flush=True)


if __name__ == '__main__':
    main()
