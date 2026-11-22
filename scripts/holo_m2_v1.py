"""Source-only V1 delivery mode for the existing M2 pilot runner."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import gzip
import re
import time

import gemmi
from scipy.spatial import cKDTree

from holo_candidate_batch import atomic_json, now, revision_from_cif, sha256
from holo_native_candidates import inspect_source, key_for
from holo_source_mapping import MappingError
from holo_v1_export import (COORD, LABEL, MODEL, choose_connected, describe_instances, export, output_label, rows, text)


def read(path):
    from json import loads
    return loads(path.read_text(encoding='utf-8-sig'))


def sequence_evidence(block, partner_entities):
    seqs = {text(r['entity_id']): ''.join(text(r['pdbx_seq_one_letter_code_can']).split())
            for r in rows(block, '_entity_poly.')}
    refs = rows(block, '_struct_ref.')
    result = []
    for entity in sorted(partner_entities & seqs.keys()):
        own_refs = [{k: text(r[k]) for k in ('db_name', 'pdbx_db_accession')} for r in refs if text(r['entity_id']) == entity]
        seq = seqs[entity]
        matches = []
        if seq:
            for r in refs:
                reference = ''.join(text(r.get('pdbx_seq_one_letter_code')).split())
                start = text(r.get('pdbx_align_begin'))
                if text(r['entity_id']) == entity or not start.isdecimal():
                    continue
                for m in re.finditer(f'(?={re.escape(seq)})', reference):
                    matches.append({'accession': text(r['pdbx_db_accession']), 'reference_entity': text(r['entity_id']),
                                    'start': int(start)+m.start(), 'end': int(start)+m.start()+len(seq)-1})
        result.append({'entity': entity, 'canonical_sequence': seq, 'deposited_references': own_refs,
                       'exact_subsequence_matches_in_other_deposited_references': matches,
                       'interpretation': 'Sequence correspondence alone does not establish the origin or biological role of this instance.'})
    return result


def projected_evidence(spec, artifact, atoms, instance_map):
    labels = {k: output_label(v) for k, v in instance_map.items()}
    partners = {labels[k] for k in spec['partners']}
    heavy = Counter(text(r['occupancy']) for r in atoms if text(r[LABEL]) in partners
                    and text(r['type_symbol']).upper() not in ('H', 'D', 'T'))
    partner_ccds = sorted({text(r['label_comp_id']) for r in atoms if text(r[LABEL]) in partners})
    residues = {k: len({(text(r['label_seq_id']), text(r['auth_seq_id']), text(r['pdbx_PDB_ins_code']), text(r['label_comp_id']))
                       for r in atoms if text(r[LABEL]) == labels[k]}) for k in spec['partners']}
    protein_labels = {labels[k] for k in spec['binding_proteins']}
    protein_alt = Counter(text(r['label_alt_id']) for r in atoms if text(r[LABEL]) in protein_labels and text(r['label_alt_id']))
    expected = spec.get('checks', {})
    for name, observed in [('partner_heavy_atom_occupancies', dict(heavy)), ('partner_ccds', partner_ccds), ('partner_residues', residues)]:
        if name in expected and expected[name] != observed:
            raise ValueError(f'Expected {name} {expected[name]!r}; observed {observed!r}')
    for check in expected.get('connections', []):
        matches = [r['output'] for r in artifact['connections'] if r['output']['id'] == check['id']]
        if len(matches) != 1 or not all(all(a.get(k) == v for k, v in b.items()) for a, b in zip(matches[0]['ends'], check['ends'])):
            raise ValueError(f'Expected connection or endpoint differs on disk: {check!r}')
    answer = {'partner_ccds': partner_ccds, 'partner_heavy_atom_occupancies': dict(heavy), 'partner_residues': residues,
              'binding_protein_named_altloc_atoms': dict(protein_alt)}
    if 'partner_heavy_atoms_by_altloc' in expected:
        alternatives = {}
        for key in spec['partners']:
            pose_counts = {}
            for r in atoms:
                if text(r[LABEL]) == labels[key] and text(r['type_symbol']).upper() not in ('H', 'D', 'T'):
                    pose_counts.setdefault(text(r['label_alt_id']), Counter())[text(r['occupancy'])] += 1
            alternatives[key] = {alt: dict(counts) for alt, counts in pose_counts.items()}
        if alternatives != expected['partner_heavy_atoms_by_altloc']:
            raise ValueError(f'Disk partner altloc/occupancy differs: {alternatives!r}')
        answer['partner_heavy_atoms_by_altloc'] = alternatives
    return answer


def disk_contact_distances(atoms, instances, partners, proteins):
    """Check contact scope from exported atoms, with the native witness's 4 A rule."""
    xyz = {k: [] for k in partners + proteins}
    labels = {output_label(instances[k]): k for k in xyz}
    for row in atoms:
        key = labels.get(text(row[LABEL]))
        if key and text(row['type_symbol']).upper() not in ('H', 'D', 'T') and float(text(row['occupancy'])) > 0:
            xyz[key].append([float(text(row[c])) for c in COORD])
    target = [p for key in partners for p in xyz[key]]
    if not target or any(not xyz[key] for key in proteins):
        raise ValueError('Contact audit is missing positive-occupancy heavy atoms on a declared member')
    tree = cKDTree(target)
    return {key: round(float(tree.query(xyz[key], k=1)[0].min()), 4) for key in proteins}


def load_source_proof(spec, source_run):
    entry = spec['pdb_id']
    if not spec.get('historical_system_manifest'):
        path = source_run/'source_proofs'/f'{entry}.json'
        return path, read(path)
    path = (source_run/spec['historical_system_manifest']).resolve()
    if not path.is_relative_to(source_run):
        raise ValueError('Historical source manifest must stay in the approved source Run')
    records = [r for r in read(path) if r['pdb_id'] == entry and r['revision'] == spec['revision']]
    if len(records) != 1:
        raise ValueError('Historical manifest must prove exactly one entry/revision')
    record = records[0]
    return path, {'identity': {k: record[k] for k in ('pdb_id', 'revision', 'source_url')},
                  'relative_path': record['source'], 'sha256': record['source_sha256'], 'bytes': record['source_bytes']}


def one(spec, source_run, run, fingerprint, verify):
    entry = spec['pdb_id']
    receipt = run/'qc'/f'{entry}.json'
    previous = read(receipt) if receipt.is_file() else {}
    proof_path, proof = load_source_proof(spec, source_run)
    ident = proof['identity']
    relative = Path(proof['relative_path'])
    if relative.parts[:1] != ('raw',) or '..' in relative.parts:
        raise ValueError('Proof must reference an original from the source Run raw directory')
    reference = source_run/relative
    source = reference.resolve()
    in_source_run = source.is_relative_to(source_run/'raw')
    # The four-case historical Run reused raw files through links to earlier isolated Runs.
    in_historical_raw = (spec.get('historical_system_manifest') and source.is_relative_to(source_run.parent)
                         and source.parent.name == 'raw')
    if not (in_source_run or in_historical_raw):
        raise ValueError('Approved source link must resolve to a raw file within the Run store')
    actual = sha256(source)
    if (ident['pdb_id'] != entry or ident['revision'] != spec['revision']
            or proof['sha256'] != spec['sha256'] or actual != spec['sha256']
            or source.stat().st_size != proof['bytes']):
        raise ValueError('Source differs from the approved fixed-revision proof')
    if proof.get('gzip_crc') != 'pass':
        with gzip.open(source, 'rb') as handle:
            while handle.read(1 << 20):
                pass
    if (not verify and previous.get('fingerprint') == fingerprint and previous.get('structure_acceptance') == 'PASS'
            and all(Path(f['path']).is_file() and sha256(Path(f['path'])) == f['sha256'] for f in previous['artifacts'].values())):
        return previous, 'reused'
    source_block, revision = revision_from_cif(source, entry)
    if revision != spec['revision']:
        raise ValueError(f'Structure model revision is {revision}, expected {spec["revision"]}')
    structure = gemmi.make_structure_from_block(source_block)
    models = [m.num for m in structure]
    aid = spec['assembly']
    representative = source_block.find_value('_pdbx_nmr_representative.conformer_id')
    model_num = spec.get('representative_model') or (int(text(representative)) if text(representative) else models[0])
    model = next(m for m in structure if m.num == model_num)
    instance_map, expand = describe_instances(structure, aid, model)
    required = set(spec['partners'] + spec['binding_proteins'] + spec.get('support_instances', []))
    if not set(spec['partners']).isdisjoint(spec['binding_proteins']) or not required <= instance_map.keys():
        raise ValueError('Manifest must identify different observed partner and receptor instances in the selected Assembly/model')
    required_labels = {output_label(instance_map[k]) for k in required}
    chosen = required_labels if expand else choose_connected(source_block, required_labels)
    labels = {output_label(v) for v in instance_map.values()}
    expanded = instance_map if expand else None
    if not chosen <= labels:
        raise MappingError('covalent_partner_outside_assembly', 'A required connected chain is outside the selected Assembly')
    observed = inspect_source(source_block, structure, entry, revision, spec.get('role_rules', {}), scope={(model_num, aid)})
    partner_set = set(spec['partners'])
    observations = [o for o in observed if {m['instance'] for m in o['members']} <= partner_set]
    if {m['instance'] for o in observations for m in o['members']} != partner_set:
        raise MappingError('candidate_identity_pending', 'Could not recover all requested members from native candidate observations')
    contacts = {k for o in observations for k in o['contact_instances_4A']}
    if not set(spec['binding_proteins']) <= contacts:
        raise MappingError('binding_chain_identity_pending', 'Declared binding chain lacks observation in the existing contact check')
    members = [m for o in observations for m in o['source_member_identity']]
    key = key_for(entry, revision, model_num, aid, members)
    art, atoms = export(run/'systems'/entry/'system.cif', source_block, entry, aid, [model_num], chosen, verify=verify, expanded=expanded)
    checks = projected_evidence(spec, art, atoms, instance_map)
    contact_checks = spec.get('checks', {})
    if 'contact_instances_4A' in contact_checks:
        if contacts != set(contact_checks['contact_instances_4A']):
            raise ValueError(f'Native contact membership differs: {sorted(contacts)}')
        result_distances = disk_contact_distances(atoms, instance_map, spec['partners'], spec['binding_proteins'])
        if {k for k, d in result_distances.items() if d <= 4.} != contacts:
            raise ValueError('System CIF contact membership differs from source')
        checks['system_contact_distances_4A'] = result_distances
    pair = spec.get('partner_group')
    if pair == 'deposited_duplex':
        label_members = {output_label(instance_map[k]) for k in partner_set}
        accepted = [r for r in art['base_pairs'] if {e[LABEL] for e in r['ends']} == label_members]
        if len(label_members) != 2 or not accepted or len(accepted) != spec['checks']['base_pair_count']:
            raise MappingError('duplex_pairing_pending', 'Requested group does not match the deposited base pairs on disk')
        checks['partner_group'] = {'type': 'deposited_duplex', 'members': sorted(partner_set),
                                   'source_and_disk_category': '_ndb_struct_na_base_pair.', 'pairs': accepted}
    artifacts = {'system': art}
    if spec.get('complete_assembly'):
        selected_models = [model_num] if spec['complete_assembly'] == 'representative' else models
        archive_models = source_block.find_value('_pdbx_nmr_ensemble.conformers_submitted_total_number')
        if spec.get('checks', {}).get('full_model_count', len(selected_models)) != len(selected_models) or (
                spec['complete_assembly'] == 'all' and text(archive_models) and int(text(archive_models)) != len(selected_models)):
            raise ValueError('Observed model count differs from the declared experiment/model expectation')
        if any(describe_instances(structure, aid, m) != (instance_map, expand) for m in structure if m.num in selected_models):
            raise MappingError('model_member_identity_changed', 'Full model scope has different Assembly identities')
        model_scope = 'all-models' if spec['complete_assembly'] == 'all' else f'model-{model_num}'
        dest = run/'assemblies'/f'{entry}-v{revision}-a{aid}-{model_scope}-all-altloc.cif'
        full, full_atoms = export(dest, source_block, entry, aid, selected_models, labels, complete=True, verify=verify, expanded=expanded)
        artifacts['complete_assembly'] = full
        if 'assembly_protein_instances' in contact_checks:
            protein_entities = {text(r['entity_id']) for r in rows(source_block, '_entity_poly.')
                                if text(r['type']).startswith('polypeptide')}
            proteins = sorted(k for k, v in instance_map.items() if v['entity'] in protein_entities)
            if proteins != sorted(contact_checks['assembly_protein_instances']):
                raise ValueError(f'Full Assembly protein membership differs: {proteins}')
            rep_atoms = [r for r in full_atoms if text(r[MODEL]) == str(model_num)]
            disk_distances = disk_contact_distances(rep_atoms, instance_map, spec['partners'], proteins)
            if {k for k, d in disk_distances.items() if d <= 4.} != contacts:
                raise ValueError('Full Assembly disk contact membership differs from source')
            checks['assembly_protein_instances'] = proteins
            checks['assembly_contact_distances_A'] = disk_distances
    alt_role = {o['role'] for o in observations}
    role = 'PASS' if alt_role == {'target_ligand'} else 'REVIEW'
    result = {'pdb_id': entry, 'revision': revision, 'assembly': aid, 'representative_source_model': model_num,
              'representative_choice': 'manifest' if spec.get('representative_model') else
              ('deposited _pdbx_nmr_representative.conformer_id' if text(representative) else 'first observed model'),
              'source': {'path': str(source), 'reference': str(reference), 'proof': str(proof_path), 'sha256': actual, 'bytes': proof['bytes'], 'gzip_crc': 'pass', 'url': ident['source_url']},
              'fingerprint': fingerprint, 'checked_utc': now(), 'structure_acceptance': 'PASS', 'coordinate_identity': 'PASS',
              'experimental_fidelity': 'PASS', 'topology': 'PASS', 'role_qc': role, 'candidate_intake': 'REVIEW' if role == 'REVIEW' else 'not_run',
              'candidate_key': key, 'assembly_instances': instance_map, 'source_member_identity': members,
              'binding_protein_instances': spec['binding_proteins'], 'partner_instances': spec['partners'],
              'support_instances': spec.get('support_instances', []), 'included_by_covalent_closure': sorted(chosen-required_labels),
              'all_native_contact_instances_4A': sorted(contacts), 'observed_roles': observations,
              'sequence_evidence_separate_from_role': sequence_evidence(source_block, {instance_map[k]['entity'] for k in partner_set}),
              'selection_rationale': spec.get('selection_rationale'), 'checks': checks, 'artifacts': artifacts}
    if not verify:
        atomic_json(receipt, result)
    return result, 'verified' if verify else 'built'


def deliver(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--source-run', type=Path, required=True)
    parser.add_argument('--v1-manifest', type=Path, required=True)
    parser.add_argument('--entries', nargs='+')
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args(argv)
    start, clock = now(), time.monotonic()
    run, src = args.run.resolve(), args.source_run.resolve()
    if run == src or run.is_relative_to(src) or src.is_relative_to(run):
        parser.error('V1 output and historical source runs must be disjoint')
    manifest = read(args.v1_manifest)
    entries = manifest['entries']
    ids = [e['pdb_id'] for e in entries]
    wanted = set(args.entries or ids)
    if not ids or len(set(ids)) != len(ids) or any(not re.fullmatch(r'[0-9A-Z]{4}', p) for p in ids) or not wanted <= set(ids) or args.jobs < 1:
        parser.error('Require unique approved IDs, a valid entry selection and positive jobs')
    code = {name: sha256(Path(__file__).with_name(name)) for name in ('holo_m2_v1.py', 'holo_v1_export.py', 'holo_m2_pilot.py',
            'holo_source_mapping.py', 'holo_native_candidates.py', 'holo_candidate_batch.py', 'holo_pilot.py')}
    fixed = {'manifest': manifest, 'manifest_sha256': sha256(args.v1_manifest), 'source_run': str(src), 'code_sha256': code, 'gemmi': gemmi.__version__}
    state = run/'v1_inputs.json'
    if state.exists() and read(state) != fixed:
        parser.error('This Run already pins different inputs or code; use a new Run directory')
    if not state.exists() and not args.verify_only:
        atomic_json(state, fixed)
    import hashlib
    import json
    fingerprint = hashlib.sha256(json.dumps(fixed, sort_keys=True).encode()).hexdigest()
    results, work = {}, {}
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        jobs = {pool.submit(one, e, src, run, fingerprint, args.verify_only): e['pdb_id'] for e in entries if e['pdb_id'] in wanted}
        for future in as_completed(jobs):
            p = jobs[future]
            try:
                result, work[p] = future.result()
            except Exception as error:
                review = isinstance(error, MappingError)
                result = {'pdb_id': p, 'structure_acceptance': 'REVIEW' if review else 'FAIL', 'coordinate_identity': 'not_confirmed',
                          'experimental_fidelity': 'not_confirmed', 'topology': 'not_confirmed', 'role_qc': 'not_run',
                          'reason': str(error), 'error_kind': getattr(error, 'kind', type(error).__name__),
                          'evidence': getattr(error, 'evidence', {}), 'checked_utc': now()}
                work[p] = 'review' if review else 'failed'
                if not args.verify_only:
                    atomic_json(run/'qc'/f'{p}.json', result)
            results[p] = result
            print(p, work[p], result['structure_acceptance'], result['role_qc'], result.get('reason', ''), flush=True)
    if not args.verify_only:
        for p in set(ids)-results.keys():
            q = run/'qc'/f'{p}.json'
            if q.is_file():
                results[p] = read(q)
    receipt = {'run': str(run), 'batch_id': manifest['batch_id'], 'planned': len(ids), 'recorded': len(results),
               'structure_states': dict(Counter(x['structure_acceptance'] for x in results.values())),
               'role_states': dict(Counter(x['role_qc'] for x in results.values())), 'invocation': work, 'jobs': args.jobs,
               'started_utc': start, 'finished_utc': now(), 'wall_seconds': round(time.monotonic()-clock, 3),
               'entries': [results[p] for p in ids if p in results]}
    if not args.verify_only:
        atomic_json(run/'run_receipt.json', receipt)
    print(json.dumps({k: v for k, v in receipt.items() if k != 'entries'}), flush=True)
    return 0 if all(results[p]['structure_acceptance'] == 'PASS' for p in wanted) else 2
