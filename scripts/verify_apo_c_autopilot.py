"""Read source mmCIF independently and audit a few automated Holo-only observations."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
from datetime import datetime, timezone

import gemmi
import numpy as np
from scipy.spatial import cKDTree


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def instances(structure, assembly_id):
    assembly = next(a for a in structure.assemblies if str(a.name) == str(assembly_id))
    offsets, result = Counter(), {}
    for generator in assembly.generators:
        for operator in generator.operators:
            for label in generator.subchains:
                offsets[label] += 1
                rows = [(chain.name, residue) for chain in structure[0] for residue in chain
                        if residue.subchain == label]
                if rows:
                    result[f'{offsets[label]}.{label}'] = (operator, rows)
    return result


def atoms(residue, operator):
    values = []
    for atom in residue:
        if atom.occ <= 0 or atom.element.is_hydrogen:
            continue
        p = operator.transform.apply(atom.pos)
        values.append([p.x, p.y, p.z])
    return values


def check_observation(sample, holo_evidence, evidence, entity, observation):
    hs = gemmi.read_structure(sample['source_uri'])
    hi = instances(hs, sample['identity']['assembly'])
    target_points = []
    for key, _, _, _, author, number, component in sample['identity']['target_residue_instances']:
        op, rows = hi[key]
        matched = [r for ch, r in rows if ch == author and str(r.seqid) == str(number) and r.name == component]
        if len(matched) != 1:
            raise ValueError('target member not unique')
        target_points.extend(atoms(matched[0], op))
    htree = cKDTree(np.array(target_points))
    hop, hrows = hi[sample['identity']['binding_instances'][0][0]]
    actual_holo_site = set()
    for _, residue in hrows:
        if residue.label_seq is None:
            continue
        points = atoms(residue, hop)
        if points and min(htree.query(np.array(points))[0]) <= 4:
            actual_holo_site.add(residue.label_seq)
    recorded_holo_site = {row['label_seq'] for row in holo_evidence['site_residues']}

    cs = gemmi.read_structure(evidence['candidate_source']['path'])
    ci = instances(cs, observation['assembly'])
    cop, crows = ci[observation['instance']]
    residues = {r.label_seq: r for _, r in crows if r.label_seq is not None}
    site_atoms, actual_observed = [], 0
    for row in observation['mapped_site']:
        residue = residues.get(row['candidate_label_seq'])
        points = atoms(residue, cop) if residue is not None else []
        if points:
            actual_observed += 1
            site_atoms.extend(points)
    if not site_atoms:
        raise ValueError('selected observation has no actual mapped site atoms')
    ctree = cKDTree(np.array(site_atoms))
    contact = {}
    for key, (op, rows) in ci.items():
        if key == observation['instance']:
            continue
        for author, residue in rows:
            if residue.is_water():
                continue
            points = atoms(residue, op)
            if not points:
                continue
            minimum = float(min(ctree.query(np.array(points))[0]))
            if minimum <= 4:
                label = key.split('.', 1)[1]
                contact[(key, label, author, str(residue.seqid), residue.name)] = minimum
    reported = {}
    for row in observation['research_flags']['nonwater_within_4A']:
        key = (row['instance'], row['label_asym'], row['author_chain'], row['author_seq'], row['component'])
        reported[key] = row['min_distance_to_mapped_site_atoms_A']
    same_contacts = set(contact) == set(reported) and all(abs(contact[k] - reported[k]) < 0.00011 for k in contact)
    return {'holo_site_matches_source': actual_holo_site == recorded_holo_site,
            'holo_site_residue_count': len(actual_holo_site),
            'candidate_site_observed_count_matches': actual_observed == observation['site_observed_residues'],
            'candidate_contact_set_and_distances_match': same_contacts,
            'actual_nonwater_within_4A': len(contact),
            'instance': observation['instance'], 'assembly': observation['assembly']}


def audit_delivery_assembly(source, original_metadata, assembly_id, model_num, output, required_connections=()):
    """Compare the final deliverable against its fixed experimental source."""
    block = gemmi.cif.read(str(source)).sole_block()
    structure = gemmi.make_structure_from_block(block)
    pdb = gemmi.cif.as_string(block.find_value('_entry.id') or '')
    revision_rows = block.find('_pdbx_audit_revision_history.', ['data_content_type', 'major_revision', 'minor_revision'])
    revisions = [(int(row[1]), int(row[2])) for row in revision_rows
                 if gemmi.cif.as_string(row[0]).casefold() == 'structure model']
    seen_revision = '.'.join(map(str, max(revisions))) if revisions else None
    assembly = next(a for a in structure.assemblies if a.name == str(assembly_id))
    model = next(m for m in structure if m.num == model_num)
    naming = gemmi.HowToNameCopiedChain.AddNumber
    expected = gemmi.make_assembly(assembly, model, naming)
    actual = gemmi.read_structure(output['path'])

    def inventory(m):
        return sorted((c.name, str(r.seqid), r.name, a.name, a.element.name, a.altloc,
                       (a.pos.x, a.pos.y, a.pos.z, a.occ))
                      for c in m for r in c for a in r)

    first, last = inventory(expected), inventory(actual[0])
    same_identifiers = len(first) == len(last) and all(x[:-1] == y[:-1] for x, y in zip(first, last))
    coordinate_delta = occupancy_delta = None
    if same_identifiers and first:
        delta = np.asarray([x[-1] for x in first]) - np.asarray([y[-1] for y in last])
        coordinate_delta = float(np.linalg.norm(delta[:, :3], axis=1).max())
        occupancy_delta = float(np.abs(delta[:, 3]).max())

    def connections(st):
        return Counter((c.name, str(c.type), str(c.partner1), str(c.partner2)) for c in st.connections)

    applicable = structure.clone()
    for i in range(len(applicable) - 1, -1, -1):
        if applicable[i].num != model_num:
            del applicable[i]
    applicable.transform_to_assembly(str(assembly_id), naming, merge_dist=0.)
    required = set(required_connections)
    same_links = (required <= {c.name for c in structure.connections}
                  and required <= {c.name for c in actual.connections}
                  and connections(applicable) == connections(actual)
                  and all(actual[0].find_cra(p).atom is not None for c in actual.connections
                          for p in (c.partner1, c.partner2)))
    source_ok = (pdb.upper() == original_metadata['pdb'].upper()
                 and seen_revision == original_metadata['revision'] and digest(source) == original_metadata['sha256'])
    return {'path': output['path'], 'fixed_original': str(source), 'source_identity_revision_and_hash_ok': source_ok,
            'entry': pdb, 'revision': seen_revision, 'assembly': str(assembly_id), 'model': model_num,
            'output_model_matches': len(actual) == 1 and actual[0].num == model_num,
            'atom_ids_and_altloc_ok': same_identifiers,
            'coordinates_ok': coordinate_delta is not None and coordinate_delta <= 0.001,
            'max_coordinate_delta_A': coordinate_delta, 'occupancy_ok': occupancy_delta is not None and occupancy_delta <= 0.0005,
            'max_occupancy_delta': occupancy_delta, 'applicable_connection_endpoints_ok': same_links,
            'applicable_connections': len(actual.connections), 'required_connections': sorted(required),
            'hash_ok': digest(output['path']) == output['qc']['sha256'], 'atoms': len(last)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', required=True)
    p.add_argument('--output-name', default='independent_autopilot_audit.json')
    p.add_argument('--pair', action='append', nargs=2, metavar=('HOLO', 'CANDIDATE'),
                   help='Check each observable Assembly copy for the selected Entry pair')
    args = p.parse_args()
    if Path(args.output_name).name != args.output_name:
        p.error('output name must be a filename')
    run = Path(args.run)
    output = run / 'receipts' / args.output_name
    if output.exists():
        raise FileExistsError(output)
    samples = {x['id']: x for x in json.loads((run / 'samples.json').read_text())}
    wanted = args.pair or [('1E2E', '1E2D'), ('1QZ5', '2PBD'), ('1FC0', '1EM6'),
                           ('1OLM', '1O6U'), ('1MJ9', '1FY7'), ('1OWZ', '2O7A')]
    casefiles = [max([p, *p.parent.glob('result.retry*.json')])
                 for p in sorted((run / 'automated').glob('*/result.json'))]
    checks = []
    for holo, candidate in wanted:
        casepath = next(p for p in casefiles if p.parent.name.startswith(f'holo_{holo}_'))
        case = json.loads(casepath.read_text())
        row = next(x for x in case['candidates'] if x['pdb'] == candidate)
        evidence = json.loads(Path(row['evidence_path']).read_text())
        sample = samples[case['system_id']]
        observable = [(e, o) for e in evidence['entities'] for o in e.get('assembly_observations', [])
                      if 'research_screen' in o]
        if not observable:
            raise ValueError(f'No physically observable site for {holo} → {candidate}')
        if not args.pair:
            observable = observable[:1]
        hhash = digest(sample['source_uri']) == sample['identity']['source_sha256']
        chash = digest(evidence['candidate_source']['path']) == evidence['candidate_source']['sha256']
        for entity, observation in observable:
            result = check_observation(sample, case['holo_evidence'], evidence, entity, observation)
            result.update({'holo': holo, 'candidate': candidate, 'evidence_path': row['evidence_path'],
                           'holo_source_hash_matches': hhash, 'candidate_source_hash_matches': chash})
            checks.append(result)
    readbacks, delivery_checks = [], []
    for casepath in casefiles:
        case = json.loads(casepath.read_text())
        if case.get('holo_system'):
            info = case['holo_system']
            s = gemmi.read_structure(info['path'])
            readbacks.append({'path': info['path'], 'hash_ok': digest(info['path']) == info['sha256'],
                              'atoms': s[0].count_atom_sites()})
        for row in case['candidates']:
            evidence = json.loads(Path(row['evidence_path']).read_text())
            for info in evidence.get('materialized_assemblies', []):
                s = gemmi.read_structure(info['path'])
                readbacks.append({'path': info['path'], 'hash_ok': digest(info['path']) == info['qc']['sha256'],
                                  'atoms': s[0].count_atom_sites(),
                                  'atom_count_ok': s[0].count_atom_sites() == info['qc']['expected_atoms']})
        for delivery in case.get('deliveries', []):
            h = delivery['holo_identity']
            c = delivery['apo_original']
            obs = delivery['apo_observation']
            targets = {e['id'] for e in case['holo_evidence']['target_connections']}
            delivery_checks.append(audit_delivery_assembly(
                delivery['holo_original'], {'pdb': h['pdb_id'], 'revision': h['revision'], 'sha256': h['source_sha256']},
                h['assembly'], case['holo_evidence']['model'], delivery['holo_assembly'], targets))
            delivery_checks.append(audit_delivery_assembly(
                c['path'], {'pdb': delivery['candidate'], 'revision': c['revision'], 'sha256': c['sha256']},
                obs['assembly'], obs['model'], delivery['apo_assembly']))
    all_pass = all(c['holo_site_matches_source'] and c['candidate_site_observed_count_matches']
                   and c['candidate_contact_set_and_distances_match'] and c['holo_source_hash_matches']
                   and c['candidate_source_hash_matches'] for c in checks)
    all_pass = all_pass and all(x['hash_ok'] and x['atoms'] > 0 and x.get('atom_count_ok', True) for x in readbacks)
    all_pass = all_pass and all(x['atoms'] > 0 and all(x[k] for k in (
        'source_identity_revision_and_hash_ok', 'output_model_matches', 'atom_ids_and_altloc_ok', 'coordinates_ok', 'occupancy_ok',
        'applicable_connection_endpoints_ok', 'hash_ok')) for x in delivery_checks)
    result = {'created_utc': datetime.now(timezone.utc).isoformat(), 'code_sha256': digest(__file__),
              'checks': checks, 'structure_readbacks': readbacks, 'delivery_assemblies': delivery_checks, 'all_pass': all_pass,
              'limits': ['Selected after seeing the technical trial; not a blind sample.',
                         'Verifies geometry and source readback, not experimental unbound status or formal Pair acceptance.']}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
    snap = run / 'code' / f'verify_apo_c_autopilot-{digest(__file__)[:12]}.py'
    if not snap.exists():
        shutil.copy2(__file__, snap)
    print(json.dumps({'output': str(output), 'source_observations': len(checks),
                      'structure_readbacks': len(readbacks), 'delivery_assemblies': len(delivery_checks), 'all_pass': all_pass}))
    if not all_pass:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
