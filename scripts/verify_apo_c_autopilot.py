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


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', required=True)
    p.add_argument('--output-name', default='independent_autopilot_audit.json')
    args = p.parse_args()
    if Path(args.output_name).name != args.output_name:
        p.error('output name must be a filename')
    run = Path(args.run)
    output = run / 'receipts' / args.output_name
    if output.exists():
        raise FileExistsError(output)
    samples = {x['id']: x for x in json.loads((run / 'samples.json').read_text())}
    wanted = [('1E2E', '1E2D'), ('1QZ5', '2PBD'), ('1FC0', '1EM6'),
              ('1OLM', '1O6U'), ('1MJ9', '1FY7'), ('1OWZ', '2O7A')]
    checks = []
    for holo, candidate in wanted:
        casepath = next((run / 'automated').glob(f'holo_{holo}_*/result.json'))
        case = json.loads(casepath.read_text())
        row = next(x for x in case['candidates'] if x['pdb'] == candidate)
        evidence = json.loads(Path(row['evidence_path']).read_text())
        entity = next(e for e in evidence['entities'] if any(
            'research_screen' in o for o in e.get('assembly_observations', [])))
        observation = next(o for o in entity['assembly_observations'] if 'research_screen' in o)
        sample = samples[case['system_id']]
        result = check_observation(sample, case['holo_evidence'], evidence, entity, observation)
        result.update({'holo': holo, 'candidate': candidate, 'evidence_path': row['evidence_path'],
                       'holo_source_hash_matches': digest(sample['source_uri']) == sample['identity']['source_sha256'],
                       'candidate_source_hash_matches': digest(evidence['candidate_source']['path']) ==
                                                        evidence['candidate_source']['sha256']})
        checks.append(result)
    readbacks = []
    for casepath in sorted((run / 'automated').glob('*/result.json')):
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
    all_pass = all(c['holo_site_matches_source'] and c['candidate_site_observed_count_matches']
                   and c['candidate_contact_set_and_distances_match'] and c['holo_source_hash_matches']
                   and c['candidate_source_hash_matches'] for c in checks)
    all_pass = all_pass and all(x['hash_ok'] and x['atoms'] > 0 and x.get('atom_count_ok', True) for x in readbacks)
    result = {'created_utc': datetime.now(timezone.utc).isoformat(), 'code_sha256': digest(__file__),
              'checks': checks, 'structure_readbacks': readbacks, 'all_pass': all_pass,
              'limits': ['Selected after seeing the technical trial; not a blind sample.',
                         'Verifies geometry and source readback, not experimental unbound status or formal Pair acceptance.']}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
    snap = run / 'code' / f'verify_apo_c_autopilot-{digest(__file__)[:12]}.py'
    if not snap.exists():
        shutil.copy2(__file__, snap)
    print(json.dumps({'output': str(output), 'source_observations': len(checks),
                      'structure_readbacks': len(readbacks), 'all_pass': all_pass}))
    if not all_pass:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
