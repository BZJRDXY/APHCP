"""Compare one observed ligand member with a candidate component in a local site frame."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import time

import gemmi
import numpy as np
from scipy.spatial import cKDTree

import apo_c_pair
import apo_c_pilot
import apo_c_triplet_regression
import holo_source_mapping
from apo_c_pair import align_coordinates, ca, get_sequence, sequence_map
from apo_c_pilot import now, read_json, sha, write_json
from apo_c_triplet_regression import positive_heavy, residue_index
from holo_source_mapping import assembly_instances


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', required=True)
    p.add_argument('--evidence', required=True)
    p.add_argument('--holo-instance', required=True)
    p.add_argument('--holo-component', required=True)
    p.add_argument('--holo-author-chain', required=True)
    p.add_argument('--holo-author-seq', required=True)
    p.add_argument('--candidate-component', required=True)
    p.add_argument('--candidate-author-chain', required=True)
    p.add_argument('--candidate-author-seq', required=True)
    args = p.parse_args()
    run, evidence_path = Path(args.run), Path(args.evidence)
    if not evidence_path.resolve().is_relative_to(run.resolve()):
        raise ValueError('evidence must be inside the Run')
    e = read_json(evidence_path)
    if e.get('execution_status') != 'MEASURED':
        raise ValueError('parent evidence is not measured')
    hp, cp = Path(e['holo_source']['path']), Path(e['candidate_source']['path'])
    if sha(hp) != e['holo_source']['sha256'] or sha(cp) != e['candidate_source']['sha256']:
        raise ValueError('fixed source changed')
    hs, cs = gemmi.read_structure(str(hp)), gemmi.read_structure(str(cp))
    hi = assembly_instances(hs, str(e['holo_identity']['assembly']), hs[0], include_water=True)
    ci = assembly_instances(cs, str(e['candidate_source']['assembly']), cs[0], include_water=True)
    hkey, ckey = e['holo_identity']['binding_instances'][0][0], e['candidate_source']['instance']
    hop, _, hrows = hi[hkey]
    cop, _, crows = ci[ckey]
    hr, cr = residue_index(hrows), residue_index(crows)
    mapping, seq = sequence_map(get_sequence(hs, hrows), get_sequence(cs, crows))
    lop, _, lrows = hi[args.holo_instance]
    found = [res for chain, res in lrows if chain == args.holo_author_chain and
             str(res.seqid) == args.holo_author_seq and res.name == args.holo_component]
    if len(found) != 1:
        raise ValueError('Holo target member is not unique in the selected Assembly')
    hlig = positive_heavy(found[0], lop)
    if not hlig:
        raise ValueError('Holo target member has no observed heavy atoms')
    hlig_xyz = np.array([x['xyz'] for x in hlig])
    htree = cKDTree(hlig_xyz)
    site, candidate_site, X, Y = [], [], [], []
    for hn, (hauth, hres) in sorted(hr.items()):
        ha = positive_heavy(hres, hop)
        if not ha or float(np.min(htree.query(np.array([a['xyz'] for a in ha]))[0])) > 4.0:
            continue
        cn = mapping.get(hn)
        row = {'holo_label_seq': hn, 'holo_author_chain': hauth, 'holo_author_seq': str(hres.seqid),
               'holo_component': hres.name, 'candidate_label_seq': cn}
        if cn in cr:
            cauth, cres = cr[cn]
            ca_atoms = positive_heavy(cres, cop)
            for atom in ca_atoms:
                candidate_site.append(atom | {'label_seq': cn, 'author_chain': cauth,
                                              'author_seq': str(cres.seqid), 'residue': cres.name})
            a, b = ca(hres, hop), ca(cres, cop)
            if a is not None and b is not None:
                X.append(a)
                Y.append(b)
            row.update({'candidate_author_chain': cauth, 'candidate_author_seq': str(cres.seqid),
                        'candidate_component': cres.name, 'candidate_has_CA': b is not None})
        else:
            row['candidate_observed'] = False
        site.append(row)
    if len(X) < 3 or not candidate_site:
        raise ValueError('insufficient mapped local site atoms')
    rotation, translation, local_rmsd = align_coordinates(np.array(X), np.array(Y))
    probe = hlig_xyz @ rotation + translation
    matches = []
    for inst, (operator, label, rows) in ci.items():
        for chain, res in rows:
            if (chain == args.candidate_author_chain and str(res.seqid) == args.candidate_author_seq
                    and res.name == args.candidate_component):
                matches.append((inst, label, operator, res))
    if len(matches) != 1:
        raise ValueError('candidate component is not unique in the selected Assembly')
    instance, label, operator, residue = matches[0]
    catoms = positive_heavy(residue, operator)
    coords = np.array([a['xyz'] for a in catoms])
    probe_dist = cKDTree(probe).query(coords)[0]
    site_dist, site_index = cKDTree(np.array([a['xyz'] for a in candidate_site])).query(coords)
    atom_evidence = []
    for i, atom in enumerate(catoms):
        atom_evidence.append({'candidate_atom': atom, 'distance_to_local_holo_member_probe_A': round(float(probe_dist[i]), 4),
                              'distance_to_actual_mapped_site_A': round(float(site_dist[i]), 4),
                              'closest_actual_site_atom': candidate_site[int(site_index[i])]})
    base = run / 'subsites' / (args.holo_component + '_' + args.candidate_component)
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
        snapshot = run / 'code' / f'{path.stem}-{digest[:12]}.py'
        if not snapshot.exists():
            shutil.copy2(path, snapshot)
        deps[path.name] = digest
    result = {'created_utc': now(), 'parent_evidence': str(evidence_path), 'parent_evidence_sha256': sha(evidence_path),
              'dependency_sha256': deps, 'request': vars(args), 'holo_fixed_sha256': sha(hp), 'candidate_fixed_sha256': sha(cp),
              'holo_member_heavy_atoms': hlig, 'holo_member_contact_residues': site, 'candidate_mapped_site_atoms': len(candidate_site),
              'local_alignment_CA_count': len(X), 'local_alignment_CA_rmsd_A': round(local_rmsd, 4),
              'candidate_component_instance': instance, 'candidate_component_label_asym': label,
              'candidate_component_atom_count': len(catoms),
              'candidate_atoms_within_3A_of_local_member_probe': int(np.sum(probe_dist <= 3.0)),
              'candidate_atoms_within_4A_of_actual_mapped_site': int(np.sum(site_dist <= 4.0)),
              'minimum_to_local_member_probe_A': round(float(np.min(probe_dist)), 4),
              'minimum_to_actual_mapped_site_A': round(float(np.min(site_dist)), 4), 'atom_evidence': atom_evidence,
              'limits': ['Projection is a diagnostic and is not a direct experimental contact.',
                         'Actual candidate component/site distances and experiment records are separate evidence.',
                         'This comparison covers the specified observed member only, not missing target regions.']}
    write_json(out / 'evidence.json', result)
    print(json.dumps({'path': str(out), 'local_rmsd_A': result['local_alignment_CA_rmsd_A'],
                      'actual_site_atoms_within_4A': result['candidate_atoms_within_4A_of_actual_mapped_site'],
                      'probe_atoms_within_3A': result['candidate_atoms_within_3A_of_local_member_probe']}), flush=True)


if __name__ == '__main__':
    main()
