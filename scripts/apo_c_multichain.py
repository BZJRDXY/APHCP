"""Inspect multi-chain assembly mappings without treating one matching chain as a pair."""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
import shutil
import sys

import gemmi
import numpy as np
from scipy.spatial import cKDTree

import apo_c_pair
import apo_c_pilot
import holo_source_mapping
from apo_c_pair import align_coordinates, ca, candidate_source, get_sequence, heavy, revision, sequence_map, write_assembly
from apo_c_pilot import now, read_json, sha, uniprot_mappings, write_json
from holo_source_mapping import assembly_instances


def accession_for(mappings, label):
    return {acc for acc, item in mappings.items()
            if any(m.get('struct_asym_id') == label for m in item.get('mappings', []))}


def residues_by_label(rr):
    out = {}
    for author, res in rr:
        if res.label_seq is None:
            continue
        if res.label_seq in out:
            raise ValueError(f'ambiguous residue at label_seq {res.label_seq}')
        out[res.label_seq] = (author, res)
    return out


def inspect(case, candidate, candidate_assembly, anchor_auth, simulate_missing, run, out):
    qbindings = [x[0] for x in case['identity']['binding_instances']]
    if len(qbindings) < 2:
        raise ValueError('this entry point requires at least two Holo binding instances')
    if len(qbindings) > 3:
        raise ValueError('more than three Holo binding instances require a separate bounded mapping plan')
    hpath = Path(case['source_uri'])
    if sha(hpath) != case['identity']['source_sha256']:
        raise ValueError('Holo source checksum differs')
    hb = gemmi.cif.read(str(hpath)).sole_block()
    hrev, hdate = revision(hb, case['pdb'])
    if hrev != case['identity']['revision']:
        raise ValueError('Holo source revision differs')
    hs = gemmi.make_structure_from_block(hb)
    cp, csource = candidate_source(run, candidate)
    cs = gemmi.read_structure(str(cp))
    hi = assembly_instances(hs, str(case['identity']['assembly']), hs[0], include_water=True)
    ci = assembly_instances(cs, candidate_assembly, cs[0], include_water=True)
    hm, hmi = uniprot_mappings(run, case['pdb'])
    cm, cmi = uniprot_mappings(run, candidate)
    if not set(qbindings) <= set(hi):
        raise ValueError('a required Holo chain instance is missing from its Assembly')
    candidates = {k: v for k, v in ci.items() if any(r.entity_type == gemmi.EntityType.Polymer for _, r in v[2])}
    if len(candidates) > 12:
        raise ValueError('more than twelve candidate chain instances require bounded preselection')
    anchor = {k for k, (_, _, rr) in candidates.items() if any(ch == anchor_auth for ch, _ in rr)}
    if not anchor:
        raise ValueError('requested candidate author chain is absent from the chosen Assembly')
    structure_evidence = {'holo_instances_required': qbindings, 'candidate_assembly_instances': sorted(candidates),
                          'anchor_author_chain': anchor_auth, 'anchor_instances': sorted(anchor),
                          'candidate_source': csource, 'holo_source': {'uri': str(hpath), 'revision': hrev,
                          'sha256': case['identity']['source_sha256'], 'first_public_model_date': hdate},
                          'sifts_requests': {'holo': hmi, 'candidate': cmi}}
    if simulate_missing:
        return structure_evidence | {'outcome': 'MISSING_HOLO_CHAIN_MAPPING',
               'mapping_attempt': {qbindings[0]: sorted(anchor)[0]}, 'unmapped_holo_instances': qbindings[1:],
               'counterfactual': 'The full experimental candidate is retained; this deliberately incomplete chain assignment tests the gate.'}

    if len(candidates) < len(qbindings):
        return structure_evidence | {'outcome': 'MISSING_CANDIDATE_CHAIN_COPIES'}
    pair_data = {}
    for hkey in qbindings:
        hop, hlabel, hrr = hi[hkey]
        hseq, hres = get_sequence(hs, hrr), residues_by_label(hrr)
        for ckey, (cop, clabel, crr) in candidates.items():
            common = accession_for(hm, hlabel) & accession_for(cm, clabel)
            if len(common) != 1:
                continue
            cseq, cres = get_sequence(cs, crr), residues_by_label(crr)
            mapping, stats = sequence_map(hseq, cseq)
            x, y = [], []
            for hn, cn in mapping.items():
                if hn in hres and cn in cres:
                    hp, cp2 = ca(hres[hn][1], hop), ca(cres[cn][1], cop)
                    if hp is not None and cp2 is not None:
                        x.append(hp)
                        y.append(cp2)
            if len(x) >= 3:
                pair_data[(hkey, ckey)] = {'x': x, 'y': y, 'map': mapping, 'stats': stats,
                                           'common': sorted(common), 'hres': hres, 'cres': cres}
    assignments = []
    for ordered in itertools.permutations(sorted(candidates), len(qbindings)):
        if not set(ordered) & anchor:
            continue
        keys = list(zip(qbindings, ordered))
        if any(k not in pair_data for k in keys):
            continue
        x = np.array([v for k in keys for v in pair_data[k]['x']])
        y = np.array([v for k in keys for v in pair_data[k]['y']])
        rotation, translation, fit = align_coordinates(x, y)
        assignments.append({'keys': keys, 'matched_ca': len(x), 'rmsd': fit,
                            'rotation': rotation, 'translation': translation})
    if not assignments:
        return structure_evidence | {'outcome': 'NO_COMPLETE_CHAIN_ASSIGNMENT', 'eligible_chain_pairs': len(pair_data)}
    max_count = max(x['matched_ca'] for x in assignments)
    assignments.sort(key=lambda x: (-x['matched_ca'], x['rmsd'], x['keys']))
    best = assignments[0]
    equivalent = [a for a in assignments if a['matched_ca'] == max_count and abs(a['rmsd'] - best['rmsd']) <= 0.01]

    atoms = []
    for instance, _, _, _, author, seq, name in case['identity']['target_residue_instances']:
        if instance not in hi:
            raise ValueError('target instance absent from Holo Assembly')
        op, _, rr = hi[instance]
        found = [r for ch, r in rr if ch == author and str(r.seqid) == str(seq) and r.name == name]
        if len(found) != 1:
            raise ValueError('target residue is not unique')
        atoms.extend(heavy(found[0], op))
    if not atoms:
        raise ValueError('target has no observed heavy atoms')
    lig = np.array([x[1] for x in atoms])
    htree = cKDTree(lig)
    site, px, py = [], [], []
    for hkey, ckey in best['keys']:
        data = pair_data[(hkey, ckey)]
        hop, _, _ = hi[hkey]
        cop, _, _ = ci[ckey]
        for hn, (author, hr) in sorted(data['hres'].items()):
            h_atoms = heavy(hr, hop)
            if not h_atoms or float(htree.query(np.array([x[1] for x in h_atoms]))[0].min()) > 4.0:
                continue
            cn = data['map'].get(hn)
            row = {'holo_instance': hkey, 'holo_label_seq': hn, 'holo_author_chain': author,
                   'holo_author_seq': str(hr.seqid), 'holo_ccd': hr.name,
                   'candidate_instance': ckey, 'candidate_label_seq': cn}
            if cn in data['cres']:
                cauth, cr = data['cres'][cn]
                hp, cp2 = ca(hr, hop), ca(cr, cop)
                row.update({'candidate_author_chain': cauth, 'candidate_author_seq': str(cr.seqid),
                            'candidate_ccd': cr.name, 'candidate_has_CA': cp2 is not None,
                            'candidate_missing_holo_heavy_atom_names': sorted(
                                {x[0] for x in h_atoms} - {x[0] for x in heavy(cr, cop)})})
                if hp is not None and cp2 is not None:
                    px.append(hp)
                    py.append(cp2)
            site.append(row)
    if len(px) < 3:
        raise ValueError('fewer than three mapped site CA atoms across all required chains')
    lr, lt, local_rmsd = align_coordinates(np.array(px), np.array(py))
    gtree = cKDTree(lig @ best['rotation'] + best['translation'])
    ltree = cKDTree(lig @ lr + lt)
    selected = {ckey for _, ckey in best['keys']}
    nearby = []
    for instance, (op, label, rr) in ci.items():
        if instance in selected:
            continue
        for author, res in rr:
            ra = heavy(res, op)
            if not ra:
                continue
            coords = np.array([x[1] for x in ra])
            gd, ld = gtree.query(coords)[0], ltree.query(coords)[0]
            if min(float(gd.min()), float(ld.min())) <= 4.0:
                nearby.append({'instance': instance, 'label_asym': label, 'author_chain': author,
                               'author_seq': str(res.seqid), 'ccd': res.name,
                               'kind': 'water' if res.is_water() else ('polymer' if res.entity_type == gemmi.EntityType.Polymer else 'nonpolymer'),
                               'global_min_A': round(float(gd.min()), 3), 'local_min_A': round(float(ld.min()), 3)})
    nearby.sort(key=lambda r: min(r['global_min_A'], r['local_min_A']))
    assembly_qc = write_assembly(cs, candidate_assembly, out / 'candidate_assembly.cif')
    shutil.copy2(case['system_uri'], out / 'holo_system.cif')
    hread = gemmi.read_structure(str(out / 'holo_system.cif'))
    if not hread or not hread[0].count_atom_sites():
        raise ValueError('Holo output did not read back')
    return structure_evidence | {
        'outcome': 'MULTICHAIN_MAPPED_PROVISIONAL',
        'best_assignment': [{'holo_instance': h, 'candidate_instance': c,
                             'common_uniprot': pair_data[(h, c)]['common'],
                             'sequence_alignment': pair_data[(h, c)]['stats']} for h, c in best['keys']],
        'assignment_count': len(assignments), 'near_equivalent_assignments': [a['keys'] for a in equivalent],
        'symmetry_or_copy_ambiguity': len(equivalent) > 1, 'global_ca_count': best['matched_ca'],
        'global_ca_rmsd_A': round(best['rmsd'], 4), 'pocket_ca_count': len(px),
        'local_pocket_ca_rmsd_A': round(local_rmsd, 4), 'site_residues': site,
        'nonselected_components_near_global_or_local_template': nearby,
        'candidate_assembly_qc': assembly_qc, 'holo_system_sha256': sha(out / 'holo_system.cif'),
        'limitations': ['The chosen homologous chain copies are a structural mapping hypothesis, not formal Pair acceptance.',
                        'All near-equivalent assignments are retained; no unique copy is inferred from a tie.',
                        'Target-template proximity is a diagnostic, especially when conformations differ.',
                        'Only the first deposited model is inspected.']}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', required=True)
    p.add_argument('--case', required=True)
    p.add_argument('--candidate', required=True)
    p.add_argument('--anchor-auth-chain', required=True)
    p.add_argument('--assembly', default='1')
    p.add_argument('--simulate-missing-chain', action='store_true')
    args = p.parse_args()
    run = Path(args.run)
    matches = [x for x in read_json(run / 'samples.json') if x['id'] == args.case]
    if len(matches) != 1:
        p.error('case must identify one selected Holo')
    base = run / 'multichain' / (args.case.replace(':', '_') + '__' + args.candidate.upper())
    base.mkdir(parents=True, exist_ok=True)
    n = 1
    while (base / f'attempt{n:03d}').exists():
        n += 1
    out = base / f'attempt{n:03d}'
    out.mkdir()
    deps = {}
    for pth in [Path(__file__), Path(apo_c_pair.__file__), Path(apo_c_pilot.__file__), Path(holo_source_mapping.__file__)]:
        digest = sha(pth)
        dest = run / 'code' / f'{pth.stem}-{digest[:12]}.py'
        if not dest.exists():
            shutil.copy2(pth, dest)
        deps[pth.name] = digest
    result = {'started_utc': now(), 'request': vars(args), 'dependency_sha256': deps}
    try:
        result.update(inspect(matches[0], args.candidate.upper(), args.assembly, args.anchor_auth_chain,
                              args.simulate_missing_chain, run, out))
        result['execution_status'] = 'INSPECTED'
    except Exception as exc:
        result.update({'execution_status': 'UNRESOLVED', 'error': f'{type(exc).__name__}: {exc}'})
    result['finished_utc'] = now()
    write_json(out / 'pair.json', result)
    print(json.dumps({'path': str(out), 'status': result['execution_status'],
                      'outcome': result.get('outcome'), 'error': result.get('error')}), flush=True)
    if result['execution_status'] != 'INSPECTED':
        sys.exit(1)


if __name__ == '__main__':
    main()
