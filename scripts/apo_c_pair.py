"""Reproducible, conservative single-chain structural inspection for the Apo pilot."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import re
import shutil
import sys

import gemmi
import numpy as np
from scipy.spatial import cKDTree

import apo_c_pilot
import holo_source_mapping
from apo_c_pilot import cache_get, now, open_catalog, read_json, relevant_uniprot, sha, uniprot_mappings, write_json
from holo_source_mapping import assembly_instances


def revision(block, pdb):
    if gemmi.cif.as_string(block.find_value('_entry.id')).upper() != pdb.upper():
        raise ValueError('coordinate entry ID differs from requested ID')
    rows = block.find('_pdbx_audit_revision_history.',
                      ['data_content_type', 'major_revision', 'minor_revision', 'revision_date'])
    model = [(int(r[1]), int(r[2]), gemmi.cif.as_string(r[3])) for r in rows
             if gemmi.cif.as_string(r[0]).casefold() == 'structure model']
    if not model:
        raise ValueError('missing Structure model revision history')
    major, minor, _ = max(model)
    return f'{major}.{minor}', min(x[2] for x in model)


def candidate_source(run, pdb):
    catalog = Path(read_json(run / 'run.json')['catalog'])
    with open_catalog(catalog) as db:
        assets = db.execute('SELECT a.uri,a.sha256,a.bytes,v.revision FROM asset_file a '
                            'JOIN entry_version v ON v.id=a.entry_version_id '
                            "WHERE v.pdb_id=? AND a.uri LIKE '%/raw/%.cif.gz'", (pdb,)).fetchall()
    assets = [a for a in assets if a[1] and Path(a[0]).is_file()]
    if assets:
        uri, digest, size, rev = max(assets, key=lambda a: tuple(map(int, a[3].split('.'))))
        path = Path(uri)
        if sha(path) != digest:
            raise ValueError(f'candidate fixed source checksum mismatch: {path}')
        block = gemmi.cif.read(str(path)).sole_block()
        observed_rev, date = revision(block, pdb)
        if observed_rev != rev:
            raise ValueError('candidate catalog revision differs from coordinate revision')
        return path, {'uri': str(path), 'sha256': digest, 'bytes': size, 'revision': rev,
                      'first_public_model_date': date, 'origin': 'sealed_catalog_asset_file'}

    metadata, info = cache_get(run, f'rcsb_entry/{pdb.lower()}',
                               f'https://data.rcsb.org/rest/v1/core/entry/{pdb}')
    accession = json.loads(metadata)['rcsb_accession_info']
    major, minor = int(accession['major_revision']), int(accession['minor_revision'])
    rev = f'{major}.{minor}'
    p = pdb.lower()
    url = (f'https://files-versioned.rcsb.org/pdb_versioned/data/entries/{p[1:3]}/'
           f'pdb_0000{p}/pdb_0000{p}_xyz_v{major}-{minor}.cif.gz')
    body, fetch = cache_get(run, f'candidate_file/{p}_v{major}-{minor}', url, timeout=90)
    if not body.startswith(b'\x1f\x8b'):
        raise ValueError('versioned coordinate response is not gzip')
    gzip.decompress(body)  # also verifies gzip integrity
    path = run / 'raw' / f'{p}-v{rev}.cif.gz'
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if sha(path) != fetch['response_sha256']:
            raise ValueError('stored fixed candidate differs from cached versioned response')
    else:
        path.write_bytes(body)
    block = gemmi.cif.read(str(path)).sole_block()
    observed_rev, date = revision(block, pdb)
    if observed_rev != rev:
        raise ValueError(f'downloaded revision {observed_rev} differs from selected {rev}')
    return path, {'uri': str(path), 'sha256': sha(path), 'bytes': path.stat().st_size,
                  'revision': rev, 'first_public_model_date': date, 'origin': 'versioned_wwpdb_download',
                  'metadata_request': info, 'coordinate_request': fetch}


def oneletter(name):
    code = gemmi.find_tabulated_residue(name).one_letter_code
    return code.upper() if code and code != ' ' else 'X'


def get_sequence(structure, residues):
    ids = {str(r.entity_id) for _, r in residues if r.entity_type == gemmi.EntityType.Polymer}
    if len(ids) != 1:
        raise ValueError(f'chain does not map to exactly one polymer entity: {sorted(ids)}')
    entity = next((e for e in structure.entities if str(e.name) == next(iter(ids))), None)
    if entity is None or not entity.full_sequence:
        raise ValueError('polymer entity has no deposited sequence')
    return [oneletter(name) for name in entity.full_sequence]


def sequence_map(query, target):
    alignment = gemmi.align_string_sequences(query, target, [], gemmi.AlignmentScoring())
    cigar = alignment.cigar_str()
    qi, ti, pairs = 0, 0, []
    for length, op in re.findall(r'(\d+)([MID])', cigar):
        n = int(length)
        if op == 'M':
            pairs.extend((qi + i + 1, ti + i + 1) for i in range(n))
            qi += n
            ti += n
        elif op == 'I':
            qi += n
        elif op == 'D':
            ti += n
    if qi != len(query) or ti != len(target):
        raise ValueError(f'unhandled alignment CIGAR {cigar}: {qi}/{len(query)}, {ti}/{len(target)}')
    identity = sum(query[q - 1] == target[t - 1] for q, t in pairs) / len(pairs) if pairs else 0.0
    mapped_query = {q for q, _ in pairs}
    mapped_target = {t for _, t in pairs}
    return dict(pairs), {'cigar': cigar, 'query_length': len(query), 'candidate_length': len(target),
                         'paired_positions': len(pairs), 'identity_on_aligned_positions': identity,
                         'query_coverage': len(pairs) / len(query), 'candidate_coverage': len(pairs) / len(target),
                         'mismatches': [{'holo_label_seq': q, 'candidate_label_seq': t,
                                         'holo': query[q - 1], 'candidate': target[t - 1]}
                                        for q, t in pairs if query[q - 1] != target[t - 1]],
                         'unaligned_holo_positions': [q for q in range(1, len(query) + 1) if q not in mapped_query],
                         'unaligned_candidate_positions': [t for t in range(1, len(target) + 1) if t not in mapped_target]}


def reference_segments(block, accession, author_chain):
    rows = block.find('_struct_ref_seq.', ['pdbx_db_accession', 'pdbx_strand_id', 'seq_align_beg', 'seq_align_end',
                                         'db_align_beg', 'db_align_end'])
    ranges = []
    for row in rows:
        acc, chains, *positions = [gemmi.cif.as_string(v) for v in row]
        if acc == accession and author_chain in (x.strip() for x in chains.split(',')) and all(p.isdecimal() for p in positions):
            ranges.append(tuple(map(int, positions)))
    return ranges


def reference_positions(block, accession, author_chain):
    """Map deposited sequence positions only through unambiguous ungapped reference spans."""
    positions = {}
    for first, last, start, end in reference_segments(block, accession, author_chain):
        if last - first != end - start or first <= 0 or start <= 0:
            continue
        for key, value in zip(range(first, last + 1), range(start, end + 1)):
            if positions.setdefault(key, value) != value:
                return {}
    return positions if len(set(positions.values())) == len(positions) else {}


def anchored_end_alignment(query, target, mapping, alignment, holo_reference, candidate_reference, site):
    """Resolve equivalent gap placements for a construct with an exact end extension."""
    if len(query) == len(target) or alignment['mismatches'] or not set(site) <= holo_reference.keys():
        return mapping, alignment
    candidate_index = {ref: seq for seq, ref in candidate_reference.items()}
    anchors = {q: candidate_index[ref] for q, ref in holo_reference.items() if ref in candidate_index}
    if not anchors or not set(site) <= anchors.keys():
        return mapping, alignment
    delta = anchors[min(anchors)] - min(anchors)
    if any(t - q != delta for q, t in anchors.items()):
        return mapping, alignment
    h_start, c_start = max(1, 1 - delta), max(1, 1 + delta)
    count = min(len(query), len(target))
    h_end, c_end = h_start + count - 1, c_start + count - 1
    if h_end > len(query) or c_end > len(target) or query[h_start - 1:h_end] != target[c_start - 1:c_end]:
        return mapping, alignment
    corrected = dict(zip(range(h_start, h_end + 1), range(c_start, c_end + 1)))
    if any(corrected.get(q) != c for q, c in anchors.items()) or corrected == mapping:
        return mapping, alignment
    operations = [(h_start - 1, 'I'), (c_start - 1, 'D'), (count, 'M'),
                  (len(query) - h_end, 'I'), (len(target) - c_end, 'D')]
    return corrected, alignment | {
        'basis': 'FIXED_UNIPROT_END_ALIGNMENT', 'reference_anchor_count': len(anchors),
        'cigar': ''.join(f'{n}{op}' for n, op in operations if n),
        'paired_positions': count, 'identity_on_aligned_positions': 1.0,
        'query_coverage': count / len(query), 'candidate_coverage': count / len(target),
        'unaligned_holo_positions': [q for q in range(1, len(query) + 1) if q not in corrected],
        'unaligned_candidate_positions': [t for t in range(1, len(target) + 1) if not c_start <= t <= c_end]}


def atom_pos(atom, operator):
    p = operator.transform.apply(atom.pos)
    return np.array([p.x, p.y, p.z])


def heavy(residue, op):
    return [(a.name, atom_pos(a, op), a.altloc if a.altloc not in ('\0', ' ', '.', '?') else '', a.occ)
            for a in residue if a.occ > 0 and not a.element.is_hydrogen]


def ca(residue, op):
    atoms = [(a.occ, a) for a in residue if a.name == 'CA' and a.occ > 0]
    return atom_pos(max(atoms, key=lambda v: v[0])[1], op) if atoms else None


def rmsd(x, y):
    return float(np.sqrt(np.mean(np.sum((x - y) ** 2, axis=1))))


def align_coordinates(holo, candidate):
    x0, y0 = holo.mean(axis=0), candidate.mean(axis=0)
    u, _, vt = np.linalg.svd((holo - x0).T @ (candidate - y0))
    flip = np.eye(3)
    flip[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ flip @ vt
    translation = y0 - x0 @ rotation
    return rotation, translation, rmsd(holo @ rotation + translation, candidate)


def write_assembly(structure, assembly_id, path, required_connections=()):
    assembly = next((a for a in structure.assemblies if a.name == assembly_id), None)
    if assembly is None:
        raise ValueError(f'candidate assembly absent: {assembly_id}')
    naming = gemmi.HowToNameCopiedChain.AddNumber
    model = gemmi.make_assembly(assembly, structure[0], naming)
    out = structure.clone()
    for i in range(len(out) - 1, 0, -1):
        del out[i]
    out.transform_to_assembly(assembly_id, naming, merge_dist=0.)
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    doc = out.make_mmcif_document()
    with path.open('x', encoding='utf8') as f:
        f.write(doc.as_string())
    reloaded = gemmi.read_structure(str(path))
    expected, actual = model.count_atom_sites(), reloaded[0].count_atom_sites()
    if expected != actual or actual == 0:
        raise ValueError(f'assembly readback atom count differs: {expected} vs {actual}')
    def inventory(m):
        return Counter((c.name, str(r.seqid), r.name, a.name, a.element.name, a.altloc, round(a.occ, 3))
                       for c in m for r in c for a in r)
    if inventory(model) != inventory(reloaded[0]):
        raise ValueError('assembly readback changed atom identity, altloc or occupancy')

    def atom_positions(m):
        return sorted((c.name, str(r.seqid), r.name, a.name, a.altloc, (a.pos.x, a.pos.y, a.pos.z))
                      for c in m for r in c for a in r)

    src, dst = atom_positions(model), atom_positions(reloaded[0])
    assert all(a[:-1] == b[:-1] for a, b in zip(src, dst))
    max_delta = max(float(np.linalg.norm(np.subtract(a[-1], b[-1]))) for a, b in zip(src, dst))
    if max_delta > 0.001:
        raise ValueError(f'assembly readback coordinates differ by {max_delta:.6f} A')

    def links(st):
        identities = []
        for link in st.connections:
            for endpoint in (link.partner1, link.partner2):
                if st[0].find_cra(endpoint).atom is None:
                    raise ValueError(f'connection endpoint absent: {link.name} {endpoint}')
            identities.append((str(link.type), str(link.partner1), str(link.partner2)))
        return Counter(identities)

    if links(out) != links(reloaded):
        raise ValueError('assembly readback changed deposited connection endpoints')
    kept_links = {link.name for link in reloaded.connections}
    if set(required_connections) - kept_links:
        raise ValueError(f'required deposited connections absent in assembly: {sorted(set(required_connections) - kept_links)}')
    # A source link can legitimately fall outside the assembly. Preserve a
    # diagnostic for the source/display difference rather than hiding it.
    return {'expected_atoms': expected, 'readback_atoms': actual, 'component_atom_inventory': 'PASS',
            'max_coordinate_delta_A': max_delta, 'source_connections': len(structure.connections),
            'readback_connection_endpoints': len(reloaded.connections), 'connection_endpoints': 'PASS',
            'required_connections': sorted(required_connections),
            'display_scope': 'first deposited model; atom coordinates transformed to selected assembly, copy chains renamed',
            'sha256': sha(path), 'bytes': path.stat().st_size}


def materialize(case, structure, assembly_id, result_dir):
    assembly_qc = write_assembly(structure, assembly_id, result_dir / 'candidate_assembly.cif')
    holo_path = result_dir / 'holo_system.cif'
    shutil.copy2(case['system_uri'], holo_path)
    holo_readback = gemmi.read_structure(str(holo_path))
    if not holo_readback or holo_readback[0].count_atom_sites() == 0:
        raise ValueError('Holo system readback contains no atoms')
    return assembly_qc, {'uri': str(holo_path), 'sha256': sha(holo_path),
                         'readback_atoms': holo_readback[0].count_atom_sites()}


def inspect_pair(case, cand_pdb, cand_auth, assembly_id, run, result_dir, candidate_instance=None):
    if len(case['identity']['binding_instances']) != 1:
        raise ValueError('multiple Holo binding chain instances require explicit multi-chain mapping')
    hpath = Path(case['source_uri'])
    if sha(hpath) != case['identity']['source_sha256']:
        raise ValueError('Holo fixed source checksum differs')
    hblock = gemmi.cif.read(str(hpath)).sole_block()
    hrev, hdate = revision(hblock, case['pdb'])
    if hrev != case['identity']['revision']:
        raise ValueError('Holo fixed source revision differs')
    hstruct = gemmi.make_structure_from_block(hblock)
    cpath, csource = candidate_source(run, cand_pdb)
    cblock = gemmi.cif.read(str(cpath)).sole_block()
    cstruct = gemmi.make_structure_from_block(cblock)

    hinstances = assembly_instances(hstruct, str(case['identity']['assembly']), hstruct[0], include_water=True)
    cinstances = assembly_instances(cstruct, assembly_id, cstruct[0], include_water=True)
    hkey = case['identity']['binding_instances'][0][0]
    if hkey not in hinstances:
        raise ValueError(f'Holo binding chain instance missing: {hkey}')
    hop, hlabel, hrr = hinstances[hkey]
    possible = [(k, v) for k, v in cinstances.items() if any(ch == cand_auth and
                r.entity_type == gemmi.EntityType.Polymer for ch, r in v[2])]
    available = [k for k, _ in possible]
    if candidate_instance is not None:
        possible = [(k, v) for k, v in possible if k == candidate_instance]
    if len(possible) != 1:
        raise ValueError(f'candidate author chain has {len(available)} assembly instances {available}; '
                         'supply one explicit candidate instance')
    ckey, (cop, clabel, crr) = possible[0]

    hm, hmi = uniprot_mappings(run, case['pdb'])
    cm, cmi = uniprot_mappings(run, cand_pdb)
    hmap = relevant_uniprot(case, hm)
    cmap = {acc: [m for m in v.get('mappings', []) if m.get('struct_asym_id') == clabel]
            for acc, v in cm.items()}
    common = sorted(set(hmap) & {k for k, v in cmap.items() if v})
    if len(common) != 1:
        raise ValueError(f'cannot establish unique shared SIFTS UniProt accession: {common}')

    qseq, cseq = get_sequence(hstruct, hrr), get_sequence(cstruct, crr)
    mapping, alignment = sequence_map(qseq, cseq)
    def unique_residues(rr):
        out = {}
        for ch, r in rr:
            if r.label_seq is None:
                continue
            if r.label_seq in out:
                raise ValueError(f'ambiguous multiple residues at label_seq {r.label_seq}')
            out[r.label_seq] = (ch, r)
        return out
    hres, cres = unique_residues(hrr), unique_residues(crr)
    targets = case['identity']['target_residue_instances']
    lig_atoms = []
    for instance, _, _, _, author, seq, name in targets:
        if instance not in hinstances:
            raise ValueError(f'Holo target instance missing: {instance}')
        op, _, rr = hinstances[instance]
        found = [r for chain, r in rr if chain == author and str(r.seqid) == str(seq) and r.name == name]
        if len(found) != 1:
            raise ValueError(f'Holo target residue not unique: {instance}/{author}/{seq}/{name}')
        lig_atoms.extend(heavy(found[0], op))
    if not lig_atoms:
        raise ValueError('Holo target has no positive-occupancy heavy atoms')
    hlig = np.array([p for _, p, _, _ in lig_atoms])
    tree = cKDTree(hlig)
    pocket = []
    for num, (author, res) in sorted(hres.items()):
        atoms = heavy(res, hop)
        if atoms and float(tree.query(np.array([a[1] for a in atoms]))[0].min()) <= 4.0:
            pocket.append(num)
    if not pocket:
        raise ValueError('no receptor residues within 4.0 A of target')

    x, y = [], []
    for hn, cn in mapping.items():
        if hn in hres and cn in cres:
            hp, cp = ca(hres[hn][1], hop), ca(cres[cn][1], cop)
            if hp is not None and cp is not None:
                x.append(hp)
                y.append(cp)
    if len(x) < 3:
        raise ValueError(f'only {len(x)} corresponding CA atoms available for alignment')
    rotation, translation, global_rmsd = align_coordinates(np.array(x), np.array(y))
    transformed_lig = hlig @ rotation + translation
    target_tree = cKDTree(transformed_lig)
    receptor_template_clashes = []
    for num, (author, res) in sorted(cres.items()):
        atoms = heavy(res, cop)
        if not atoms:
            continue
        distances = target_tree.query(np.array([a[1] for a in atoms]))[0]
        if float(np.min(distances)) <= 2.0:
            receptor_template_clashes.append({'label_seq': num, 'auth_chain': author, 'auth_seq': str(res.seqid),
                                              'ccd': res.name, 'min_distance_A': round(float(np.min(distances)), 3),
                                              'atoms_within_2A': int(np.sum(distances <= 2.0))})

    residues, pocket_x, pocket_holo, pocket_y = [], [], [], []
    for hn in pocket:
        hauth, hr = hres[hn]
        cn = mapping.get(hn)
        row = {'holo_label_seq': hn, 'holo_auth_chain': hauth, 'holo_auth_seq': str(hr.seqid),
               'holo_ccd': hr.name, 'candidate_label_seq': cn}
        if cn in cres:
            cauth, cr = cres[cn]
            hnms = {name for name, _, _, _ in heavy(hr, hop)}
            cnms = {name for name, _, _, _ in heavy(cr, cop)}
            hp, cp = ca(hr, hop), ca(cr, cop)
            row.update({'candidate_auth_chain': cauth, 'candidate_auth_seq': str(cr.seqid),
                        'candidate_ccd': cr.name, 'same_residue': oneletter(hr.name) == oneletter(cr.name),
                        'holo_heavy_atoms_absent_in_candidate': sorted(hnms - cnms), 'candidate_has_CA': cp is not None})
            if hp is not None and cp is not None:
                pocket_x.append(hp @ rotation + translation)
                pocket_holo.append(hp)
                pocket_y.append(cp)
        else:
            row.update({'candidate_observed': False})
        residues.append(row)

    if len(pocket_holo) < 3:
        assembly_qc, holo_qc = materialize(case, cstruct, assembly_id, result_dir)
        return {'sample': case['id'], 'holo_status': case['holo_status'], 'holo_source': {
                    'uri': str(hpath), 'revision': hrev, 'sha256': case['identity']['source_sha256'],
                    'first_public_model_date': hdate, 'assembly': case['identity']['assembly'], 'model': hstruct[0].num},
                'candidate_pdb': cand_pdb, 'candidate_source': csource, 'candidate_assembly': assembly_id,
                'holo_binding_instance': hkey, 'candidate_binding_instance': ckey,
                'candidate_author_chain': cand_auth, 'candidate_label_asym': clabel,
                'sifts_common_uniprot': common, 'holo_sifts_mappings': hmap, 'candidate_sifts_mappings': cmap,
                'sifts_request': {'holo': hmi, 'candidate': cmi}, 'sequence_alignment': alignment,
                'target_residue_instances': targets, 'target_heavy_atoms': len(lig_atoms),
                'site_definition': 'Holo receptor residues with a positive-occupancy heavy atom within 4.0 A of observed target',
                'site_residues': residues, 'geometric_screen': 'SITE_OBSERVATION_OR_MAPPING_INSUFFICIENT',
                'site_observation_status': 'INSUFFICIENT_OR_UNMAPPED',
                'reason': 'Fewer than three mapped pocket CA atoms; local alignment and occupancy are not determined.',
                'alignment_ca_count': len(x), 'global_ca_rmsd_A': round(global_rmsd, 4),
                'pocket_ca_count': len(pocket_holo), 'candidate_assembly_qc': assembly_qc, 'holo_system_qc': holo_qc,
                'limits': ['SIFTS is a live chain/segment mapping; residue mapping uses deposited-sequence alignment.',
                           'Only model 1 and one selected receptor chain are inspected.',
                           'No occupancy or formal Pair conclusion is possible for this site.']}
    local_rotation, local_translation, local_rmsd = align_coordinates(np.array(pocket_holo), np.array(pocket_y))
    local_tree = cKDTree(hlig @ local_rotation + local_translation)
    local_receptor_clashes = []
    for num, (author, res) in sorted(cres.items()):
        atoms = heavy(res, cop)
        if not atoms:
            continue
        distances = local_tree.query(np.array([a[1] for a in atoms]))[0]
        if float(np.min(distances)) <= 2.0:
            local_receptor_clashes.append({'label_seq': num, 'auth_chain': author, 'auth_seq': str(res.seqid),
                                            'ccd': res.name, 'min_distance_A': round(float(np.min(distances)), 3),
                                            'atoms_within_2A': int(np.sum(distances <= 2.0))})
    nearby = []
    for instance, (op, label, rr) in cinstances.items():
        if instance == ckey:
            continue
        for author, res in rr:
            atoms = heavy(res, op)
            if not atoms:
                continue
            distances = target_tree.query(np.array([a[1] for a in atoms]))[0]
            local_distances = local_tree.query(np.array([a[1] for a in atoms]))[0]
            minimum = float(np.min(distances))
            local_minimum = float(np.min(local_distances))
            if min(minimum, local_minimum) <= 4.0:
                kind = ('water' if res.is_water() else 'polymer' if res.entity_type == gemmi.EntityType.Polymer
                        else 'single_atom' if len(atoms) == 1 else 'nonpolymer')
                nearby.append({'instance': instance, 'label_asym': label, 'author_chain': author,
                               'author_seq': str(res.seqid), 'ccd': res.name, 'kind': kind,
                               'min_distance_to_transformed_holo_ligand_A': round(minimum, 3),
                               'atoms_within_3A': int(np.sum(distances <= 3.0)),
                               'local_probe_min_distance_A': round(local_minimum, 3),
                               'local_probe_atoms_within_3A': int(np.sum(local_distances <= 3.0)),
                               'altlocs': sorted({a[2] for a in atoms})})
    nearby.sort(key=lambda x: x['min_distance_to_transformed_holo_ligand_A'])
    overlap = [v for v in nearby if v['kind'] != 'water' and v['atoms_within_3A'] > 0]
    local_overlap = [v for v in nearby if v['kind'] != 'water' and v['local_probe_atoms_within_3A'] > 0]
    global_components = {(v['instance'], v['kind']) for v in overlap}
    local_components = {(v['instance'], v['kind']) for v in local_overlap}
    complete = all(r.get('candidate_has_CA') and r.get('same_residue') and
                   not r.get('holo_heavy_atoms_absent_in_candidate') for r in residues)
    if local_rmsd > 2.0:
        state = 'POCKET_ALIGNMENT_UNCERTAIN_REVIEW'
    elif global_components != local_components:
        state = 'ALIGNMENT_DEPENDENT_OCCUPANCY_REVIEW'
    elif overlap:
        state = 'OTHER_OR_TARGET_COMPONENT_NEAR_TEMPLATE'
    elif receptor_template_clashes or local_receptor_clashes:
        state = 'RECEPTOR_TEMPLATE_OVERLAP_REQUIRES_REVIEW'
    elif complete:
        state = 'NO_NONWATER_COMPONENT_NEAR_TEMPLATE_PROVISIONAL'
    else:
        state = 'SITE_OBSERVATION_OR_MAPPING_INSUFFICIENT'
    assembly_qc, holo_qc = materialize(case, cstruct, assembly_id, result_dir)
    return {'sample': case['id'], 'holo_status': case['holo_status'], 'holo_source': {
                'uri': str(hpath), 'revision': hrev, 'sha256': case['identity']['source_sha256'],
                'first_public_model_date': hdate, 'assembly': case['identity']['assembly'], 'model': hstruct[0].num},
            'candidate_pdb': cand_pdb, 'candidate_source': csource, 'candidate_assembly': assembly_id,
            'holo_binding_instance': hkey, 'candidate_binding_instance': ckey,
            'candidate_author_chain': cand_auth, 'candidate_label_asym': clabel,
            'sifts_common_uniprot': common, 'holo_sifts_mappings': hmap, 'candidate_sifts_mappings': cmap,
            'sifts_request': {'holo': hmi, 'candidate': cmi}, 'sequence_alignment': alignment,
            'target_residue_instances': targets, 'target_heavy_atoms': len(lig_atoms),
            'site_definition': 'Holo receptor residues with a positive-occupancy heavy atom within 4.0 A of observed target',
            'site_residues': residues, 'candidate_components_within_4A_of_transformed_holo_target': nearby,
            'candidate_receptor_atoms_within_2A_of_transformed_holo_target': receptor_template_clashes,
            'candidate_receptor_atoms_within_2A_of_local_probe': local_receptor_clashes,
            'geometric_screen': state, 'pocket_all_holo_observed_heavy_atoms_present': complete,
            'site_observation_status': ('COMPLETE_RELATIVE_TO_HOLO_OBSERVED_HEAVY_ATOMS' if complete
                                        else 'INSUFFICIENT_OR_UNMAPPED'),
            'alignment_ca_count': len(x), 'global_ca_rmsd_A': round(global_rmsd, 4),
            'pocket_ca_rmsd_after_global_alignment_A': round(rmsd(np.array(pocket_x), np.array(pocket_y)), 4)
                if pocket_x else None,
            'pocket_ca_rmsd_after_local_alignment_A': round(local_rmsd, 4),
            'local_rmsd_review_threshold_A': 2.0,
            'pocket_ca_count': len(pocket_x), 'candidate_assembly_qc': assembly_qc,
            'holo_system_qc': holo_qc,
            'limits': ['SIFTS API supplies live chain/segment mapping; residue mapping is the explicitly recorded sequence alignment.',
                       'Only model 1 and one explicitly selected receptor chain are inspected.',
                       'The superposed target is a geometric probe, not a ligand fit or proof of an empty site.',
                       '2.0 A local RMSD is a provisional review trigger, not a pair acceptance or exclusion threshold.',
                       'Alternative conformers and local rearrangements require structural review.',
                       'Holo target-role and final APHCP acceptance are separate decisions.']}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', required=True)
    p.add_argument('--case', required=True)
    p.add_argument('--candidate', required=True)
    p.add_argument('--author-chain', required=True)
    p.add_argument('--candidate-instance', help='explicit Assembly copy, e.g. 1.A')
    p.add_argument('--assembly', default='1')
    p.add_argument('--source', required=True, help='candidate provenance, e.g. AHoJ or PLINDER')
    args = p.parse_args()
    run = Path(args.run)
    cases = [c for c in read_json(run / 'samples.json') if c['id'] == args.case]
    if len(cases) != 1:
        p.error('case must identify exactly one selected sample')
    name = args.case.replace(':', '_') + '__' + args.candidate.upper() + '_' + args.author_chain
    if args.candidate_instance:
        name += '__copy_' + args.candidate_instance
    base = run / 'pairs' / name
    base.mkdir(parents=True, exist_ok=True)
    n = 1
    while (base / f'attempt{n:03d}').exists():
        n += 1
    out = base / f'attempt{n:03d}'
    out.mkdir()
    dependencies = {}
    for code_path in [Path(__file__), Path(apo_c_pilot.__file__), Path(holo_source_mapping.__file__)]:
        digest = sha(code_path)
        code = run / 'code' / f'{code_path.stem}-{digest[:12]}.py'
        if not code.exists():
            shutil.copy2(code_path, code)
        dependencies[code_path.name] = digest
    result = {'started_utc': now(), 'source_route': args.source, 'code_sha256': sha(Path(__file__)),
              'dependency_sha256': dependencies,
              'request': vars(args)}
    try:
        result.update(inspect_pair(cases[0], args.candidate.upper(), args.author_chain, args.assembly, run, out,
                                   args.candidate_instance))
        result['execution_status'] = 'INSPECTED'
    except Exception as exc:
        result.update({'execution_status': 'UNRESOLVED', 'error': f'{type(exc).__name__}: {exc}'})
    result['finished_utc'] = now()
    write_json(out / 'pair.json', result)
    print(json.dumps({'path': str(out), 'status': result['execution_status'],
                      'screen': result.get('geometric_screen'), 'error': result.get('error')}), flush=True)
    if result['execution_status'] != 'INSPECTED':
        sys.exit(1)


if __name__ == '__main__':
    main()
