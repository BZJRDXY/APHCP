"""Read existing fixed originals for a frozen candidate union; do not infer Apo from proximity."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sqlite3
import time

import gemmi

import apo_c_pair
import apo_c_pilot
import apo_c_triplet_regression
import holo_source_mapping
from apo_c_pair import anchored_end_alignment, get_sequence, reference_positions, revision, sequence_map
from apo_c_pilot import now, read_json, sha, write_json
from apo_c_triplet_regression import (component_contacts, deposited_metadata, positive_heavy,
                                      residue_index, unobserved_records)
from holo_source_mapping import assembly_instances


HOLO_CASES = {'3ABX': '3ABX_3A64', '2SRC': '2SRC_1FMK', '1AKE': '1AKE_4AKE'}


def accessions(block, entity):
    cat = block.get_mmcif_category('_struct_ref.')
    return sorted({str(cat['pdbx_db_accession'][i]) for i, db in enumerate(cat.get('db_name', []))
                   if str(db).upper() in {'UNP', 'UNIPROT'} and cat['pdbx_db_accession'][i]
                   and str(cat['entity_id'][i]) == str(entity)})


def existing_original(db, pdb, required_revision=None):
    rows = db.execute('SELECT a.uri,a.sha256,a.bytes,v.revision FROM asset_file a '
                      'JOIN entry_version v ON a.entry_version_id=v.id '
                      "WHERE v.pdb_id=? AND a.uri LIKE '%/raw/%.cif.gz'", (pdb,)).fetchall()
    rows = [x for x in rows if x[1] and Path(x[0]).is_file() and
            (required_revision is None or x[3] == required_revision)]
    if not rows:
        return None
    chosen = max(rows, key=lambda x: (tuple(map(int, x[3].split('.'))), x[0]))
    if sha(Path(chosen[0])) != chosen[1]:
        raise ValueError(f'candidate checksum differs: {chosen[0]}')
    return chosen


def candidate_evidence(run, holo, pdb, source_rows, q, qblock, qstruct, qseq, qsite, db,
                       fixed_source=None, fixed_origin=None):
    source = fixed_source or existing_original(db, pdb)
    source_origin = fixed_origin or 'sealed_catalog_asset_file'
    supplemental = run / 'rank_comparison' / 'supplemental_raw' / f'{pdb}.json'
    if not source and supplemental.is_file():
        manifest = read_json(supplemental)
        if sha(Path(manifest['path'])) != manifest['sha256']:
            raise ValueError(f'supplemental source checksum differs: {pdb}')
        source = (manifest['path'], manifest['sha256'], manifest['bytes'], manifest['revision'])
        source_origin = 'wwpdb_fixed_version_download'
    if not source:
        return {'status': 'SOURCE_NOT_IN_SEALED_ASSETS', 'pdb': pdb,
                'reason': 'Existing fixed original not found; no download attempted in this pass.'}
    path, digest, size, declared_revision = source
    block = gemmi.cif.read(path).sole_block()
    actual_revision, first_date = revision(block, pdb)
    if actual_revision != declared_revision:
        raise ValueError(f'candidate {pdb} revision differs from catalog')
    structure = gemmi.make_structure_from_block(block)
    entity_ids = {row['entity'].split('_', 1)[1] for row in source_rows if row.get('entity')}
    if not entity_ids:
        author_chains = {row['author_chain'] for row in source_rows if row.get('author_chain')}
        entity_ids = {str(res.entity_id) for ch in structure[0] for res in ch
                      if ch.name in author_chains and res.entity_type == gemmi.EntityType.Polymer}
    entity_ids = sorted(entity_ids)
    if not entity_ids:
        return {'status': 'SOURCE_ENTITY_NOT_RESOLVED', 'pdb': pdb, 'source': path}
    output = {'status': 'SOURCE_AND_ENTITY_READ', 'pdb': pdb,
              'candidate_source': {'path': path, 'sha256': digest, 'bytes': size,
                                   'revision': actual_revision, 'first_public_model_date': first_date,
                                   'origin': source_origin},
              'candidate_title': gemmi.cif.as_string(block.find_value('_struct.title')),
              'candidate_deposited_metadata': deposited_metadata(block),
              'model_inspected': structure[0].num, 'model_count': len(structure), 'entities': []}
    for entity_id in entity_ids:
        entity = next((x for x in structure.entities if str(x.name) == entity_id), None)
        if entity is None or not entity.full_sequence:
            output['entities'].append({'entity_id': entity_id, 'status': 'DEPOSITED_SEQUENCE_MISSING'})
            continue
        chain_match = [(ch, res) for ch in structure[0] for res in ch
                       if str(res.entity_id) == entity_id and res.entity_type == gemmi.EntityType.Polymer]
        if not chain_match:
            output['entities'].append({'entity_id': entity_id, 'status': 'ENTITY_COORDINATES_MISSING'})
            continue
        cseq = get_sequence(structure, chain_match)
        cacc = accessions(block, entity_id)
        common = sorted(set(q['sifts_common_uniprot'] or []) & set(cacc))
        mapping, alignment = sequence_map(qseq, cseq)
        if len(common) == 1 and len(qseq) != len(cseq) and not alignment['mismatches']:
            h_authors = {row['author_chain'] for row in qsite.values()}
            c_authors = {chain.name for chain, _ in chain_match}
            if len(h_authors) == 1:
                h_ref = reference_positions(qblock, common[0], next(iter(h_authors)))
                refs = [reference_positions(block, common[0], chain) for chain in sorted(c_authors)]
                if refs and refs[0] and all(ref == refs[0] for ref in refs):
                    mapping, alignment = anchored_end_alignment(qseq, cseq, mapping, alignment, h_ref, refs[0], qsite)
        covered = {x: mapping[x] for x in qsite if x in mapping}
        record = {'entity_id': entity_id, 'source_accessions': cacc, 'common_uniprot': common,
                  'sequence_alignment': alignment, 'holo_site_residues': len(qsite),
                  'site_residues_covered_by_construct_sequence': len(covered), 'assembly_observations': []}
        if not common:
            record['status'] = 'PROTEIN_IDENTITY_UNRESOLVED'
        elif not covered:
            record['status'] = 'NO_TARGET_SITE_SEQUENCE_COVERAGE'
        else:
            record['status'] = 'SITE_SEQUENCE_COVERAGE_REQUIRES_OBSERVATION_CHECK'
        for assembly in structure.assemblies:
            try:
                instances = assembly_instances(structure, assembly.name, structure[0], include_water=True)
            except Exception as exc:
                record['assembly_observations'].append({'assembly': assembly.name, 'status': 'TECHNICAL_FAILURE',
                                                         'error': f'{type(exc).__name__}: {exc}'})
                continue
            for key, (op, label, rows) in instances.items():
                if not any(str(res.entity_id) == entity_id and res.entity_type == gemmi.EntityType.Polymer
                           for _, res in rows):
                    continue
                try:
                    residues = residue_index(rows)
                    unobs = unobserved_records(block, '_pdbx_unobs_or_zero_occ_residues.', label, structure[0].num)
                    site_rows, site_atoms = [], []
                    for hn in sorted(qsite):
                        cn = mapping.get(hn)
                        row = {'holo_label_seq': hn, 'holo_author_seq': qsite[hn]['author_seq'],
                               'holo_component': qsite[hn]['component'], 'candidate_label_seq': cn,
                               'candidate_unobserved_records': unobs.get(cn, [])}
                        if cn in residues:
                            author, residue = residues[cn]
                            atoms = positive_heavy(residue, op)
                            row.update({'author_chain': author, 'author_seq': str(residue.seqid),
                                        'component': residue.name, 'same_residue_component': residue.name == qsite[hn]['component'],
                                        'observed_positive_heavy_atoms': len(atoms),
                                        'missing_relative_holo_heavy_atoms': sorted(qsite[hn]['atom_names'] - {a['name'] for a in atoms})})
                            site_atoms.extend(a | {'label_seq': cn, 'author_chain': author,
                                                   'author_seq': str(residue.seqid), 'residue': residue.name}
                                              for a in atoms)
                        else:
                            row['observed_positive_heavy_atoms'] = 0
                        site_rows.append(row)
                    contacts = component_contacts(instances, {key}, site_atoms, None)
                    record['assembly_observations'].append({'assembly': assembly.name, 'instance': key,
                        'label_asym': label, 'author_chain': rows[0][0], 'model': structure[0].num, 'mapped_site': site_rows,
                        'site_observed_residues': sum(x['observed_positive_heavy_atoms'] > 0 for x in site_rows),
                        'site_complete_relative_holo_residues': sum(x['observed_positive_heavy_atoms'] > 0 and
                            x.get('same_residue_component') and not x.get('missing_relative_holo_heavy_atoms') for x in site_rows),
                        'actual_component_contacts': contacts,
                        'status': 'SITE_UNMAPPED_OR_UNOBSERVED' if not site_atoms else 'ACTUAL_CONTACTS_MEASURED'})
                except Exception as exc:
                    record['assembly_observations'].append({'assembly': assembly.name, 'instance': key,
                                                             'status': 'TECHNICAL_FAILURE',
                                                             'error': f'{type(exc).__name__}: {exc}'})
        if not record['assembly_observations']:
            record['status'] = 'NO_ASSEMBLY_CHAIN_OBSERVATION'
        output['entities'].append(record)
    output['limits'] = ['First deposited model only.', 'Sequence mapping is not a residue-level SIFTS file.',
                        'Actual component proximity alone does not prove site occupancy or unbound state.',
                        'Different Assembly and chain-copy observations remain separate.',
                        'This is material triage, not an independently adjudicated formal Pair.']
    return output


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', required=True)
    p.add_argument('--catalog', required=True)
    p.add_argument('--output-name', default='material_triage')
    p.add_argument('--receipt-name', default='material_triage_receipt.json')
    args = p.parse_args()
    run, catalog = Path(args.run), Path(args.catalog).resolve()
    if Path(args.output_name).name != args.output_name or Path(args.receipt_name).name != args.receipt_name:
        raise ValueError('output and receipt names must be simple names')
    frozen = read_json(run / 'rank_comparison' / 'frozen_candidates.json')
    refs = read_json(run / 'rank_comparison' / 'known_references.json')
    known = {(x['holo'], x['candidate']) for x in refs['references']}
    inputs = {x['id']: x for x in read_json(run / 'inputs.json')}
    output_dir = run / 'rank_comparison' / args.output_name
    output_dir.mkdir(parents=True, exist_ok=True)
    code_paths = [Path(__file__), Path(apo_c_pair.__file__), Path(apo_c_pilot.__file__),
                  Path(apo_c_triplet_regression.__file__), Path(holo_source_mapping.__file__)]
    for path in code_paths:
        dest = run / 'code' / f'{path.stem}-{sha(path)[:12]}.py'
        if not dest.exists():
            dest.write_bytes(path.read_bytes())
    start = time.monotonic()
    with sqlite3.connect(catalog.as_uri() + '?mode=ro', uri=True) as db:
        for case in frozen['cases']:
            holo = case['holo']
            entry = inputs[HOLO_CASES[holo]]
            qpath = Path(entry['sample']['source_uri'])
            if sha(qpath) != entry['sample']['identity']['source_sha256']:
                raise ValueError(f'Holo source changed: {holo}')
            qblock = gemmi.cif.read(str(qpath)).sole_block()
            qstruct = gemmi.make_structure_from_block(qblock)
            hi = assembly_instances(qstruct, str(entry['sample']['identity']['assembly']), qstruct[0])
            hkey = entry['sample']['identity']['binding_instances'][0][0]
            qop, _, qrows = hi[hkey]
            qseq = get_sequence(qstruct, qrows)
            qres = residue_index(qrows)
            qev_path = sorted((run / 'pairs' / HOLO_CASES[holo]).glob('attempt*/evidence.json'))[-1]
            qev = read_json(qev_path)
            qsite = {x['holo_label_seq']: {'author_seq': x['holo_author_seq'],
                     'component': x['holo_component'],
                     'atom_names': {a['name'] for a in positive_heavy(qres[x['holo_label_seq']][1], qop)}}
                     for x in qev['mapped_site']}
            for pdb, origins in case['top10_union'].items():
                if (holo, pdb) in known:
                    continue
                out = output_dir / f'{holo}_{pdb}.json'
                if out.exists():
                    print('TRIAGE_REUSE', holo, pdb, flush=True)
                    continue
                started = now()
                try:
                    rows = [row for source in ('RCSB', 'PLINDER', 'AHoJ')
                            for row in origins.get(source, {}).get('source_rows', [])]
                    d = candidate_evidence(run, holo, pdb, rows, qev, qblock, qstruct, qseq, qsite, db)
                except Exception as exc:
                    d = {'status': 'TECHNICAL_FAILURE', 'pdb': pdb, 'error': f'{type(exc).__name__}: {exc}'}
                d.update({'holo': holo, 'holo_system_id': entry['sample']['id'], 'frozen_origins': origins,
                          'holo_evidence': str(qev_path), 'code_sha256': sha(Path(__file__)),
                          'started_utc': started, 'finished_utc': now()})
                write_json(out, d)
                print('TRIAGE', holo, pdb, d['status'], flush=True)
    summary = {'finished_utc': now(), 'elapsed_s': round(time.monotonic() - start, 3),
               'files': len(list(output_dir.glob('*.json'))), 'new_coordinate_downloads': 0,
               'scope': 'Sealed assets and explicitly pinned supplemental originals; no downloads in this triage step or acceptance inferred from proximity.'}
    write_json(run / 'rank_comparison' / args.receipt_name, summary)
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
