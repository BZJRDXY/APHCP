"""Prepare a bounded, metadata-only shortlist for a future independent Apo/Holo audit."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random
import shutil
import sqlite3

import gemmi

from apo_c_pilot import now, read_json, sha, write_json


def source_metadata(path, binding_labels=None):
    block = gemmi.cif.read(str(path)).sole_block()
    asym = block.get_mmcif_category('_struct_asym.')
    label_entity = dict(zip(asym.get('id', []), asym.get('entity_id', [])))
    entities = {label_entity[x] for x in binding_labels or [] if x in label_entity}
    refs = block.get_mmcif_category('_struct_ref.')
    accessions = set()
    for i, accession in enumerate(refs.get('pdbx_db_accession', [])):
        database = refs.get('db_name', [])
        if (accession and i < len(database) and str(database[i]).upper() in {'UNP', 'UNIPROT'}
                and (not entities or str(refs.get('entity_id', [])[i]) in entities)):
            accessions.add(str(accession))
    poly = block.get_mmcif_category('_entity_poly.')
    sequences = set()
    for i, entity in enumerate(poly.get('entity_id', [])):
        if entities and str(entity) not in entities:
            continue
        raw = poly.get('pdbx_seq_one_letter_code_can', [])[i]
        if raw:
            seq = ''.join(str(raw).split()).upper()
            sequences.add(hashlib.sha256(seq.encode('ascii')).hexdigest())
    citations = block.get_mmcif_category('_citation.')
    literature = set()
    for i, cid in enumerate(citations.get('id', [])):
        if str(cid).lower() != 'primary':
            continue
        for key, prefix in [('pdbx_database_id_PubMed', 'pmid:'), ('pdbx_database_id_DOI', 'doi:')]:
            values = citations.get(key, [])
            if i < len(values) and values[i] not in (None, False, '', '.', '?'):
                literature.add(prefix + str(values[i]).lower())
    return {'accessions': sorted(accessions), 'sequence_sha256': sorted(sequences), 'primary_citation_ids': sorted(literature)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', required=True)
    p.add_argument('--catalog', required=True)
    p.add_argument('--seed', type=int, default=20261125)
    p.add_argument('--scan-limit', type=int, default=4000)
    p.add_argument('--target-count', type=int, default=30)
    p.add_argument('--multi-target', type=int, default=5)
    p.add_argument('--output-name', default='provisional_nominees.json')
    args = p.parse_args()
    run, catalog = Path(args.run), Path(args.catalog).resolve()
    if Path(args.output_name).name != args.output_name:
        raise ValueError('output name must be a simple filename')
    out = run / 'holdout_preparation' / args.output_name
    if out.exists():
        raise FileExistsError(f'shortlist already exists: {out}')
    baseline = Path(read_json(run / 'run.json')['baseline_run'])
    prior_samples = read_json(baseline / 'samples.json')
    prior_paths = {Path(x['source_uri']) for x in prior_samples}
    prior_pdbs = {x['pdb'] for x in prior_samples}
    for folder in [baseline / 'pairs', baseline / 'multichain']:
        for pth in folder.glob('*/attempt*/pair.json'):
            d = read_json(pth)
            source = d.get('candidate_source') or {}
            if source.get('uri') and Path(source['uri']).is_file():
                prior_paths.add(Path(source['uri']))
                candidate = d.get('candidate_pdb') or d.get('request', {}).get('candidate')
                if candidate:
                    prior_pdbs.add(str(candidate).upper())
    for pth in (baseline / 'raw').glob('*.cif.gz'):
        prior_paths.add(pth)
        prior_pdbs.add(pth.name.split('-v')[0].upper())
    exclude_acc, exclude_seq, exclude_citation = set(), set(), set()
    metadata_failures = []
    for path in sorted(prior_paths):
        try:
            m = source_metadata(path)
            exclude_acc.update(m['accessions'])
            exclude_seq.update(m['sequence_sha256'])
            exclude_citation.update(m['primary_citation_ids'])
        except Exception as exc:
            metadata_failures.append({'path': str(path), 'error': f'{type(exc).__name__}: {exc}'})
    if metadata_failures:
        raise ValueError(f'prior-source metadata incomplete: {metadata_failures[:3]}')
    uri = catalog.as_uri() + '?mode=ro'
    with sqlite3.connect(uri, uri=True) as db:
        rows = db.execute("SELECT object_id,pdb_id,error_message,output_uri FROM run_item "
                          "WHERE stage='STRUCTURE_OBSERVATION' AND status='REGISTERED' AND output_uri IS NOT NULL "
                          "ORDER BY pdb_id,object_id LIMIT ?", (args.scan_limit,)).fetchall()
        dedup = {}
        for oid, pdb, message, system_uri in rows:
            evidence = json.loads(message)
            identity = evidence['identity']
            if pdb in prior_pdbs or len(identity['target_residue_instances']) != 1:
                continue
            multi = len(identity['binding_instances']) > 1
            dedup.setdefault((pdb, multi), (oid, pdb, system_uri, evidence))
        pools = {True: [], False: []}
        for (_, multi), row in dedup.items():
            pools[multi].append(row)
        rng = random.Random(args.seed)
        rng.shuffle(pools[True])
        rng.shuffle(pools[False])
        candidates = [(True, x) for x in pools[True]] + [(False, x) for x in pools[False]]
        selected, seen_pdb, seen_acc, seen_seq, seen_citation = [], set(), set(), set(), set()
        skipped = Counter()
        skipped_details = []
        for multi, (oid, pdb, system_uri, evidence) in candidates:
            if len(selected) >= args.target_count:
                break
            if multi and sum(x['multi_chain'] for x in selected) >= args.multi_target:
                continue
            if pdb in seen_pdb:
                skipped['duplicate_pdb'] += 1
                continue
            ident = evidence['identity']
            assets = db.execute('SELECT a.uri,a.sha256 FROM asset_file a JOIN entry_version v ON v.id=a.entry_version_id '
                                "WHERE v.pdb_id=? AND v.revision=? AND a.sha256=? AND a.uri LIKE '%.cif.gz'",
                                (pdb, ident['revision'], ident['source_sha256'])).fetchall()
            paths = [Path(x[0]) for x in assets if Path(x[0]).is_file()]
            if not paths or not Path(system_uri).is_file():
                skipped['fixed_source_or_system_missing'] += 1
                continue
            labels = {x[1] for x in ident['binding_instances']}
            try:
                m = source_metadata(paths[0], labels)
            except Exception as exc:
                skipped['source_metadata_unreadable'] += 1
                skipped_details.append({'pdb': pdb, 'source_uri': str(paths[0]),
                                        'error': f'{type(exc).__name__}: {exc}'})
                continue
            acc, seq, cit = set(m['accessions']), set(m['sequence_sha256']), set(m['primary_citation_ids'])
            if not acc or not seq or not cit:
                skipped['missing_protein_sequence_or_primary_citation'] += 1
                continue
            if acc & (exclude_acc | seen_acc) or seq & (exclude_seq | seen_seq):
                skipped['protein_or_exact_sequence_overlap'] += 1
                continue
            if cit & (exclude_citation | seen_citation):
                skipped['primary_citation_overlap'] += 1
                continue
            if sha(paths[0]) != ident['source_sha256']:
                raise ValueError(f'fixed source checksum differs: {paths[0]}')
            selected.append({'system_id': oid, 'pdb': pdb, 'multi_chain': multi,
                             'holo_identity': ident, 'source_uri': str(paths[0]), 'system_uri': system_uri,
                             'protein_and_citation_metadata': m, 'holo_status': 'REGISTERED'})
            seen_pdb.add(pdb)
            seen_acc.update(acc)
            seen_seq.update(seq)
            seen_citation.update(cit)
    result = {'created_utc': now(), 'status': 'PROVISIONAL_NOT_FROZEN', 'seed': args.seed,
              'catalog': str(catalog), 'selection_scope': 'First ordered bounded catalog rows, shuffled with the recorded seed; not a representative sample.',
              'scan_limit': args.scan_limit, 'rows_read': len(rows), 'candidate_pools': {str(k): len(v) for k, v in pools.items()},
              'target_count': args.target_count, 'requested_multi_chain': args.multi_target,
              'prior_entries_excluded': len(prior_pdbs), 'prior_accessions_excluded': len(exclude_acc),
              'prior_sequence_hashes_excluded': len(exclude_seq), 'prior_primary_citation_identifiers_excluded': len(exclude_citation),
              'skipped': dict(skipped), 'skipped_details': skipped_details, 'selected': selected,
              'limits': ['No Apo candidates or labels were searched or assessed for these nominees.',
                         'Only exact accessions, exact sequences and recorded primary DOI/PMID overlap were excluded.',
                         'Publication aliases, related experiments and structural-family relationships still require independent review.',
                         'Pair-level coverage of conformational motion, missing loops and competing ligands is not established.',
                         'Freeze sample, candidate budget, rules and independent audit plan before any validation run.']}
    write_json(out, result)
    code = run / 'code' / f'{Path(__file__).stem}-{sha(Path(__file__))[:12]}.py'
    if not code.exists():
        shutil.copy2(Path(__file__), code)
    print(json.dumps({'path': str(out), 'selected': len(selected),
                      'multi_chain': sum(x['multi_chain'] for x in selected), 'skipped': dict(skipped)}), flush=True)


if __name__ == '__main__':
    main()
