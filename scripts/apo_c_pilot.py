"""Small, versioned M6/M7 pilot. Inputs and data remain outside the source tree."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import html
import io
import json
from pathlib import Path
import re
import shutil
import sqlite3
import sys
import time

import gemmi
import pyarrow.parquet as pq
import requests


VERSION = 'apo-c-pilot-0.1'


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f'refusing to overwrite evidence: {path}')
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf8')


def read_json(path):
    return json.loads(path.read_text(encoding='utf8'))


def cache_get(run, name, url, params=None, timeout=25, accepted=(200, 204), validate=None):
    """Cache successful responses and their request evidence; failures remain explicit."""
    base = run / 'http' / name
    body, meta = base.with_suffix('.gz'), base.with_suffix('.json')
    if body.exists() and meta.exists():
        with gzip.open(body, 'rb') as f:
            payload, info = f.read(), read_json(meta)
        if info['status'] not in accepted:
            raise ValueError(f'unexpected cached HTTP status: {base}')
        if validate:
            validate(payload, info['status'])
        return payload, info
    if body.exists() or meta.exists():
        raise RuntimeError(f'incomplete cache entry: {base}')
    started, clock = now(), time.monotonic()
    r = requests.get(url, params=params, timeout=(8, timeout), headers={'User-Agent': 'APHCP-apo-research-pilot/0.1'})
    elapsed = time.monotonic() - clock
    if r.status_code not in accepted:
        raise RuntimeError(f'HTTP {r.status_code} for {r.url}: {r.text[:180]}')
    if validate:
        validate(r.content, r.status_code)
    base.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(body, 'xb') as f:
        f.write(r.content)
    info = {'url': r.url, 'status': r.status_code, 'started_utc': started, 'elapsed_s': round(elapsed, 3),
            'response_bytes': len(r.content), 'response_sha256': hashlib.sha256(r.content).hexdigest(),
            'content_type': r.headers.get('content-type')}
    write_json(meta, info)
    return r.content, info


def open_catalog(path):
    return sqlite3.connect(path.resolve().as_uri() + '?mode=ro', uri=True)


def initial_samples(args):
    sample_path, run = Path(args.samples), Path(args.run)
    specification = read_json(sample_path)
    catalog = Path(args.catalog).resolve()
    if run.exists():
        raise FileExistsError(f'run already exists: {run}')
    rows = []
    with open_catalog(catalog) as db:
        for spec in specification['samples']:
            found = db.execute("SELECT pdb_id,status,output_uri,error_message FROM run_item "
                               "WHERE object_id=? AND stage='STRUCTURE_OBSERVATION' "
                               "AND output_uri IS NOT NULL ORDER BY updated_at DESC,run_id DESC LIMIT 1",
                               (spec['id'],)).fetchone()
            if not found:
                raise ValueError(f'sample absent from sealed catalog: {spec["id"]}')
            pdb, status, system_path, error = found
            evidence = json.loads(error)
            ident = evidence['identity']
            assets = db.execute('SELECT a.uri,a.sha256,a.bytes FROM asset_file a '
                                'JOIN entry_version v ON v.id=a.entry_version_id '
                                "WHERE v.pdb_id=? AND v.revision=? AND a.sha256=? AND a.uri LIKE '%.cif.gz'",
                                (pdb, ident['revision'], ident['source_sha256'])).fetchall()
            assets = [a for a in assets if Path(a[0]).is_file()]
            source_record = 'sealed_catalog_asset_file'
            if not assets:
                # Some inherited observations retain a source proof in their original Run
                # even though asset_file has no matching row. Require that proof and hash.
                original_run = Path(system_path).parents[3]
                proof_path = original_run / 'source_proofs' / f'{pdb}.json'
                if proof_path.is_file():
                    proof = read_json(proof_path)
                    pident = proof.get('identity', {})
                    candidate = (original_run / proof.get('relative_path', '')).resolve()
                    if (proof.get('status') == 'PASS' and proof.get('gzip_crc') == 'PASS'
                            and pident.get('pdb_id') == pdb and pident.get('revision') == ident['revision']
                            and proof.get('sha256') == ident['source_sha256']
                            and candidate.is_relative_to(original_run.resolve()) and candidate.is_file()):
                        assets = [(str(candidate), proof['sha256'], candidate.stat().st_size)]
                        source_record = str(proof_path)
            if not assets or not Path(system_path).is_file():
                raise ValueError(f'missing fixed source or system for {spec["id"]}')
            raw = assets[0]
            if sha(Path(raw[0])) != ident['source_sha256']:
                raise ValueError(f'raw source checksum differs for {spec["id"]}')
            rows.append({'id': spec['id'], 'split': spec['split'], 'question': spec['question'],
                         'pdb': pdb, 'holo_status': status, 'identity': ident, 'qc': evidence.get('qc'),
                         'system_uri': system_path, 'source_uri': raw[0], 'source_bytes': raw[2],
                         'source_record': source_record})
    run.mkdir(parents=True)
    (run / 'code').mkdir()
    shutil.copy2(Path(__file__), run / 'code' / f'{Path(__file__).stem}-{sha(Path(__file__))[:12]}.py')
    write_json(run / 'sample_spec.json', specification)
    write_json(run / 'samples.json', rows)
    write_json(run / 'run.json', {'version': VERSION, 'started_utc': now(), 'catalog': str(catalog),
                                 'index': str(Path(args.index).resolve()), 'samples': len(rows),
                                 'holdout': sum(r['split'] == 'holdout' for r in rows),
                                 'code_sha256': sha(Path(__file__)), 'sample_spec_sha256': sha(sample_path),
                                 'principle': 'Research candidate evidence only; no formal APHCP Pair acceptance.'})
    print(json.dumps({'run': str(run), 'samples': len(rows), 'holdout': sum(r['split'] == 'holdout' for r in rows)}))


def plinder(case, index):
    pid = case['pdb'].lower()
    revisions = pq.read_table(index / 'entry_sources.parquet', filters=[('entry_pdb_id', '=', pid)]).to_pylist()
    columns = ['entry_pdb_id', 'system_id', 'system_biounit_id', 'ligand_id', 'ligand_ccd_code',
               'ligand_instance_chains', 'ligand_protein_chains_asym_id', 'ligand_member_asym_ids',
               'ligand_auth_id', 'ligand_residue_numbers']
    ann = pq.read_table(index / 'annotation_table.parquet', columns=columns,
                        filters=[('entry_pdb_id', '=', pid)]).to_pylist()
    ident = case['identity']
    targets = {x[0] for x in ident['target_residue_instances']}
    binding = {x[0] for x in ident['binding_instances']}
    rev = (f"{revisions[0]['source_mmcif_major_revision']}.{revisions[0]['source_mmcif_minor_revision']}"
           if revisions else None)
    matched, chain_hints = [], []
    if ann:
        import holo_source_mapping as mapping
        block = gemmi.cif.read(case['source_uri']).sole_block()
        structure = gemmi.make_structure_from_block(block)
        instances = mapping.assembly_instances(structure, str(ident['assembly']), structure[0])
        for row in ann:
            if (str(row['system_biounit_id']) != str(ident['assembly'])
                    or set(row['ligand_protein_chains_asym_id'] or []) != binding):
                continue
            chains = set(row['ligand_instance_chains'] or [])
            if chains == targets:
                granularity = plinder_instance_match(row, ident, instances, rev, block)
                if granularity == 'INSTANCE':
                    matched.append(row | {'match_granularity': granularity})
                elif granularity == 'CHAIN_HINT':
                    chain_hints.append(row | {'match_granularity': granularity})
    selected = matched or chain_hints
    systems = sorted({r['system_id'] for r in selected})
    links = (pq.read_table(index / 'linked_apo_structures.parquet',
                           filters=[('reference_system_id', 'in', systems)]).to_pylist() if systems else [])
    for link in links:
        link['target_match_granularity'] = 'INSTANCE' if matched else 'CHAIN_HINT'
    status = ('no_entry_coverage' if not revisions else 'entry_without_system_annotation' if not ann
              else 'no_matching_target_instance' if not selected
              else 'chain_level_candidate_hints' if not matched and links
              else 'candidates' if links else 'matched_no_linked_candidates')
    return {'status': status, 'release': 'PLINDER 2026-09', 'holo_revision': ident['revision'],
            'source_revision': rev, 'same_revision': rev == ident['revision'] if rev else None,
            'matched_annotations': matched, 'chain_level_hints': chain_hints, 'links': links}


def plinder_instance_match(row, identity, instances, source_revision, block):
    """Resolve author residue identity from the matching PLINDER assembly copy."""
    targets = identity['target_residue_instances']
    if set(row['ligand_member_asym_ids'] or []) != {t[1] for t in targets}:
        return None
    if len(targets) != 1:
        return 'CHAIN_HINT'
    instance, label, _, entity, author, native, ccd = targets[0]
    if (str(row['ligand_ccd_code']).upper() != ccd.upper() or row.get('ligand_auth_id') != author):
        return None
    if source_revision != identity['revision'] or instance not in instances:
        return 'CHAIN_HINT'
    _, observed_label, residues = instances[instance]
    if observed_label != label:
        return None
    matching = [r for ach, r in residues if ach == author and r.name == ccd and
                str(r.entity_id) == str(entity) and not r.is_water()]
    verified = [r for r in matching if str(r.seqid) == str(native)]
    if len(verified) != 1:
        return None
    scheme = block.get_mmcif_category('_pdbx_nonpoly_scheme.')
    keys = ['asym_id', 'entity_id', 'mon_id', 'pdb_strand_id', 'pdb_seq_num', 'pdb_ins_code', 'ndb_seq_num']
    if not all(k in scheme for k in keys):
        return 'CHAIN_HINT'
    residue = verified[0]
    native_rows = [ord_ for lab, ent, comp, ach, number, ins, ord_ in zip(*(scheme[k] for k in keys))
                   if (lab, str(ent), comp, ach, str(number), str(ins or '').strip()) ==
                   (label, str(entity), ccd, author, str(residue.seqid.num), residue.seqid.icode.strip())]
    if len(native_rows) == 1 and row['ligand_residue_numbers']:
        return 'INSTANCE' if {str(v) for v in row['ligand_residue_numbers']} == set(native_rows) else None
    return 'CHAIN_HINT'


def uniprot_mappings(run, pdb):
    body, info = cache_get(run, f'sifts/{pdb.lower()}',
                           f'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb.lower()}',
                           accepted=(200, 204, 404))
    if info['status'] == 404:
        if json.loads(body).get('message') != 'Requested endpoint does not contain any data':
            raise RuntimeError(f'unexpected PDBe 404 response for {pdb}')
        return {}, info
    if info['status'] == 204:
        return {}, info
    data = json.loads(body).get(pdb.lower(), {})
    return data.get('UniProt', {}), info


def relevant_uniprot(case, mappings):
    labels = {x[1] for x in case['identity']['binding_instances']}
    return {acc: {'name': v.get('name'), 'mappings': [m for m in v.get('mappings', [])
                                                       if m.get('struct_asym_id') in labels]}
            for acc, v in mappings.items()
            if any(m.get('struct_asym_id') in labels for m in v.get('mappings', []))}


class RcsbPageError(RuntimeError):
    def __init__(self, name, offset, expected, hits, pages, exc):
        super().__init__(f'{name} page {offset}: {exc}')
        self.progress = {'query_id': name, 'missing_offset': offset, 'total_count': expected,
                         'retained_hits': len(hits), 'pages': pages}


def rcsb_query(run, name, node):
    all_hits, pages = [], []
    start, total = 0, None
    seen = set()
    while True:
        query = {'query': node,
            'return_type': 'polymer_entity', 'request_options': {
                'results_content_type': ['experimental'], 'paginate': {'start': start, 'rows': 1000}}}
        def validate(body, status):
            if status == 204:
                if start:
                    raise ValueError('204 on a later page')
                return
            record = json.loads(body)
            count, items = record['total_count'], record.get('result_set', [])
            if not isinstance(count, int) or count < start or (total is not None and total != count):
                raise ValueError('result total changed or is invalid')
            if len(items) != min(1000, count - start):
                raise ValueError('missing or excess entries in page')
            identifiers = {x['identifier'] for x in items}
            if len(identifiers) != len(items) or identifiers & seen:
                raise ValueError('duplicate entities within or across pages')

        try:
            body, info = cache_get(run, f'rcsb/{name}_{start}', 'https://search.rcsb.org/rcsbsearch/v2/query',
                                   params={'json': json.dumps(query, separators=(',', ':'))}, validate=validate)
        except Exception as exc:
            raise RcsbPageError(name, start, total, all_hits, pages, exc) from exc
        pages.append(info)
        if info['status'] == 204:
            return {'query_id': name, 'total_count': 0, 'hits': [], 'pages': pages, 'query': query}
        result = json.loads(body)
        hits = result.get('result_set', [])
        all_hits.extend(hits)
        total = result['total_count']
        seen.update(x['identifier'] for x in hits)
        if len(all_hits) == total:
            break
        if not hits:
            raise RuntimeError(f'RCSB pagination ended early at {len(all_hits)}/{total}')
        start += len(hits)
    identifiers = [h['identifier'] for h in all_hits]
    if len(set(identifiers)) != len(identifiers):
        raise RuntimeError(f'RCSB pagination returned duplicate entities for {name}')
    return {'query_id': name, 'total_count': total, 'hits': all_hits, 'pages': pages, 'query': query}


def rcsb(run, accession):
    node = {'type': 'terminal', 'service': 'text', 'parameters': {
        'attribute': 'rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession',
        'operator': 'exact_match', 'value': accession}}
    return {'accession': accession, **rcsb_query(run, accession, node)}


def rcsb_sequence(run, sequence):
    digest = hashlib.sha256(sequence.encode('ascii')).hexdigest()
    node = {'type': 'terminal', 'service': 'sequence', 'parameters': {
        'sequence_type': 'protein', 'value': sequence, 'identity_cutoff': 0.9, 'evalue_cutoff': 1}}
    return {'sequence_sha256': digest, 'identity_cutoff': 0.9, 'evalue_cutoff': 1,
            **rcsb_query(run, 'seq_' + digest, node)}


def source_sequences(case, labels):
    structure = gemmi.read_structure(case['source_uri'])
    entities = {str(e.name): e for e in structure.entities}
    found = {}
    for label in labels:
        entity_ids = {str(r.entity_id) for c in structure[0] for r in c
                      if r.subchain == label and r.entity_type == gemmi.EntityType.Polymer}
        if len(entity_ids) != 1:
            raise ValueError(f'no unique source polymer for {label}: {entity_ids}')
        entity = entities[next(iter(entity_ids))]
        sequence = ''.join((gemmi.find_tabulated_residue(r).one_letter_code or 'X').upper()
                           for r in entity.full_sequence).replace(' ', 'X')
        if not sequence:
            raise ValueError(f'no deposited sequence for {label}')
        found.setdefault(sequence, []).append(label)
    return found


def parse_ahoj_table(page, kind):
    m = re.search(r'<table\s+id="' + re.escape(kind) + r'-results"[^>]*>(.*?)</table>', page, re.S)
    if not m:
        return []
    out = []
    for row in re.findall(r'<tr\s+data-row-id="[^"]+"[^>]*>.*?</tr>', m.group(1), re.S):
        p = re.search(r'/structure/([a-zA-Z0-9]{4,})', row)
        chain = re.search(r'class="chain-id"[^>]*>(.*?)</span>', row, re.S)
        if not p or not chain:
            raise ValueError('AHoJ result row missing entry or chain')
        cols = [html.unescape(re.sub('<[^>]+>', '', c)).strip() for c in re.findall(r'<td(?:\s[^>]*)?>(.*?)</td>', row, re.S)]
        out.append({'pdb': p.group(1).upper(), 'auth_chain': html.unescape(chain.group(1)).strip(),
                    'source_order': len(out) + 1, 'columns': cols})
    return out


def ahoj(case, run):
    pdb = case['pdb'].lower()
    body, info = cache_get(run, f'ahoj/{pdb}_summary', 'https://apoholo.cz/db/search-results/export-csv',
                           params={'pdb_ids': pdb})
    records = list(csv.DictReader(io.StringIO(body.decode('utf-8-sig')), skipinitialspace=True))
    t = case['identity']['target_residue_instances']
    if len(t) != 1:
        return {'status': 'complex_target_requires_instance_mapping', 'summary': info, 'records': records}
    author_chain, author_seq, ccd = str(t[0][4]), str(t[0][5]), str(t[0][6])
    suffix = f'-{ccd}-{author_seq}'.lower()
    matches = [r for r in records if (r.get('entry_key') or '').lower().endswith(suffix)
               and (author_chain in re.split(r'[ ,;]+', r.get('target_chains') or ''))]
    details = []
    for r in matches:
        key = r['entry_key']
        page, page_info = cache_get(run, f'ahoj/{key}', f'https://apoholo.cz/db/entry/{key}', timeout=40)
        text = page.decode('utf8')
        apo_rows, holo_rows = parse_ahoj_table(text, 'apo'), parse_ahoj_table(text, 'holo')
        if (len(apo_rows) != int(r['num_apo_sites']) or len(holo_rows) != int(r['num_holo_sites'])):
            raise ValueError(f'AHoJ rendered rows do not match summary counts for {key}')
        details.append({'entry': r, 'page': page_info, 'apo': apo_rows, 'holo': holo_rows,
                        'unobserved': parse_ahoj_table(text, 'unobserved')})
    return {'status': 'matched_target' if details else 'no_matching_target_instance',
            'version': 'live AHoJ-DB; exact correspondence to local v2c archive not yet verified',
            'summary': info, 'matches': details}


def case_file(run, case):
    return run / 'search' / (case['id'].replace(':', '_') + '.json')


def attempts(path):
    return ([path] if path.exists() else []) + sorted(path.parent.glob(path.stem + '.retry*.json'))


def next_attempt(path):
    prior = attempts(path)
    return path if not prior else path.with_name(path.stem + f'.retry{len(prior):06d}.json')


def search_cases(args):
    run = Path(args.run)
    cases = read_json(run / 'samples.json')
    selected = [x for x in cases if not args.case or x['id'] in args.case]
    if args.split:
        selected = [x for x in selected if x['split'] == args.split]
    index = Path(read_json(run / 'run.json')['index'])
    code_dir = run / 'code'
    code_snapshot = code_dir / f'{Path(__file__).stem}-{sha(Path(__file__))[:12]}.py'
    if not code_snapshot.exists():
        shutil.copy2(Path(__file__), code_snapshot)
    for case in selected:
        path = case_file(run, case)
        previous = attempts(path)
        latest = read_json(previous[-1]) if previous else {}
        if latest:
            if not args.retry_failed or not any(latest[k]['status'] == 'technical_failure'
                                                for k in ('plinder', 'sifts_rcsb', 'ahoj')):
                print('REUSE', case['id'], flush=True)
                continue
            path = next_attempt(path)
        result = {'sample': case['id'], 'started_utc': now(), 'code_sha256': sha(Path(__file__))}
        for name in ('plinder', 'sifts_rcsb', 'ahoj'):
            if name in latest and latest[name]['status'] != 'technical_failure':
                result[name] = latest[name]
                continue
            try:
                if name == 'plinder':
                    result[name] = plinder(case, index)
                elif name == 'sifts_rcsb':
                    m, info = uniprot_mappings(run, case['pdb'])
                    relevant = relevant_uniprot(case, m)
                    labels = {x[1] for x in case['identity']['binding_instances']}
                    covered = {row['struct_asym_id'] for v in relevant.values() for row in v['mappings']}
                    missing = labels - covered
                    seqs = source_sequences(case, sorted(missing)) if missing else {}
                    searches = [rcsb(run, acc) for acc in sorted(relevant)]
                    searches += [{'source_chain_labels': chain_labels, **rcsb_sequence(run, seq)}
                                 for seq, chain_labels in seqs.items()]
                    result[name] = {'status': 'searched_with_sequence_for_unmapped_chains' if missing else 'searched',
                                    'sifts_info': info, 'mappings': relevant, 'unmapped_chain_labels': sorted(missing),
                                    'searches': searches}
                else:
                    result[name] = ahoj(case, run)
            except RcsbPageError as exc:
                result[name] = {'status': 'technical_failure', 'error': str(exc), 'pagination': exc.progress}
            except Exception as exc:
                result[name] = {'status': 'technical_failure', 'error': f'{type(exc).__name__}: {exc}'}
        sets = {
            'plinder': {x['source_entry_id'].upper() for x in result['plinder'].get('links', [])},
            'ahoj': {x['pdb'] for m in result['ahoj'].get('matches', []) for x in m['apo']},
            'rcsb': {h['identifier'].split('_')[0].upper() for s in result['sifts_rcsb'].get('searches', [])
                     for h in s['hits']}}
        result['entry_overlap'] = {k: len(v) for k, v in sets.items()}
        result['entry_overlap']['all_three'] = len(sets['plinder'] & sets['ahoj'] & sets['rcsb'])
        result['entry_overlap']['union'] = len(set.union(*sets.values()))
        result['finished_utc'] = now()
        write_json(path, result)
        print('SEARCH', case['id'], {k: result[k]['status'] for k in ('plinder', 'sifts_rcsb', 'ahoj')},
              result['entry_overlap'], flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('stage', choices=['init', 'search'])
    p.add_argument('--run', required=True)
    p.add_argument('--catalog')
    p.add_argument('--index')
    p.add_argument('--samples')
    p.add_argument('--case', action='append')
    p.add_argument('--split', choices=['development', 'holdout'])
    p.add_argument('--retry-failed', action='store_true')
    args = p.parse_args()
    if args.stage == 'init':
        if not args.catalog or not args.index or not args.samples:
            p.error('init requires --catalog, --index and --samples')
        initial_samples(args)
    else:
        search_cases(args)


if __name__ == '__main__':
    main()
