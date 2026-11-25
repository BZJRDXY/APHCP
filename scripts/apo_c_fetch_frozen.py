"""Fetch only missing fixed-revision originals already present in a frozen union."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil

import gemmi
import pyarrow.parquet as pq

from apo_c_pair import revision
from apo_c_pilot import cache_get, now, read_json, sha, write_json


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', required=True)
    p.add_argument('--index', required=True)
    p.add_argument('--pdb', action='append', required=True)
    p.add_argument('--receipt-name', default='supplemental_fetch_receipt.json')
    args = p.parse_args()
    run, index = Path(args.run), Path(args.index)
    if Path(args.receipt_name).name != args.receipt_name:
        raise ValueError('receipt name must be a simple filename')
    frozen = read_json(run / 'rank_comparison' / 'frozen_candidates.json')
    allowed = {pdb for case in frozen['cases'] for pdb in case['top10_union']}
    triage_dir = run / 'rank_comparison' / 'material_triage'
    selected = sorted(set(x.upper() for x in args.pdb))
    for pdb in selected:
        if pdb not in allowed:
            raise ValueError(f'candidate is outside the frozen union: {pdb}')
        prior = [read_json(x) for x in triage_dir.glob('*_' + pdb + '.json')]
        if not prior or not any(x['status'] == 'SOURCE_NOT_IN_SEALED_ASSETS' for x in prior):
            raise ValueError(f'candidate was not recorded as missing from sealed assets: {pdb}')
    rows = pq.read_table(index / 'entry_sources.parquet',
                         filters=[('entry_pdb_id', 'in', [x.lower() for x in selected])]).to_pylist()
    revisions = {r['entry_pdb_id'].upper(): (int(r['source_mmcif_major_revision']),
                                           int(r['source_mmcif_minor_revision'])) for r in rows}
    raw = run / 'rank_comparison' / 'supplemental_raw'
    raw.mkdir(parents=True, exist_ok=True)
    results = []
    for pdb in selected:
        manifest = raw / f'{pdb}.json'
        if manifest.exists():
            d = read_json(manifest)
            if sha(Path(d['path'])) != d['sha256']:
                raise ValueError(f'cached fixed original changed: {pdb}')
            results.append(d | {'reused': True})
            continue
        started = now()
        try:
            major, minor = revisions[pdb]
            lower = pdb.lower()
            url = (f'https://files-versioned.rcsb.org/pdb_versioned/data/entries/{lower[1:3]}/'
                   f'pdb_0000{lower}/pdb_0000{lower}_xyz_v{major}-{minor}.cif.gz')
            body, response = cache_get(run, f'frozen_candidate/{lower}_v{major}-{minor}', url, timeout=90)
            if not body.startswith(b'\x1f\x8b'):
                raise ValueError('response is not gzip')
            text = gzip.decompress(body).decode('utf8')
            block = gemmi.cif.read_string(text).sole_block()
            actual, first_date = revision(block, pdb)
            expected = f'{major}.{minor}'
            if actual != expected:
                raise ValueError(f'revision mismatch: expected {expected}, found {actual}')
            target = raw / f'{lower}-v{expected}.cif.gz'
            with target.open('xb') as f:
                f.write(body)
            d = {'pdb': pdb, 'path': str(target), 'sha256': hashlib.sha256(body).hexdigest(),
                 'bytes': len(body), 'revision': actual, 'first_public_model_date': first_date,
                 'source_version_basis': 'PLINDER 2026-09 entry_sources', 'request': response,
                 'started_utc': started, 'finished_utc': now(),
                 'limit': 'AHoJ online source revision is not thereby proven identical.'}
            write_json(manifest, d)
            results.append(d | {'reused': False})
            print('FETCH', pdb, actual, len(body), flush=True)
        except Exception as exc:
            failure = {'pdb': pdb, 'status': 'TECHNICAL_FAILURE', 'started_utc': started,
                       'finished_utc': now(), 'error': f'{type(exc).__name__}: {exc}'}
            n = 1
            while (raw / f'{pdb}.failure{n:03d}.json').exists():
                n += 1
            write_json(raw / f'{pdb}.failure{n:03d}.json', failure)
            results.append(failure)
            print('FETCH_ERROR', pdb, failure['error'], flush=True)
    code = run / 'code' / f'{Path(__file__).stem}-{sha(Path(__file__))[:12]}.py'
    if not code.exists():
        shutil.copy2(Path(__file__), code)
    receipt = {'created_utc': now(), 'requested': selected, 'results': results,
               'new_files': sum(x.get('reused') is False for x in results),
               'failures': sum(x.get('status') == 'TECHNICAL_FAILURE' for x in results)}
    write_json(run / 'rank_comparison' / args.receipt_name, receipt)
    print(json.dumps({'receipt': args.receipt_name, 'new_files': receipt['new_files'],
                      'failures': receipt['failures']}), flush=True)


if __name__ == '__main__':
    main()
