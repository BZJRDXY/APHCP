"""Regression tests for paginated search, deterministic source selection and retry."""
import gzip
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from apo_c_autopilot import fetch_candidate, run_in_download_order, saved_original
from apo_c_pilot import RcsbPageError, rcsb_query, sha
from apo_c_union_triage import existing_original
from apo_c_pair import revision
import gemmi


def response(status, ids=(), total=1001):
    body = json.dumps({'total_count': total, 'result_set': [{'identifier': i} for i in ids]}).encode()
    return SimpleNamespace(status_code=status, content=body, url='https://example.invalid/search',
                           headers={}, text=body.decode())


def fixed_cif(version):
    return gzip.compress(f"""data_1ABC
_entry.id 1ABC
loop_
_pdbx_audit_revision_history.data_content_type
_pdbx_audit_revision_history.major_revision
_pdbx_audit_revision_history.minor_revision
_pdbx_audit_revision_history.revision_date
'Structure model' 1 {version} 2020-01-01
""".encode())


class RecoveryTests(unittest.TestCase):
    def test_later_204_and_short_page_recover_without_refetching_first_page(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            first = response(200, [f'x{i}_1' for i in range(1000)])
            with patch('apo_c_pilot.requests.get', side_effect=[first, response(204), response(200),
                                                                response(200, ['last_1'])]) as get:
                for _ in range(2):
                    with self.assertRaises(RcsbPageError) as caught:
                        rcsb_query(run, 'query', {})
                    self.assertEqual(caught.exception.progress['missing_offset'], 1000)
                    self.assertEqual(caught.exception.progress['retained_hits'], 1000)
                    self.assertFalse((run / 'http/rcsb/query_1000.json').exists())
                out = rcsb_query(run, 'query', {})
                self.assertEqual(len({i['identifier'] for i in out['hits']}), 1001)
                self.assertEqual(out, rcsb_query(run, 'query', {}))
                offsets = [json.loads(c.kwargs['params']['json'])['request_options']['paginate']['start']
                           for c in get.call_args_list]
                self.assertEqual(offsets, [0, 1000, 1000, 1000])

    def test_first_page_204_is_empty_and_cached(self):
        with tempfile.TemporaryDirectory() as tmp, patch('apo_c_pilot.requests.get', return_value=response(204)) as get:
            for _ in range(2):
                out = rcsb_query(Path(tmp), 'zero', {})
                self.assertEqual((out['total_count'], out['hits']), (0, []))
            get.assert_called_once()

    def test_two_catalog_revisions_can_choose_the_required_older_original(self):
        with tempfile.TemporaryDirectory() as tmp, sqlite3.connect(':memory:') as db:
            db.executescript('CREATE TABLE entry_version(id TEXT,pdb_id TEXT,revision TEXT);'
                             'CREATE TABLE asset_file(entry_version_id TEXT,uri TEXT,sha256 TEXT,bytes INT);')
            folder = Path(tmp) / 'raw'
            folder.mkdir()
            for minor in (1, 2):
                path = folder / f'1abc-v1.{minor}.cif.gz'
                path.write_bytes(fixed_cif(minor))
                db.execute('INSERT INTO entry_version VALUES (?, ?, ?)', (str(minor), '1ABC', f'1.{minor}'))
                db.execute('INSERT INTO asset_file VALUES (?, ?, ?, ?)',
                           (str(minor), str(path), sha(path), path.stat().st_size))
            self.assertEqual(existing_original(db, '1ABC')[3], '1.2')
            older = existing_original(db, '1ABC', '1.1')
            self.assertEqual(older[3], '1.1')
            self.assertEqual(revision(gemmi.cif.read(older[0]).sole_block(), '1ABC')[0], '1.1')

    def test_failed_fetch_resume_only_once_and_reuses_validated_original(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            cfg = {'max_compressed_bytes': 20000, 'download_deadline_seconds': 5,
                   'max_candidate_download_attempts': 2}
            with patch('apo_c_autopilot.download_bounded',
                       side_effect=[TimeoutError('injected'), (fixed_cif(1), {'status': 200})]) as download:
                self.assertIsNone(fetch_candidate(run, '1ABC', {}, cfg, '1.1'))
                self.assertIsNone(fetch_candidate(run, '1ABC', {}, cfg, '1.1'))
                second = {k: v for k, v in cfg.items() if k != 'download_session'} | {'retry_failed': True}
                source = fetch_candidate(run, '1ABC', {}, second, '1.1')
                self.assertEqual(source['sha256'], sha(Path(source['path'])))
                self.assertEqual(saved_original(run, '1ABC', '1.1')[0][3], '1.1')
                self.assertEqual(fetch_candidate(run, '1ABC', {}, second, '1.1'), source)
                self.assertEqual(download.call_count, 2)
                self.assertEqual(len(list((run / 'fetch').glob('*.json'))), 2)
                self.assertEqual(len(list((run / 'raw').glob('*.gz'))), 1)

    def test_waiting_holo_get_first_download_before_earlier_holo_get_second(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            config = {'max_compressed_bytes': 20000, 'download_deadline_seconds': 5,
                      'max_candidate_download_attempts': 3}

            def requests(entries):
                for pdb in entries:
                    yield
                    fetch_candidate(run, pdb, {}, config, '1.1')
                return len(entries)

            with patch('apo_c_autopilot.download_bounded', side_effect=TimeoutError('injected')) as download:
                results = run_in_download_order([requests(['1AAA', '1AAB', '1AAC']),
                                                 requests([]), requests(['2AAA', '2AAB'])])
            self.assertEqual(results, [3, 0, 2])
            self.assertEqual([call.args[0].rsplit('/', 1)[1] for call in download.call_args_list],
                             ['pdb_00001aaa_xyz_v1-1.cif.gz', 'pdb_00002aaa_xyz_v1-1.cif.gz',
                              'pdb_00001aab_xyz_v1-1.cif.gz'])


if __name__ == '__main__':
    unittest.main()
