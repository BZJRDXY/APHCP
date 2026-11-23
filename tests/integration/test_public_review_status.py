"""Synthetic receipt-to-publication example: no production data or SSH required."""
from hashlib import sha256
import json
from pathlib import Path
import sys
import tempfile
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from public_review_status import SOURCES, render, runtime_fingerprint


class PublicReviewStatusTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.runtime = runtime_fingerprint({'src/example.py': b'pass\n'})
        rollup = self.put('rollup.json', {
            'counts': {'entries': 3, 'plinder_annotation_entries': 2,
                       'plinder_no_annotation_entries': 1, 'plinder_candidate_systems': 5},
            'window_count': 1, 'denominator': {'total_count': 3},
            'tables': {'internal-path': 'PRIVATE_EXAMPLE_PATH'},
        })
        self.put('verify.json', {
            'entries': 3, 'windows': 1, 'cross_window_status': 'INVENTORY_COMPLETE',
            'merged_summary_sha256': sha256(rollup.read_bytes()).hexdigest(),
            'candidate_systems': 5, 'plinder_no_annotation_entries': 1,
            'verification_finished_utc': '2026-01-02T09:00:00+08:00',
        })
        result = {'REGISTERED': 1, 'REVIEW': 1, 'OUT_OF_SCOPE': 0}
        self.first = {'manifest_sha256': 'a'*64, 'results_sha256': 'b'*64, 'dispositions': result,
                      'input_candidate_count': 2, 'status': 'M3_STAGING_UNPUBLISHED',
                      'registry_path': 'PRIVATE_EXAMPLE_PATH'}
        self.put('first.json', self.first)
        self.put('rerun.json', self.first | {'database_registration_reused': True})
        self.put('catalog.json', {
            'sqlite_integrity': 'PASS', 'protected_and_rerun_sha256': 'PASS', 'missing_files': 0,
            'release_count': 0, 'previous_m3_decisions': result, 'source_rows': 4,
            'outcomes': result | {'OBSERVED': 1, 'OUT_OF_SCOPE': 1}, 'queryable_systems': 3,
            'unknown_primary_class': 1, 'checked': {'example': 'PRIVATE_EXAMPLE_ID'},
        })
        self.sources = {
            'm1': {'rollup': 'rollup.json', 'verification': 'verify.json',
                   'verified_public_runtime_sha256': self.runtime},
            'm3': {'catalog_verification': 'catalog.json', 'verified_on': '2026-01-04',
                   'registration_runs': [{'initial': 'first.json', 'rerun': 'rerun.json'}]},
            'open_review': ['m3_new_batch'],
        }
        self.write(self.root / SOURCES, self.sources)

    @staticmethod
    def write(path, data):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data), encoding='utf-8')
        return path

    def put(self, name, data):
        return self.write(self.root / 'reports/inventory' / name, data)

    def test_public_view_separates_registration_and_observations(self):
        page = render(self.root, self.runtime)
        self.assertIn('**2** proposals', page)
        self.assertIn('**4** records, **3** retrievable', page)
        self.assertIn('2026-01-02T01:00:00+00:00', page)
        self.assertIn('accepted evidence covers this public runtime fingerprint', page)
        self.assertIn('historical Run is not pinned', page)
        self.assertNotIn('PRIVATE_EXAMPLE', page)
        self.assertNotIn('rollup.json', page)
        self.assertEqual(page, render(self.root, self.runtime))

    def test_code_change_preserves_evidence_date_but_invalidates_match(self):
        same = runtime_fingerprint({'src/example.py': b'pass\r\n', 'docs/public/review.md': b'new'})
        self.assertEqual(same, self.runtime)
        changed = runtime_fingerprint({'src/example.py': b'x = 1\n'})
        page = render(self.root, changed)
        self.assertIn('Current public runtime is not covered by this evidence', page)
        self.assertIn('2026-01-02T01:00:00+00:00', page)

    def test_continuous_summary_requires_replay_and_safe_aggregate(self):
        values = {'status': 'STAGED_COMPLETE', 'pending_technical': [], 'window': ['2018-02-07', '2018-02-15'],
                  'verified_utc': '2026-01-04T09:00:00+08:00', 'pipeline_seconds': 123.4,
                  'frozen_entries': 3, 'originals_verified': 3, 'entry_statuses': {'REVIEW': 3},
                  'candidates': 5, 'structural_states': {'EXPORTED_READBACK': 4, 'REVIEW': 1},
                  'review_dispositions': {'REGISTERED': 2, 'REVIEW': 1}, 'archive_network_bytes': 12345,
                  'new_catalog_queryable': 4, 'catalog_queryable_total': 10, 'same_input_reused': True,
                  'unified_lookup_verified': True, 'targeted_retry_verified': True,
                  'examples': ['PRIVATE_EXAMPLE_ID']}
        self.put('continuous.json', values)
        self.write(self.root / SOURCES, self.sources | {'m3_continuous_verification': 'continuous.json'})
        page = render(self.root, self.runtime)
        self.assertIn('**3** exact-revision originals verified', page)
        self.assertIn('**4** new retrievable systems', page)
        self.assertNotIn('PRIVATE_EXAMPLE', page)
        self.put('continuous.json', values | {'same_input_reused': False})
        with self.assertRaisesRegex(ValueError, 'Continuous M3 verification or replay is incomplete'):
            render(self.root, self.runtime)

    def test_rolling_gate_separates_residual_queue_and_preserves_scope(self):
        batch = {'frozen': 5, 'pending_final': [], 'wall_from_previous_gate_s': 45.5}
        values = {'status': 'PAUSED_OBSERVATION_LIMIT', 'verified_utc': '2026-01-04T01:00:00Z',
                  'main_new_batches': {'b0005': batch, 'b0006': batch | {'pending_final': ['PRIVATE_EXAMPLE_ID']}},
                  'pending': {'b0006': ['PRIVATE_EXAMPLE_ID']}, 'settled': 19, 'attempted': 20,
                  'unattempted': 60, 'frozen_rolling_entries': 80, 'b4_reused': 4, 'b4_recovered': 1,
                  'confirmed': 2, 'queryable': 24, 'releases': 0,
                  'old_hashes_unchanged': True, 'observed_no_batch_7': True,
                  'writer_locks_free': True, 'seed_records_preserved': True}
        self.put('rolling-gate.json', values)
        self.write(self.root / SOURCES, self.sources | {'m3_rolling_production_gate': 'rolling-gate.json'})
        page = render(self.root, self.runtime)
        self.assertIn('two new production batches attempted **10** entries', page)
        self.assertIn('**19 / 80** frozen rolling entries', page)
        self.assertNotIn('PRIVATE_EXAMPLE', page)
        self.put('rolling-gate.json', values | {'observed_no_batch_7': False})
        with self.assertRaisesRegex(ValueError, 'Production recovery gate evidence is incomplete'):
            render(self.root, self.runtime)

    def test_changed_rollup_or_replay_cannot_be_published(self):
        self.put('rerun.json', self.first | {'database_registration_reused': True, 'results_sha256': 'e'*64})
        with self.assertRaisesRegex(ValueError, 'registration and replay receipts disagree'):
            render(self.root, self.runtime)
        self.put('rollup.json', {'different': True})
        with self.assertRaisesRegex(ValueError, 'not the independently verified'):
            render(self.root, self.runtime)


if __name__ == '__main__':
    unittest.main()
