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

    def test_three_main_separates_live_recovery_and_shutdown(self):
        batch = {'frozen': 5, 'first_pending': 1, 'pending_at_stop': 0, 'settled_at_stop': 5}
        values = {
            'status': 'PASS_BOUNDED_EXECUTION', 'verified_utc': '2026-01-04T01:00:00Z',
            'execution_commit': 'a'*40, 'verified_public_runtime_sha256': self.runtime,
            'batch': {'b0011': batch, 'b0012': batch, 'b0013': batch | {'pending_at_stop': 1, 'settled_at_stop': 4}},
            'progress': {'status': 'PAUSED_OBSERVATION_LIMIT', 'partitions_completed': 13,
                         'new_settled_entries': 64, 'attempted_pending_count': 1, 'unattempted_entries': 15},
            'catalog': {'REGISTERED': 2, 'REVIEW': 8, 'NOT_EXPORTED': 1}, 'frozen_rolling_entries': 80,
            'mainline': {'new_gaps': 3, 'merged_before_or_at_third_commit': 1,
                         'merged_after_third_commit': 1, 'backlog_at_third_commit': 2},
            'partition_14_preflight_or_dispatch': False, 'sqlite_integrity_ok': True,
            'sqlite_foreign_keys_ok': True, 'controller_exit_clean': True, 'staging_release_count': 0,
            'query_checks': [{'path': 'PRIVATE_ASSET_PATH'}],
        }
        self.put('three-main.json', values)
        self.write(self.root / SOURCES, self.sources | {'m3_three_main': 'three-main.json'})
        page = render(self.root, self.runtime)
        self.assertIn('**1** were recovered and merged by the last main seal and **1** more at shutdown', page)
        self.assertIn('**10** retrievable systems', page)
        self.assertNotIn('PRIVATE_ASSET_PATH', page)
        self.put('three-main.json', values | {'mainline': values['mainline'] | {'backlog_at_third_commit': 1}})
        with self.assertRaisesRegex(ValueError, 'Three-batch production totals'):
            render(self.root, self.runtime)

    def test_archive_admission_keeps_backlog_and_release_separate(self):
        base = {'status': 'PASS_BOUNDED_EXECUTION', 'no_next_preflight': True,
                'staging_releases': 0, 'code_commit': 'a'*40}
        retry = base | {'retry_input': 9, 'merged_total': 7, 'final': {'partitions': 13, 'pending': 2}}
        main = base | {'verified_utc': '2026-01-04T01:00:00Z', 'retry_input': 2, 'new_gaps': 3,
                       'merged_total': 2, 'merged_during_main': 1, 'main_input': 30, 'source_reconcile': 100,
                       'catalog_count': 95, 'final': {'partitions': 16, 'pending': 3, 'settled': 60, 'unattempted': 37},
                       'main_batches': {'b0014': {'head': {'local_queue_timeout': 0, 'sent': 10},
                                                   'get': {'partial_body_timeout': 2}}},
                       'catalog_run': 'PRIVATE_CATALOG_RUN'}
        self.write(self.root / SOURCES, self.sources | {'m3_archive_admission': 'admission.json'})
        self.put('admission.json', {'existing_gaps': retry, 'mixed_load': main, 'gate': 'HOLD_BACKLOG'})
        page = render(self.root, self.runtime)
        self.assertIn('leaving **3** pending overall', page)
        self.assertIn('Expansion was held because retries did not keep up', page)
        self.assertNotIn('PRIVATE_CATALOG_RUN', page)
        self.put('admission.json', {'existing_gaps': retry, 'mixed_load': main | {'retry_input': 1}, 'gate': 'HOLD_BACKLOG'})
        with self.assertRaisesRegex(ValueError, 'Archive admission production evidence does not reconcile'):
            render(self.root, self.runtime)

    def test_changed_rollup_or_replay_cannot_be_published(self):
        self.put('rerun.json', self.first | {'database_registration_reused': True, 'results_sha256': 'e'*64})
        with self.assertRaisesRegex(ValueError, 'registration and replay receipts disagree'):
            render(self.root, self.runtime)
        self.put('rollup.json', {'different': True})
        with self.assertRaisesRegex(ValueError, 'not the independently verified'):
            render(self.root, self.runtime)

    def test_apo_aggregate_keeps_zero_hidden_enrollment_distinct_from_recall(self):
        code = self.root / 'scripts/apo_c_autopilot.py'
        code.parent.mkdir()
        code.write_bytes(b'pass\n')
        counts = dict(natural_targets=980, natural_targets_with_terminal_results=980,
                      natural_first_readback_passed=35, natural_delivered_pairs=39, natural_chain_copies=43,
                      natural_pairs_sampled=20, natural_pairs_author_confirmed=19, natural_readback_errors=0,
                      natural_sampled_author_unresolved=1, natural_sampled_readback_failed=0,
                      hidden_expected=20, hidden_actual=0, hidden_first_verified=0)
        row = dict(counts=counts, created_utc='2026-01-04T09:00:00+08:00', run_elapsed_s=3600, rule='trial-rule',
                   formal_pair_acceptances=0, fixed_controls_all_pass=True, final_output_seal_sha256='a'*64,
                   autopilot_sha256=sha256(code.read_bytes()).hexdigest(), run='PRIVATE_RUN_PATH')
        self.write(self.root / SOURCES, self.sources | {'apo_trial': 'apo-trial.json'})
        self.put('apo-trial.json', row)
        page = render(self.root, self.runtime)
        self.assertIn('**35/980**', page)
        self.assertIn('**1** author checks remained unresolved.', page)
        self.assertIn('Zero independently verified hidden positives enrolled. Recall unavailable.', page)
        self.assertNotIn('PRIVATE_RUN_PATH', page)
        self.put('apo-trial.json', row | {'counts': counts | {'hidden_first_verified': 1}})
        with self.assertRaisesRegex(ValueError, 'Apo trial aggregates do not reconcile'):
            render(self.root, self.runtime)

        self.put('apo-trial.json', row)
        self.write(self.root / SOURCES, self.sources | {'apo_trial': 'apo-trial.json', 'apo_hidden_trial': 'hidden.json'})
        supplemental = dict(denominator=3, reference_delivered=1, terminal_targets=3,
                            reference_searched=3, reference_selected=3, reference_inspected=3,
                            scientific_rule=row['rule'], autopilot_sha256=row['autopilot_sha256'],
                            parent_final_output_seal_sha256=row['final_output_seal_sha256'],
                            production_final_seal_sha256='b'*64, private_evidence='PRIVATE_REFERENCE_KEY')
        self.put('hidden.json', supplemental)
        page = render(self.root, self.runtime)
        self.assertIn('separately sealed blind run recovered **1/3** positives', page)
        self.assertNotIn('PRIVATE_REFERENCE_KEY', page)
        review = dict(automated_confirmed=19, adjudicated_additional_confirmed=1, misleading_excerpt_references=3,
                      evidence_sha256='d'*64, failure_sample_count=10, verified_missed_targets=2)
        reviewed = row | {'counts': counts | {'natural_pairs_author_confirmed': 20, 'natural_sampled_author_unresolved': 0},
                          'independent_author_review': review, 'verified_missed_delivery_lower_bound': 2}
        self.put('apo-trial.json', reviewed)
        page = render(self.root, self.runtime)
        self.assertIn('**20/20** sampled relations', page)
        self.assertIn('**3** incorrect excerpts', page)
        self.assertIn('**2** missed targets among **10**', page)
        self.put('apo-trial.json', reviewed | {'verified_missed_delivery_lower_bound': 3})
        with self.assertRaisesRegex(ValueError, 'Apo independent author adjudication does not reconcile'):
            render(self.root, self.runtime)
        self.put('apo-trial.json', reviewed)
        self.put('hidden.json', supplemental | {'parent_final_output_seal_sha256': 'c'*64})
        with self.assertRaisesRegex(ValueError, 'Supplemental Apo trial aggregates do not reconcile'):
            render(self.root, self.runtime)


    def test_batch_200_requires_verified_stop_and_separates_the_two_denominators(self):
        checks = dict(accepted=True, same_input_replay=True, provenance_and_foreign_keys_ok=True)
        gate = dict(checks, window_start=101, sealed=105, plan_entries_this_window=9,
                    retry={'exhausted': 1}, entries_per_hour=90, next_source_feed=None,
                    observed_recent_primary_interval_baseline_entries_per_hour=95)
        progress = dict(active_partition=None, partitions_completed=200, new_settled_entries=70,
                        unattempted_entries=24, attempted_pending_count=6, m1_entries=120,
                        exclusions={'HISTORICAL_SAME_REVISION': 5, 'OUT_OF_SCOPE': 15})
        retry = dict(queued=6, exhausted=6, running=0, retryable_not_ready=0,
                     verified_waiting_for_catalog=0)
        final = dict(checks, window_start=106, sealed=200, plan_entries_this_window=31,
                     verified_utc='2026-01-04T01:00:00Z', progress=progress, retry=retry,
                     queryable=130, confirmed=14, published=0, entries_per_hour=105,
                     head='PRIVATE_BATCH_200_CATALOG')
        data = dict(stage_final=final, gate=gate, execution_commit='a' * 40,
                    frozen_101_200=40, rolling_entries=100, no_batch_201=True,
                    writer_locks_free=True, watcher={'state': 'FINISHED'},
                    supervisor={'state': 'NORMAL_EXIT'})
        self.write(self.root / SOURCES, self.sources | {'m3_through200': 'through200.json'})
        self.put('through200.json', {'data': data})
        page = render(self.root, self.runtime)
        self.assertIn('**35/40** frozen PDB entries settled', page)
        self.assertIn('**70/100** settled, **6** exhausted, **24** not yet attempted', page)
        self.assertNotIn('PRIVATE_BATCH_200_CATALOG', page)
        for change in ({'no_batch_201': False},
                       {'stage_final': final | {'same_input_replay': False}},
                       {'stage_final': final | {'retry': retry | {'verified_waiting_for_catalog': 1}}}):
            self.put('through200.json', {'data': data | change})
            with self.assertRaisesRegex(ValueError, 'Invalid batch-200 closure'):
                render(self.root, self.runtime)


if __name__ == '__main__':
    unittest.main()
