"""Build the public APHCP status from selected, accepted local receipts."""
from datetime import date, datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re


STATUS = 'docs/public/status.md'
SOURCES = 'config/public_review_sources.local.json'
REGISTRATION = ('REGISTERED', 'REVIEW', 'OUT_OF_SCOPE')
OBSERVATIONS = (*REGISTRATION, 'OBSERVED', 'PREPARATION_USE_RECORDED')
PENDING = {
    'm3_new_batch': 'Use a genuinely new input batch to verify registration and chemistry annotation end to end.',
    'm3_effective_decision': 'Verify that lookup selects the effective scientific decision while retaining its history.',
    'm3_long_tail': 'Reduce repeat Assembly preparation during independent readback for entries with many systems; measure throughput on the next authorized window.',
    'm3_rolling_measurement': 'Measure sustained whole-catalog throughput and allocated storage; the bounded scale gate is not a completion date.',
    'm3_scale': 'Measure sustained throughput on a separately approved production increment; the bounded pilot is not a scale estimate.',
    'm3_retry_balance': 'Verify that sustained recovery keeps pace with mainline technical gaps before resuming large-scale rolling production.',
    'm3_connectivity': 'Keep partner-connection and role exceptions under review; the bounded pilot does not broaden scientific acceptance.'
}
SHA = re.compile(r'[0-9a-f]{64}')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def runtime_fingerprint(payload):
    """Hash the public executable inputs, independent of review docs and tests."""
    prefixes = ('src/', 'scripts/', 'migrations/', 'config/', 'workflow/')
    dependencies = {'pyproject.toml', 'requirements-local.txt', 'requirements-server.txt', '.streamlit/config.toml'}
    review_tools = {'scripts/public_review_status.py', 'scripts/sync_github.py'}
    rows = [[name, sha256(data.replace(b"\r\n", b"\n")).hexdigest()] for name, data in sorted(payload.items())
            if (name.startswith(prefixes) or name in dependencies) and name not in review_tools]
    return sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()


def receipt(root, name):
    base = (root / 'reports/inventory').resolve()
    path = (base / name).resolve()
    require(path.is_relative_to(base) and path.suffix == '.json', 'Receipt must be JSON inside reports/inventory')
    raw = path.read_bytes()
    return json.loads(raw.decode('utf-8-sig')), sha256(raw).hexdigest()


def count(value):
    require(type(value) is int and value >= 0, 'Expected a nonnegative count in receipt')
    return value


def counts(section, names):
    return {key: count(section.get(key, 0)) for key in names}


def utc(value):
    dt = datetime.fromisoformat(value.replace('Z', '+00:00'))
    require(dt.tzinfo is not None, 'Verification timestamp must carry a time zone')
    return dt.astimezone(timezone.utc).isoformat(timespec='seconds')


def code_link(verified, current):
    if verified is None:
        return 'Unknown: the historical Run is not pinned to a public runtime fingerprint.'
    require(isinstance(verified, str) and SHA.fullmatch(verified), 'Invalid verified public runtime fingerprint')
    if verified == current:
        return 'The accepted evidence covers this public runtime fingerprint.'
    return f'Current public runtime is not covered by this evidence (verified runtime: `{verified}`).'


def pilot_section(root, source):
    pilot, _ = receipt(root, source)
    transition = pilot['transition']
    keys = ('INDEX_MATCH', 'INDEX_VERSION_MISMATCH', 'NO_LIGAND_ANNOTATION')
    routes = counts(pilot['strata'], keys)
    exported = counts(pilot['entry_routes_with_readback'], keys)
    states = counts(pilot['m3_outcomes'], REGISTRATION)
    systems = pilot['system_statuses']
    require(pilot['status'] == 'PASS' and pilot['formal_release_count'] == 0 and
            pilot['protected_registries_unchanged'] is True and pilot['replay_preserved_files_and_results'] is True and
            all(transition[k] is True for k in ('test_only','reverse_order_pass','same_system','same_files','catalog_pass')) and
            sum(routes.values()) == pilot['entries'] and all(0 < exported[k] <= routes[k] for k in keys) and
            sum(map(count,systems.values())) == pilot['candidate_count'] and
            sum(states.values()) == pilot['registration_count'] == systems['EXPORTED_READBACK'],
            'Incremental M3 pilot receipt is incomplete or inconsistent')
    return f"""## M3 · New-batch pilot (unpublished)

- Independently verified: **{utc(pilot['verified_utc'])}**. **{sum(routes.values()):,}** new entries: **{routes[keys[0]]:,}** same-revision, **{routes[keys[1]]:,}** different-revision, and **{routes[keys[2]]:,}** without PLINDER ligand annotations.
- **{count(pilot['candidate_count']):,}** candidate systems, **{systems['EXPORTED_READBACK']:,}** structures read back. **{states['REGISTERED']:,}** newly confirmed in this isolated staging batch and **{states['REVIEW']:,}** held with files retained. These numbers are separate from the historical totals above.
- A test-only `REVIEW` → `REGISTERED` transition preserved the same system and files, retained both decisions, and was independent of run order. Replaying the real batch reused its outputs.

"""


def production_section(root, source, heading="First continuous production batch"):
    batch, _ = receipt(root, source)
    frozen = count(batch['frozen_entries']); originals = count(batch['originals_verified'])
    structures = batch['structural_states']; decisions = batch['review_dispositions']
    fetched = count(batch['archive_network_bytes']); systems = count(batch['new_catalog_queryable'])
    require(batch['status'] == 'STAGED_COMPLETE' and not batch['pending_technical'] and
            frozen == originals == sum(map(count, batch['entry_statuses'].values())) and
            count(batch['candidates']) == sum(map(count, structures.values())) and
            sum(map(count, decisions.values())) <= count(structures['EXPORTED_READBACK']) and
            systems <= structures['EXPORTED_READBACK'] and
            count(batch['catalog_queryable_total']) >= systems and batch['same_input_reused'] is True and
            batch['unified_lookup_verified'] is True and batch['targeted_retry_verified'] is True,
            'Continuous M3 verification or replay is incomplete')
    left, right = map(date.fromisoformat, batch['window'])
    elapsed = batch['pipeline_seconds']
    require(left < right and type(elapsed) in (int, float) and 0 < elapsed < 86400 * 30,
            'Invalid continuous production window or measurement')
    return f'''## M3 · {heading} (unpublished)

- Independently verified: **{utc(batch['verified_utc'])}**. First release dates from **{left}** (included) to **{right}** (excluded).
- **{frozen:,}** exact-revision originals verified; **{count(batch['candidates']):,}** candidates, **{structures['EXPORTED_READBACK']:,}** read back and **{systems:,}** new retrievable systems. **{count(decisions.get('REGISTERED',0)):,}** confirmed in staging; **{count(decisions.get('REVIEW',0)):,}** held for review. Unified catalog: **{batch['catalog_queryable_total']:,}** retrievable systems.
- **{fetched:,} bytes** of archival response bodies; **{elapsed:.1f} s** for the production pipeline in this window. Replaying the same input reused committed phases. In a separate control, the resume command repaired one failed Entry without modifying prior results.

'''


def recovery_section(root, name):
    r, _ = receipt(root, name)
    first, later = r['fourth_batch'], r['isolated']
    require(r['scope'] == 'ISOLATED_ONLY' and r['live_head_unchanged'] is True and
            count(first['entries']) == count(first['reused']) + count(first['recovered']) and
            count(later['settled']) == count(later['entries']) and
            later['next_primary_before_retry'] is True and later['replay_unchanged'] is True and
            count(later['retry_receipts']) == 1 and count(first['releases']) == 0,
            'Missing completed isolated rolling recovery evidence')
    return f'''## M3 · Rolling recovery control (isolated, unpublished)

- Verified: **{utc(r['verified_utc'])}**. The interrupted batch reconciled **{first['entries']:,}** entries in an isolated catalog: **{first['reused']:,}** reused, **{first['recovered']:,}** recovered. The live catalog pointer was left in place.
- A separate **{later['entries']:,}**-entry control sealed its next primary batch before the slow-file retry merged. A bounded Range retry recovered the previously partial archive bodies. The retry wave merged once, and replay did not change committed receipts.
- This is engineering verification on already scoped entries. Sustained production throughput and the live controller handoff are outside this control.

'''


def rolling_gate_section(root, name):
    r, _ = receipt(root, name)
    batches = r['main_new_batches']
    pending = sum(len(v) for v in r['pending'].values())
    require(r['status'] == 'PAUSED_OBSERVATION_LIMIT' and len(batches) == 2 and
            all(r[key] is True for key in ('old_hashes_unchanged', 'observed_no_batch_7',
                                           'writer_locks_free', 'seed_records_preserved')) and
            count(r['settled']) + pending == count(r['attempted']) and
            r['attempted'] + count(r['unattempted']) == count(r['frozen_rolling_entries']) and
            count(r['confirmed']) <= count(r['queryable']) and count(r['releases']) == 0,
            'Production recovery gate evidence is incomplete')
    attempts = sum(count(b['frozen']) for b in batches.values())
    pending_main = sum(len(b['pending_final']) for b in batches.values())
    attempt_times = [b['wall_from_previous_gate_s'] for _, b in sorted(batches.items())]
    require(attempts > 0 and 0 <= pending_main <= attempts and
            all(type(t) in (float, int) and t > 0 for t in attempt_times), 'Invalid production recovery samples')
    seal = ' / '.join(f'{t:.1f} s' for t in attempt_times)
    return f'''## M3 · Bounded production recovery gate (unpublished)

- Independently verified: **{utc(r['verified_utc'])}**. The previous controller was closed, its sealed history was preserved, and the interrupted batch was completed by reusing **{count(r['b4_reused']):,}** settled entries and recovering **{count(r['b4_recovered']):,}** missing ones.
- Exactly two new production batches attempted **{attempts:,}** entries; **{pending_main:,}** were still pending in those batches at the observation limit. Successive batch seals took **{seal}**, including handoff/preflight where applicable. These observations do not establish a full-catalog ETA.
- The staging catalog contained **{r['queryable']:,}** retrievable structures, including **{r['confirmed']:,}** confirmed in staging. Overall, **{r['settled']:,} / {r['frozen_rolling_entries']:,}** frozen rolling entries were settled; **{pending:,}** attempted entries remained pending. Production paused at the gate; no dataset release was created.

'''


def after_ten_section(root, name):
    r, _ = receipt(root, name)
    retry=r['retry']; paired=r['paired']; timing=r['timing']; control=r['control']
    require(r['scope']=='ISOLATED_EXISTING_GAPS' and retry['input']==160 and
            count(retry['recovered'])+count(retry['remaining'])==retry['input'] and
            control['main_unchanged'] is True and count(control['main_batches'])==10 and
            control['restart_verified'] is True and control['same_structural_evidence'] is True and
            count(control['releases'])==0 and paired['whole']['n']==paired['ranges']['n']==24,
            'Isolated after-ten evidence is incomplete or inconsistent')
    registry=r.get('registry_speedup'); replay=r.get('rerun'); detail=''
    if registry:
        require(all(registry[k] is True for k in ('equivalent_review_sha256','equal_registry_counts',
                                                   'quick_check_and_foreign_keys_passed')),
                'Registry speed control did not retain the same checks and decisions')
        detail+=f"- With the same {count(registry['candidates'])}-system staging input, registration took **{registry['before_full_s']:.1f} → {registry['after_full_s']:.1f} s**. The database integrity and foreign-key checks still ran, and the review result hash and registry counts matched.\n"
    if replay:
        require(replay['unchanged'] is True and replay['remaining']==retry['remaining'] and replay['formal_unchanged'] is True,
                'Same-input isolate replay does not match the accepted state')
        detail+=f"- Same-input replay preserved **{count(replay['evidence_files'])}** checked evidence files, including the catalog, receipts, and source proofs.\n"
    return f'''## M3 · Engineering verification after ten main batches (isolated)

- Evidence snapshot **{utc(r['verified_utc'])}**: **{retry['recovered']:,}** of **{retry['input']:,}** existing technical gaps recovered; **{retry['remaining']:,}** remained at the snapshot. The main catalog pointer was unchanged and no new main batch started.
- Verified that retries continue during a main-stage wait. In an isolated process interruption, recovery waited for the inherited writer lease and reused the committed receipt and source proof. Intentional exit did not restart.
- Under the same 30-second transfer budget, whole GET verified **{count(paired['whole']['verified'])} of 24** and Range verified **{count(paired['ranges']['verified'])} of 24**, with **{count(paired['whole']['requests'])}** and **{count(paired['ranges']['requests'])}** requests respectively. The small-file Range threshold was retained.
- On the same existing heavy entry, structure build took **{timing['build_before_s']:.1f} → {timing['build_after_s']:.1f} s**; the separate independent review control took **{timing['review_before_s']:.1f} → {timing['review_after_s']:.1f} s**, with equivalent candidate decisions and structure evidence. These are isolated timings, not whole-pipeline throughput.
{detail}
'''



def three_main_section(root, name, runtime):
    r, _ = receipt(root, name)
    progress, main, catalog = r['progress'], r['mainline'], r['catalog']
    batches = r['batch']
    require(sorted(batches) == ['b0011', 'b0012', 'b0013'] and
            r['status'] == 'PASS_BOUNDED_EXECUTION' and progress['status'] == 'PAUSED_OBSERVATION_LIMIT' and
            progress['partitions_completed'] == 13 and r['partition_14_preflight_or_dispatch'] is False and
            all(r[k] is True for k in ('sqlite_integrity_ok', 'sqlite_foreign_keys_ok', 'controller_exit_clean')) and
            sum(count(b['first_pending']) for b in batches.values()) == count(main['new_gaps']) and
            sum(count(b['pending_at_stop']) for b in batches.values()) == count(progress['attempted_pending_count']) and
            count(r['staging_release_count']) == 0,
            'Three-batch production boundary or per-batch evidence does not reconcile')
    attempted = sum(count(b['frozen']) for b in batches.values())
    during, after = count(main['merged_before_or_at_third_commit']), count(main['merged_after_third_commit'])
    settled, pending, unattempted = (count(progress[k]) for k in
                                     ('new_settled_entries', 'attempted_pending_count', 'unattempted_entries'))
    frozen = count(r['frozen_rolling_entries'])
    q = sum(count(catalog.get(k, 0)) for k in ('OBSERVED', 'PREPARATION_USE_RECORDED', 'REGISTERED', 'REVIEW'))
    commit = r['execution_commit']
    require(isinstance(commit, str) and re.fullmatch(r'[0-9a-f]{40}', commit) and
            settled + pending + unattempted == frozen and
            sum(count(b['settled_at_stop']) for b in batches.values()) + pending == attempted and
            during + count(main['backlog_at_third_commit']) == main['new_gaps'] and
            after + pending == main['backlog_at_third_commit'],
            'Three-batch production totals or source commit do not reconcile')
    return f'''## M3 · Three consecutive production batches (unpublished)

- Verified **{utc(r['verified_utc'])}** against execution commit `{commit}`. {code_link(r.get('verified_public_runtime_sha256'), runtime)}
- Batches 11–13 attempted **{attempted:,}** entries. **{main['new_gaps']:,}** technical gaps appeared; **{during:,}** were recovered and merged by the last main seal and **{after:,}** more at shutdown, leaving **{pending:,}**. Batch 14 was neither preflighted nor dispatched. Three batches do not establish sustained recovery capacity or a full-catalog ETA.
- **{settled:,} / {frozen:,}** rolling entries were settled and **{unattempted:,}** had not started. The staging catalog contained **{q:,}** retrievable systems, including **{catalog.get('REGISTERED', 0):,}** confirmed in staging. Database integrity, foreign keys and file retrieval were checked. No dataset release was created.

'''


def archive_admission_section(root, name):
    r, _ = receipt(root, name)
    retry, main = r['existing_gaps'], r['mixed_load']
    commit = retry['code_commit']
    require(re.fullmatch(r'[0-9a-f]{40}', commit) and commit == main['code_commit'] and
            retry['status'] == main['status'] == 'PASS_BOUNDED_EXECUTION' and
            retry['final']['partitions'] == 13 and main['final']['partitions'] == 16 and
            retry['retry_input'] - retry['merged_total'] == retry['final']['pending'] == main['retry_input'] and
            main['retry_input'] + main['new_gaps'] - main['merged_total'] == main['final']['pending'] and
            all(x['staging_releases'] == 0 and x['no_next_preflight'] for x in (retry, main)),
            'Archive admission production evidence does not reconcile')
    after = count(main['merged_total']) - count(main['merged_during_main'])
    final = main['final']
    head = sum(count(b['head']['local_queue_timeout']) for b in main['main_batches'].values())
    sent = sum(count(b['head']['sent']) for b in main['main_batches'].values())
    get = sum(count(b['get']['partial_body_timeout']) for b in main['main_batches'].values())
    require(after >= 0 and count(final['settled']) + count(final['pending']) + count(final['unattempted']) ==
            count(main['source_reconcile']), 'Archive admission final totals do not reconcile')
    go = r['gate']
    require(go in ('HOLD_BACKLOG', 'HOLD_OTHER', 'PASS_BOUNDED_SCALE'), 'Unknown archive admission decision')
    phrase = {'HOLD_BACKLOG': 'Expansion was held because retries did not keep up.',
              'HOLD_OTHER': 'Expansion was held pending further engineering evidence.',
              'PASS_BOUNDED_SCALE': 'The measured window met the bounded scale gate.'}[go]
    return f'''## M3 · Archive admission and retry capacity (unpublished)

- Verified **{utc(main['verified_utc'])}** against execution commit `{commit}`. From **{retry['retry_input']:,}** existing gaps, **{retry['merged_total']:,}** were recovered and merged before the next main batch.
- Batches 14–16 attempted **{main['main_input']:,}** entries. **{sent:,}** HEAD requests reached HTTP; **{head:,}** timed out waiting locally. Main GET had **{get:,}** partial-body timeouts. The shared archive slot limit was unchanged.
- **{main['new_gaps']:,}** new technical gaps appeared. **{main['merged_during_main']:,}** were merged during production and **{after:,}** more at shutdown, leaving **{final['pending']:,}** pending overall. {phrase}
- **{final['settled']:,} / {main['source_reconcile']:,}** rolling entries were settled. The staging catalog contained **{count(main['catalog_count']):,}** retrievable systems; no dataset release was created.

'''


def render(root: Path, current_runtime: str):
    require(SHA.fullmatch(current_runtime), 'Invalid current public runtime fingerprint')
    src = json.loads((root / SOURCES).read_text(encoding='utf-8-sig'))
    m1, m3 = src['m1'], src['m3']
    rollup, rollup_hash = receipt(root, m1['rollup'])
    verified, _ = receipt(root, m1['verification'])
    require(verified['cross_window_status'] == 'INVENTORY_COMPLETE' and
            verified['merged_summary_sha256'] == rollup_hash,
            'M1 rollup is not the independently verified complete inventory')
    m1_counts = rollup['counts']
    entries = count(m1_counts['entries'])
    windows = count(rollup['window_count'])
    annotated = count(m1_counts['plinder_annotation_entries'])
    unannotated = count(m1_counts['plinder_no_annotation_entries'])
    candidates = count(m1_counts['plinder_candidate_systems'])
    require(entries == annotated + unannotated == rollup['denominator']['total_count'] == verified['entries'] and
            windows == verified['windows'] and candidates == verified['candidate_systems'] and
            unannotated == verified['plinder_no_annotation_entries'], 'M1 totals disagree')

    catalog, _ = receipt(root, m3['catalog_verification'])
    require(catalog['sqlite_integrity'] == catalog['protected_and_rerun_sha256'] == 'PASS' and
            count(catalog['missing_files']) == count(catalog['release_count']) == 0,
            'M3 catalog must be independently verified, retrievable and unpublished')
    review = dict.fromkeys(REGISTRATION, 0)
    registration_runs = m3['registration_runs']
    require(registration_runs, 'No accepted M3 registration receipts selected')
    seen = set()
    for run in registration_runs:
        first, _ = receipt(root, run['initial'])
        rerun, _ = receipt(root, run['rerun'])
        manifest, result = first['manifest_sha256'], first['results_sha256']
        require(SHA.fullmatch(manifest) and SHA.fullmatch(result) and manifest not in seen,
                'Invalid or duplicate M3 registration manifest')
        seen.add(manifest)
        run_counts = counts(first['dispositions'], REGISTRATION)
        require(first['status'] == rerun['status'] == 'M3_STAGING_UNPUBLISHED' and
                rerun['database_registration_reused'] is True and rerun['manifest_sha256'] == manifest and
                rerun['results_sha256'] == result and rerun['dispositions'] == first['dispositions'] and
                sum(run_counts.values()) == first['input_candidate_count'] == rerun['input_candidate_count'],
                'M3 registration and replay receipts disagree')
        for key in review:
            review[key] += run_counts[key]
    outcomes = counts(catalog['outcomes'], OBSERVATIONS)
    source_rows = count(catalog['source_rows'])
    retrievable = count(catalog['queryable_systems'])
    unknown_type = count(catalog['unknown_primary_class'])
    require(review == counts(catalog['previous_m3_decisions'], REGISTRATION) and
            all(review[k] == outcomes[k] for k in ('REGISTERED', 'REVIEW')) and
            sum(outcomes.values()) == source_rows == retrievable + outcomes['OUT_OF_SCOPE'] and
            unknown_type <= retrievable, 'M3 catalog and registration totals disagree')

    m1_at = utc(verified['verification_finished_utc'])
    m3_on = date.fromisoformat(m3['verified_on']).isoformat()
    m1_link = code_link(m1.get('verified_public_runtime_sha256'), current_runtime)
    m3_link = code_link(m3.get('verified_public_runtime_sha256'), current_runtime)
    pilot = pilot_section(root, src['m3_pilot_verification']) if src.get('m3_pilot_verification') else ''
    production = production_section(root, src['m3_continuous_verification']) if src.get('m3_continuous_verification') else ''
    next_window = src.get('m3_scale_gate_verification')
    if next_window:
        gate, _ = receipt(root, next_window)
        require(gate['accepted'] is True and 300 <= count(gate['frozen_entries']) <= 500,
                'Missing bounded production scale gate')
        production += production_section(root, next_window, 'Scale gate: second continuous production batch')
    recovery = recovery_section(root, src['m3_rolling_recovery']) if src.get('m3_rolling_recovery') else ''
    production_gate = rolling_gate_section(root, src['m3_rolling_production_gate']) if src.get('m3_rolling_production_gate') else ''
    after_ten = after_ten_section(root, src['m3_after_ten']) if src.get('m3_after_ten') else ''
    three_main = three_main_section(root, src['m3_three_main'], current_runtime) if src.get('m3_three_main') else ''
    admission = archive_admission_section(root, src['m3_archive_admission']) if src.get('m3_archive_admission') else ''
    pending = '\n'.join(f'- {PENDING[k]}' for k in src['open_review']) or '- None listed in this snapshot.'
    return f'''# APHCP public verification snapshot

Generated during GitHub synchronization from accepted local receipts. The dates below are **evidence dates**, not the time this page was generated. This snapshot describes inventory and unpublished staging; it is not a dataset release or a live service status.

Public runtime fingerprint (SHA-256): `{current_runtime}`

## M1 · Metadata inventory

- Independently verified: **{m1_at}**.
- **{entries:,} unique experimental protein PDB entries** across **{windows:,} date windows**.
- PLINDER ligand annotations: **{annotated:,}** entries with annotations; **{unannotated:,}** without. **{candidates:,}** indexed candidate systems, before structure QC.
- Verification and current public code: {m1_link}

## M3 · Staging registration and structure observation catalog

- Catalog independently verified: **{m3_on} (Asia/Shanghai)**.
- Registration review: **{sum(review.values()):,}** proposals; **{review['REGISTERED']:,}** confirmed in staging, **{review['REVIEW']:,}** held for review with files retained, **{review['OUT_OF_SCOPE']:,}** out of scope.
- Separate observation catalog: **{source_rows:,}** records, **{retrievable:,}** retrievable, including the **{review['REGISTERED']:,}** confirmed staging systems. The catalog has **{outcomes['OUT_OF_SCOPE']:,}** out-of-scope records without supplied assets. Its denominator differs from the registration review.
- **{unknown_type:,}** retrievable observations have an unknown primary chemical type. Unknown type does not by itself block eligible structures.
- Verification and current public code: {m3_link}

{pilot}{production}{recovery}{production_gate}{after_ten}{three_main}{admission}## Open review

{pending}

See the [review guide](review.md) for scope, evidence definitions and review instructions.
'''
