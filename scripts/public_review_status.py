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
    'm3_scale': 'Measure sustained throughput on a separately approved production increment; the bounded pilot is not a scale estimate.',
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

{pilot}## Open review

{pending}

See the [review guide](review.md) for scope, evidence definitions and review instructions.
'''
