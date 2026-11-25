"""Read a sealed Holo catalogue and compare it with a fixed PLINDER index.

Outputs are research-only and never write to the source catalogue. See --help.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import codecs
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3
import time

import pyarrow.parquet as pq


def compact(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def utc():
    return datetime.now(timezone.utc).isoformat()


def array_stream(path, expected_hash):
    """Stream a large JSON array, checking its already-recorded production hash."""
    decoder = json.JSONDecoder()
    utf8 = codecs.getincrementaldecoder('utf-8')()
    hasher = hashlib.sha256()
    buffer = ''
    started = False
    ended = False
    eof = False
    with path.open('rb') as stream:
        while True:
            if not eof and len(buffer) < 1024 * 1024:
                chunk = stream.read(1024 * 1024)
                if chunk:
                    hasher.update(chunk)
                    buffer += utf8.decode(chunk)
                else:
                    buffer += utf8.decode(b'', final=True)
                    eof = True
            buffer = buffer.lstrip()
            if not started:
                if not buffer and not eof:
                    continue
                if not buffer.startswith('['):
                    raise ValueError(f'Expected JSON array: {path}')
                buffer = buffer[1:]
                started = True
                continue
            buffer = buffer.lstrip()
            if buffer.startswith(','):
                buffer = buffer[1:].lstrip()
            if buffer.startswith(']'):
                ended = True
                buffer = buffer[1:]
                if not eof:
                    for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                        hasher.update(chunk)
                        buffer += utf8.decode(chunk)
                    buffer += utf8.decode(b'', final=True)
                    eof = True
                break
            if not buffer and not eof:
                continue
            try:
                row, end = decoder.raw_decode(buffer)
            except json.JSONDecodeError:
                if eof:
                    raise
                chunk = stream.read(1024 * 1024)
                if chunk:
                    hasher.update(chunk)
                    buffer += utf8.decode(chunk)
                else:
                    buffer += utf8.decode(b'', final=True)
                    eof = True
                continue
            yield row
            buffer = buffer[end:]
    if not ended or buffer.strip():
        raise ValueError(f'Invalid trailing JSON in {path}')
    if hasher.hexdigest() != expected_hash:
        raise ValueError(f'Production receipt hash mismatch: {path}')


def open_source(path):
    return sqlite3.connect(path.resolve().as_uri() + '?mode=ro', uri=True)


def create(out):
    out.mkdir(parents=True, exist_ok=False)
    db = sqlite3.connect(out / 'research.sqlite3')
    db.executescript('''
        PRAGMA journal_mode=WAL;
        CREATE TABLE observation (system_id TEXT PRIMARY KEY, run_id TEXT,
            updated_at TEXT, document TEXT NOT NULL);
        CREATE TABLE assessment (run_id TEXT, system_id TEXT, registered_status TEXT,
            document TEXT, candidate_key TEXT, PRIMARY KEY (run_id, system_id));
        CREATE INDEX assessment_system ON assessment(system_id);
        CREATE TABLE result (system_id TEXT PRIMARY KEY, pdb_id TEXT, revision TEXT,
            disposition TEXT, reasons TEXT, primary_class TEXT, instance_key TEXT,
            evidence_key TEXT, prior_role TEXT, current_role TEXT, conflicts TEXT,
            partner_count INTEGER, target_residue_count INTEGER, old_qc TEXT,
            new_qc TEXT, source_evidence TEXT, decision_state TEXT, match_level TEXT,
            annotation_interpretation TEXT, candidate_count INTEGER, annotations TEXT);
        CREATE INDEX result_grouping ON result(disposition, match_level, annotation_interpretation);
    ''')
    return db


def resolve(assessment_rows):
    if not assessment_rows:
        return None, 'NO_M3_ASSESSMENT'
    by_id = {run + '/' + row['system_id']: row for run, row in assessment_rows}
    parents = set()
    for row in by_id.values():
        for parent in row.get('supersedes') or []:
            if parent not in by_id:
                raise ValueError('A superseded assessment is missing')
            parents.add(parent)
    heads = sorted(set(by_id) - parents)
    if not heads:
        raise ValueError('Assessment graph has no head')
    if len({by_id[key]['disposition'] for key in heads}) > 1:
        return None, 'CONFLICT'
    return by_id[heads[0]], 'EFFECTIVE' if len(heads) == 1 else 'CONCORDANT'


def identifiers(identity):
    targets = identity.get('target_residue_instances') or []
    inst = (identity['pdb_id'], identity['revision'], identity['assembly'],
            identity.get('model_scope'), targets)
    # Evidence reuse proxy: ignore assembly operator copies, not source chains or residues.
    # Distinct PDB records can still belong to the same publication/experiment.
    source_targets = sorted({(str(t[1]), str(t[3]), str(t[4]), str(t[5]), str(t[6]))
                             for t in targets if len(t) >= 7})
    source_proteins = sorted({str(x[1]) for x in identity.get('binding_instances', []) if len(x) >= 2})
    evidence = (identity['pdb_id'], identity['revision'], source_targets, source_proteins)
    return compact(inst), compact(evidence)


def assess_match(identity, rows, revisions, chain_members):
    """Return conservative linkage; matching does not establish relevance.

    Strongest level is single-residue chain/assembly/receptor correspondence. The
    source table does not carry a model field, and the site residues must still be
    independently checked; neither property is asserted by this linkage alone.
    """
    pdb = identity['pdb_id'].lower()
    if pdb not in revisions:
        return 'NO_ENTRY_IN_PLINDER', 'UNKNOWN', []
    if not rows:
        return 'NO_ANNOTATION', 'UNKNOWN', []
    targets = identity.get('target_residue_instances') or []
    if not targets:
        return 'TARGET_IDENTITY_UNAVAILABLE', 'UNKNOWN', []
    codes = {str(t[6]).upper() for t in targets if len(t) >= 7}
    same_ccd = [r for r in rows if str(r.get('ligand_ccd_code') or '').upper() in codes]
    if not same_ccd:
        return 'NO_CCD_MATCH', 'UNKNOWN', []
    single = len(targets) == 1 and len(targets[0]) >= 7
    if not single:
        return 'CCD_ONLY_MULTIMEMBER', 'UNKNOWN', []
    target = targets[0]
    exact = [r for r in same_ccd if str(r.get('system_biounit_id')) == str(identity['assembly'])
             and r.get('ligand_instance_chain') == target[0]
             and r.get('ligand_asym_id') == target[1]
             and r.get('ligand_auth_id') == target[4]]
    if not exact:
        return 'CCD_ONLY', 'UNKNOWN', []
    annotations = []
    for r in exact:
        annotations.append({k: r.get(k) for k in ('system_id','ligand_id','ligand_ccd_code',
            'ligand_is_proper','ligand_is_artifact','ligand_in_artifact_list','ligand_is_ion',
            'ligand_is_cofactor_like','ligand_is_covalent','ligand__is_multi_residue',
            'ligand_member_asym_ids','ligand_instance_chains','ligand_residue_numbers',
            'ligand_protein_chains_asym_id','ligand_neighboring_residues')})
    if revisions[pdb] != str(identity['revision']):
        return 'DIFFERENT_REVISION_UNVERIFIED', 'UNKNOWN', annotations
    if identity.get('model_scope') != [1]:
        return 'MODEL_NOT_COMPARABLE', 'UNKNOWN', annotations
    proteins = {x[0] for x in identity.get('binding_instances', [])}
    checked = []
    for r in exact:
        ligand_chains = set(r.get('ligand_instance_chains') or [])
        labels = set(r.get('ligand_member_asym_ids') or [])
        receptor_chains = set(r.get('ligand_protein_chains_asym_id') or [])
        chainkeys = {(pdb, str(identity['assembly']), x) for x in ligand_chains | proteins}
        if (r.get('ligand__is_multi_residue') is False and labels == {target[1]}
            and ligand_chains == {target[0]} and len(r.get('ligand_residue_numbers') or []) == 1
            and proteins and receptor_chains == proteins and chainkeys <= chain_members):
            checked.append(r)
    if not checked:
        return 'SAME_REVISION_SCOPE_OR_SITE_AMBIGUOUS', 'UNKNOWN', annotations
    unique_ids = {r['ligand_id'] for r in checked}
    if len(unique_ids) != 1:
        return 'SAME_REVISION_MULTIPLE_LIGAND_IDS', 'CONFLICT', annotations
    labels = {(r.get('ligand_is_proper'), r.get('ligand_is_artifact'), r.get('ligand_is_ion')) for r in checked}
    if len(labels) != 1:
        return 'SAME_REVISION_INSTANCE_CHAIN_MATCH', 'CONFLICT', annotations
    proper, artifact, ion = next(iter(labels))
    if proper is True and artifact is False and ion is False:
        interpretation = 'PLINDER_PROPER_RULE_SUPPORT'
    elif artifact is True:
        interpretation = 'PLINDER_ARTIFACT_FLAG'
    elif ion is True:
        interpretation = 'PLINDER_ION_FLAG'
    else:
        interpretation = 'UNKNOWN'
    return 'SAME_REVISION_INSTANCE_CHAIN_MATCH', interpretation, annotations


def ingest_source(source, db, max_reports=None):
    c = open_source(source)
    sql = '''SELECT object_id, run_id, updated_at, error_message FROM run_item
             WHERE stage='STRUCTURE_OBSERVATION' AND object_type='observation'
             AND output_uri IS NOT NULL ORDER BY updated_at,run_id'''
    n = 0
    for sid, run, updated, raw in c.execute(sql):
        obj = json.loads(raw)
        if obj.get('system_id') != sid:
            raise ValueError('Observation system identity mismatch')
        db.execute('INSERT OR REPLACE INTO observation VALUES (?,?,?,?)', (sid, run, updated, raw))
        n += 1
        if n % 20000 == 0:
            db.commit()
    db.commit()
    print('observation_rows', n, 'unique', db.execute('select count(*) from observation').fetchone()[0], flush=True)
    db.executemany('INSERT INTO assessment(run_id,system_id,registered_status) VALUES (?,?,?)',
        c.execute("SELECT run_id,object_id,status FROM run_item WHERE stage='M3_REVIEW' AND object_type='system'"))
    db.commit()
    reports = c.execute("""SELECT DISTINCT run.id,run.report_path,run.output_manifest_hash FROM run
        JOIN run_item ON run.id=run_item.run_id WHERE run_item.stage='M3_REVIEW' ORDER BY run.id""").fetchall()
    if max_reports is not None:
        reports = reports[:max_reports]
    for number, (run, filename, expected) in enumerate(reports, 1):
        if not filename or not expected:
            raise ValueError(f'Missing sealed report path or hash: {run}')
        for result in array_stream(Path(filename), expected):
            sid = result.get('system_id')
            if not sid:
                continue
            candidate = str(result.get('candidate_key', ''))
            compact_row = {k: result.get(k) for k in ('system_id','candidate_key','disposition',
                'reasons','prior_role','role','qc','previous_qc','partner_instances','supersedes')}
            prev = db.execute('select registered_status,candidate_key from assessment where run_id=? and system_id=?',
                              (run, sid)).fetchone()
            if prev is None or prev[0] != result['disposition']:
                raise ValueError(f'Assessment disagrees with sealed registry: {run} {sid}')
            if prev[1] is None or candidate < prev[1]:
                db.execute('update assessment set document=?,candidate_key=? where run_id=? and system_id=?',
                           (compact(compact_row), candidate, run, sid))
        db.commit()
        if number % 20 == 0 or number == len(reports):
            print('sealed_reports', number, '/', len(reports), flush=True)
    if max_reports is None:
        missing = db.execute('select count(*) from assessment where document is null').fetchone()[0]
        if missing:
            raise ValueError(f'{missing} registered assessments missing from sealed reports')
    return {'observation_rows': n, 'reports': len(reports)}


def load_plinder(index):
    paths = {name: index / (name + '.parquet') for name in
             ('entry_sources','entry_biounit_chains','annotation_table')}
    revisions = {}
    for r in pq.read_table(paths['entry_sources']).to_pylist():
        revisions[r['entry_pdb_id']] = f"{r['source_mmcif_major_revision']}.{r['source_mmcif_minor_revision']}"
    cols = ['entry_pdb_id','biounit_id','chain_instance']
    chains = set()
    for batch in pq.ParquetFile(paths['entry_biounit_chains']).iter_batches(columns=cols):
        chains.update((r['entry_pdb_id'], str(r['biounit_id']), r['chain_instance']) for r in batch.to_pylist())
    cols = ['entry_pdb_id','system_biounit_id','system_id','ligand_id','ligand_ccd_code',
            'ligand_asym_id','ligand_auth_id','ligand_instance_chain','ligand_instance_chains',
            'ligand_member_asym_ids','ligand_residue_numbers','ligand__is_multi_residue',
            'ligand_protein_chains_asym_id','ligand_neighboring_residues','ligand_is_proper',
            'ligand_is_artifact','ligand_in_artifact_list','ligand_is_ion','ligand_is_cofactor_like',
            'ligand_is_covalent']
    annotations = defaultdict(list)
    total = 0
    for batch in pq.ParquetFile(paths['annotation_table']).iter_batches(columns=cols):
        for row in batch.to_pylist():
            annotations[row['entry_pdb_id']].append(row)
            total += 1
    print('plinder_rows', total, 'entries', len(annotations), 'assembly_chains', len(chains), flush=True)
    return revisions, chains, annotations, {k: {'bytes': v.stat().st_size,
        'rows': pq.ParquetFile(v).metadata.num_rows} for k, v in paths.items()}


def derive(db, index):
    revisions, chains, annotations, index_info = load_plinder(index)
    count = 0
    cursor = db.execute('select system_id,document from observation order by system_id')
    for sid, raw in cursor:
        obs = json.loads(raw)
        identity = obs['identity']
        assessments = [(r, json.loads(d)) for r, d in db.execute(
            'select run_id,document from assessment where system_id=? and document is not null', (sid,))]
        selected, state = resolve(assessments)
        if state == 'CONFLICT':
            disposition, reasons, role, prior, conflicts, qc, members = 'CONFLICT', ['incompatible_unsuperseded_decisions'], {}, None, None, None, []
        elif selected:
            disposition, reasons = selected['disposition'], selected.get('reasons') or []
            role = selected.get('role') or {}
            prior, conflicts, qc = selected.get('prior_role'), role.get('conflicts'), selected.get('qc')
            members = selected.get('partner_instances') or []
        else:
            disposition, reasons = obs['status'], obs.get('m3_reasons') or []
            role, prior, conflicts, qc, members = {}, None, None, None, []
        inst, evidence = identifiers(identity)
        level, annotation, hits = assess_match(identity, annotations.get(identity['pdb_id'].lower(), []), revisions, chains)
        args = (sid, identity['pdb_id'], str(identity['revision']), disposition, compact(reasons),
                (obs.get('type') or {}).get('primary_class'), inst, evidence, prior, role.get('status'),
                compact(conflicts), len(members) if selected else None,
                len(identity.get('target_residue_instances') or []), compact(obs.get('qc')),
                compact(qc), compact(role.get('source_evidence')), state, level, annotation,
                len(hits), compact(hits))
        db.execute('insert into result values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', args)
        count += 1
        if count % 20000 == 0:
            db.commit()
            print('derived', count, flush=True)
    db.commit()
    return index_info


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--plinder-index', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-reports', type=int, help='diagnostic only; never use for full results')
    parser.add_argument('--resume-derive', action='store_true', help='resume an interrupted derive after complete ingest')
    args = parser.parse_args()
    started = utc(); clock = time.monotonic()
    if args.resume_derive:
        if args.max_reports is not None or (args.output / 'run.json').exists():
            raise ValueError('Cannot resume a diagnostic or completed research run')
        db = sqlite3.connect(args.output / 'research.sqlite3')
        if (db.execute('select count(*) from result').fetchone()[0] or
            db.execute('select count(*) from assessment where document is null').fetchone()[0]):
            raise ValueError('Resume requires a completed ingest and no derived rows')
    else:
        db = create(args.output)
    meta = {'started_utc': started, 'snapshot': str(args.snapshot.resolve()),
            'plinder_index': str(args.plinder_index.resolve()), 'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'partial_max_reports': args.max_reports}
    try:
        if args.resume_derive:
            meta['ingest'] = {'unique_observations': db.execute('select count(*) from observation').fetchone()[0],
                              'assessments': db.execute('select count(*) from assessment').fetchone()[0],
                              'resumed_after_complete_ingest': True}
        else:
            meta['ingest'] = ingest_source(args.snapshot, db, args.max_reports)
        if args.max_reports is None:
            meta['plinder'] = derive(db, args.plinder_index)
        meta['status'] = 'COMPLETE' if args.max_reports is None else 'PARTIAL_DIAGNOSTIC'
        meta['ended_utc'] = utc(); meta['elapsed_s'] = round(time.monotonic()-clock, 2)
        (args.output / 'run.json').write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding='utf-8')
        print(json.dumps(meta, ensure_ascii=False), flush=True)
    finally:
        db.close()


if __name__ == '__main__':
    main()
