"""Index verified structure files in the existing APHCP staging registry."""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import sys
import time

from holo_m3_register import (
    canonical_id, connect, context_gate, database, digest, file_ref,
    frozen_registry_root, jsonl, read, read_only_path, recorded_assessments, require, seed_registry,
    sha, utc_now, write,
)
from aphcp_control.registry import backup, recover_empty_registry
from holo_m3_decisions import decision_view, resolve as resolve_decisions


STAGE = 'STRUCTURE_OBSERVATION'
VERSION = 'holo-m3-observations-v1'
INPUT_TABLES = ('baseline_system_labels.jsonl', 'supplemental_system_labels.jsonl')


def key(row):
    return str(row['batch']), row['pdb_id'], row['candidate_key']


def input_observations(inputs, selected):
    replacements = {key(row): row for row in jsonl(inputs[1])}
    for baseline in jsonl(inputs[0]):
        row = replacements.pop(key(baseline), baseline)
        if row['prior_status'] == 'EXPORTED_READBACK' and (not selected or key(row) in selected):
            yield row
    for row in replacements.values():
        if row['prior_status'] == 'EXPORTED_READBACK' and (not selected or key(row) in selected):
            yield row


def review_records(db, runs):
    by_identity = defaultdict(list)
    candidate_identities = defaultdict(set)
    for rid, path, result in recorded_assessments(db, runs):
        k = key(result)
        identity = result.get('system_id') or k
        candidate_identities[k].add(identity)
        by_identity[identity].append((rid, path, result))
    resolved = {identity: resolve_decisions(group) for identity, group in by_identity.items()}
    # A changed partner composition is a different identity; keep both rather than pick by run order.
    return {k: {identity: resolved[identity] for identity in ids} for k, ids in candidate_identities.items()}


def observation(row, entry, review, original_role_evidence, roots):
    if (review.get(key(row), {}).get('disposition') == 'OUT_OF_SCOPE' or
            context_gate({'role_evidence': {'matches': [{'text': entry['source_title']}]}}) or
            context_gate({'role_evidence': original_role_evidence})):
        return None, {'key': key(row), 'status': 'OUT_OF_SCOPE'}

    mapping_path = read_only_path(row['mapping_path'], roots)
    mapping_bytes = mapping_path.read_bytes()
    mapping = json.loads(mapping_bytes)
    artifact = mapping['artifact']
    require(all(str(mapping[k]) == str(row[k]) for k in
                ('candidate_key', 'pdb_id', 'revision', 'assembly', 'source_sha256')) and
            str(mapping['source_model']) == str(row['model']) and
            artifact['sha256'] == row['system_sha256'], 'observation_mapping_conflicts_with_ledger')
    cif = read_only_path(row['system_path'], roots)
    raw = read_only_path(row['source_path'], roots)
    require(cif.stat().st_size == artifact['bytes'], 'observation_cif_size_changed')
    instances = mapping['instances']
    bound = sorted(mapping['binding_protein_instances'])
    identity = {
        'pdb_id': row['pdb_id'], 'revision': row['revision'], 'source_sha256': row['source_sha256'],
        'assembly': row['assembly'], 'model_scope': [mapping['source_model']],
        'binding_instances': [(label, instances[label]['label_asym'], instances[label]['operator']) for label in bound],
        'target_residue_instances': sorted(mapping['partner_member_identities']),
        'altloc_scope': artifact['altloc_policy'],
    }
    sid = canonical_id('holo:' + row['pdb_id'], identity)
    review = review.get(sid, {})
    cls = row['classification']
    status = review.get('disposition') or ('PREPARATION_USE_RECORDED' if row['intake_after'] == 'EXCLUDED' else 'OBSERVED')
    files = {
        'source_mmcif': {'path': str(raw), 'sha256': row['source_sha256']},
        'system_cif': {'path': str(cif), 'sha256': row['system_sha256']},
        'source_mapping': {'path': str(mapping_path), 'sha256': hashlib.sha256(mapping_bytes).hexdigest()},
    }
    if review.get('files', {}).get('complete_assembly'):
        files['complete_assembly'] = review['files']['complete_assembly']
    role_review = review.get('role', {})
    return sid, {
        'system_id': sid, 'identity': identity, 'key': key(row), 'status': status,
        'evidence_layer': row['evidence_layer'], 'type': cls, 'role': role_review.get('status', row['role_after']),
        'role_evidence': role_review.get('source_evidence', original_role_evidence),
        'm3_reasons': review.get('reasons', []),
        'current_decision_id': review.get('current_decision_id'), 'decision_state': review.get('decision_state'),
        'm3_readback_qc': review.get('qc'),
        'qc': {fld: row.get(fld) for fld in ('coordinate_identity', 'experimental_fidelity', 'topology', 'partner_scope_qc')},
        'files': files,
    }


def add_batch(db, run, seal, started):
    label = seal.get("label", run.name)
    db.execute('INSERT INTO config_version VALUES (?,?,?)', (label, digest(seal), str(run/'frozen_catalog.json')))
    db.execute('''INSERT OR IGNORE INTO batch(id,module_id,status,source_release,input_manifest,config_version,candidate_count)
                  VALUES (?,'holo','OBSERVATIONS_STAGED','wwPDB-fixed-mmCIF',?,?,?)''',
               (run.name, str(run/'frozen_catalog.json'), label, seal['expected']))
    db.execute('INSERT INTO task VALUES (?,?,?,?,?)', (label, run.name, 'Holo 结构观察目录', label, 'OBSERVATIONS_STAGED'))
    db.execute('''INSERT INTO run(id,batch_id,task_id,status,git_commit,code_version,config_version,config_hash,
                  source_release,input_manifest_hash,environment_lock_hash,host,start_time_utc,report_path,input_records)
                  VALUES (?,?,?,'OBSERVATIONS_STAGED',?,?,?,?,?,?,?,?,?,?,?)''',
               (label, run.name, label, seal['git_head'], VERSION, label, digest(seal), 'wwPDB-fixed-mmCIF', digest(seal),
                seal['runtime_hash'], 'Pro5000', started, str(run/'catalog_receipt.json'), seal['expected']))


def add_observation(db, label, sid, row, entry, stamp):
    _, pdb, ck = row['key']
    db.execute('INSERT OR IGNORE INTO pdb_entry VALUES (?,?,?)', (pdb, entry['source_title'] if sid else None, None))
    if not sid:
        db.execute('INSERT INTO run_item VALUES (?,?,?,?,?,?,?,?,?,?)',
                   (label, 'observation', 'candidate:'+ck, pdb, STAGE, row['status'], None,
                    json.dumps(row), None, stamp))
        return False

    ident = row['identity']
    version = f"wwpdb:{pdb}:{ident['revision']}"
    db.execute('INSERT OR IGNORE INTO entry_version VALUES (?,?,?,?,?,?)',
               (version, pdb, 'wwPDB', 'fixed-archive', pdb, ident['revision']))
    db.execute('INSERT OR IGNORE INTO work_item(id,entry_version_id,kind,label) VALUES (?,?,?,?)',
               (sid, version, 'system', ck))
    previous = db.execute('SELECT error_message FROM run_item WHERE run_id=? AND object_type=? AND object_id=? AND stage=?',
                          (label, 'observation', sid, STAGE)).fetchone()
    if previous:
        saved = json.loads(previous['error_message'])
        require(saved['identity'] == json.loads(json.dumps(ident)), 'observation_key_collision')
        rank = {'OBSERVED': 0, 'PREPARATION_USE_RECORDED': 1, 'REVIEW': 2, 'REGISTERED': 3, 'CONFLICT': 4, 'TECHNICAL_FAILURE': -1}
        if rank[row['status']] > rank[saved['status']]:
            saved, row = row | {'aliases': saved.get('aliases', [])}, saved
        saved.setdefault('aliases', []).append({k: row[k] for k in ('key', 'evidence_layer', 'status', 'files')})
        db.execute('''UPDATE run_item SET status=?,error_message=?,output_uri=?
                      WHERE run_id=? AND object_type=? AND object_id=? AND stage=?''',
                   (saved['status'], json.dumps(saved, ensure_ascii=False), saved['files']['system_cif']['path'],
                    label, 'observation', sid, STAGE))
        return False
    db.execute('INSERT INTO run_item VALUES (?,?,?,?,?,?,?,?,?,?)',
               (label, 'observation', sid, pdb, STAGE, row['status'], None,
                json.dumps(row, ensure_ascii=False), row['files']['system_cif']['path'], stamp))
    return True


def build(args):
    root = args.root.resolve()
    run = args.run.resolve()
    dbpath = database(frozen_registry_root(root, run, create=True))
    seed_path, seed_hash = seed_registry(root, run, args)
    inputs = tuple(args.inputs.resolve()/name for name in INPUT_TABLES)
    selected = {tuple(k) for k in read(args.only)} if args.only else set()
    entries_file = args.entries.resolve()
    roles_file = args.roles.resolve()
    seal = {
        'version': VERSION, 'inputs': {str(p): sha(p) for p in (*inputs, entries_file, roles_file)},
        'seed': {'path': str(seed_path), 'sha256': seed_hash}, 'code_sha256': sha(__file__),
        'git_head': args.git_head, 'runtime_hash': digest(sys.version),
        'selection': [list(k) for k in sorted(selected)], 'expected': args.expected,
    }
    frozen = run/'frozen_catalog.json'
    if frozen.exists():
        require(read(frozen) == seal, 'catalog_inputs_changed_use_new_run')
    else:
        write(frozen, seal)
    recover_empty_registry(dbpath,(run/'catalog_receipt.json',run/'phase_receipt.json'))
    if not dbpath.exists():
        backup(seed_path, dbpath)
    entries = {row['pdb_id']: row for row in jsonl(entries_file)}
    stamp = utc_now()
    t0 = time.monotonic()
    labels = Counter()
    ingested = set()
    canonical_count = 0
    with connect(dbpath) as db:
        db.execute('BEGIN IMMEDIATE')
        prior = db.execute('SELECT config_hash,error_summary FROM run WHERE id=?', (run.name,)).fetchone()
        if prior:
            require(prior['config_hash'] == digest(seal), 'catalog_registry_seal_changed')
            receipt = json.loads(prior['error_summary'])
            if not (run/'catalog_receipt.json').is_file():
                write(run/'catalog_receipt.json', receipt)
            print(json.dumps(receipt | {'reused': True}))
            return
        known = review_records(db, root/'runs')
        role_evidence = {key(row): row['role_evidence_after'] for row in jsonl(roles_file)
                         if not selected or key(row) in selected}
        add_batch(db, run, seal, stamp)
        for row in input_observations(inputs, selected):
            k = key(row)
            require(k not in ingested, 'duplicate_input_observation')
            ingested.add(k)
            e = entries[k[1]]
            # A scoped structural observation keeps the output that was actually verified.
            sid, record = observation(row, e, known.get(k, {}), role_evidence.get(k), (root/'runs',))
            labels[record['status']] += 1
            canonical_count += add_observation(db, run.name, sid, record, e, stamp)
        require(len(ingested) == args.expected and (not selected or selected == ingested), 'catalog_observation_count_differs')
        require(not db.execute('PRAGMA foreign_key_check').fetchall() and
                db.execute('PRAGMA integrity_check').fetchone()[0] == 'ok', 'catalog_integrity_failed')
        require(db.execute('SELECT COUNT(*) FROM release').fetchone()[0] == 0, 'catalog_contains_release')
        row_count = db.execute('SELECT COUNT(*) FROM run_item WHERE run_id=? AND stage=?', (run.name, STAGE)).fetchone()[0]
        require(row_count == canonical_count + labels['OUT_OF_SCOPE'], 'catalog_registry_count_differs')
        duration = round(time.monotonic() - t0, 3)
        receipt = {'status': 'OBSERVATIONS_STAGED', 'source_rows': len(ingested), 'queryable_systems': canonical_count,
                   'outcomes': dict(labels), 'registry': str(dbpath), 'seed_sha256': seed_hash, 'elapsed_s': duration}
        db.execute('''UPDATE run SET end_time_utc=?,elapsed_seconds=?,input_records=?,success_records=?,
                      uncertain_records=?,error_summary=? WHERE id=?''',
                   (utc_now(), duration, len(ingested), canonical_count, labels['OBSERVED']+labels['REVIEW']+labels['CONFLICT'],
                    json.dumps(receipt), run.name))
    write(run/'catalog_receipt.json', receipt)
    print(json.dumps(receipt))


def append_increment(root, run, seed, evidence_path, chemistry_path, git_head):
    """Append one increment to its isolated copy of the previous catalog."""
    from types import SimpleNamespace
    from holo_m3_register import clone_sealed_registry

    started = time.monotonic()
    dbpath = database(frozen_registry_root(root, run, create=True))
    seed_path, seed_hash = seed_registry(root, run, SimpleNamespace(seed_run=seed))
    recover_empty_registry(dbpath,(run/'phase_receipt.json',run/'catalog_receipt.json',run/'registration_receipt.json'))
    if not dbpath.exists():
        clone_sealed_registry(seed_path, dbpath)
    label = run.name + '-observations'
    data = read(evidence_path)
    reviews_path = run/'review_results.json'
    reviews = {key(r): r for r in read(reviews_path)} if reviews_path.exists() else {}
    classified = {key(r): r['classification'] for r in jsonl(chemistry_path)}
    require(set(classified) == {key(r) for r in data['candidates']}, 'increment_catalog_chemistry_keys_differ')
    seal = {'version': VERSION, 'label': label, 'expected': len(data['candidates']), 'git_head': git_head,
            'runtime_hash': digest(sys.version), 'seed_sha256': seed_hash, 'code_sha256': sha(__file__),
            'inputs': {str(p): sha(p) for p in (evidence_path, chemistry_path, *([reviews_path] if reviews_path.exists() else []))}}
    f = run/'frozen_catalog.json'
    if f.exists():
        require(read(f) == seal, 'catalog_inputs_changed_use_new_run')
    else:
        write(f, seal)
    with connect(dbpath) as db:
        db.execute('BEGIN IMMEDIATE')
        prior = db.execute('SELECT config_hash,error_summary FROM run WHERE id=?', (label,)).fetchone()
        if prior:
            require(prior[0] == digest(seal), 'catalog_registry_seal_changed')
            receipt = json.loads(prior[1])
        else:
            add_batch(db, run, seal, utc_now())
            labels = Counter(); canonical_count = 0; entries = {}
            stamp = utc_now()
            for row in data['candidates']:
                k = key(row); p = k[1]; review = reviews.get(k, {})
                if p not in entries:
                    discovery = file_ref(row['discovery_path'], (root/'runs',), row['discovery_sha256'])
                    entries[p] = read(discovery['path'])
                entry = entries[p]
                if row['prior_status'] == 'EXPORTED_READBACK':
                    known = {k: review, review.get('system_id'): review}
                    sid, record = observation(row | {'classification': classified[k]}, entry, known,
                                              row['role_evidence_after'], (root/'runs',))
                    if sid:
                        require(not review.get('system_id') or review['system_id'] == sid, 'review_and_observation_identity_conflict')
                        if review.get('disposition') == 'TECHNICAL_FAILURE':
                            record.update(status='TECHNICAL_FAILURE', m3_reasons=review['reasons'])
                else:
                    status = 'TECHNICAL_FAILURE' if row['prior_status'] == 'EXPORT_TECHNICAL_FAILURE' else 'NOT_EXPORTED'
                    sid, record = None, {'key': k, 'status': status, 'reason': row.get('export_reason'), 'type': classified[k]}
                labels[record['status']] += 1
                canonical_count += add_observation(db, label, sid, record, entry, stamp)
            for p, evidence in data['entry_evidence'].items():
                db.execute('INSERT OR IGNORE INTO pdb_entry VALUES (?,?,?)', (p, None, None))
                db.execute('INSERT INTO run_item VALUES (?,?,?,?,?,?,?,?,?,?)',
                           (label, 'entry', p, p, STAGE, evidence['state'], None,
                            json.dumps({'pdb_id': p, 'status': evidence['state'], 'entry_evidence': evidence,
                                        'source_run': data['source_run']}), None, stamp))
            require(not db.execute('PRAGMA foreign_key_check').fetchone(), 'catalog_foreign_key_error')
            duration = round(time.monotonic()-started, 3)
            receipt = {'status': 'OBSERVATIONS_STAGED', 'source_rows': len(data['candidates']),
                       'input_entries': len(data['entries']), 'queryable_systems': canonical_count, 'outcomes': dict(labels),
                       'registry': str(dbpath), 'seed_sha256': seed_hash, 'elapsed_s': duration}
            db.execute("UPDATE run SET end_time_utc=?,elapsed_seconds=?,success_records=?,uncertain_records=?,error_summary=? WHERE id=?",
                       (utc_now(), duration, canonical_count, labels['REVIEW']+labels['OBSERVED'], json.dumps(receipt), label))
    if not (run/'catalog_receipt.json').exists():
        write(run/'catalog_receipt.json', receipt)
    return receipt


def lookup_records(run, *, pdb=None, system=None, verify_files=False, limit=None):
    run = Path(run).resolve()
    root = run.parent.parent
    require(bool(pdb) != bool(system), 'specify_pdb_or_system')
    with connect(database(frozen_registry_root(root, run)), readonly=True) as db:
        column, value = ('pdb_id', pdb.upper()) if pdb else ('object_id', system)
        grouped = defaultdict(list)
        for item in db.execute(f'''SELECT run_id,object_id,object_type,error_message FROM run_item
            WHERE stage=? AND {column}=? ORDER BY object_type DESC, object_id, updated_at, run_id''', (STAGE, value)):
            grouped[item['object_type'], item['object_id']].append((item['run_id'], json.loads(item['error_message'])))
        decisions = defaultdict(list)
        if any('system_id' in records[0][1] for records in grouped.values()):
            for rid, path, result in recorded_assessments(db, root/'runs', pdb=pdb.upper() if pdb else None, system=system):
                if result.get('system_id'):
                    decisions[result['system_id']].append((rid, path, result))
        rows = []
        for records in grouped.values():
            if limit is not None and len(rows) >= limit:
                break
            rid, row = records[-1]
            if row.get('system_id') in decisions:
                row = decision_view(row, decisions[row['system_id']])
            if len(records) > 1 or row.get('entry_evidence'):
                row = row | {'observations': [{'run': run_id, 'status': r['status'], 'key': r.get('key'),
                                               'entry_evidence': r.get('entry_evidence')} for run_id, r in records]}
            if verify_files:
                for asset in row.get('files', {}).values():
                    file_ref(asset['path'], (root/'runs',), asset.get('sha256'))
            rows.append(row)
    return rows


def lookup(args):
    print(json.dumps(lookup_records(args.run, pdb=args.pdb, system=args.system,
                                    verify_files=args.verify_files), ensure_ascii=False))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    build_parser = commands.add_parser('build')
    for name in ('root', 'run', 'seed-run', 'inputs', 'entries', 'roles'):
        build_parser.add_argument('--'+name, required=True, type=Path)
    build_parser.add_argument('--only', type=Path)
    build_parser.add_argument('--expected', type=int, required=True)
    build_parser.add_argument('--git-head', required=True)
    lookup_parser = commands.add_parser('lookup')
    lookup_parser.add_argument('--run', type=Path, required=True)
    query = lookup_parser.add_mutually_exclusive_group(required=True)
    query.add_argument('--pdb')
    query.add_argument('--system')
    lookup_parser.add_argument('--verify-files', action='store_true')
    args = parser.parse_args()
    (build if args.command == 'build' else lookup)(args)
