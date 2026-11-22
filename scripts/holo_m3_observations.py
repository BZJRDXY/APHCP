"""Index verified structure files in the existing APHCP staging registry."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import time

from holo_m3_register import (
    canonical_id, connect, context_gate, database, digest, file_ref,
    frozen_registry_root, jsonl, read, read_only_path, require, seed_registry,
    sha, utc_now, write,
)
from aphcp_control.registry import backup


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
    statuses = {}
    for run in db.execute('SELECT DISTINCT r.* FROM run r JOIN run_item i ON i.run_id=r.id WHERE i.stage=?',
                          ('M3_REVIEW',)):
        report = read_only_path(run['report_path'], (runs,))
        require(sha(report) == run['output_manifest_hash'], 'staged_review_hash_differs')
        for result in read(report):
            k = key(result)
            previous = statuses.setdefault(k, result)
            require(previous['disposition'] == result['disposition'] and
                    previous.get('system_id') == result.get('system_id'), 'conflicting_staged_review')
    return statuses


def observation(row, entry, review, original_role_evidence, roots):
    if (review.get('disposition') == 'OUT_OF_SCOPE' or
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
    require(not review.get('system_id') or review['system_id'] == sid, 'observation_and_staged_identity_differ')
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
        'qc': {fld: row.get(fld) for fld in ('coordinate_identity', 'experimental_fidelity', 'topology', 'partner_scope_qc')},
        'files': files,
    }


def add_batch(db, run, seal, started):
    label = run.name
    db.execute('INSERT INTO config_version VALUES (?,?,?)', (label, digest(seal), str(run/'frozen_catalog.json')))
    db.execute('''INSERT INTO batch(id,module_id,status,source_release,input_manifest,config_version,candidate_count)
                  VALUES (?,'holo','OBSERVATIONS_STAGED','wwPDB-fixed-mmCIF',?,?,?)''',
               (label, str(run/'frozen_catalog.json'), label, seal['expected']))
    db.execute('INSERT INTO task VALUES (?,?,?,?,?)', (label, label, 'Holo 结构观察目录', label, 'OBSERVATIONS_STAGED'))
    db.execute('''INSERT INTO run(id,batch_id,task_id,status,git_commit,code_version,config_version,config_hash,
                  source_release,input_manifest_hash,environment_lock_hash,host,start_time_utc,report_path,input_records)
                  VALUES (?,?,?,'OBSERVATIONS_STAGED',?,?,?,?,?,?,?,?,?,?,?)''',
               (label, label, label, seal['git_head'], VERSION, label, digest(seal), 'wwPDB-fixed-mmCIF', digest(seal),
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
        rank = {'OBSERVED': 0, 'PREPARATION_USE_RECORDED': 1, 'REVIEW': 2, 'REGISTERED': 3}
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
                   (utc_now(), duration, len(ingested), canonical_count, labels['OBSERVED']+labels['REVIEW'],
                    json.dumps(receipt), run.name))
    write(run/'catalog_receipt.json', receipt)
    print(json.dumps(receipt))


def lookup(args):
    run = args.run.resolve()
    root = run.parent.parent
    with connect(database(frozen_registry_root(root, run)), readonly=True) as db:
        if args.pdb:
            matches = db.execute('SELECT error_message FROM run_item WHERE run_id=? AND stage=? AND pdb_id=? ORDER BY object_id',
                                 (run.name, STAGE, args.pdb.upper()))
        else:
            matches = db.execute('SELECT error_message FROM run_item WHERE run_id=? AND stage=? AND object_id=?',
                                 (run.name, STAGE, args.system))
        rows = [json.loads(x[0]) for x in matches]
    if args.verify_files:
        for row in rows:
            for asset in row.get('files', {}).values():
                file_ref(asset['path'], (root/'runs',), asset.get('sha256'))
    print(json.dumps(rows, ensure_ascii=False))


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
