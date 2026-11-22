"""Verify and register the completed M1 metadata run without writing structure tables."""
from collections import Counter
import csv, hashlib, json, shutil, platform, socket, sqlite3, sys
from datetime import datetime, timezone
from pathlib import Path
import pyarrow as pa

from m1_holo_config import RELEASE_FIELD, metadata_complete, read_config, window

run=Path(sys.argv[1]).resolve(); inv=run/'inventory'; budget=run/'budget'; root=run.parent.parent
read=lambda p:json.loads((run/p).read_text(encoding='utf8'))
jsonl=lambda p:[json.loads(l) for l in (run/p).read_text(encoding='utf8').splitlines()]
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def sha_json(x):return hashlib.sha256(json.dumps(x,sort_keys=True).encode()).hexdigest()
def write(name,payload):
    p=run/name; tmp=p.with_suffix('.tmp'); tmp.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n',encoding='utf8'); tmp.replace(p)
def loadcsv(p):
    with (run/p).open(encoding='utf8',newline='') as f:return list(csv.DictReader(f))
def total(paths):return sum(p.stat().st_size for p in paths if p.is_file())

Q = read_config(run)
start_day, end_day, cutoff = window(Q)
start, end = start_day.isoformat(), end_day.isoformat()
S=read('search/summary.json'); G=read('data/summary.json')
if not metadata_complete(S, G):
    raise SystemExit('Unresolved Search or Entry metadata gaps; no completion receipt or registry write')
I=read('inventory/inventory_summary.json'); B=read('budget/budget_summary.json'); C=I['counts']
E=jsonl('inventory/candidate_entries.jsonl'); L=jsonl('inventory/plinder_ligand_rows.jsonl'); N=read('budget/m2_nominations.json'); extra=read('inventory/plinder_period_only_reconcile.json')
entry_ids = {row['pdb_id'] for row in E}
assert [r['pdb_id'] for r in E]==S['entry_ids'] and S['concordant'] and S['total_count']==len(E)==G['requested']
expected_unavailable = G.get('confirmed_unavailable', {})
assert expected_unavailable.keys() <= entry_ids
assert I['entry_metadata_sha256'] == sha(run/'data'/'entries.json') == G.get('entries_sha256', I['entry_metadata_sha256'])
assert not any(row.get('_technical_missing') for row in extra), 'Unresolved PLINDER-only Data API requests'
assert len(L)==len({(x['system_id'],x['ligand_id']) for x in L})==C['plinder_annotation_rows']
assert len({x['system_id'] for x in L})==C['plinder_candidate_systems'] and len({x['pdb_id'] for x in L})==C['plinder_annotation_entries']
assert C['plinder_annotation_entries'] + C['plinder_no_annotation_entries'] == len(E)
shared_period_entries = sum(
    row['plinder_index_presence'] and start <= (row['plinder_metadata_release_day'] or '') < end
    for row in E
)
assert C['plinder_same_period_annotated_entries'] == shared_period_entries + C['plinder_only_calendar_subset']
assert len(extra) == C['plinder_only_calendar_metadata']
extra_statuses = Counter(row['rcsb_protein_eligibility'] for row in extra)
if cutoff:
    cutoff_rows = loadcsv(f'inventory/boundary_{cutoff.isoformat()}.csv')
    assert len(cutoff_rows) == C['cutoff_entry_count']
    assert all(row['first_release_date'][:10] == cutoff.isoformat() for row in cutoff_rows)
files = loadcsv('budget/monthly_file_plan.csv')
assert len(files) + B['files_without_current_revision'] == len(E)
assert len(files) == len({(x['pdb_id'], x['revision'], x['file_type']) for x in files})
expected_files = {(x['pdb_id'], x['rcsb_latest_revision']) for x in E if x['rcsb_latest_revision']}
assert {(x['pdb_id'], x['revision']) for x in files} == expected_files
assert len(N) == B['m2_nomination_count'] == len({row['pdb_id'] for row in N})
assert not any(x['selected_revision'] for x in E) and all(x['structure_qc']=='NOT_RUN' for x in L)
A=loadcsv('inventory/rcsb_assemblies.csv'); entity=loadcsv('inventory/rcsb_entities.csv'); ids={x['pdb_id'] for x in E}; H=read('inventory/index_review.json'); aset={x['pdb_id'] for x in H}
assert len({(x['pdb_id'],x['assembly_id']) for x in A})==C['rcsb_assembly_pairs']
assert all(r['pdb_id'] in ids for r in entity)
source_revision_diff=[x['pdb_id'] for x in E if x['rcsb_latest_revision']!=x['plinder_source_revision_if_available']]
metadata_notices={r['pdb_id'] for r in E if r['metadata_review']}; index_review=set(source_revision_diff)|aset
raw_chain_totals=Counter()
for x in entity:
    k=x['entity_kind']; kind=('protein' if x['polymer_type']=='Protein' else 'other_polymer') if k=='polymer' else k
    raw_chain_totals[kind]+=len(json.loads(x['label_asym_ids']))
assembly_instances={}
for field in ['polymer_entity_instance_count_protein','nonpolymer_entity_instance_count','branched_entity_instance_count']:
    vals=[json.loads(x['size']).get(field) for x in A if x['size']]
    assembly_instances[field]={'assemblies_reporting':sum(v is not None for v in vals),'reported_instance_count_sum':sum(v for v in vals if v is not None)}
assembly_generator_coverage=sum(bool(json.loads(x['generators'])) for x in A)

def phase(id_suffix,name,started,finished,entries,uncertain,sources,inputs,paths,meta_transfer=None,memory=None):
    config_path=run/'inputs'/'query.json'; stage_paths=sorted(p for p in paths if p.is_file()); stage_outputs={str(p.relative_to(run)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in stage_paths}; start=datetime.fromisoformat(started); stop=datetime.fromisoformat(finished)
    return {'id':f'{Q["run_id"]}-{id_suffix}{recovery_suffix}','task':f'{Q["batch_id"]}-{id_suffix}','name':name,'status':'SUCCEEDED','git_commit':Q['git_head'],'code_version':'sha256:'+sha_json({p.name:sha(p) for p in code_dir.glob('*.py')}),'config_version':Q['config_version'],'config_hash':sha(config_path),'source_release':sources,'input_manifest_hash':sha_json({str(p.relative_to(run)):sha(p) for p in inputs}),'output_manifest_hash':sha_json(stage_outputs),'outputs':stage_outputs,'host':socket.gethostname(),'start_time_utc':started,'end_time_utc':finished,'elapsed_seconds':round((stop-start).total_seconds(),3),'input_records':entries,'success_records':entries-uncertain,'failed_records':0,'uncertain_records':uncertain,'bytes_written':sum(x['bytes'] for x in stage_outputs.values()),'downloaded_bytes':meta_transfer,'peak_memory_bytes':memory}

def response_bytes(paths):
    return sum((read(p).get('attempts') or [{}])[-1].get('bytes',0) for p in paths)

input_dir=run/'inputs'; code_dir=Path(__file__).resolve().parent
recovery_suffix='-recovery' if (run/'superseded_run_receipt.json').is_file() else ''
read_paths=list((run/'search').glob('page*.json'))+list((run/'search').glob('count*.json'))+list((run/'data').glob('page*.json'))+list((run/'data').glob('retry*.json'))
period_response = inv / 'plinder_period_only_data_api.json'
period_requests = [period_response, *inv.glob('period_retry_*.json'), *inv.glob('period_holdings_*.json')] if extra else []
index_finished = max((read(p.relative_to(run))['attempts'][-1]['finished_utc'] for p in period_requests), default=I['finished_utc'])
period_response_bytes = response_bytes(period_requests)
R=[phase('rcsb','M1.2 RCSB inventory',S['started_utc'],G['finished_utc'],len(E),len(metadata_notices),'RCSB archive (live snapshot)',[input_dir/'query.json',code_dir/'m1_holo_fetch.py'],list((run/'search').glob('*.json'))+list((run/'data').glob('*.json')),response_bytes(read_paths)),
   phase('index','M1.3 PLINDER cross-reference and coverage', I['started_utc'], index_finished, len(E), len(index_review), f'PLINDER {Q["plinder_release"]} + RCSB', [run/'data'/'entries.json',code_dir/'m1_holo_join.py'], list(inv.glob('*')), period_response_bytes, I['peak_rss_kib']*1024),
   phase('budget','M1.4–M1.5 取数预算与M2提名' if N else 'M1.4 取数预算（独立盘点）',B['started_utc'],B['finished_utc'],len(files),0,'Fixed archive HEAD + existing Holo pilot receipts',[inv/'candidate_entries.csv',inv/'plinder_ligand_rows.csv',code_dir/'m1_holo_budget.py'] + ([input_dir/'m2_selection.json'] if (input_dir/'m2_selection.json').is_file() else []),list(budget.glob('*')),0)]
for st in R:write(st['id']+'.json',st)
lock={'python':platform.python_version(),'pyarrow':pa.__version__,'executable':sys.executable,'mode':'existing pilot environment'}; lock_hash=sha_json(lock)

receipt = {
    'batch_id': Q['batch_id'],
    'run_directory': str(run),
    'status': 'INVENTORY_COMPLETE',
    'window': {'start_inclusive': start, 'end_exclusive': end, 'field': RELEASE_FIELD},
    'cutoff_date': cutoff.isoformat() if cutoff else None,
    'input_query': str(input_dir / 'query.json'),
    'environment': lock,
    'environment_lock_hash': lock_hash,
    'stage_receipts': [stage['id'] + '.json' for stage in R],
    'stage_ids': [stage['id'] for stage in R],
    'source_revision_diff_ids': source_revision_diff,
    'entry_metadata': {'requested':G['requested'],'received':G['received'],'missing_ids':G['missing_ids'],'confirmed_unavailable':expected_unavailable},
    'assembly_or_instance_review_ids': sorted(aset),
    'plinder_period_only': {
        'total': len(extra),
        'annotated': sum(int(row['annotation_rows']) > 0 for row in extra),
        'eligibility_counts': dict(extra_statuses),
    },
    'rcsb_asu_label_counts_by_kind': raw_chain_totals,
    'rcsb_assembly_declared_instances': assembly_instances,
    'rcsb_assemblies_with_generator_rows': assembly_generator_coverage,
    'source_metadata_snapshot_bytes': sum(stage['bytes_written'] for stage in R),
    'data_filesystem_available_bytes_at_receipt': shutil.disk_usage(run).free,
    'coordinate_downloads': 0,
    'structure_tables_written': [],
    'next_step': 'Wait for separate M2 authorization and check any unresolved archive locators.',
}

if '--no-register' not in sys.argv[2:]:
    db_file=root/'registry'/'aphcp_registry.sqlite3'
    with sqlite3.connect(db_file.as_uri() + '?mode=rw', uri=True) as db:
        db.execute('PRAGMA foreign_keys=ON'); db.execute('BEGIN IMMEDIATE')
        batch = Q['batch_id']
        prior = db.execute(
            'SELECT id,module_id,window_start,window_end FROM batch '
            'WHERE id=? OR (module_id=? AND window_start=? AND window_end=?)',
            (batch, 'holo', start, end),
        ).fetchall()
        if prior and prior != [(batch, 'holo', start, end)]:
            raise RuntimeError(f'Existing batch or window conflicts with this run: {prior}')
        protected=['pdb_entry','entry_version','asset_file','work_item','system_member','run_item','release']
        before={t:db.execute(f'SELECT count(*) FROM {t}').fetchone()[0] for t in protected}
        db.execute('''INSERT INTO batch(id,module_id,status,source_release,window_start,window_end,input_manifest,config_version,git_commit,candidate_count,next_step)
            VALUES (?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(id) DO UPDATE SET status=excluded.status,candidate_count=excluded.candidate_count,next_step=excluded.next_step''', (batch,'holo','INVENTORY_COMPLETE',f'RCSB live; PLINDER {Q["plinder_release"]}',start,end,str(inv/'candidate_entries.csv'),Q['config_version'],Q['git_head'],len(E),'等待单独授权 M2'))
        for r in R:
            task_key=sha_json([batch,r['task'],Q['config_version']]); db.execute('INSERT INTO task(id,batch_id,name,task_key,status) VALUES (?,?,?,?,?) ON CONFLICT(id) DO NOTHING',(r['task'],batch,r['name'],task_key,'SUCCEEDED'))
            mapping={**r,'batch_id':batch,'task_id':r['task'],'environment_lock_hash':lock_hash,'report_path':str(run/(r['id']+'.json')),'peak_memory':r['peak_memory_bytes']}
            columns='id batch_id task_id status git_commit code_version config_version config_hash source_release input_manifest_hash output_manifest_hash environment_lock_hash host start_time_utc end_time_utc elapsed_seconds input_records success_records failed_records uncertain_records bytes_written downloaded_bytes peak_memory report_path'.split()
            db.execute('INSERT INTO run('+','.join(columns)+') VALUES ('+','.join(':'+k for k in columns)+') ON CONFLICT(id) DO NOTHING',mapping)
        run_ids = [stage['id'] for stage in R]
        placeholders = ','.join('?' for _ in run_ids)
        check = db.execute(
            f'SELECT id,status,output_manifest_hash FROM run WHERE batch_id=? AND id IN ({placeholders})',
            [batch, *run_ids],
        ).fetchall()
        assert sorted(check)==sorted((r['id'],r['status'],r['output_manifest_hash']) for r in R)
        assert before=={t:db.execute(f'SELECT count(*) FROM {t}').fetchone()[0] for t in protected}
        receipt['registry']={'location':str(db_file),'batch_count':len(db.execute('SELECT id FROM batch WHERE id=?',(batch,)).fetchall()),'runs':[x[0] for x in check],'structure_table_counts_unchanged':before}
receipt['finished_utc']=datetime.now(timezone.utc).isoformat(); write('run_receipt.json',receipt); print(json.dumps(receipt,ensure_ascii=False),flush=True)
