"""Continue M3 staging over the accepted M1 windows using the existing increment runner."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from contextlib import ExitStack
from csv import DictReader
import datetime as dt
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from holo_m2_blind import raw_path, fetch
from holo_archive_runtime import Archive, sync_path
from holo_m2_final_selection import SENSITIVE
from holo_m3_increment import once
from holo_m3_register import (digest, file_ref, frozen_registry_root, jsonl,
                              read, require, sha, utc_now, write)

GIB = 1 << 30
SPEC = ('pdb_id', 'revision', 'source_url', 'file_type')
PDB = re.compile(r'[0-9][A-Z0-9]{3}')
REV = re.compile(r'\d+\.\d+')


def identity(row):
    return row['pdb_id'], row.get('revision')


def archived(row):
    p, v = identity(row)
    if not PDB.fullmatch(p) or not REV.fullmatch(v or ''):
        return False
    low = p.lower()
    base = 'pdb_0000' + low
    expected = f'https://files-versioned.rcsb.org/pdb_versioned/data/entries/{low[1:3]}/{base}/{base}_xyz_v{v.replace(".","-")}.cif.gz'
    return row['file_type'] == 'pdbx_mmcif_gzip' and row['source_url'] == expected


def previous_entries(sources):
    completed = {}
    for r in sources:
        merged = r/'handoff/entry_results.csv'
        if merged.is_file():
            with merged.open(encoding='utf-8-sig', newline='') as stream:
                records = list(DictReader(stream))
        else:
            records = [read(p) for p in (r/'entry_results').glob('*.json')]
        planned = {identity(row) for row in read(r/'frozen.json')['files']}
        require(records and len({identity(row) for row in records}) == len(records)
                and all(identity(row) in planned for row in records), 'historical_entries_differ_from_frozen_input')
        for row in records:
            if row['status'] and row['status'] != 'TECHNICAL_FAILURE' and row.get('source_sha256'):
                completed[identity(row)] = r.name
    return completed


def build_plan(root, run, m1, prior, cap):
    plan = read(m1/'plan.json')
    summary = read(m1/'full_inventory/summary.json')
    plan_sha = sha(m1/'plan.json')
    scope='/'.join(plan['scope'][k] for k in ('start_inclusive','end_exclusive'))
    cross=read(m1/'step2_cross_window.json')
    require(cross['status']=='INVENTORY_COMPLETE' and cross['unique_ids']==plan['denominator']['total_count'] and
            summary['plan_sha256']==cross['plan_sha256']==plan_sha and summary['scope']==cross['scope']==scope and
            summary['denominator']==cross['denominator']==plan['denominator'], 'M1_rollup_not_linked_to_frozen_plan')
    windows = plan['windows']
    bounds = [(w['start_inclusive'], w['end_exclusive']) for w in windows]
    require(bounds and all(a < b for a,b in bounds) and all(b == c for (_,b),(c,_) in zip(bounds,bounds[1:])) and
            '/'.join((bounds[0][0],bounds[-1][1])) == scope, 'M1_window_boundaries_differ')
    source = m1/'full_inventory'; entries = source/'candidate_entries.jsonl'; file_plan = source/'monthly_file_plan.csv'
    files = {}
    with file_plan.open(encoding='utf8', newline='') as stream:
        for x in DictReader(stream):
            p = x['pdb_id']; require(p not in files,'M1_file_plan_duplicate')
            files[p] = ({k:x[k] for k in SPEC}, int(x['compressed_bytes_p90'] or 0))
    done = previous_entries(prior)
    groups = [[] for _ in windows]; seen = set(); reasons = []
    ended = [b for _,b in bounds]
    from bisect import bisect_right
    for x in jsonl(entries):
        p = x['pdb_id']; day = (x.get('first_release_date') or '')[:10]
        i = bisect_right(ended,day)
        require(p not in seen and i < len(bounds) and bounds[i][0] <= day < bounds[i][1], 'M1_entry_duplicate_or_out_of_window')
        seen.add(p)
        pin, high = files[p]
        require(pin['revision'] == (x.get('rcsb_latest_revision') or ''), 'M1_entry_and_archive_revision_disagree')
        route = x['plinder_revision_status'] if x['plinder_annotation_rows'] else 'NO_LIGAND_ANNOTATION'
        k = identity(pin)
        limited = ('HISTORICAL_SAME_REVISION' if k in done else
                   'OUT_OF_SCOPE' if SENSITIVE.search(x.get('title') or '') else
                   'UNRESOLVED_M1_ARCHIVE_IDENTITY' if not archived(pin) else None)
        record = {'source':pin, 'date':day, 'route':route, 'window':i, 'raw_p90_estimate':high}
        if limited:
            reasons.append({'pdb_id':p, 'revision':pin['revision'], 'window':i, 'status':limited,
                            **({'source_run':done[k]} if k in done else {})})
        groups[i].append((day,p,None if limited else record))
    require(seen == files.keys() == set(cross['entry_ids']) and len(seen) == plan['denominator']['total_count'] and
            len(windows) == summary['window_count'], 'M1_full_denominator_differs')
    path = run/'partitions'; path.mkdir(exist_ok=True)
    batch_plan = []; serial = 0
    for i, (spec,group) in enumerate(zip(windows,groups)):
        require(len(group) == spec['total_count'], 'M1_actual_window_count_differs')
        items=[]; raw_budget=0
        def emit():
            nonlocal serial, items, raw_budget
            if not items:
                return
            serial += 1
            name = f'b{serial:04d}'
            f = path/f'{name}.json'
            once(f,items)
            batch_plan.append({'id':name, 'window':i, 'first_release':items[0]['date'], 'last_release':items[-1]['date'],
                               'entries':len(items),'raw_p90_estimate':raw_budget,'sha256':sha(f)})
            items=[]; raw_budget=0
        for _,_,row in sorted(group):
            if row is None:
                continue
            if len(items) >= cap or raw_budget + row['raw_p90_estimate'] > 512 << 20:
                emit()
            items.append(row); raw_budget += row['raw_p90_estimate']
        emit()
    once(run/'excluded.json',sorted(reasons,key=lambda r:r['pdb_id']))
    require(len(reasons) + sum(s['entries'] for s in batch_plan) == len(seen), 'M1_to_production_count_differs')
    return {'schema':1, 'scope':scope, 'denominator':len(seen), 'window_count':len(windows),
            'm1_plan_sha256':plan_sha, 'm1_entries_sha256':sha(entries), 'm1_file_plan_sha256':sha(file_plan),
            'exclusions':dict(Counter(x['status'] for x in reasons)), 'excluded_sha256':sha(run/'excluded.json'),
            'partitions':batch_plan}


def cache_index(before):
    result={}
    for r in before:
        for p in (r/'source_proofs').glob('*.json'):
            proof = read(p); pin = proof.get('identity', {})
            if proof.get('status') != 'PASS' or proof.get('gzip_crc') != 'PASS' or not all(k in pin for k in SPEC):
                continue
            if proof.get('observed_revision') == pin['revision'] and proof.get('bytes',0) > 0 and proof.get('relative_path'):
                result[identity(pin)] = (r,proof)
    return result


def prepare(run, label, rows, cache, jobs, archive):
    folder = run/'preflight'; folder.mkdir(exist_ok=True)
    path = folder/f'{label}.json'
    if path.is_file():
        result=read(path)
        require(result['input_digest'] == digest(rows), 'rolling_preflight_input_changed')
        return result
    start = time.monotonic()
    def probe(record):
        pin=record['source']; cached=cache.get(identity(pin))
        if cached:
            origin,proof=cached
            if all(proof['identity'][k]==pin[k] for k in SPEC):
                try:
                    cached_file=file_ref(origin/proof['relative_path'],(origin,),proof['sha256'])
                    if cached_file['bytes']==proof['bytes']:
                        return pin | {'status':200,'content_length':proof['bytes'],'availability':'VERIFIED',
                                      'verified_cached_file':cached_file,'checked_utc':utc_now()}
                except (OSError,ValueError):
                    pass
        return archive.head(pin)
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        checks=list(pool.map(probe,rows))
    result={'input_digest':digest(rows),'checks':checks,'elapsed_s':round(time.monotonic()-start,3)}
    once(path,result)
    return result


def stage_batch(a, run, label, rows, seed, preflight):
    t=time.monotonic()
    runs=a.root/'runs'; m3=runs/f'{run.name}-{label}'; m2=runs/f'holo-m2-{run.name.removeprefix("holo-m3-")}-{label}'
    checks=preflight['checks']; linked={identity(r['source']):r['source'] for r in rows}
    require(len(linked)==len(checks) and set(linked)=={identity(c) for c in checks}, 'rolling_preflight_does_not_cover_batch')
    allowed=[c for c in checks if c['availability']=='VERIFIED']
    deferred=sorted(c['pdb_id'] for c in checks if c['availability']!='VERIFIED')
    if not allowed:
        return {'status':'ARCHIVE_PENDING','catalog_run':seed.name,'total':len(rows),'pending':deferred,'settled':0,'network_bytes':0}
    m2.mkdir(exist_ok=True); m3.mkdir(exist_ok=True)
    inputs=m2/'inputs'; inputs.mkdir(exist_ok=True)
    once(inputs/'blind.json', [linked[identity(c)] for c in allowed])
    once(inputs/'archive_head.json', allowed)
    current = m3/'increment_receipt.json'
    if not current.is_file():
        for c in allowed:
            if proof:=c.get('verified_cached_file'):
                dst=raw_path(m2,c); dst.parent.mkdir(exist_ok=True)
                if not dst.exists():
                    src=Path(proof['path']); require(src.is_relative_to(runs) and sha(src)==proof['sha256'], 'cached_raw_changed_after_preflight')
                    os.link(src,dst)
        pipeline=Path(__file__).resolve().with_name('holo_m3_increment.py')
        subprocess.run([sys.executable,'-B',str(pipeline),'--root',str(a.root),'--m2-run',str(m2),'--m3-run',str(m3),
                        '--batch',m3.name,'--seed-run',str(seed),'--index',str(a.index),'--reference-run',str(a.reference_run),
                        '--expected-entries',str(len(allowed)),'--jobs',str(a.jobs),'--fetch-jobs',str(a.fetch_jobs),
                        '--git-head',a.git_head,'--defer-pending','--archive-slots',str(a.archive.folder),
                        '--archive-limit',str(a.archive.limit),'--archive-ancillary',str(a.archive.ancillary),
                        '--primary-budget',str(a.primary_budget),
                         *(v for fd in a.owner_fds for v in ('--controller-lease-fd',str(fd)))],
                       check=True,pass_fds=a.owner_fds)
    summary=read(current)
    require(read(m3/'increment_plan.json')['seed_run']==str(seed), 'rolling_seed_branch_disagrees')
    return batch_outcome(runs,rows,m2,m3,summary,deferred,preflight,t)


def batch_outcome(runs, rows, m2, m3, summary, deferred, preflight, started):
    phase=Path(summary['catalog_run'])
    require((phase/'phase_receipt.json').is_file() and phase.parent==runs,
            'increment_did_not_seal_catalog')
    cases=Counter(); targets=Counter(); successful={}; phases=summary['phases']; network=0
    for phase_path in map(Path,phases):
        r=read(phase_path/'phase_receipt.json'); src=Path(r['m2_run'])
        targets.update(r['catalog']['outcomes']); cases.update(r['review_dispositions'])
        for p in r['selected']:
            proof=src/'source_proofs'/f'{p}.json'
            if proof.is_file():
                item=read(proof); network+=sum(x.get('body_bytes',0) for x in item.get('attempts',[]))
        for p in set(r['selected'])-set(r['m2_pending'])-set(r['m3_pending']):
            successful[p]=phase_path.name
    pending=sorted(set(summary['pending_technical']) | set(deferred))
    require(set(successful).isdisjoint(pending) and len(successful)+len(pending)==len(rows), 'rolling_batch_entry_accounting_differs')
    return {'status':summary['status'] if not deferred else 'ARCHIVE_PENDING', 'catalog_run':phase.name, 'total':len(rows),
            'pending':pending,'archive_deferred':deferred,'temporary_failures':summary['pending_technical'], 'settled':len(successful),
            'catalog_outcomes':dict(targets),'review_dispositions':dict(cases),'network_bytes':network,
            'cache_hits':sum(bool(r.get('verified_cached_file')) for r in preflight['checks'] if r['availability']=='VERIFIED'), 'elapsed_s':round(time.monotonic()-started,3),
            'pipeline_seconds':summary['elapsed_s'], 'preflight_seconds':preflight['elapsed_s'], 'm2_run':m2.name,'m3_run':m3.name}


def legacy_phase(root, legacy, label, rows, seed):
    """Rebuild a primary receipt from an intact old phase, without adopting its partial retry."""
    from aphcp_control.registry import connect
    runs=root/'runs'; parent=runs/f'{legacy.name}-{label}'; phase=parent.with_name(parent.name+'-phase-001')
    saved=phase/'phase_receipt.json'
    if not saved.is_file(): return None
    require(not (legacy/'receipts'/f'{label}.json').exists(),'legacy_batch_already_sealed')
    r=read(saved); require(Path(r['previous_catalog_run'])==seed and Path(r['registry_run'])==phase,
                          'legacy_phase_seed_differs')
    m2=Path(r['m2_run']); require(m2.parent==runs and m2.name.startswith('holo-m2-'),'legacy_phase_m2_outside_runs')
    preflight=read(legacy/'preflight'/f'{label}.json'); require(preflight['input_digest']==digest(rows),'legacy_preflight_differs')
    checks=preflight['checks']; require(len(checks)==len(rows) and {tuple(c[k] for k in SPEC) for c in checks}==
        {tuple(row['source'][k] for k in SPEC) for row in rows},'legacy_partition_and_preflight_differ')
    allowed={c['pdb_id'] for c in checks if c['availability']=='VERIFIED'}; require(allowed==set(r['selected']), 'legacy_selected_differs')
    frozen=m2/'frozen.json'; specs=read(frozen)['files']; pinned={x['source']['pdb_id']:x['source'] for x in rows}
    require(allowed=={x['pdb_id'] for x in specs} and all(all(x[k]==pinned[x['pdb_id']][k] for k in SPEC) for x in specs),
            'legacy_m2_source_differs')
    require(read(parent/'increment_plan.json')['m2_manifest_sha256']==sha(frozen),'legacy_m2_manifest_changed')
    catalog=read(phase/'catalog_receipt.json'); require(catalog==r['catalog'],'legacy_phase_catalog_differs')
    db=Path(catalog['registry']).resolve(strict=True)
    require(db.is_relative_to(root/'batches'/phase.name),'legacy_database_outside_batch')
    with connect(db,readonly=True) as con:
        require(con.execute('PRAGMA quick_check').fetchone()[0]=='ok' and not con.execute('PRAGMA foreign_key_check').fetchone(),
                'legacy_phase_sqlite_integrity')
        proof=con.execute('SELECT config_hash,error_summary FROM run WHERE id=? AND status=?',
                          (phase.name+'-observations','OBSERVATIONS_STAGED')).fetchone()
        require(proof and read(phase/'frozen_catalog.json') and digest(read(phase/'frozen_catalog.json'))==proof[0]
                and json.loads(proof[1])==catalog,'legacy_catalog_transaction_differs')
    pending=set(r['m2_pending'])|set(r['m3_pending']); require(pending<=allowed,'legacy_pending_outside_input')
    out=batch_outcome(runs,rows,m2,parent,{'catalog_run':str(phase),'phases':[str(phase)],
        'pending_technical':sorted(pending),'status':'STAGED_PENDING_RETRY' if pending else 'STAGED_COMPLETE',
        'elapsed_s':r['elapsed_s']},sorted({c['pdb_id'] for c in checks}-allowed),preflight,time.monotonic())
    return out | {'adopted_phase_sha256':sha(saved),'adopted_frozen_sha256':sha(frozen)}


def retry_download(run, row, archive, previous, known_head, max_attempts, consumed=0, prefer_ranges=False):
    """One bounded attempt writes only its own cache and source proof."""
    pin=row['source']; key=tuple(pin[k] for k in SPEC); item=run/'retry_cache'/digest(key)
    item.mkdir(parents=True,exist_ok=True)
    with (item/'claim.lock').open('a') as claim:
        fcntl.flock(claim,fcntl.LOCK_EX|fcntl.LOCK_NB)
        status=item/'attempt.json'; state=read(status) if status.is_file() else {'attempts':0}
        interrupted=state.get('state')=='STARTED'
        require(not state.get('pin') or state['pin']==pin,'retry_pin_changed')
        if state.get('state')=='READY' and state['attempts']>consumed: return state
        if state['attempts']>=max_attempts and not interrupted: return None
        # An interrupted attempt already consumed its slot; retry the checkpoint, not the counter.
        attempt=state['attempts'] if interrupted else state['attempts']+1
        started_epoch=time.time()
        proof_path=item/'source_proofs'/f"{pin['pdb_id']}.json"
        old=read(proof_path) if proof_path.is_file() else {}
        old_attempts=state['prior_proof_attempts'] if interrupted else len(old.get('attempts',[]))
        write(status,{'attempts':attempt,'next_epoch':started_epoch+15*2**(attempt-1),
                      'state':'STARTED','pin':pin,'prior_proof_attempts':old_attempts})
        check=state['check'] | {'retry_network_bytes':0,'retry_requests':[]} if state.get('state')=='READY' else None
        if check is None and previous:
            origin,proof=previous
            if all(proof['identity'].get(k)==pin[k] for k in SPEC):
                try:
                    verified=file_ref(origin/proof['relative_path'],(origin,),proof['sha256'])
                    if verified['bytes']==proof['bytes']:
                        check=pin | {'status':200,'content_length':proof['bytes'],'availability':'VERIFIED',
                                     'verified_cached_file':verified}
                except (OSError,ValueError):
                    pass
        if check is None:
            check=state.get('head') if interrupted else None
            if not check:
                check=known_head if known_head and known_head['availability']=='VERIFIED' else archive.head(pin)
            if check['availability']=='VERIFIED':
                write(status,read(status) | {'head':check})
            if check['availability']=='VERIFIED':
                expected=pin | {'expected_bytes':check['content_length']}
                if interrupted and len(old.get('attempts',[]))>old_attempts:
                    require(old['identity']==expected,'retry_saved_proof_identity_changed')
                    result='failed' if old['status']=='FAIL' else 'recovered_verified'
                    proof=old
                else:
                    _, result=fetch(item,expected,archive,'retry',prefer_ranges=prefer_ranges)
                    proof=read(proof_path)
                added=proof.get('attempts',[])[old_attempts:]
                check=check | {'retry_network_bytes':sum(a.get('body_bytes',0) for a in added),
                               'retry_requests':[x for a in added for x in a.get('requests',[])]}
                if result!='failed':
                    verified=file_ref(item/proof['relative_path'],(item,),proof['sha256'])
                    check=check | {'verified_cached_file':verified}
                else:
                    check=check | {'availability':'DEFERRED_BEFORE_FREEZE','error':proof.get('error')}
        result={'pin':pin,'attempts':attempt,'next_epoch':time.time()+15*2**(attempt-1),
                'state':'READY' if check['availability']=='VERIFIED' else 'FAILED','check':check}
        write(status,result)
        return result


def rolling(a):
    root=a.root.resolve(strict=True); a.root=root; runs=root/'runs'; run=a.run.resolve(strict=True)
    frozen_registry_root(root,run)
    for r in [a.seed_run,*a.prior]:
        require(r.resolve(strict=True).parent==runs, 'rolling_run_outside_data_root')
    code=Path(__file__).resolve().parent
    config={'seed_run':str(a.seed_run), 'historical_runs':[str(x) for x in a.prior], 'm1':str(a.m1),
            'index':str(a.index),'reference_run':str(a.reference_run), 'git_head':a.git_head,
            'code_sha256':sha(__file__), 'execution_files':{f.name:sha(f) for f in sorted(code.glob('holo*.py'))},
            'batch_size':a.batch_size}
    with ExitStack() as stack:
        lock=stack.enter_context((run/'rolling.lock').open('a'))
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        # An orphaned writer or downloader keeps the corresponding lease until it exits.
        a.owner_fds=(lock.fileno(),)
        legacy=Path(a.continue_from).resolve(strict=True) if getattr(a,'continue_from',None) else None
        if legacy:
            require(legacy.parent==runs and legacy!=run,'invalid_legacy_rolling_run')
            old_lock=stack.enter_context((legacy/'rolling.lock').open('r'))
            fcntl.flock(old_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            a.owner_fds+=(old_lock.fileno(),)
            old_config=read(legacy/'rolling_config.json')
            require(all(old_config[k]==config[k] for k in ('seed_run','m1','index','reference_run','batch_size')),
                    'legacy_rolling_scope_differs')
        config.update(archive_limit=a.fetch_jobs,archive_ancillary=a.archive_ancillary,primary_budget=a.primary_budget,
                      max_retries=a.max_retries,retry_jobs=a.retry_jobs,max_pending=a.max_pending,continued_from=str(legacy) if legacy else None)
        if legacy:
            config['legacy_plan_sha256']=sha(legacy/'rolling_plan.json')
        once(run/'rolling_config.json',config)
        a.archive=Archive(run/'archive_slots',a.fetch_jobs,a.archive_ancillary,lease_fds=a.owner_fds)
        path=run/'rolling_plan.json'
        if not path.is_file():
            once(path,read(legacy/'rolling_plan.json') if legacy else build_plan(root,run,a.m1,a.prior,a.batch_size))
        plan=read(path)
        print({k:v for k,v in plan.items() if k not in ('partitions',)} | {'batches':len(plan['partitions'])},flush=True)
        if a.plan_only:
            return
        proof=read(a.acceptance)
        require(proof['accepted'] and proof['catalog_run']==a.seed_run.name and proof['pending_technical']==[],
                'production_acceptance_gate_failed')
        cache=cache_index(tuple(r for r in runs.glob('holo-m2-*') if (r/'source_proofs').is_dir()))
        receipts=run/'receipts'; receipts.mkdir(exist_ok=True)
        seed=a.seed_run; done=0; primary=set(); pending={}; retries=Counter(); bytes_network=0; started=time.monotonic(); sequence=0
        by_id={p['id']:p for p in plan['partitions']}; consumed_attempts={}
        def progress(status, active=None):
            head={'catalog_run':seed.name}; pointer=run/'catalog_head.json'
            if not pointer.exists() or read(pointer)!=head: write(pointer,head)
            write(run/'rolling_progress.json',{'status':status,'updated_utc':utc_now(),'scope':plan['scope'], 'm1_entries':plan['denominator'],'active_partition':active,
                'exclusions':plan['exclusions'],'partitions_completed':len(primary),'partitions_total':len(plan['partitions']),
                'new_settled_entries':done,'attempted_pending_count':sum(map(len,pending.values())),
                'unattempted_entries':sum(x['entries'] for x in plan['partitions'] if x['id'] not in primary), 'catalog_run':seed.name,
                'unresolved_by_partition':{k:sorted(v) for k,v in pending.items() if v},'network_bytes':bytes_network,
                'invocation_elapsed_s':round(time.monotonic()-started,3),'free_disk_bytes':shutil.disk_usage(root).free})
        def rows_for(part):
            f=(legacy or run)/'partitions'/f"{part['id']}.json"
            require(sha(f)==part['sha256'],'rolling_partition_manifest_modified')
            return read(f)
        def load_receipt(label,saved):
            nonlocal seed,done,bytes_network,sequence
            bid,sep,attempt=label.partition('-r'); require(bid in by_id, 'rolling_receipt_unknown_partition')
            require(saved['sequence']==sequence+1 and saved['previous_catalog']==seed.name,'rolling_receipt_seed_chain_changed')
            if sep:
                chosen=set(saved.get('retry_selected',pending.get(bid,())))
                require(bid in primary and int(attempt)==retries[bid]+1 and len(chosen)==saved['total'] and
                        chosen<=pending[bid] and set(saved['pending'])<=chosen, 'rolling_retry_not_a_subset_of_pending_keys')
                checkpoints=saved.get('retry_attempts',{})
                require(not checkpoints or set(checkpoints)==chosen,'retry_attempt_checkpoint_differs')
                for pdb,n in checkpoints.items():
                    require(n>consumed_attempts.get((bid,pdb),0),'retry_attempt_already_consumed')
                    consumed_attempts[bid,pdb]=n
                retries[bid]+=1
            else:
                require(bid not in primary and saved['total']==by_id[bid]['entries'],'rolling_primary_receipt_conflict')
                primary.add(bid)
            remaining=set(saved['pending'])
            pending[bid]=((pending[bid]-chosen)|remaining) if sep else remaining
            require(saved.get('settled',0)+len(remaining)==saved['total'],'rolling_receipt_entry_count_differs')
            seed=runs/saved['catalog_run']; bytes_network+=saved.get('network_bytes',0)
            require((seed/'phase_receipt.json').is_file(), 'rolling_receipt_missing_catalog')
            done+=saved.get('settled',0); sequence+=1
        # Reconstruct the true commit order, including interleaved retries. A restart
        # never points clients back at an earlier catalog while replaying receipts.
        inherited=[(f.stem,read(f)) for f in (legacy/'receipts').glob('*.json')] if legacy else []
        for label,saved in inherited: once(receipts/f'{label}.json',saved)
        journal=[(f.stem,read(f)) for f in receipts.glob('*.json')]
        for label,saved in sorted(journal,key=lambda item:item[1]['sequence']):
            load_receipt(label,saved)
        dispatch=run/'active_dispatch.json'; request=run/'pause.request.json'
        if getattr(a,'resume_paused',False) and request.is_file():
            request.unlink(); sync_path(request.parent)

        def process(label, rows, prep=None, adopt=False):
            f=receipts/f'{label}.json'
            require(not f.is_file(), 'rolling_already_committed')
            if shutil.disk_usage(root).free < a.reserve_gib*GIB:
                progress('PAUSED_DISK_RESERVE')
                raise RuntimeError('rolling_disk_reserve_reached')
            current={'label':label,'pdb_ids':sorted(row['source']['pdb_id'] for row in rows),
                     'previous_catalog':seed.name,'adopt':adopt}
            (run/'preflight').mkdir(exist_ok=True)
            if dispatch.is_file() and read(dispatch):
                require(read(dispatch)==current,'unfinished_dispatch_must_resume_first')
            elif prep is not None and not hasattr(prep,'result'):
                once(run/'preflight'/f'{label}.json',prep)
            write(dispatch,current)
            progress('ROLLING',label)
            free_before=shutil.disk_usage(root).free
            try:
                if adopt:
                    saved=legacy_phase(root,legacy,label,rows,seed)
                    require(saved is not None,'legacy_adoption_phase_missing')
                else:
                    preflight=prep.result() if hasattr(prep,'result') else (prep or prepare(run,label,rows,cache,a.archive.ancillary,a.archive))
                    saved=stage_batch(a,run,label,rows,seed,preflight)
                saved |= {'previous_catalog':seed.name,'sequence':sequence+1,'committed_utc':utc_now(),
                          'filesystem_used_delta_bytes':free_before-shutil.disk_usage(root).free}
                if '-r' in label:
                    saved['retry_selected']=current['pdb_ids']
                    saved['retry_attempts']={x['pdb_id']:x['retry_attempt'] for x in preflight['checks']}
                    saved['network_bytes']+=sum(x.get('retry_network_bytes',0) for x in preflight['checks'])
            except (OSError,ValueError,RuntimeError,subprocess.SubprocessError):
                progress('PAUSED_UNSEALED_PARTITION',label)
                raise
            once(f,saved)
            load_receipt(label,saved)
            write(dispatch,None)
            print(label,saved,flush=True)
            progress('ROLLING')
            return saved

        if request.is_file():
            progress('PAUSED_AT_BOUNDARY'); return
        progress('RESUMING' if journal else 'ROLLING')
        if dispatch.is_file() and (active:=read(dispatch)):
            label=active['label']; bid=label.split('-r')[0]
            if (receipts/f'{label}.json').is_file():
                write(dispatch,None)  # Commit was durable before the process died.
            else:
                wanted=set(active['pdb_ids'])
                subset=[row for row in rows_for(by_id[bid]) if row['source']['pdb_id'] in wanted]
                require(len(subset)==len(wanted),'dispatch_source_plan_differs')
                process(label,subset,adopt=active['adopt'])

        (run/'preflight').mkdir(exist_ok=True)
        active_jobs={}; partition_cache={}; known_checks={}
        def retry_rows(bid):
            if bid not in partition_cache:
                partition_cache[bid]={row['source']['pdb_id']:row for row in rows_for(by_id[bid])}
            return partition_cache[bid]

        def head_for(bid,pdb):
            if bid not in known_checks:
                paths=[run/'preflight'/f'{bid}.json']+([legacy/'preflight'/f'{bid}.json'] if legacy else [])
                check_file=next((f for f in paths if f.is_file()),None)
                known_checks[bid]={c['pdb_id']:c for c in read(check_file)['checks']} if check_file else {}
            return known_checks[bid].get(pdb)

        def job_state(row):
            key=tuple(row['source'][k] for k in SPEC)
            path=run/'retry_cache'/digest(key)/'attempt.json'
            return read(path) if path.is_file() else {'attempts':0,'next_epoch':0}

        def schedule(pool):
            if request.is_file(): return
            ready=time.time()
            for bid,remaining in pending.items():
                for pdb in sorted(remaining):
                    if len(active_jobs)>=a.retry_jobs: return
                    if (bid,pdb) in active_jobs: continue
                    row=retry_rows(bid)[pdb]; state=job_state(row)
                    if unmerged(bid,pdb,state) or state['next_epoch']>ready: continue
                    if state['attempts']>=a.max_retries and state.get('state')!='STARTED': continue
                    pin=row['source']; primary_receipt=read(receipts/f'{bid}.json')
                    origin=primary_receipt.get('m2_run'); partial=False
                    if origin:
                        proof=runs/origin/'source_proofs'/f'{pdb}.json'
                        if proof.is_file():
                            partial=any(x.get('exit_code')==28 and x.get('http_status')==200 and x.get('body_bytes',0)>0
                                        for attempt in read(proof).get('attempts',[]) for x in attempt.get('requests',[]))
                            if not partial:  # Legacy transfer receipts predate per-request timing.
                                partial=any(x.get('http_status')==200 and x.get('body_bytes',0)>0 and 'Timeout' in x.get('error','')
                                            for x in read(proof).get('attempts',[]))
                    active_jobs[(bid,pdb)]=pool.submit(retry_download,run,row,a.archive,cache.get(identity(pin)),
                                                       head_for(bid,pdb),a.max_retries,consumed_attempts.get((bid,pdb),0),partial)

        def unmerged(bid,pdb,state):
            return 'check' in state and state['attempts']>consumed_attempts.get((bid,pdb),0)

        def collect():
            for key,future in list(active_jobs.items()):
                if future.done():
                    future.result(); del active_jobs[key]
            if request.is_file(): return
            in_flight={bid for bid,_ in active_jobs}
            for bid,remaining in list(pending.items()):
                # A bounded fetch wave writes one catalog increment, not one clone per file.
                # Main partitions still advance while this wave finishes.
                if bid in in_flight: continue
                rows=[]; checks=[]
                for pdb in sorted(remaining):
                    if (bid,pdb) in active_jobs: continue
                    row=retry_rows(bid)[pdb]; state=job_state(row)
                    if unmerged(bid,pdb,state):
                        rows.append(row); checks.append(state['check'] | {'retry_attempt':state['attempts']})
                if rows:
                    label=f'{bid}-r{retries[bid]+1:02d}'
                    process(label,rows,{'input_digest':digest(rows),'checks':checks,'elapsed_s':0})

        def next_retry():
            epochs=[]; exhausted=0
            for bid,remaining in pending.items():
                for pdb in remaining:
                    state=job_state(retry_rows(bid)[pdb]); ep=state['next_epoch']
                    if unmerged(bid,pdb,state): epochs.append(0)
                    elif state['attempts']<a.max_retries or state.get('state')=='STARTED': epochs.append(ep)
                    else: exhausted+=1
            return min(epochs,default=None),exhausted

        def retry_status():
            count=sum(map(len,pending.values())); _,exhausted=next_retry()
            oldest=min(((receipts/f'{bid}.json').stat().st_mtime for bid,keys in pending.items() if keys),default=None)
            return {'queued':count,'running':sum(not f.done() for f in active_jobs.values()),'exhausted':exhausted,
                    'oldest_wait_s':round(max(0,time.time()-oldest),1) if oldest is not None else None,
                    'settled_by_retry':sum(read(f).get('settled',0) for f in receipts.glob('b*-r*.json'))}

        unstarted=[part for part in plan['partitions'] if part['id'] not in primary]
        with ThreadPoolExecutor(max_workers=1) as ahead, ThreadPoolExecutor(max_workers=a.retry_jobs) as pool:
            future=None; consecutive_failures=0; completed_here=0; pause_reason=None
            for i,part in enumerate(unstarted):
                if request.is_file(): break
                collect()
                if request.is_file(): break
                if sum(map(len,pending.values()))>=a.max_pending:
                    pause_reason='PAUSED_RETRY_BACKLOG'; break
                schedule(pool)
                if i+1<len(unstarted) and (a.max_main_batches is None or completed_here+1<a.max_main_batches):
                    nxt=unstarted[i+1]
                    lookahead=ahead.submit(prepare,run,nxt['id'],rows_for(nxt),cache,a.archive.ancillary,a.archive)
                else:
                    lookahead=None
                label=part['id']; adopt=bool(legacy and i==0 and
                    (runs/f'{legacy.name}-{label}-phase-001'/'phase_receipt.json').is_file())
                result=process(label,rows_for(part),future,adopt=adopt); future=lookahead
                completed_here+=1
                if request.is_file(): break
                collect(); schedule(pool)
                write(run/'retry_progress.json',retry_status())
                consecutive_failures=consecutive_failures+1 if len(result['pending'])*4>result['total'] else 0
                if consecutive_failures>=3:
                    pause_reason='PAUSED_SOURCE_FAILURE'; break
                if a.max_main_batches is not None and completed_here>=a.max_main_batches:
                    pause_reason='PAUSED_OBSERVATION_LIMIT'; break
            if not pause_reason and not request.is_file():
                # Drain bounded retries at the end; don't keep retrying forever.
                while not request.is_file():
                    collect(); schedule(pool)
                    wake,_=next_retry()
                    if not active_jobs and wake is None: break
                    if active_jobs: wait(active_jobs.values(),timeout=.25,return_when=FIRST_COMPLETED)
                    else: time.sleep(min(.25,max(.01,wake-time.time())))
        # No asynchronous writer is left when publishing an actual paused state.
        if request.is_file(): pause_reason='PAUSED_AT_BOUNDARY'
        if pause_reason:
            if not request.is_file(): collect()
            write(run/'retry_progress.json',retry_status())
            progress(pause_reason); return
        write(run/'retry_progress.json',retry_status())
        exclusions=sum(plan['exclusions'].values())
        require(done+sum(map(len,pending.values()))+exclusions==plan['denominator'],'rolling_final_denominator_differs')
        progress('STAGED_PENDING_RETRY' if any(pending.values()) or plan['exclusions'].get('UNRESOLVED_M1_ARCHIVE_IDENTITY') else 'STAGED_COMPLETE')


def main():
    if len(sys.argv)==3 and sys.argv[1]=='--request-pause':
        run=Path(sys.argv[2]).resolve(strict=True)
        require((run/'rolling_config.json').is_file(),'not_a_rolling_run')
        write(run/'pause.request.json',{'requested_utc':utc_now()})
        return
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','run','m1','seed-run','index','reference-run','acceptance'):
        p.add_argument('--'+k,required=True,type=Path)
    p.add_argument('--prior',nargs='+',required=True,type=Path)
    p.add_argument('--git-head',required=True)
    p.add_argument('--batch-size',type=int,default=550)
    p.add_argument('--jobs',type=int,default=48)
    p.add_argument('--fetch-jobs',type=int,default=24)
    p.add_argument('--reserve-gib',type=int,default=1024)
    p.add_argument('--continue-from',type=Path,help='Adopt a stopped controller history from a distinct Run')
    p.add_argument('--primary-budget',type=float,default=30)
    p.add_argument('--archive-ancillary',type=int,default=6)
    p.add_argument('--retry-jobs',type=int,default=3)
    p.add_argument('--max-retries',type=int,default=3)
    p.add_argument('--max-pending',type=int,default=550)
    p.add_argument('--max-main-batches',type=int)
    p.add_argument('--resume-paused',action='store_true')
    p.add_argument('--plan-only',action='store_true')
    a=p.parse_args()
    if (min(a.jobs,a.fetch_jobs,a.batch_size,a.reserve_gib,a.primary_budget,a.max_pending,a.max_retries)<1 or
        a.max_retries>6 or not 1<=a.retry_jobs<=a.archive_ancillary<a.fetch_jobs or
        (a.max_main_batches is not None and a.max_main_batches<1)):
        p.error('invalid scheduling limits')
    rolling(a)


if __name__=='__main__':
    main()