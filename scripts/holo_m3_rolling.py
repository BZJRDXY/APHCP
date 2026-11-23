"""Continue M3 staging over the accepted M1 windows using the existing increment runner."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from csv import DictReader
import datetime as dt
import fcntl
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from holo_m2_blind import raw_path
from holo_m2_final_selection import SENSITIVE, head
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


def prepare(run, label, rows, cache, jobs):
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
        return head(pin)
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
        return {'status':'ARCHIVE_PENDING','catalog_run':seed.name,'total':len(rows),'pending':deferred,'published':0}
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
                        '--git-head',a.git_head,'--passes','2'],check=True)
    summary=read(current)
    require(read(m3/'increment_plan.json')['seed_run']==str(seed), 'rolling_seed_branch_disagrees')
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
            'cache_hits':sum(bool(r.get('verified_cached_file')) for r in allowed), 'elapsed_s':round(time.monotonic()-t,3),
            'pipeline_seconds':summary['elapsed_s'], 'preflight_seconds':preflight['elapsed_s'], 'm2_run':m2.name,'m3_run':m3.name}


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
    once(run/'rolling_config.json',config)
    with (run/'rolling.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        path=run/'rolling_plan.json'
        if not path.is_file():
            once(path,build_plan(root,run,a.m1,a.prior,a.batch_size))
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
        by_id={p['id']:p for p in plan['partitions']}
        def progress(status, active=None):
            write(run/'catalog_head.json',{'catalog_run':seed.name})
            write(run/'rolling_progress.json',{'status':status,'updated_utc':utc_now(),'scope':plan['scope'], 'm1_entries':plan['denominator'],'active_partition':active,
                'exclusions':plan['exclusions'],'partitions_completed':len(primary),'partitions_total':len(plan['partitions']),
                'new_settled_entries':done,'attempted_pending_count':sum(map(len,pending.values())),
                'unattempted_entries':sum(x['entries'] for x in plan['partitions'] if x['id'] not in primary), 'catalog_run':seed.name,
                'unresolved_by_partition':{k:sorted(v) for k,v in pending.items() if v},'network_bytes':bytes_network,
                'invocation_elapsed_s':round(time.monotonic()-started,3),'free_disk_bytes':shutil.disk_usage(root).free})
        def rows_for(part):
            f=run/'partitions'/f"{part['id']}.json"
            require(sha(f)==part['sha256'],'rolling_partition_manifest_modified')
            return read(f)
        def load_receipt(label,saved):
            nonlocal seed,done,bytes_network,sequence
            bid,sep,attempt=label.partition('-r'); require(bid in by_id, 'rolling_receipt_unknown_partition')
            require(saved['sequence']==sequence+1 and saved['previous_catalog']==seed.name,'rolling_receipt_seed_chain_changed')
            if sep:
                require(bid in primary and int(attempt)==retries[bid]+1 and len(pending[bid])==saved['total'] and
                        set(saved['pending'])<=pending[bid], 'rolling_retry_not_a_subset_of_pending_keys')
                retries[bid]+=1
            else:
                require(bid not in primary and saved['total']==by_id[bid]['entries'],'rolling_primary_receipt_conflict')
                primary.add(bid)
            pending[bid]=set(saved['pending'])
            require(saved.get('settled',0)+len(pending[bid])==saved['total'],'rolling_receipt_entry_count_differs')
            seed=runs/saved['catalog_run']; bytes_network+=saved.get('network_bytes',0)
            require((seed/'phase_receipt.json').is_file(), 'rolling_receipt_missing_catalog')
            done+=saved.get('settled',0); sequence+=1
        # Reconstruct the true commit order, including interleaved retries. A restart
        # never points clients back at an earlier catalog while replaying receipts.
        journal=[(f.stem,read(f)) for f in receipts.glob('*.json')]
        for label,saved in sorted(journal,key=lambda item:item[1]['sequence']):
            load_receipt(label,saved)
        progress('ROLLING')
        def process(label, rows, prep=None):
            f=receipts/f'{label}.json'
            require(not f.is_file(), 'rolling_already_committed')
            if shutil.disk_usage(root).free < a.reserve_gib*GIB:
                progress('PAUSED_DISK_RESERVE')
                raise RuntimeError('rolling_disk_reserve_reached')
            progress('ROLLING',label)
            free_before=shutil.disk_usage(root).free
            try:
                preflight=prep.result() if prep else prepare(run,label,rows,cache,a.fetch_jobs)
                saved=stage_batch(a,run,label,rows,seed,preflight) | {'previous_catalog':seed.name,'sequence':sequence+1,
                                                                   'filesystem_used_delta_bytes':free_before-shutil.disk_usage(root).free}
            except (OSError,ValueError,RuntimeError,subprocess.SubprocessError):
                progress('PAUSED_UNSEALED_PARTITION',label)
                raise
            once(f,saved)
            load_receipt(label,saved)
            print(label,saved,flush=True)
            progress('ROLLING')
            return saved
        unstarted=[part for part in plan['partitions'] if part['id'] not in primary]
        with ThreadPoolExecutor(max_workers=1) as ahead:
            future=None; consecutive_failures=0
            for i,part in enumerate(unstarted):
                if i+1<len(unstarted):
                    nxt=unstarted[i+1]
                    lookahead=ahead.submit(prepare,run,nxt['id'],rows_for(nxt),cache,a.fetch_jobs)
                else:
                    lookahead=None
                result=process(part['id'],rows_for(part),future); future=lookahead
                consecutive_failures=consecutive_failures+1 if len(result['pending'])*4>result['total'] else 0
                if consecutive_failures>=3:
                    progress('PAUSED_SOURCE_FAILURE'); return
        # Rebase delayed retries on the latest sealed catalog. Do not rewind the head.
        for part in plan['partitions']:
            bid=part['id']
            if not pending.get(bid) or retries[bid]>=99:
                continue
            rows=[r for r in rows_for(part) if r['source']['pdb_id'] in pending[bid]]
            # Reuse verified raw from a partially successful previous pass when present.
            for label in (bid,f"{bid}-r{retries[bid]:02d}"):
                record_path=receipts/f'{label}.json'
                if record_path.is_file():
                    m2_name=read(record_path).get('m2_run')
                    if m2_name:
                        cache.update(cache_index((runs/m2_name,)))
            process(f'{bid}-r{retries[bid]+1:02d}',rows)
        exclusions=sum(plan['exclusions'].values())
        require(done+sum(map(len,pending.values()))+exclusions==plan['denominator'],'rolling_final_denominator_differs')
        progress('STAGED_PENDING_RETRY' if any(pending.values()) or plan['exclusions'].get('UNRESOLVED_M1_ARCHIVE_IDENTITY') else 'STAGED_COMPLETE')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','run','m1','seed-run','index','reference-run','acceptance'):
        p.add_argument('--'+k,required=True,type=Path)
    p.add_argument('--prior',nargs='+',required=True,type=Path)
    p.add_argument('--git-head',required=True)
    p.add_argument('--batch-size',type=int,default=550)
    p.add_argument('--jobs',type=int,default=48)
    p.add_argument('--fetch-jobs',type=int,default=24)
    p.add_argument('--reserve-gib',type=int,default=1024)
    p.add_argument('--plan-only',action='store_true')
    a=p.parse_args()
    if min(a.jobs,a.fetch_jobs,a.batch_size,a.reserve_gib)<1:
        p.error('positive limits are required')
    rolling(a)


if __name__=='__main__':
    main()