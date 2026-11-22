"""Dispatch incomplete M1 windows and reconcile their frozen on-disk results."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from copy import deepcopy
from datetime import date, datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from m1_holo_config import (RELEASE_FIELD, date_filters, metadata_complete, read_json,
                           reusable_window, search_complete, sha256, window, write_json)


def population(request):
    """Compare selection semantics, ignoring only first-release bounds and presentation."""
    key = deepcopy(request)
    query = key['query']
    if query.get('type') != 'group' or query.get('logical_operator') != 'and':
        raise ValueError('Expected a group of ANDed M1 predicates')
    query['nodes'] = sorted((n for n in query['nodes']
        if n.get('parameters', {}).get('attribute') != RELEASE_FIELD),
        key=lambda n: json.dumps(n, sort_keys=True))
    key.pop('request_info', None)
    for field in ('paginate','return_counts','sort'):
        key.get('request_options', {}).pop(field, None)
    return key


def span(spec):
    start, end = (date.fromisoformat(spec[k]) for k in ('start_inclusive','end_exclusive'))
    if start >= end:
        raise ValueError(f'Invalid half-open interval: {spec}')
    return start, end


def label(start, end):
    return f'{start}/{end}'


def plan_windows(plan_path):
    plan = read_json(plan_path)
    lower, upper = span(plan['scope']); denominator = plan['denominator']
    sample = denominator['request']; population_key = population(sample)
    if sample.get('return_type') != 'entry' or sample['request_options'].get('return_counts') is not True:
        raise ValueError('Denominator must be a Search Entry count request')
    if type(denominator['total_count']) is not int or denominator['total_count'] < 0:
        raise ValueError('Denominator total_count must be a nonnegative integer')
    if not datetime.fromisoformat(denominator['queried_utc']).tzinfo:
        raise ValueError('Denominator needs a timestamp with timezone')
    if any(n.get('parameters', {}).get('attribute') == RELEASE_FIELD for n in sample['query']['nodes']):
        bounds = date_filters(sample)
        if (bounds['greater_or_equal']['value'],bounds['less']['value']) != (str(lower),str(upper)):
            raise ValueError('Denominator date bounds differ from the declared scope')

    windows = []
    planned_counts = {}
    for item in plan['windows']:
        config_path = (plan_path.parent/item['config']).resolve()
        run = (plan_path.parent/item['run_directory']).resolve()
        cfg = read_json(config_path); req = cfg['m1_1_request']
        if population(req) != population_key:
            raise ValueError(f'Count and window predicates differ: {config_path}')
        start, end, _ = window(cfg)
        frozen = read_json(run/'inputs'/'query.json')
        if frozen is not None and frozen != cfg:
            raise ValueError(f'Window differs from frozen input: {run}')
        windows.append((start,end,config_path,run,cfg))
        if 'total_count' in item:
            planned = item['total_count']
            if type(planned) is not int or planned < 0:
                raise ValueError(f'Invalid planned count: {config_path}')
            planned_counts[cfg['run_id']] = planned
    if len({w[3] for w in windows}) != len(windows) or len({w[4]['run_id'] for w in windows}) != len(windows):
        raise ValueError('Window run paths and run IDs must be unique')
    if planned_counts and len(planned_counts) == len(windows):
        if sum(planned_counts.values()) != denominator['total_count']:
            raise ValueError('Planned window counts do not reconcile with the denominator')
    return lower,upper,denominator,sorted(windows),planned_counts


def coverage_gaps(lower, upper, windows):
    gaps=[]; overlaps=[]; out=[]; ceiling=lower
    for start,end,*_ in windows:
        if start<lower or end>upper:
            out.append(label(start,end))
        a,b=max(start,lower),min(end,upper)
        if a>=b:
            continue
        if a>ceiling:
            gaps.append(label(ceiling,a))
        if a<ceiling:
            overlaps.append(label(a,min(ceiling,b)))
        ceiling=max(ceiling,b)
    if ceiling<upper:
        gaps.append(label(ceiling,upper))
    return gaps,overlaps,out


def collect(lower, upper, denominator, windows, recounts, planned_counts=None):
    gaps, overlaps, out = coverage_gaps(lower,upper,windows)
    planned_counts = planned_counts or {}
    planned_differences = []
    owners={}; duplicates={}; statuses=[]; missing=[]; unresolved=[]; off_window=[]; no_dates=[]; integrity=[]; drift=[]; unlocalized=[]
    for start,end,_,run,cfg in windows:
        id_=cfg['run_id']; search=read_json(run/'search'/'summary.json') or {}; data=read_json(run/'data'/'summary.json') or {}
        progress=read_json(run/'progress.json') or {}; receipt=read_json(run/'run_receipt.json') or {}
        good_search=search_complete(search); good=reusable_window(run)
        ids=search.get('entry_ids',[]) if good_search else []
        planned = planned_counts.get(id_)
        if good_search and planned is not None and len(ids) != planned:
            planned_differences.append({'window':id_, 'planned':planned, 'observed':len(ids)})
        technical=sorted(data.get('missing_ids',[])); period=sorted(progress.get('period_missing_ids',[]))
        if not good:
            (unresolved if search or progress or receipt else missing).append(id_)
        dates_path=run/'inventory'/'candidate_entries.jsonl'; seen=set()
        if ids and dates_path.is_file():
            with dates_path.open(encoding='utf-8') as rows:
                for line in rows:
                    row=json.loads(line); pdb=row['pdb_id']; day=row.get('first_release_date')
                    if pdb in seen:
                        integrity.append({'window':id_,'duplicate_candidate_id':pdb})
                    seen.add(pdb)
                    if day:
                        try:
                            day_key=date.fromisoformat(day[:10])
                        except (ValueError,TypeError):
                            day_key=None
                        if day_key is None or not start <= day_key < end:
                            off_window.append({'window':id_,'pdb_id':pdb,'first_release_date':day})
                    elif pdb not in data.get('confirmed_unavailable',{}) and pdb not in technical:
                        no_dates.append({'window':id_,'pdb_id':pdb})
            if seen != set(ids):
                integrity.append({'window':id_,'candidate_rows_missing_ids':sorted(set(ids)-seen),
                                  'candidate_rows_extra_ids':sorted(seen-set(ids))})
        if receipt:
            actual=receipt.get('window',{})
            if (actual.get('start_inclusive'),actual.get('end_exclusive'),receipt.get('batch_id')) != (str(start),str(end),cfg['batch_id']):
                integrity.append({'window':id_,'receipt_window_identity_mismatch':True})
        for pdb in ids:
            key=pdb.upper()
            if key in owners:
                duplicates.setdefault(key,[owners[key]]).append(id_)
            else:
                owners[key]=id_
        count=recounts.get(id_,{})
        if count and count.get('total_count') is not None and good_search and count['total_count'] != len(ids):
            drift.append({'window':id_,'observed_count':len(ids),'recount':count['total_count'],
                          'search_started_utc':search.get('started_utc'),'recount_utc':count['queried_utc']})
        if not count or count.get('total_count') is None:
            unlocalized.append(id_)
        statuses.append({'window':id_,'bounds':label(start,end),'run_directory':str(run),'complete':good,
                         'planned_count':planned,
                         'search_utc':search.get('started_utc'),'search_total_count':search.get('total_count'),
                         'search_refresh_offsets':search.get('refresh_offsets',[]), 'metadata_missing_ids':technical,
                         'period_missing_ids':period,'confirmed_unavailable':data.get('confirmed_unavailable',{}),
                         'missing_metadata_summary':bool(ids) and not data, 'receipt_status':receipt.get('status'),
                         'metadata_complete':not ids or metadata_complete(search,data)})
    delta=denominator['total_count']-len(owners)
    complete=not any((gaps,overlaps,out,missing,unresolved,off_window,no_dates,integrity,duplicates,delta,drift,planned_differences))
    return {'status':'INVENTORY_COMPLETE' if complete else 'PARTIAL','scope':label(lower,upper),
            'denominator':denominator,'unique_ids':len(owners),'entry_ids':sorted(owners),'denominator_minus_unique_ids':delta,
            'missing_date_ranges':gaps,'overlapping_date_ranges':overlaps,'windows_outside_scope':out,
            'missing_windows':missing,'unfinished_windows':unresolved,'duplicate_ids':duplicates,
            'off_window':off_window,'unverified_release_dates':no_dates,'integrity_issues':integrity,
            'planned_count_differences':planned_differences,'count_drift_windows':drift,'recount_evidence':recounts,
            'windows_without_recount_when_total_differs':unlocalized if delta else [],
            'denominator_note':('The count and window receipts are separate snapshots. Recount suspected windows; '
                                'if in-range counts agree, check planned date bounds and denominator timestamp.') if delta else None,
            'windows':statuses}


def recount(cfg):
    req=deepcopy(cfg['m1_1_request']); options=req['request_options']
    for name in ('paginate','sort'):
        options.pop(name,None)
    options['return_counts']=True
    url='https://search.rcsb.org/rcsbsearch/v2/query'; now=datetime.now(timezone.utc).isoformat()
    try:
        with urlopen(Request(url,json.dumps(req).encode(),{'Content-Type':'application/json'}),timeout=38) as response:
            total=0 if response.status==204 else json.load(response)['total_count']
            if type(total) is not int or total<0:
                raise ValueError('Invalid Search total_count')
            return {'request':req,'queried_utc':now,'total_count':total,'http_status':response.status}
    except (OSError,HTTPError,ValueError,KeyError) as error:
        return {'request':req,'queried_utc':now,'total_count':None,'error':str(error)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan',type=Path)
    parser.add_argument('--output',required=True,type=Path,help='Cross-window JSON receipt, separate from each Run')
    parser.add_argument('--inventory-dir',type=Path,help='Write the merged tables and budget after reconciliation')
    parser.add_argument('--jobs',type=int,default=1)
    parser.add_argument('--check-only',action='store_true')
    parser.add_argument('--register',action='store_true')
    parser.add_argument('--recount',action='append',default=[],metavar='RUN_ID',help='Explicitly request a fresh count for a suspect window')
    args=parser.parse_args()
    if args.jobs<1:
        parser.error('--jobs must be positive')
    lower,upper,denominator,windows,planned_counts=plan_windows(args.plan.resolve()); plan_hash=sha256(args.plan)
    known={w[4]['run_id']:w for w in windows}
    if set(args.recount)-known.keys():
        parser.error('--recount requires a selected run ID')
    if args.output.resolve() == args.plan.resolve() or any(args.output.resolve().is_relative_to(w[3]) for w in windows):
        parser.error('--output must be separate from the plan and outside individual window Run directories')
    if args.inventory_dir:
        inventory_dir = args.inventory_dir.resolve()
        if args.check_only or any(inventory_dir.is_relative_to(w[3]) for w in windows):
            parser.error('--inventory-dir requires a run and must be outside individual window directories')
        if args.output.resolve().is_relative_to(inventory_dir) or args.plan.resolve().is_relative_to(inventory_dir):
            parser.error('--inventory-dir must be separate from the plan and cross-window receipt')
    future_errors={}; launched=[]
    old=read_json(args.output) or {}
    counts=old.get('recount_evidence',{}) if old.get('plan_sha256')==plan_hash else {}
    if not args.check_only and not any(coverage_gaps(lower,upper,windows)):
        def dispatch(w):
            _,_,config,run,cfg=w
            cmd=[sys.executable,'-B',str(Path(__file__).with_name('m1_holo_inventory.py')),str(run),'--config',str(config)]
            if args.register: cmd.append('--register')
            return subprocess.run(cmd).returncode
        todo=[w for w in windows if not reusable_window(w[3]) or (args.register and
            (read_json(w[3]/'search'/'summary.json') or {}).get('unique_ids') and
            not (read_json(w[3]/'run_receipt.json') or {}).get('registry'))]
        launched=[w[4]['run_id'] for w in todo]
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            tasks={pool.submit(dispatch,w):w[4]['run_id'] for w in todo}
            for task in as_completed(tasks):
                try:
                    result=task.result()
                    if result: future_errors[tasks[task]]=result
                except OSError as error:
                    future_errors[tasks[task]]=str(error)
    for key in args.recount:
        counts[key]=recount(known[key][4])
    summary=collect(lower,upper,denominator,windows,counts,planned_counts)
    summary.update(plan_sha256=plan_hash, dispatched=launched, child_errors=future_errors,
                   checked_utc=datetime.now(timezone.utc).isoformat())
    if future_errors:
        summary['status']='PARTIAL'
    if args.inventory_dir and summary['status'] == 'INVENTORY_COMPLETE':
        from m1_holo_rollup import write_inventory
        try:
            summary['inventory'] = write_inventory(inventory_dir, summary)
        except (OSError, ValueError, KeyError) as error:
            summary['inventory_error'] = str(error)
            summary['status'] = 'PARTIAL'
    write_json(args.output,summary)
    print(f"{summary['status']}: {summary['unique_ids']} / {denominator['total_count']} unique Entry, {len(launched)} windows dispatched; {args.output}")
    return 0 if summary['status']=='INVENTORY_COMPLETE' else 2


if __name__=='__main__':
    raise SystemExit(main())
