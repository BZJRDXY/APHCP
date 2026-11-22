"""Freeze a count-only M1 plan for m1_holo_windows.py; fetch just two boundary IDs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
from datetime import date, datetime, timedelta, timezone
import json
from pathlib import Path, PurePosixPath
import re
import time
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from m1_holo_config import RELEASE_FIELD, date_filters, metadata_complete, read_json, search_complete, sha256, window, write_json
from m1_holo_windows import coverage_gaps, population

SEARCH = 'https://search.rcsb.org/rcsbsearch/v2/query'
DATA = 'https://data.rcsb.org/rest/v1/core/entry/'
PAGE_LIMIT = 10_000  # RCSB Search API documentation, Pagination section.


def now():
    return datetime.now(timezone.utc).isoformat()


def total(receipt):
    if receipt.get('http_status') == 204 and receipt.get('response') is None:
        return 0
    obj = receipt.get('response')
    if receipt.get('http_status') == 200 and isinstance(obj, dict):
        count = obj.get('total_count')
        if type(count) is int and count >= 0 and obj.get('result_type') == 'entry' and not obj.get('errors'):
            return count
    return None


def count_body(query, start=None, end=None):
    q = population(query)
    if start is not None:
        nodes = deepcopy(query['query']['nodes'])
        bounded = {'query': {'nodes': nodes}}
        filters = date_filters(bounded)
        filters['greater_or_equal']['value'], filters['less']['value'] = str(start), str(end)
        q['query']['nodes'] = nodes
    else:  # Keep the same predicate order as the template.
        q['query']['nodes'] = [deepcopy(n) for n in query['query']['nodes'] if n['parameters']['attribute'] != RELEASE_FIELD]
    q['request_options']['return_counts'] = True
    return q


def search(path, body, valid):
    prior = read_json(path) or {}
    if prior.get('url') == SEARCH and prior.get('request') == body and valid(prior):
        return prior
    attempts = prior.get('attempts', []) if prior.get('request') == body else []
    for wait in (2, 4, 8, None):
        start = now(); clock = time.monotonic(); status = None; headers = {}; response = None
        try:
            request = Request(SEARCH, json.dumps(body).encode(), {'Content-Type': 'application/json', 'Accept': 'application/json'})
            with urlopen(request, timeout=38) as r:
                status, headers, raw = r.status, dict(r.headers), r.read()
                response = json.loads(raw) if raw else None
        except HTTPError as error:
            status, headers, response = error.code, dict(error.headers), error.read(4000).decode(errors='replace')
        except (OSError, ValueError) as error:
            response = str(error)
        attempt = dict(queried_utc=start, responded_utc=now(), seconds=round(time.monotonic()-clock, 3),
                       http_status=status, response_headers=headers, response=response)
        attempts.append(attempt)
        receipt = dict(url=SEARCH, method='POST', request=body, **attempt, attempts=attempts)
        write_json(path, receipt)
        if valid(receipt):
            return receipt
        if wait is None or status in (400, 401, 403, 404, 422):
            break
        delay = next((v for k, v in headers.items() if k.lower() == 'retry-after'), '')
        time.sleep(max(wait, int(delay)) if delay.isdecimal() else wait)
    raise RuntimeError(f'Unresolved Search request: {path}')


def boundary(out, query, label):
    q = count_body(query)
    q['request_options'].pop('return_counts')
    q['request_options']['paginate'] = {'start': 0, 'rows': 1}
    q['request_options']['sort'] = deepcopy(query['request_options']['sort'])
    q['request_options']['sort'][0]['direction'] = 'asc' if label == 'earliest' else 'desc'
    def valid(x):
        obj = x.get('response')
        rows = obj.get('result_set', []) if isinstance(obj, dict) else []
        return total(x) and len(rows) == 1 and re.fullmatch(r'[0-9A-Z]+', rows[0].get('identifier', ''))
    sr = search(out / f'{label}.json', q, valid)
    pdb = sr['response']['result_set'][0]['identifier']
    path = out / f'boundary-{pdb}.json'; url = DATA + pdb
    rec = read_json(path) or {}
    if rec.get('url') != url or rec.get('http_status') != 200:
        started = now()
        with urlopen(url, timeout=38) as r:
            status, data = r.status, json.load(r)
        rec = dict(url=url, method='GET', queried_utc=started, http_status=status,
                   response_subset={k:data.get(k) for k in ('rcsb_id', 'rcsb_accession_info', 'exptl')})
        write_json(path, rec)
    result = rec['response_subset']
    if result['rcsb_id'] != pdb:
        raise ValueError(f'Data API identity mismatch: {path}')
    day = date.fromisoformat(result['rcsb_accession_info']['initial_release_date'][:10])
    return day, {'pdb_id':pdb, 'date':str(day), 'search_total':total(sr), 'search_receipt':f'{label}.json', 'data_receipt':path.name}


def partition(start, end, months):
    boundaries = [start] + [date(y, m, 1) for y in range(start.year, end.year + 1) for m in months
                            if start < date(y, m, 1) < end] + [end]
    return list(zip(boundaries, boundaries[1:]))


def merge_sparse(atoms, locked, limit):
    """Coalesce consecutive low-volume segments; never cross either edge of a locked window."""
    groups = []
    for a,b,n in atoms:
        if (a,b) not in locked and groups:
            prev_a,prev_b,prev_n = groups[-1]
            if (prev_a,prev_b) not in locked and prev_b == a and prev_n + n <= limit:
                groups[-1] = prev_a,b,prev_n+n
                continue
        groups.append((a,b,n))
    return groups


def build(args):
    started = now(); clock = time.monotonic()
    out = args.output; out.mkdir(parents=True, exist_ok=True)
    template = read_json(args.template); query = template['m1_1_request']
    # A resumed plan keeps the original observation time; start in a new directory for a fresh baseline.
    def take(path, body):
        return search(out / path, body, lambda rec: total(rec) is not None)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        f = pool.submit(take, 'full.json', count_body(query))
        l = pool.submit(boundary, out, query, 'earliest'); h = pool.submit(boundary, out, query, 'latest')
        full = f.result(); lower, first = l.result(); newest, last = h.result()
        upper = datetime.fromisoformat(full['queried_utc']).date() + timedelta(days=1)
        if not lower <= newest < upper:
            raise ValueError(f'Unexpected release dates: earliest {lower}, latest {newest}, exclusive end {upper}')
        scope_future = pool.submit(take, 'in-scope.json', count_body(query, lower, upper))
        missing_filter = count_body(query)
        missing_filter['query']['nodes'].append({'type': 'terminal', 'service': 'text', 'parameters': {
            'attribute': RELEASE_FIELD, 'operator': 'exists', 'negation': True}})
        missing_future = pool.submit(take, 'without-release-date.json', missing_filter)
        reuse = {}
        for config, run in args.reuse:
            path = Path(config); cfg = read_json(path); a,b,_ = window(cfg); key = a,b
            if (key in reuse or not lower <= a < b <= upper or population(cfg['m1_1_request']) != population(query)
                    or PurePosixPath(run).name != cfg['run_id'] or not run.startswith('/')):
                raise ValueError(f'Invalid reuse window: {config}, {run}')
            reuse[key] = path, run
        base = sorted({d for pair in partition(lower, upper, (1,)) for d in pair}
                      | {d for pair in reuse for d in pair})
        counts = {}
        def count(pair):
            a,b = pair; path = f'counts/{a}__{b}.json'
            rec = take(path, count_body(query, a, b))
            return total(rec), path, rec['queried_utc']
        def get(pairs):
            todo = sorted(set(pairs) - counts.keys())
            counts.update(zip(todo, pool.map(count, todo)))
        remaining = list(zip(base, base[1:])); atoms = []
        while remaining:
            get(remaining)
            split = []
            for a,b in remaining:
                n = counts[a,b][0]
                if n <= args.target or (a,b) in reuse or (b-a).days == 1:
                    atoms.append((a,b,n))
                    continue
                pieces = partition(a,b,(1,4,7,10))
                if len(pieces) == 1:
                    pieces = partition(a,b,tuple(range(1,13)))
                if len(pieces) == 1:
                    midpoint = a + (b-a)//2
                    pieces = [(a,midpoint),(midpoint,b)]
                split.extend(pieces)
            remaining = split
        groups = merge_sparse(sorted(atoms), set(reuse), args.target)
        get((a,b) for a,b,_ in groups)
        try:
            missing = total(missing_future.result())
        except RuntimeError:
            missing = None  # Test recorded; bounded vs. unbounded counts still quantify any gap.
        in_scope = total(scope_future.result())
        tail = total(take('full-after.json', count_body(query)))

    planned=[]; reuse_rows=[]; drift=[]; cap_issues=[]
    for a,b,expected in groups:
        n, receipt, observed = counts[a,b]; key = a,b
        if n != expected:
            drift.append({'start':str(a),'end':str(b),'component_sum':expected,'direct_count':n})
        if n > PAGE_LIMIT:
            cap_issues.append({'start':str(a),'end':str(b),'count':n})
        cfg = deepcopy(template)
        old = reuse.get(key); run = None; status = 'NEW_WINDOW'
        config = f'windows/{a:%Y%m%d}-{b:%Y%m%d}/query.json'
        if old is not None:
            old_cfg, old_run = old
            ss = read_json(old_cfg.with_name('search_summary.json')) or {}
            ds = read_json(old_cfg.with_name('data_summary.json')) or {}
            rr_path = old_cfg.with_name('run_receipt.json')
            rr = read_json(rr_path) or {}
            reusable = search_complete(ss) and metadata_complete(ss,ds) and rr.get('status') == 'INVENTORY_COMPLETE'
            status = 'REUSE_ACCEPTED' if reusable and ss['unique_ids'] == n else 'PRIOR_WINDOW_REVIEW'
            prior = {'status':status,'old_run':old_run,'old_unique_ids':ss.get('unique_ids'),
                     'old_query_sha256':sha256(old_cfg),'old_receipt_sha256':sha256(rr_path) if rr else None,
                     'old_missing_ids':ds.get('missing_ids'),'new_count':n}
            reuse_rows.append(prior)
            if status == 'REUSE_ACCEPTED':
                cfg, run = read_json(old_cfg), old_run
        if run is None:
            filt = date_filters(cfg['m1_1_request'])
            filt['greater_or_equal']['value'], filt['less']['value'] = str(a), str(b)
            batch = f'holo-m1-{a:%Y%m%d}-{b:%Y%m%d}'
            cfg.update(batch_id=batch,run_id=f'{batch}-{args.tag}', config_version=f'm1-full-{args.tag}-v1')
            run = str(PurePosixPath(args.run_root) / cfg['run_id'])
        write_json(out / config, cfg)
        row = dict(config=config,run_directory=run,start_inclusive=str(a),end_exclusive=str(b),total_count=n,
                   queried_utc=observed,count_receipt=receipt,reuse_status=status)
        if old is not None:
            row['prior_window'] = prior
        planned.append(row)
    gaps, overlaps, outside = coverage_gaps(lower,upper,[(a,b) for a,b,_ in groups])
    n_full = total(full); n_windows = sum(w['total_count'] for w in planned)
    cross_checks = dict(full_count=n_full,in_scope_count=in_scope,full_count_after=tail,without_release_date_count=missing,
                        full_minus_in_scope=n_full-in_scope,window_count_sum=n_windows,full_minus_windows=n_full-n_windows,
                        in_scope_minus_windows=in_scope-n_windows,window_count=len(planned),missing_date_ranges=gaps,
                        overlapping_date_ranges=overlaps,outside_scope=outside,pagination_cap_issues=cap_issues,
                        component_count_drift=drift,historical_reuse=reuse_rows)
    pass_ = (not any((gaps, overlaps, outside, cap_issues, drift)) and n_full == n_windows == in_scope == tail
             and first['search_total'] == last['search_total'] == n_full and missing in (0,None))
    scope = dict(start_inclusive=str(lower),end_exclusive=str(upper),field=RELEASE_FIELD,
                 end_rule='UTC calendar day after initial denominator request')
    denominator = dict(total_count=n_full, queried_utc=full['queried_utc'], request=full['request'],
                       response=full['response'], http_status=full['http_status'], receipt='full.json')
    planning = dict(status='READY_FOR_M1_STEP_2' if pass_ else 'REVIEW', started_utc=started, finished_utc=now(),
                    seconds=round(time.monotonic()-clock,3), count_worker_limit=args.jobs, range_count_requests=len(counts),
                    target_entries_per_window=args.target, documented_pagination_max=PAGE_LIMIT,
                    boundary_entries={'earliest':first,'latest':last}, denominator_concordance=cross_checks)
    write_json(out/'plan.json',dict(batch_id=f'holo-m1-full-{args.tag}',scope=scope,denominator=denominator,windows=planned,planning=planning))
    print(json.dumps(planning,ensure_ascii=False,indent=2))
    return 0 if pass_ else 2


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--template', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--run-root', required=True)
    p.add_argument('--tag', required=True)
    p.add_argument('--reuse', nargs=2, action='append', default=[], metavar=('QUERY_JSON','RUN_DIRECTORY'))
    p.add_argument('--jobs', type=int, default=8)
    p.add_argument('--target', type=int, default=8000)
    args = p.parse_args()
    if not args.run_root.startswith('/') or not re.fullmatch(r'[0-9a-z]+',args.tag) or not 3 <= args.jobs or not 1 <= args.target <= PAGE_LIMIT:
        p.error('Use an absolute Pro5000 run root, a lowercase alphanumeric tag, at least three workers, and target <= 10000')
    return build(args)


if __name__ == '__main__':
    raise SystemExit(main())
