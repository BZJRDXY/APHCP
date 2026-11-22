"""Plan and preflight one continuous M3 date window from the accepted M1 snapshot."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from csv import DictReader
from datetime import date, timedelta
from pathlib import Path
import json
import re
import time

from holo_candidate_batch import atomic_json, now
from holo_m2_final_selection import SENSITIVE, head


def preflight(root, run, start, target, jobs, historical_ids):
    t = time.monotonic()
    inventory = root/'batches/holo-m1-full-20261121a/planning/full_inventory'
    exclusions = set()
    for p in (root/'runs').glob('holo-m2-*/source_proofs/*.json'):
        with p.open(encoding='utf8') as stream:
            proof = json.load(stream)
        if proof.get('status') == 'PASS':
            ident = proof['identity']; exclusions.add((ident['pdb_id'],ident['revision']))
    candidates = []
    end_limit = start + timedelta(days=60)
    with (inventory/'candidate_entries.jsonl').open(encoding='utf8') as stream:
        for line in stream:
            row = json.loads(line)
            day = (row.get('first_release_date') or '')[:10]
            if start.isoformat() <= day < end_limit.isoformat():
                candidates.append(row)
    def exclude(row):
        p, rev = row['pdb_id'], row.get('rcsb_latest_revision')
        if p in historical_ids or (p,rev) in exclusions:
            return 'ALREADY_PROCESSED'
        if SENSITIVE.search(row.get('title') or ''):
            return 'OUT_OF_SCOPE'
        if not re.fullmatch(r'\d+\.\d+',rev or '') or not re.fullmatch(r'[0-9][A-Z0-9]{3}',p):
            return 'UNRESOLVED_M1_ARCHIVE_IDENTITY'
        return None
    days = Counter(r['first_release_date'][:10] for r in candidates if not exclude(r))
    options = [(abs(sum(n for d,n in days.items() if d<end)-target),end) for end in
               [(date.fromisoformat(d)+timedelta(days=1)).isoformat() for d in days]
               if 300 <= sum(n for d,n in days.items() if d<end) <= 500]
    if not options:
        raise ValueError('No whole-day window with 300–500 eligible Entries')
    end = min(options)[1]
    window = [r for r in candidates if r['first_release_date'][:10]<end]
    wanted = {r['pdb_id']:r for r in window if not exclude(r)}
    pins = {}
    with (inventory/'monthly_file_plan.csv').open(encoding='utf8',newline='') as stream:
        for r in DictReader(stream):
            if r['pdb_id'] in wanted:
                p = r['pdb_id']
                if p in pins or r['revision'] != wanted[p]['rcsb_latest_revision']:
                    raise ValueError(f'M1 file plan identity conflict: {p}')
                pins[p] = {k:r[k] for k in ('pdb_id','revision','source_url','file_type')}
    if set(pins) != set(wanted):
        raise ValueError('Not all eligible Entries have an exact M1 file plan')
    inputs = run/'inputs'; inputs.mkdir(parents=True,exist_ok=True)
    if (run/'frozen.json').exists():
        raise ValueError('Preflight cannot replace an immutable M2 freeze')
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        proofs = list(pool.map(head, sorted(pins.values(),key=lambda x:x['pdb_id'])))
    available = {x['pdb_id'] for x in proofs if x['availability']=='VERIFIED'}
    deferred = sorted(set(pins)-available)
    atomic_json(inputs/'archive_head_all.json',proofs)
    atomic_json(inputs/'window_scope.json', [{'pdb_id':x['pdb_id'],'first_release_date':x['first_release_date'],
        'revision':x['rcsb_latest_revision'],'outcome':exclude(x) or ('PRE_FREEZE_ARCHIVE_UNAVAILABLE' if x['pdb_id'] in deferred else 'FROZEN')}
        for x in window])
    routes = Counter(('NO_LIGAND_ANNOTATION' if not wanted[p]['plinder_annotation_rows'] else
                      wanted[p]['plinder_revision_status']) for p in available)
    receipt = {'window':[start.isoformat(),end],'window_entries':len(window), 'exclusions':dict(Counter(
        x for r in window if (x:=exclude(r)))), 'requested_archive':len(pins),'archive_deferred':len(deferred),
        'frozen_entries':len(available),'planned_archive_bytes':sum(x['content_length'] for x in proofs if x['pdb_id'] in available),
        'm1_routes':dict(routes),'completed_utc':now(),'elapsed_s':round(time.monotonic()-t,3)}
    atomic_json(inputs/'preflight_receipt.json',receipt)
    if not 300<=len(available)<=500:
        raise ValueError('Archive preflight fell outside the authorized freeze size')
    atomic_json(inputs/'blind.json',[pins[p] for p in sorted(available)])
    atomic_json(inputs/'archive_head.json',[x for x in proofs if x['pdb_id'] in available])
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True,type=Path)
    parser.add_argument('--run',required=True,type=Path)
    parser.add_argument('--start',required=True,type=date.fromisoformat)
    parser.add_argument('--target',default=380,type=int)
    parser.add_argument('--jobs',default=16,type=int)
    parser.add_argument('--historical-ids',type=Path)
    args=parser.parse_args()
    old=set(args.historical_ids.read_text(encoding='utf-8-sig').split()) if args.historical_ids else set()
    preflight(args.root,args.run,args.start,args.target,args.jobs,old)
