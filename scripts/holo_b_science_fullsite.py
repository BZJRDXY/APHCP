"""Read-only pocket-residue comparison for same-revision Holo B matches."""
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import argparse
import json
from pathlib import Path
import sqlite3
import time

from holo_b_science_sitecheck import check


def utc():
    return datetime.now(timezone.utc).isoformat()


def work(item):
    sid, raw, annotations, qc = item
    try:
        result = check(json.loads(raw), json.loads(annotations))
    except (OSError, ValueError, KeyError, RuntimeError, IndexError, TypeError) as exc:
        result = {'state':'CHECK_FAILED','reason':str(exc)}
    if result['state'] == 'COMPARED':
        a, b, inter = (result['mapped_near_residues'],result['plinder_neighbor_residues'],result['intersection'])
        result['pocket_agreement'] = ('EXACT_NONEMPTY' if a == b == inter and a > 0 else
            'EMPTY_BOTH' if a == b == 0 else 'DIFFERENT')
    else:
        result['pocket_agreement'] = 'NOT_COMPARED'
    data = json.loads(qc) if qc else None
    source_models = data.get('all_deposited_models') if isinstance(data,dict) else None
    model = 'SOURCE_SINGLE_MODEL_1' if source_models == [1] else 'MULTIPLE_OR_NOT_RECORDED'
    return sid, model, result['state'], result['pocket_agreement'], json.dumps(result,sort_keys=True,separators=(',',':'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--database', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--workers', type=int, default=6)
    p.add_argument('--limit', type=int, help='pilot only')
    args = p.parse_args()
    if args.workers < 1:
        raise ValueError('workers must be positive')
    if args.output.exists():
        raise ValueError('output already exists')
    start = time.monotonic()
    con = sqlite3.connect(args.database.resolve().as_uri()+'?mode=ro',uri=True)
    dst = sqlite3.connect(args.output)
    dst.execute('create table site_check (system_id text primary key,model_evidence text,state text,pocket_agreement text,details text)')
    sql = """select result.system_id,observation.document,result.annotations,result.new_qc
             from result join observation using(system_id)
             where result.match_level='SAME_REVISION_INSTANCE_CHAIN_MATCH'
             order by result.system_id"""
    if args.limit is not None:
        sql += ' limit '+str(args.limit)
    started = utc()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for n, row in enumerate(pool.map(work, con.execute(sql), chunksize=16),1):
            dst.execute('insert into site_check values (?,?,?,?,?)',row)
            if n % 1000 == 0:
                dst.commit()
            if n % 10000 == 0:
                print('checked',n,'elapsed_s',round(time.monotonic()-start,1),flush=True)
    dst.commit()
    counts = dst.execute('select model_evidence,state,pocket_agreement,count(*) from site_check group by 1,2,3').fetchall()
    info = {'started_utc':started,'ended_utc':utc(),'elapsed_s':round(time.monotonic()-start,2),
            'workers':args.workers,'limit':args.limit,'groups':counts}
    dst.execute('create table run (document text not null)')
    dst.execute('insert into run values (?)',(json.dumps(info),));dst.commit();dst.close();con.close()
    print(json.dumps(info),flush=True)


if __name__ == '__main__':
    main()
