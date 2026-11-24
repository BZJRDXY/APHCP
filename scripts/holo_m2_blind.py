"""Isolated, resumable APHCP blind M2: fetch pinned sources, build and audit."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter
import csv
import fcntl
import gzip
from http.client import HTTPException
import json
from pathlib import Path
import re
import sys
from threading import BoundedSemaphore
import time
from urllib.request import HTTPRedirectHandler, Request, build_opener

from holo_archive_runtime import Archive, durable_json, sync_path
from holo_candidate_batch import atomic_json, now, revision_from_cif, sha256


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def freeze(run):
    a, b = [run/'inputs'/x for x in ('blind.json', 'archive_head.json')]
    manifest, evidence = read(a), read(b)
    def key(r):
        return tuple(r[k] for k in ('pdb_id', 'revision', 'source_url', 'file_type'))
    heads = {key(r): r for r in evidence}
    if not manifest or len(heads) != len(evidence) or set(heads) != {key(r) for r in manifest}:
        raise ValueError('Frozen manifest and archive proof disagree')
    plan = []
    for r in manifest:
        p, ver, url, kind = key(r)
        low = p.lower(); base = 'pdb_0000'+low
        if (set(r) != {'pdb_id','revision','source_url','file_type'} or not re.fullmatch(r'[0-9][A-Z0-9]{3}',p)
                or not re.fullmatch(r'\d+\.\d+', ver) or kind != 'pdbx_mmcif_gzip'
                or url != f'https://files-versioned.rcsb.org/pdb_versioned/data/entries/{low[1:3]}/{base}/{base}_xyz_v{ver.replace(".","-")}.cif.gz'):
            raise ValueError(f'Unrecognized exact source {p}')
        h = heads[key(r)]
        if h['status'] != 200 or not isinstance(h['content_length'], int) or h['content_length'] <= 0:
            raise ValueError(f'Missing archive size proof: {p}')
        plan.append(r | {'expected_bytes': h['content_length']})
    if len({r['pdb_id'] for r in plan}) != len(plan):
        raise ValueError('The frozen list contains repeated Entry IDs')
    frozen = run/'frozen.json'
    details = {'run_id':run.name,'input_sha256':sha256(a),'head_sha256':sha256(b),
               'plan_count':len(plan),'expected_bytes':sum(r['expected_bytes'] for r in plan),'files':plan}
    if frozen.exists():
        if any(read(frozen).get(k) != v for k,v in details.items()):
            raise ValueError('Frozen input changed')
    else:
        atomic_json(frozen, details | {'frozen_utc':now()})
    return plan


class NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def raw_path(run, r):
    return run/'raw'/f"{r['pdb_id'].lower()}-v{r['revision']}.cif.gz"


RANGE_SIZE = 512 << 10
RANGE_MINIMUM = 8 << 20
RANGE_SLOTS = BoundedSemaphore(24)


def fetch_ranges(pin, temp, attempt, archive=None, deadline=None):
    """Resume fixed archive bodies by validated byte ranges, when supported."""
    url = pin['source_url']; total = pin['expected_bytes']
    headers = {'Accept-Encoding':'identity','User-Agent':'APHCP-M2-fixed-source/1'}
    if archive:
        trace={}; attempt.setdefault('requests', []).append(trace)
        fields=archive.request(url, 'range', '/dev/null', 15, trace, headers=('Range: bytes=0-0',), deadline=deadline)
        status=trace.get('http_status'); etag=fields.get('etag',''); span=fields.get('content-range')
        attempt['body_bytes']+=trace.get('body_bytes',0)
    else:
        with build_opener(NoRedirect()).open(Request(url,headers=headers | {'Range':'bytes=0-0'}),timeout=15) as sample:
            status=sample.status if sample.geturl()==url else 0; etag=sample.headers.get('ETag',''); span=sample.headers.get('Content-Range')
    if status != 206 or span != f'bytes 0-0/{total}' or not re.fullmatch(r'"[^"\r\n]+"',etag):
        return None
    folder=temp.with_suffix('.ranges'); folder.mkdir(exist_ok=True)
    identity=folder/'identity.json'; seal={'source_url':url,'bytes':total,'etag':etag,'range_size':RANGE_SIZE}
    if identity.exists() and read(identity)!=seal:
        raise ValueError('Archive ETag or range plan changed for the same fixed revision')
    if not identity.exists():
        (durable_json if archive else atomic_json)(identity,seal)
    attempt.update(http_status=206, etag=etag, transfer='ranges')

    def receive(start):
        end=min(start+RANGE_SIZE,total)-1; n=end-start+1
        out=folder/f'{start:012d}.range'; partial=out.with_suffix('.tmp')
        if out.is_file() and out.stat().st_size==n:
            return (out,0,None)
        if archive:
            trace={}; attempt.setdefault('requests', []).append(trace)
            try:
                fields=archive.request(url,'range',partial,65,trace,
                                       headers=(f'Range: bytes={start}-{end}',f'If-Match: {etag}'),deadline=deadline)
                if (trace['http_status']!=206 or fields.get('etag')!=etag or
                        fields.get('content-range')!=f'bytes {start}-{end}/{total}' or
                        fields.get('content-length')!=str(n) or partial.stat().st_size!=n):
                    raise ValueError('Archive range identity, ETag or length differs from the pin')
                sync_path(partial); partial.replace(out); sync_path(out.parent)
                return out,trace['body_bytes'],None
            except (OSError,ValueError) as e:
                partial.unlink(missing_ok=True)
                return out,trace.get('body_bytes',0),f'{type(e).__name__}: {e}'
        with RANGE_SLOTS:
            req=Request(url,headers=headers | {'Range':f'bytes={start}-{end}','If-Match':etag})
            received=0; socket_deadline=time.monotonic()+65
            try:
                with build_opener(NoRedirect()).open(req,timeout=12) as resp, partial.open('wb') as f:
                    if (resp.status!=206 or resp.geturl()!=url or resp.headers.get('ETag')!=etag or
                            resp.headers.get('Content-Range')!=f'bytes {start}-{end}/{total}' or
                            resp.headers.get('Content-Length')!=str(n)):
                        raise ValueError('Archive range identity, ETag or length differs from the pin')
                    while received<n:
                        if time.monotonic()>socket_deadline:
                            raise TimeoutError('Archive range exceeded its time budget')
                        chunk=resp.read1(min(64<<10,n-received))
                        if not chunk:
                            break
                        f.write(chunk); received+=len(chunk)
                if received!=n:
                    raise ValueError('Incomplete archive range')
                partial.replace(out)
                return (out,received,None)
            except (OSError,ValueError,HTTPException) as exc:
                partial.unlink(missing_ok=True)
                return (out,received,f'{start}-{end}: {type(exc).__name__}: {exc}')
    offsets=range(0,total,RANGE_SIZE)
    with ThreadPoolExecutor(max_workers=min(archive.ancillary if archive else 12,len(offsets))) as workers:
        results=list(workers.map(receive,offsets))
    attempt['body_bytes']+=sum(n for _,n,_ in results)
    attempt['reused_ranges']=sum(n==0 and err is None for _,n,err in results)
    errors=[err for _,_,err in results if err]
    if errors:
        attempt['range_errors']=errors
        raise ValueError(f'{len(errors)} of {len(results)} archive ranges need retry')
    with temp.open('wb') as f:
        for path,_,_ in results:
            f.write(path.read_bytes())
    return folder



def fetch(run, r, archive=None, policy="historical", budget=None, *, prefer_ranges=False):
    p = r['pdb_id']; raw = raw_path(run, r); proof = run/'source_proofs'/f'{p}.json'
    raw.parent.mkdir(parents=True, exist_ok=True)
    proof.parent.mkdir(parents=True, exist_ok=True)
    previous = read(proof) if proof.exists() else {}
    if (previous.get('status') == 'PASS' and previous['identity'] == r and raw.is_file()
            and raw.stat().st_size == r['expected_bytes'] and previous['sha256'] == sha256(raw)):
        return p, 'reused_verified'
    attempts = previous.get('attempts', [])
    temp = raw.with_suffix('.part.gz')
    count = 1 if archive else 3
    save = durable_json if archive else atomic_json
    for i in range(count):
        started=time.monotonic()
        limit=budget if budget is not None else (90 if r['expected_bytes']>=RANGE_MINIMUM or policy=='retry' else 30)
        attempt = {'started_utc':now(),'body_bytes':0}; attempts.append(attempt)
        try:
            range_storage=None
            if i == 0 and raw.is_file() and raw.stat().st_size==r['expected_bytes']:
                held = raw
            else:
                held = temp
                range_threshold=(1<<20) if archive and policy=='retry' and prefer_ranges else RANGE_MINIMUM
                if r['expected_bytes']>=range_threshold and (not archive or policy=='retry'):
                    range_storage=fetch_ranges(r,temp,attempt,archive,started+limit if archive else None)
                if range_storage is None and archive:
                    trace={}; attempt.setdefault('requests', []).append(trace)
                    try:
                        fields=archive.request(r['source_url'], 'retry' if policy=='retry' else 'get', temp, limit, trace,
                                               deadline=started+limit)
                        if trace['http_status']!=200 or int(fields.get('content-length',r['expected_bytes']))!=r['expected_bytes']:
                            raise ValueError('Exact archive response or length differs from HEAD proof')
                    finally:
                        attempt['http_status']=trace.get('http_status'); attempt['body_bytes']+=trace.get('body_bytes',0)
                elif range_storage is None:
                    request = Request(r['source_url'],headers={'Accept-Encoding':'identity','User-Agent':'APHCP-M2-fixed-source/1'})
                    deadline=time.monotonic()+90
                    with build_opener(NoRedirect()).open(request,timeout=15) as resp, temp.open('wb') as f:
                        attempt['http_status'] = resp.status
                        if resp.status!=200 or resp.geturl()!=r['source_url'] or int(resp.headers.get('Content-Length',r['expected_bytes']))!=r['expected_bytes']:
                            raise ValueError('Exact archive response or length differs from HEAD proof')
                        while attempt['body_bytes'] < r['expected_bytes']:
                            if time.monotonic()>deadline:
                                raise TimeoutError('Archive response exceeded its time budget')
                            chunk=resp.read1(min(1<<20,r['expected_bytes']-attempt['body_bytes']))
                            if not chunk:
                                break
                            attempt['body_bytes']+=len(chunk); f.write(chunk)
            if held.stat().st_size!=r['expected_bytes']:
                raise ValueError('Incomplete source')
            with gzip.open(held,'rb') as gz:
                while gz.read(1<<20):
                    pass
            _, rev = revision_from_cif(held,p)
            if rev!=r['revision']:
                raise ValueError(f'Archive Structure model {rev} differs from frozen {r["revision"]}')
            digest=sha256(held)
            if archive: sync_path(held)
            if held !=raw:
                held.replace(raw)
            if archive: sync_path(raw.parent)
            attempt.update(finished_utc=now(),status='PASS',elapsed_s=round(time.monotonic()-started,4))
            save(proof,{'status':'PASS','identity':r,'sha256':digest,'gzip_crc':'PASS','bytes':raw.stat().st_size,
                               'relative_path':str(raw.relative_to(run)),'attempts':attempts,'observed_revision':rev})
            if range_storage:
                for part in range_storage.glob('*.range'):
                    part.unlink()
                (range_storage/'identity.json').unlink()
                range_storage.rmdir()
            return p, 'downloaded' if attempt['body_bytes'] else 'recovered_verified'
        except Exception as exc:
            temp.unlink(missing_ok=True)
            attempt.update(finished_utc=now(),status='FAIL',error=f'{type(exc).__name__}: {exc}',elapsed_s=round(time.monotonic()-started,4))
            save(proof,{'status':'FAIL','identity':r,'attempts':attempts,'error':attempt['error']})
            if i+1<count:
                time.sleep(i+1)
    return p, 'failed'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('run',type=Path); ap.add_argument('--index',type=Path)
    ap.add_argument('--jobs',type=int,default=12); ap.add_argument('--fetch-jobs',type=int,default=16)
    ap.add_argument('--entries',nargs='+'); ap.add_argument('--mode',choices=('all','fetch','build','audit'),default='all')
    ap.add_argument('--archive-slots',type=Path); ap.add_argument('--archive-limit',type=int,default=24)
    ap.add_argument('--archive-ancillary',type=int,default=6); ap.add_argument('--primary-budget',type=float)
    ap.add_argument('--stream-build',action='store_true',help='Build each Entry as soon as its fixed source has been checked')
    args=ap.parse_args(); run=args.run.resolve(); start=time.monotonic()
    archive=Archive(args.archive_slots,args.archive_limit,args.archive_ancillary) if args.archive_slots else None
    if args.primary_budget is not None and (not archive or args.primary_budget<=0): ap.error('positive budget requires archive slots')
    if args.mode=='all' and not args.index: ap.error('--index is required')
    if args.stream_build and (args.mode!='all' or not archive or args.entries is not None):
        ap.error('stream build requires all fixed Entries and shared archive slots')
    with (run/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        plan=freeze(run); selected=set(args.entries or (x['pdb_id'] for x in plan))
        if not selected<=set(r['pdb_id'] for r in plan):
            ap.error('selected Entry is outside the freeze')
        if args.mode=='audit':
            from holo_m2_blind_audit import audit
            audit(run,plan); return
        if args.mode=='build':
            from holo_m2_blind_structure import build
            build(run,plan,selected,args.index,args.jobs,{},start); return
        checkpoint=run/'stream_receipt.json'
        summary=run/'summary.json'
        if args.stream_build and checkpoint.is_file():
            from holo_m2_blind_structure import index_stamp
            previous=read(checkpoint)
            if (previous['frozen_sha256']!=sha256(run/'frozen.json') or previous['summary_sha256']!=sha256(summary) or
                previous['index']!=str(args.index.resolve()) or
                read(run/'index_snapshot.json')['source_tables']!=index_stamp(args.index) or
                any(sha256(Path(__file__).parent/f)!=h for f,h in read(run/'build_inputs.json')['code_files'].items())):
                raise ValueError('Completed M2 stream changed; use a new Run')
            print('M2 stream checkpoint reused',run.name,flush=True)
            return
        fetch_counts=Counter(); total=len(selected)
        def record(p,result,last):
            fetch_counts[result]+=1; n=sum(fetch_counts.values())
            if last or n%20==0 or result=='failed':
                print('fetch',n,'/',total,p,dict(fetch_counts),round(time.monotonic()-start,1),flush=True)
            if last:
                durable_json(run/'last_fetch_invocation.json',{'finished_utc':now(),'counts':dict(fetch_counts),
                                                                'elapsed_s':round(time.monotonic()-start,2)})
        with ThreadPoolExecutor(max_workers=args.fetch_jobs) as pool:
            def downloads():
                return [pool.submit(fetch,run,r,archive,"primary",
                                    args.primary_budget if r["expected_bytes"]<RANGE_MINIMUM else None)
                        for r in plan if r['pdb_id'] in selected]
            if args.stream_build:
                from holo_m2_blind_structure import build
                build(run,plan,selected,args.index,args.jobs,fetch_counts,start,downloads=downloads,on_fetch=record)
                durable_json(checkpoint,{'frozen_sha256':sha256(run/'frozen.json'),'summary_sha256':sha256(summary),
                                          'index':str(args.index.resolve())})
                return
            for n,f in enumerate(as_completed(downloads()),1):
                p,result=f.result(); record(p,result,n==total)
        if args.mode=='all':
            from holo_m2_blind_structure import build
            build(run,plan,selected,args.index,args.jobs,fetch_counts,start)


if __name__=='__main__':
    main()



