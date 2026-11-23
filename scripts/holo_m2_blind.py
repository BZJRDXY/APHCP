"""Isolated, resumable APHCP blind M2: fetch pinned sources, build and audit."""
import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
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


def fetch_ranges(pin, temp, attempt):
    """Resume large fixed archive bodies by validated byte ranges, when supported."""
    url = pin['source_url']; total = pin['expected_bytes']
    headers = {'Accept-Encoding':'identity','User-Agent':'APHCP-M2-fixed-source/1'}
    with build_opener(NoRedirect()).open(Request(url,headers=headers | {'Range':'bytes=0-0'}),timeout=15) as sample:
        etag = sample.headers.get('ETag','')
        if (sample.status != 206 or sample.geturl()!=url or sample.headers.get('Content-Range')!=f'bytes 0-0/{total}'
                or not re.fullmatch(r'"[^"\r\n]+"',etag)):
            return None
    folder=temp.with_suffix('.ranges'); folder.mkdir(exist_ok=True)
    identity=folder/'identity.json'; seal={'source_url':url,'bytes':total,'etag':etag,'range_size':RANGE_SIZE}
    if identity.exists() and read(identity)!=seal:
        raise ValueError('Archive ETag or range plan changed for the same fixed revision')
    if not identity.exists():
        atomic_json(identity,seal)
    attempt.update(http_status=206, etag=etag, transfer='ranges')

    def receive(start):
        end=min(start+RANGE_SIZE,total)-1; n=end-start+1
        out=folder/f'{start:012d}.range'; partial=out.with_suffix('.tmp')
        if out.is_file() and out.stat().st_size==n:
            return (out,0,None)
        with RANGE_SLOTS:
            req=Request(url,headers=headers | {'Range':f'bytes={start}-{end}','If-Match':etag})
            received=0; deadline=time.monotonic()+65
            try:
                with build_opener(NoRedirect()).open(req,timeout=12) as resp, partial.open('wb') as f:
                    if (resp.status!=206 or resp.geturl()!=url or resp.headers.get('ETag')!=etag or
                            resp.headers.get('Content-Range')!=f'bytes {start}-{end}/{total}' or
                            resp.headers.get('Content-Length')!=str(n)):
                        raise ValueError('Archive range identity, ETag or length differs from the pin')
                    while received<n:
                        if time.monotonic()>deadline:
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
    with ThreadPoolExecutor(max_workers=min(12,len(offsets))) as workers:
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



def fetch(run, r):
    p = r['pdb_id']; raw = raw_path(run, r); proof = run/'source_proofs'/f'{p}.json'
    raw.parent.mkdir(parents=True, exist_ok=True)
    previous = read(proof) if proof.exists() else {}
    if (previous.get('status') == 'PASS' and previous['identity'] == r and raw.is_file()
            and raw.stat().st_size == r['expected_bytes'] and previous['sha256'] == sha256(raw)):
        return p, 'reused_verified'
    attempts = previous.get('attempts', [])
    temp = raw.with_suffix('.part.gz')
    for i in range(3):
        attempt = {'started_utc':now(),'body_bytes':0}; attempts.append(attempt)
        try:
            range_storage=None
            if i == 0 and raw.is_file() and raw.stat().st_size==r['expected_bytes']:
                held = raw
            else:
                held = temp
                if r['expected_bytes']>=RANGE_MINIMUM:
                    range_storage=fetch_ranges(r,temp,attempt)
                if range_storage is None:
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
            if held !=raw:
                held.replace(raw)
            attempt.update(finished_utc=now(),status='PASS')
            atomic_json(proof,{'status':'PASS','identity':r,'sha256':digest,'gzip_crc':'PASS','bytes':raw.stat().st_size,
                               'relative_path':str(raw.relative_to(run)),'attempts':attempts,'observed_revision':rev})
            if range_storage:
                for part in range_storage.glob('*.range'):
                    part.unlink()
                (range_storage/'identity.json').unlink()
                range_storage.rmdir()
            return p, 'downloaded' if attempt['body_bytes'] else 'recovered_verified'
        except Exception as exc:
            temp.unlink(missing_ok=True)
            attempt.update(finished_utc=now(),status='FAIL',error=f'{type(exc).__name__}: {exc}')
            atomic_json(proof,{'status':'FAIL','identity':r,'attempts':attempts,'error':attempt['error']})
            if i<2:
                time.sleep(i+1)
    return p, 'failed'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('run',type=Path); ap.add_argument('--index',type=Path)
    ap.add_argument('--jobs',type=int,default=12); ap.add_argument('--fetch-jobs',type=int,default=16)
    ap.add_argument('--entries',nargs='+'); ap.add_argument('--mode',choices=('all','fetch','build','audit'),default='all')
    args=ap.parse_args(); run=args.run.resolve(); start=time.monotonic()
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
        fetch_counts=Counter()
        with ThreadPoolExecutor(max_workers=args.fetch_jobs) as pool:
            tasks=[pool.submit(fetch,run,r) for r in plan if r['pdb_id'] in selected]
            for n,f in enumerate(as_completed(tasks),1):
                p,result=f.result(); fetch_counts[result]+=1
                if n==len(tasks) or n%20==0 or result=='failed':
                    print('fetch',n,'/',len(tasks),p,dict(fetch_counts),round(time.monotonic()-start,1),flush=True)
        atomic_json(run/'last_fetch_invocation.json',{'finished_utc':now(),'counts':dict(fetch_counts),'elapsed_s':round(time.monotonic()-start,2)})
        if args.mode=='all':
            if not args.index:
                ap.error('--index is required')
            from holo_m2_blind_structure import build
            build(run,plan,selected,args.index,args.jobs,dict(fetch_counts),start)


if __name__=='__main__':
    main()



