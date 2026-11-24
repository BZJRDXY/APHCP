"""Feed verified fixed archive files to the rolling M3 source cache over SSH."""
import argparse
from contextlib import contextmanager
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
from threading import Timer
import time

NAME=re.compile(r'holo-m2-prefetch-[a-z0-9-]{1,60}')
PDB=re.compile(r'[0-9][A-Z0-9]{3}')
CONTROLLER=re.compile(r'holo-m3-[a-z0-9-]{1,80}')
REV=re.compile(r'\d+\.\d+')
MAX_FILE=256 << 20
KEYS=('pdb_id','revision','source_url','file_type')


def pin(pdb, rev):
    if not PDB.fullmatch(pdb) or not REV.fullmatch(rev):
        raise ValueError('Invalid fixed archive identity')
    low=pdb.lower(); base='pdb_0000'+low
    return {'pdb_id':pdb,'revision':rev,'file_type':'pdbx_mmcif_gzip',
            'source_url':f'https://files-versioned.rcsb.org/pdb_versioned/data/entries/{low[1:3]}/{base}/{base}_xyz_v{rev.replace(".","-")}.cif.gz'}


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


class FeedConnectionError(RuntimeError):
    pass


def remote_command(args, *command):
    from aphcp_control.transport import load_settings, ssh_command
    settings=load_settings(Path(__file__).resolve().parents[1])
    if any(not p.is_absolute() or '..' in p.parts for p in (args.python,args.receiver)):
        raise ValueError('Use absolute server runtime paths')
    ssh=ssh_command(settings)
    return [ssh[0],'-o','ServerAliveInterval=5','-o','ServerAliveCountMax=2',*ssh[1:],
            shlex.join([str(args.python),'-B',str(args.receiver),command[0],settings['server']['data_root'],*command[1:]])]


def ssh_json(command, stdin=b'', timeout=30):
    try:
        inputs={'input':stdin} if isinstance(stdin,bytes) else {'stdin':stdin}
        child=subprocess.run(command,**inputs,capture_output=True,timeout=timeout)
    except (OSError,subprocess.TimeoutExpired) as exc:
        raise FeedConnectionError('alternate SSH timed out or was unavailable') from exc
    if child.returncode:
        raise FeedConnectionError(f'alternate SSH command failed ({child.returncode})')
    return json.loads(child.stdout)


def deliver(command, directory, files):
    packet=directory/'feed.tar'
    with tarfile.open(packet,'w') as archive:
        for path in sorted(files):
            entry=archive.gettarinfo(str(path),arcname=path.name)
            entry.pax_headers={'APHCP.sha256':sha(path)}
            with path.open('rb') as stream:
                archive.addfile(entry,stream)
    with packet.open('rb') as stream:
        return ssh_json(command,stream,timeout=120)


def collect(args, rows=None, *, available=None, partition=None):
    if rows is None:
        rows=json.loads(args.input.read_bytes())
    seen={}
    for item in rows:
        src=item.get('source',item); expected=pin(src['pdb_id'],src['revision'])
        if any(src[k]!=expected[k] for k in KEYS) or expected['pdb_id'] in seen:
            raise ValueError('Manifest must contain distinct fixed archive pins')
        seen[expected['pdb_id']]=expected
    if not seen or len(seen)>2000:
        raise ValueError('Source feed accepts 1–2000 Entries per invocation')
    ssh=remote_command(args,'receive',args.cache)
    cached=ssh_json(ssh)['available'] if available is None else available
    curl=shutil.which('curl.exe') or shutil.which('curl')
    if not curl:
        raise RuntimeError('curl is required on the alternate download machine')
    start=time.monotonic(); errors={}; loaded={}; send={}; replayed=0; transfer_s=0; first_delivery_s=None
    targets={p:r for p,r in seen.items() if cached.get(p)!=r['revision']}
    with tempfile.TemporaryDirectory(prefix='aphcp-source-feed-') as temp:
        def get(row):
            p=row['pdb_id']; target=Path(temp)/(p+'-v'+row['revision']+'.cif.gz')
            try:
                result=subprocess.run([curl,'-q','--proto','=https','--fail','--silent','--show-error',
                          '--max-time',str(args.budget),'--connect-timeout','12','--max-filesize',str(MAX_FILE),
                          '--header','Accept-Encoding: identity','--output',str(target),'--write-out','%{http_code}',row['source_url']],
                          timeout=args.budget+2,capture_output=True)
                if result.returncode or result.stdout!=b'200' or not 0<target.stat().st_size<=MAX_FILE:
                    return p,None,'archive_transfer_failed'
                return p,target,None
            except (OSError,subprocess.TimeoutExpired):
                return p,None,'archive_transport_error'

        if targets:
            with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                pending={pool.submit(get,row) for row in targets.values()}
                flush_by=None
                while pending or send:
                    timeout=max(0,flush_by-time.monotonic()) if flush_by is not None else None
                    finished,pending=wait(pending,timeout=timeout,return_when=FIRST_COMPLETED)
                    for work in finished:
                        p,path,error=work.result()
                        if error: errors[p]=error
                        else: send[p]=path
                    if send and flush_by is None:
                        flush_by=time.monotonic()+4
                    if not send or (len(send)<4 and pending and time.monotonic()<flush_by):
                        continue
                    begin=time.monotonic()
                    try:
                        received=deliver(ssh,Path(temp),send.values())
                    except FeedConnectionError:
                        print(json.dumps({'cache':args.cache,'partition':partition,'delivery':'ACK_LOST',
                                          'sent_files':len(send),'sent_bytes':sum(f.stat().st_size for f in send.values())}),flush=True)
                        for work in pending: work.cancel()
                        raise
                    transfer_s+=time.monotonic()-begin
                    if first_delivery_s is None: first_delivery_s=time.monotonic()-start
                    if not set(received['loaded'])<=send.keys() or not set(received['errors'])<=send.keys():
                        raise ValueError('Receiver acknowledgement did not match sent pins')
                    loaded.update(received['loaded']); errors.update(received['errors']); cached.update(received['available'])
                    replayed+=received.get('replayed',0)
                    for p,path in send.items():
                        if cached.get(p)!=targets[p]['revision']:
                            errors.setdefault(p,'archive_delivery_unverified')
                        path.unlink()
                    send.clear(); flush_by=None
    result={'cache':args.cache,'partition':partition,'input':len(seen),'reused':len(seen)-len(targets),'received':len(loaded),
            'received_bytes':sum(loaded.values()),'replayed_delivery':replayed,'errors':errors,
            'first_delivery_s':round(first_delivery_s,3) if first_delivery_s is not None else None,
            'transfer_verify_s':round(transfer_s,3),'elapsed_s':round(time.monotonic()-start,3)}
    print(json.dumps(result,sort_keys=True),flush=True)
    return result


def cached_revision(cache, pdb):
    proof=cache/'source_proofs'/f'{pdb}.json'
    if proof.is_file():
        saved=json.loads(proof.read_bytes()); pin_data=saved['identity']
        body=cache/saved['relative_path']
        if (pin_data['pdb_id']==pdb and saved['status']==saved.get('gzip_crc')=='PASS'
                and pin_data['revision']==saved.get('observed_revision') and body.is_file()
                and body.stat().st_size==saved['bytes']):
            return pin_data['revision']
    return None


def bounded_queue(root, controller, cache, last):
    runs=(Path(root).resolve(strict=True)/'runs').resolve(strict=True); run=runs/controller
    if run.resolve()!=run or not run.is_dir(): raise ValueError('Unknown rolling Run')
    read=lambda p:json.loads(p.read_bytes())
    cfg=read(run/'rolling_config.json'); pg=read(run/'rolling_progress.json'); sup=read(run/'supervisor_status.json')
    stop=cfg.get('source_feed_stop')
    if cfg.get('source_feed')!=cache or not stop or last>stop:
        raise ValueError('The controller must use this feed and approve the requested stop')
    sealed=pg['partitions_completed']
    if sealed>=last or sup.get('state')!='RUNNING' or (run/'pause.request.json').is_file():
        return {'active':False,'sealed':sealed,'partitions':{}}
    pid=sup.get('controller_pid')
    try:
        argv=(Path('/proc')/str(int(pid))/'cmdline').read_bytes()
        if controller.encode() not in argv or b'holo_m3_rolling.py' not in argv:
            return {'active':False,'sealed':sealed,'partitions':{}}
    except (OSError,TypeError,ValueError): return {'active':False,'sealed':sealed,'partitions':{}}
    parts=read(run/'rolling_plan.json')['partitions']; output={}
    for part in parts[sealed:min(sealed+2,last)]:
        owner=run; visited=set()
        while True:
            if owner.parent!=runs or owner in visited: raise ValueError('Invalid rolling Run history')
            visited.add(owner); path=owner/'partitions'/(part['id']+'.json')
            if path.is_file(): break
            owner=Path(read(owner/'rolling_config.json')['continued_from']).resolve(strict=True)
        if sha(path)!=part['sha256']: raise ValueError('Frozen production partition changed')
        output[part['id']]=[r['source'] for r in read(path)]
    verified={row['pdb_id']:row['revision'] for items in output.values() for row in items
              if cached_revision(runs/cache,row['pdb_id'])==row['revision']}
    return {'active':True,'sealed':sealed,'partitions':output,'available':verified,'feed_reserve':cfg['feed_reserve']}


@contextmanager
def client_lease(args):
    cmd=remote_command(args,'hold',args.cache)
    try:
        with subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE) as child:
            guard=Timer(30,child.kill); guard.start()
            try:
                if child.stdout.readline()!=b'locked\n':
                    raise FeedConnectionError('Could not acquire the alternate source feed lease')
            finally:
                guard.cancel()
            try:
                yield
            finally:
                child.stdin.close()
                try: child.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    child.kill(); child.wait()
    except OSError as exc:
        raise FeedConnectionError('alternate source feed lease connection lost') from exc


def hold(args):
    import fcntl
    runs=(Path(args.root).resolve(strict=True)/'runs').resolve(strict=True); run=runs/args.cache
    if run.resolve()!=run: raise ValueError('Source feed directory must be inside runs')
    run.mkdir(exist_ok=True)
    with (run/'client.lock').open('ab') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        print('locked',flush=True)
        while sys.stdin.buffer.read(4096): pass


def follow(args):
    probe=remote_command(args,'queue',args.controller,args.cache,str(args.stop_after))
    failures=Counter(); consecutive=0
    while True:
        try:
            with client_lease(args):
                while True:
                    state=ssh_json(probe)
                    if not state['active']: return
                    available=state['available']
                    if args.jobs>state['feed_reserve']:
                        raise ValueError('Alternate download jobs exceed the reserved archive slots')
                    pins={(r['pdb_id'],r['revision']) for group in state['partitions'].values() for r in group}
                    failures=Counter({key:count for key,count in failures.items() if key in pins})
                    for label,items in state['partitions'].items():
                        outstanding=[r for r in items if available.get(r['pdb_id'])!=r['revision']
                                     and failures[r['pdb_id'],r['revision']]<2]
                        if not outstanding: continue
                        group=outstanding[:3*args.jobs]
                        result=collect(args,group,available=available,partition=label)
                        for row in group:
                            key=row['pdb_id'],row['revision']
                            if key[0] in result['errors']: failures[key]+=1
                        consecutive=0
                        break
                    else:
                        consecutive=0; time.sleep(5)
        except FeedConnectionError:
            consecutive+=1
            if consecutive>5: raise
            delay=min(2**consecutive,20)
            print(json.dumps({'cache':args.cache,'state':'RECONNECTING','attempt':consecutive,'delay_s':delay}),flush=True)
            time.sleep(delay)


def receive(args):
    import fcntl
    from holo_archive_runtime import durable_json,sync_path
    from holo_candidate_batch import revision_from_cif
    runs=(Path(args.root).resolve(strict=True)/'runs').resolve(strict=True); run=runs/args.cache
    if run.resolve()!=run: raise ValueError('Source feed directory must be inside runs')
    run.mkdir(exist_ok=True); raw=run/'raw'; proofs=run/'source_proofs'
    for path in (raw,proofs):
        if path.resolve()!=path: raise ValueError('Invalid source feed directory')
        path.mkdir(exist_ok=True)
    loaded={}; errors={}; available={}; replayed=0
    with (run/'receive.lock').open('ab') as lease:
        fcntl.flock(lease,fcntl.LOCK_EX)
        if not sys.stdin.buffer.peek(1):
            available={f.stem: rev for f in proofs.glob('*.json') if (rev:=cached_revision(run,f.stem))}
        else:
            with tarfile.open(fileobj=sys.stdin.buffer,mode='r|') as packet:
                for member in packet:
                    match=re.fullmatch(r'([0-9][A-Z0-9]{3})-v(\d+\.\d+)\.cif\.gz',member.name)
                    checksum=member.pax_headers.get('APHCP.sha256','')
                    if not member.isfile() or not match or not re.fullmatch(r'[0-9a-f]{64}',checksum) or not 0<member.size<=MAX_FILE:
                        raise ValueError('Unexpected source feed member')
                    pdb,rev=match.groups(); source=pin(pdb,rev)|{'expected_bytes':member.size}
                    body=raw/(member.name.lower()); proof=proofs/(pdb+'.json')
                    if proof.exists() and cached_revision(run,pdb)!=rev:
                        errors[pdb]='revision_conflicts_with_cache'; continue
                    temp=None
                    try:
                        with tempfile.NamedTemporaryFile(suffix='.cif.gz',dir=raw,delete=False) as f:
                            temp=Path(f.name)
                            with packet.extractfile(member) as obj: shutil.copyfileobj(obj,f)
                        if temp.stat().st_size!=member.size or sha(temp)!=checksum:
                            raise ValueError('transport_sha256_or_length_mismatch')
                        with gzip.open(temp,'rb') as stream:
                            while stream.read(1<<20): pass
                        _, observed=revision_from_cif(temp,pdb)
                        if observed!=rev: raise ValueError('pinned_revision_mismatch')
                        if proof.exists():
                            previous=json.loads(proof.read_bytes())
                            if previous['sha256']!=checksum or not body.is_file() or sha(body)!=checksum:
                                raise ValueError('fixed_archive_differs_from_existing_cache')
                            replayed+=1
                        else:
                            sync_path(temp); temp.replace(body); sync_path(raw)
                            durable_json(proof,{'status':'PASS','gzip_crc':'PASS','identity':source,'bytes':member.size,
                                'sha256':checksum,'relative_path':str(body.relative_to(run)),'observed_revision':observed,
                                'attempts':[],'evidence':'alternate_route_ssh_feed'})
                            loaded[pdb]=member.size
                        available[pdb]=rev
                    except (OSError,ValueError,RuntimeError,EOFError) as exc: errors[pdb]=str(exc)
                    finally:
                        if temp: temp.unlink(missing_ok=True)
    print(json.dumps({'available':available,'loaded':loaded,'errors':errors,'replayed':replayed},sort_keys=True))


def main():
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='mode',required=True)
    a=sub.add_parser('fetch'); a.add_argument('input',type=Path); a.add_argument('cache')
    b=sub.add_parser('follow'); b.add_argument('controller'); b.add_argument('cache'); b.add_argument('--stop-after',type=int,required=True)
    for q in (a,b):
        q.add_argument('--python',required=True,type=PurePosixPath); q.add_argument('--receiver',required=True,type=PurePosixPath)
        q.add_argument('--jobs',type=int,default=6); q.add_argument('--budget',type=int,default=60)
    for mode in ('receive','hold'):
        a=sub.add_parser(mode); a.add_argument('root'); a.add_argument('cache')
    a=sub.add_parser('queue'); a.add_argument('root'); a.add_argument('controller'); a.add_argument('cache'); a.add_argument('stop_after',type=int)
    args=p.parse_args()
    if not NAME.fullmatch(args.cache) or (args.mode in ('fetch','follow') and (not 1<=args.jobs<=64 or not 10<=args.budget<=180)):
        p.error('Invalid cache name or transfer bounds')
    if args.mode in ('queue','follow') and (not CONTROLLER.fullmatch(args.controller) or args.stop_after<1):
        p.error('Invalid rolling Run or absolute stop')
    if args.mode=='queue': print(json.dumps(bounded_queue(args.root,args.controller,args.cache,args.stop_after)))
    elif args.mode=='follow': follow(args)
    elif args.mode=='fetch':
        with client_lease(args): collect(args)
    else: (receive if args.mode=='receive' else hold)(args)


if __name__=='__main__': main()
