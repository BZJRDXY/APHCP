"""Feed verified fixed archive files to the rolling M3 source cache over SSH."""
import argparse
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor, as_completed
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


def collect(args, rows=None):
    from aphcp_control.transport import load_settings,ssh_command
    project=Path(__file__).resolve().parents[1]; settings=load_settings(project)
    if rows is None: rows=json.loads(args.input.read_bytes())
    seen={}
    for item in rows:
        src=item.get('source',item); expected=pin(src['pdb_id'],src['revision'])
        if any(src[k]!=expected[k] for k in KEYS) or expected['pdb_id'] in seen:
            raise ValueError('Manifest must contain distinct fixed archive pins')
        seen[expected['pdb_id']]=expected
    if not seen or len(seen)>2000: raise ValueError('Source feed accepts 1–2000 Entries per invocation')
    for path in (args.python,args.receiver):
        if not path.is_absolute() or '..' in path.parts: raise ValueError('Use absolute server runtime paths')
    ssh=ssh_command(settings)+[shlex.join([str(args.python),'-B',str(args.receiver),'receive',
                                           settings['server']['data_root'],args.cache])]
    check=subprocess.run(ssh,input=b'',capture_output=True,timeout=30)
    if check.returncode: raise RuntimeError(check.stderr.decode(errors='replace')[:500])
    cached=json.loads(check.stdout)['available']; curl=shutil.which('curl.exe') or shutil.which('curl')
    if not curl: raise RuntimeError('curl is required on the alternate download machine')
    t=time.monotonic(); errors={}; loaded={}; send={}
    targets={p:r for p,r in seen.items() if cached.get(p)!=r['revision']}
    with tempfile.TemporaryDirectory(prefix='aphcp-source-feed-') as temp:
        def get(row):
            p=row['pdb_id']; name=p+'-v'+row['revision']+'.cif.gz'; target=Path(temp)/name
            try:
                result=subprocess.run([curl,'-q','--proto','=https','--fail','--silent','--show-error',
                          '--max-time',str(args.budget),'--connect-timeout','12','--max-filesize',str(MAX_FILE),
                          '--header','Accept-Encoding: identity','--output',str(target),'--write-out','%{http_code}',row['source_url']],
                          timeout=args.budget+2,capture_output=True)
                if result.returncode or result.stdout!=b'200' or not 0<target.stat().st_size<=MAX_FILE:
                    return p,None,'archive_transfer_failed'
                return p,target,None
            except (OSError,subprocess.TimeoutExpired): return p,None,'archive_transport_error'
        if targets:
            with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                for work in as_completed([pool.submit(get,r) for r in targets.values()]):
                    p,f,error=work.result()
                    if error: errors[p]=error
                    else: send[p]=f
        fetch_s=round(time.monotonic()-t,3); xfer=time.monotonic()
        if send:
            packet=Path(temp)/'feed.tar'
            with tarfile.open(packet,'w') as tf:
                for f in sorted(send.values()):
                    entry=tf.gettarinfo(str(f),arcname=f.name); entry.pax_headers={'APHCP.sha256':sha(f)}
                    with f.open('rb') as s: tf.addfile(entry,s)
            with packet.open('rb') as stream:
                child=subprocess.run(ssh,stdin=stream,capture_output=True,timeout=1800)
            if child.returncode: raise RuntimeError(child.stderr.decode(errors='replace')[:700])
            received=json.loads(child.stdout); loaded=received['loaded']; errors.update(received['errors'])
    result={'cache':args.cache,'input':len(seen),'reused':len(seen)-len(targets),'received':len(loaded),
            'received_bytes':sum(loaded.values()),'errors':errors,'fetch_s':fetch_s,
            'transfer_verify_s':round(time.monotonic()-xfer,3),'elapsed_s':round(time.monotonic()-t,3)}
    print(json.dumps(result,sort_keys=True),flush=True)


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
    return {'active':True,'sealed':sealed,'partitions':output,'feed_reserve':cfg['feed_reserve']}


@contextmanager
def client_lease(args):
    from aphcp_control.transport import load_settings,ssh_command
    settings=load_settings(Path(__file__).resolve().parents[1])
    if any(not p.is_absolute() or '..' in p.parts for p in (args.python,args.receiver)):
        raise ValueError('Use absolute server runtime paths')
    cmd=[*ssh_command(settings),shlex.join([str(args.python),'-B',str(args.receiver),'hold',
                 settings['server']['data_root'],args.cache])]
    with subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE) as child:
        if child.stdout.readline()!=b'locked\n': raise RuntimeError('Could not acquire the alternate source feed lease')
        try: yield
        finally:
            child.stdin.close(); child.wait(timeout=20)


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
    from aphcp_control.transport import load_settings,ssh_command
    settings=load_settings(Path(__file__).resolve().parents[1]); done=set()
    if any(not p.is_absolute() or '..' in p.parts for p in (args.python,args.receiver)):
        raise ValueError('Use absolute server runtime paths')
    probe=[*ssh_command(settings),shlex.join([str(args.python),'-B',str(args.receiver),'queue',
               settings['server']['data_root'],args.controller,args.cache,str(args.stop_after)])]
    while True:
        state=json.loads(subprocess.check_output(probe,timeout=30))
        if not state['active']: return
        if args.jobs>state['feed_reserve']: raise ValueError('Alternate download jobs exceed the reserved archive slots')
        for label, items in state['partitions'].items():
            if label not in done:
                collect(args,items); done.add(label)
                # Refresh the approved window after each feed to respect pause/stop.
                break
        else: time.sleep(5)


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
    loaded={}; errors={}; available={}
    with (run/'receive.lock').open('ab') as lease:
        fcntl.flock(lease,fcntl.LOCK_EX)
        for f in proofs.glob('*.json'):
            saved=json.loads(f.read_bytes()); p=saved['identity']
            body=run/saved['relative_path']
            if saved['status']=='PASS' and body.is_file() and body.stat().st_size==saved['bytes']:
                available[p['pdb_id']]=p['revision']
        if sys.stdin.buffer.peek(1):
            with tarfile.open(fileobj=sys.stdin.buffer,mode='r|') as packet:
                for member in packet:
                    match=re.fullmatch(r'([0-9][A-Z0-9]{3})-v(\d+\.\d+)\.cif\.gz',member.name)
                    checksum=member.pax_headers.get('APHCP.sha256','')
                    if not member.isfile() or not match or not re.fullmatch(r'[0-9a-f]{64}',checksum) or not 0<member.size<=MAX_FILE:
                        raise ValueError('Unexpected source feed member')
                    pdb,rev=match.groups(); source=pin(pdb,rev)|{'expected_bytes':member.size}
                    body=raw/(member.name.lower()); proof=proofs/(pdb+'.json')
                    if proof.exists() and available.get(pdb)!=rev:
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
                        else:
                            sync_path(temp); temp.replace(body); sync_path(raw)
                            durable_json(proof,{'status':'PASS','gzip_crc':'PASS','identity':source,'bytes':member.size,
                                'sha256':checksum,'relative_path':str(body.relative_to(run)),'observed_revision':observed,
                                'attempts':[],'evidence':'alternate_route_ssh_feed'})
                        available[pdb]=rev; loaded[pdb]=member.size
                    except (OSError,ValueError,RuntimeError,EOFError) as exc: errors[pdb]=str(exc)
                    finally:
                        if temp: temp.unlink(missing_ok=True)
    print(json.dumps({'available':available,'loaded':loaded,'errors':errors},sort_keys=True))


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
    elif args.mode in ('fetch','follow'):
        with client_lease(args): (collect if args.mode=='fetch' else follow)(args)
    else: (receive if args.mode=='receive' else hold)(args)


if __name__=='__main__': main()
