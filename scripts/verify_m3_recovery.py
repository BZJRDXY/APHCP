"""Run isolated Linux controls for the rolling downloader and its commit journal."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack, closing
import fcntl
import gzip
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import multiprocessing as mp
import os
import signal
import subprocess
import sys
from pathlib import Path
import sqlite3
import tempfile
from threading import Thread
import time
from types import SimpleNamespace
from unittest.mock import patch

from holo_archive_runtime import Archive, durable_json
from holo_m2_blind import fetch
import holo_m3_rolling as rolling
from holo_m3_register import digest, read, sha
from holo_m3_supervise import supervise

PDB='1ABC'
CIF=b'''data_1ABC
_entry.id 1ABC
loop_
_pdbx_audit_revision_history.data_content_type
_pdbx_audit_revision_history.major_revision
_pdbx_audit_revision_history.minor_revision
'Structure model' 1 2
'''
BODY=gzip.compress(CIF,mtime=0)


class Server(BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def respond(self, get):
        path=self.path.partition('?')[0]
        if path=='/headers': time.sleep(.8)
        if path=='/cap': time.sleep(.26)
        data=BODY if not self.headers.get('Range') else BODY[:1]
        try:
            self.send_response(206 if len(data)==1 else 200)
            self.send_header('Content-Length',str(len(data)))
            self.end_headers()
            if get and path=='/body':
                self.wfile.write(data[:1]); self.wfile.flush(); time.sleep(.8); data=data[1:]
            if get: self.wfile.write(data)
        except (BrokenPipeError,ConnectionResetError):
            pass

    def do_HEAD(self): self.respond(False)
    def do_GET(self): self.respond(True)


def take(root, url, kind, index, queue):
    trace={}; a=Archive(root/'slots',4,2,local_test=True)
    headers=['Range: bytes=0-0'] if kind=='range' else []
    fields=a.request(url+'/cap',kind,root/f'http-{index}',2,trace,headers=headers)
    assert int(fields['content-length'])==(1 if headers else len(BODY))
    queue.put(trace)


def transport(root,url):
    a=Archive(root/'slots',4,2,local_test=True)
    deadlines={}
    for target in ('/headers','/body'):
        trace={}; t=time.monotonic()
        try:
            a.request(url+target,'get',root/'test-http.part',.25,trace,deadline=t+.25)
            raise AssertionError('unbounded request')
        except TimeoutError:
            wall=time.monotonic()-t; assert .18<=wall<=.65 and trace['exit_code']==28, (wall,trace)
            deadlines[target]=(round(wall,3),trace['body_bytes'])
    q=mp.get_context('spawn').Queue(); ps=[]
    for i,kind in enumerate(['get','head','retry','range','get','head','get','retry','get','range']):
        ps.append(mp.get_context('spawn').Process(target=take,args=(root,url,kind,i,q)))
    for p in ps: p.start()
    traces=[q.get(timeout=10) for _ in ps]
    for p in ps: p.join(timeout=2); assert p.exitcode==0
    def high(xs):
        return max(sum(t['started_monotonic']<=x['started_monotonic']<t['finished_monotonic'] for t in xs) for x in xs)
    max_parallel=high(traces); max_ancillary=high([t for t in traces if t['kind']!='get'])
    assert max_parallel==4 and max_ancillary==2, traces
    pin=dict(pdb_id=PDB,revision='1.2',file_type='pdbx_mmcif_gzip',source_url=url+'/body')
    expected=pin | {'expected_bytes':len(BODY)}; case=root/'slow'
    assert fetch(case,expected,a,'primary',.25)[1]=='failed'
    assert len(read(case/'source_proofs'/f'{PDB}.json')['attempts'])==1
    pin['source_url']=url+'/ok'; row={'source':pin}
    check=a.head(pin); retry=root/'recover'; retry.mkdir()
    # Kill after the written native download proof, before READY publication.
    def crash_on_ready(path, value):
        if path.name=='attempt.json' and value['state']=='READY': raise InterruptedError('after proof')
        durable_json(path,value)
    with patch.object(rolling,'write',crash_on_ready):
        try: rolling.retry_download(retry,row,a,None,check,1)
        except InterruptedError: pass
        else: raise AssertionError('missing crash')
    attempt_file=next((retry/'retry_cache').glob('*/attempt.json'))
    rec=rolling.retry_download(retry,row,a,None,None,1)
    assert rec['state']=='READY' and rec['attempts']==1 and rec['check']['retry_network_bytes']==len(BODY)
    proof=read(attempt_file.parent/'source_proofs'/f'{PDB}.json')
    assert len(proof['attempts'])==1 and proof['observed_revision']=='1.2' and proof['gzip_crc']=='PASS'
    assert rolling.retry_download(retry,row,a,None,None,1)==rec
    rec=rolling.retry_download(retry,row,a,None,None,2,consumed=1)
    assert rec['state']=='READY' and rec['attempts']==2 and rec['check']['retry_network_bytes']==0
    assert rolling.retry_download(retry,row,a,None,None,2,consumed=2) is None
    # Even if shared and main slots are occupied, HEAD and retry can make progress.
    with ExitStack() as locks:
        for slot in (2,3):
            f=locks.enter_context((a.folder/f'{slot:03d}.lock').open('ab'))
            fcntl.flock(f,fcntl.LOCK_EX | fcntl.LOCK_NB)
        with a.slot('retry') as (retry_slot,_):
            assert retry_slot==1 and a.head(pin,.5)['availability']=='VERIFIED'
        with a.slot('head') as (head_slot,_):
            assert head_slot==0
            check=a.head(pin,.18,queue_seconds=.18); assert check['availability']!='VERIFIED' and check['request']['phase']=='QUEUE'
            assert .17<=check['request']['queued_s']<=.28 and 'started_monotonic' not in check['request']
        with a.slot('range',time.monotonic()+.1) as (retry_slot,_): assert retry_slot==1
    return {'timeouts':deadlines,'max_parallel_across_processes':max_parallel,'max_ancillary':max_ancillary,'reserved_auxiliary_lanes_verified':True,'head_queue_budget_verified':True,
            'single_primary_attempt':True,'recovered_committed_source_proof':True,'retry_limit_enforced':True}


def journal(root, url, cut):
    root=root/cut; runs=root/'runs'; runs.mkdir(parents=True)
    seed=runs/'holo-m3-test-seed'; seed.mkdir(); durable_json(seed/'phase_receipt.json',{})
    run=runs/'holo-m3-test-recovery'; run.mkdir(); (run/'partitions').mkdir()
    accept=root/'acceptance.json'; durable_json(accept,{'accepted':True,'catalog_run':seed.name,'pending_technical':[]})
    partitions=[]
    for i,p in enumerate((PDB,'1ABD'),1):
        bid=f'b{i:04d}'; f=run/'partitions'/f'{bid}.json'
        durable_json(f,[{'source':dict(pdb_id=p,revision='1.2',file_type='pdbx_mmcif_gzip',source_url=url+'/ok')}]); partitions.append({'id':bid,'entries':1,'sha256':sha(f)})
    durable_json(run/'rolling_plan.json',{'denominator':2,'exclusions':{},'scope':'test-only','partitions':partitions})
    a=SimpleNamespace(root=root,run=run,seed_run=seed,prior=[],m1=root,index=root,reference_run=root,git_head='test',
                      acceptance=accept,batch_size=1,jobs=1,fetch_jobs=4,reserve_gib=0,primary_budget=1,
                      archive_ancillary=2,retry_jobs=1,max_retries=1,max_pending=4,max_main_batches=None,plan_only=False)
    tripped=False
    def fault(kind):
        nonlocal tripped
        if kind==cut and not tripped:
            tripped=True
            raise InterruptedError(kind)
    def stage(a, run, label, rows, previous, preflight):
        assert preflight['checks'][0]['availability']=='VERIFIED'
        target=runs/f'{run.name}-{label}'; target.mkdir(exist_ok=True)
        with closing(sqlite3.connect(target/'catalog.db')) as db:
            with db:
                db.execute('CREATE TABLE IF NOT EXISTS obs(key TEXT PRIMARY KEY)')
                db.execute('INSERT OR IGNORE INTO obs VALUES(?)',(rows[0]['source']['pdb_id'],))
        fault('sqlite-before-phase')
        durable_json(target/'phase_receipt.json',{'seed':previous.name})
        fault('phase-before-batch')
        return {'catalog_run':target.name,'total':1,'pending':[],'settled':1,'network_bytes':0}
    def once(path, value):
        rolling.holo_original_once(path,value)
        if path.name=='b0001.json' and path.parent.name=='receipts': fault('batch-before-head')
    def local_archive(*args,**kwargs): return Archive(*args,**kwargs,local_test=True)
    rolling.holo_original_once=rolling.once
    with patch.multiple(rolling,stage_batch=stage,once=once,Archive=local_archive):
        try: rolling.rolling(a)
        except InterruptedError: pass
        else: raise AssertionError('fault was not injected')
        assert read(run/'catalog_head.json')['catalog_run']==seed.name
        if cut=='batch-before-head': assert (run/'receipts'/'b0001.json').is_file()
        else: assert not (run/'receipts'/'b0001.json').exists()
        rolling.rolling(a)
        ledger={p.name:sha(p) for p in (run/'receipts').glob('*.json')}
        rolling.rolling(a)
    assert len(ledger)==2 and ledger=={p.name:sha(p) for p in (run/'receipts').glob('*.json')}
    assert read(run/'catalog_head.json')['catalog_run']==run.name+'-b0002'
    assert read(run/'rolling_progress.json')['new_settled_entries']==2
    for part in partitions:
        with closing(sqlite3.connect(runs/f"{run.name}-{part['id']}"/'catalog.db')) as db:
            assert db.execute('SELECT COUNT(*) FROM obs').fetchone()[0]==1
    # A second controller is refused while the old writer's lease is held.
    newer=runs/'holo-m3-test-next'; newer.mkdir()
    a.run=newer; a.continue_from=run
    with (run/'rolling.lock').open('r') as lease:
        fcntl.flock(lease,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try: rolling.rolling(a)
        except BlockingIOError: pass
        else: raise AssertionError('competing controller accepted')
    return {'fault':cut,'same_receipt_hashes_after_replay':True,'unique_records':2,'competing_writer_refused':True}


def scheduler(root, url, *, child_mode=None):
    """Use the real rolling journal and waits, with bounded substitute stage work."""
    runs=root/'runs'; run=runs/'holo-m3-test-scheduler'; seed=runs/'holo-m3-test-seed'
    pin=lambda p: dict(pdb_id=p,revision='1.2',file_type='pdbx_mmcif_gzip',source_url=url+'/ok')
    inputs=(('1ABC','1ABD'),('1ABE','1ABF'),('1ABG',))
    acc=root/'acceptance.json'
    if not (run/'rolling_plan.json').is_file():
        (run/'partitions').mkdir(parents=True); seed.mkdir(); durable_json(seed/'phase_receipt.json',{})
        durable_json(acc,{'accepted':True,'catalog_run':seed.name,'pending_technical':[]})
        parts=[]
        for i,ids in enumerate(inputs,1):
            bid=f'b{i:04d}'; f=run/'partitions'/f'{bid}.json'
            durable_json(f,[{'source':pin(p)} for p in ids]); parts.append({'id':bid,'entries':len(ids),'sha256':sha(f)})
        durable_json(run/'rolling_plan.json',{'denominator':sum(map(len,inputs)),'exclusions':{},'scope':'test-only','partitions':parts})
    a=SimpleNamespace(root=root,run=run,seed_run=seed,prior=[],m1=root,index=root,reference_run=root,git_head='test',
                      acceptance=acc,batch_size=2,jobs=1,fetch_jobs=4,reserve_gib=0,primary_budget=1,
                      archive_ancillary=2,retry_jobs=1,max_retries=1,max_pending=6,max_main_batches=None,
                      stop_after_partition=(1 if child_mode=='stop' else 3),plan_only=False)
    if child_mode in ('balance-kill','balance-timeout'):
        a.balance_retries=True; a.balance_budget=1 if child_mode=='balance-timeout' else 12
        a.max_pending=2; a.retry_only=False
    wave=[]
    waiting=False
    def retry(run,row,archive,previous,known,n,consumed=0,prefer_ranges=False):
        t=time.monotonic(); identity=row['source']; p=identity['pdb_id']; wave.append((p,waiting,t))
        time.sleep(.14)
        state=({'attempts':0,'next_epoch':time.time()+60,'state':'FAILED'} if child_mode=='balance-timeout' else
               {'attempts':1,'next_epoch':0,'state':'READY','check':identity | {'retry_network_bytes':0,'availability':'VERIFIED'}})
        key=digest(tuple(identity[k] for k in rolling.SPEC)); queue=run/'retry_cache'/key
        queue.mkdir(exist_ok=True,parents=True); durable_json(queue/'attempt.json',state)
        return state
    def stage(a, run, label, rows, previous, preflight):
        nonlocal waiting
        keys=[r['source']['pdb_id'] for r in rows]; done=runs/f'{run.name}-{label}'; done.mkdir(exist_ok=True)
        if label=='b0003' and child_mode is None:
            waiting=True
            with subprocess.Popen([sys.executable,'-c','import time; time.sleep(2.6)']) as child:
                assert rolling.wait_stage(child,a.pump_retry)==0
            waiting=False
        if child_mode=='balance-kill' and label=='b0002':
            saved=read(run/'retry_balance'/'b0002.json'); marker=read(root/'balance_killed.json')
            assert saved['state']=='DRAINED' and saved['final']['queued'] < 2
            assert saved['started_epoch']==marker['started_epoch']
            durable_json(root/'balance_resumed.json',saved)
        pending=keys if (child_mode is None and label in ('b0001','b0002') or
                         child_mode in ('balance-kill','balance-timeout') and label=='b0001') else []
        if child_mode=='download' and label=='b0001':
            frozen=pin(PDB) | {'expected_bytes':len(BODY)}
            proof=done/'source_proofs'/f'{PDB}.json'; marker=root/'killed_download.json'
            archive=Archive(run/'archive_slots',4,2,local_test=True,lease_fds=a.owner_fds)
            if not marker.is_file():
                assert fetch(done,frozen,archive,'primary',1)[1]=='downloaded'
                # Leave a real archive transfer alive while killing its owner. Its
                # inherited lease blocks the new writer until the request exits.
                curl=subprocess.Popen(['curl','-q','-sS','--max-time','2','--output',str(done/'abandoned.bin'),url+'/body'],
                                      pass_fds=a.owner_fds)
                durable_json(marker,{'proof_sha':sha(proof),'http_pid':curl.pid})
                os.kill(os.getpid(),signal.SIGKILL)
            else:
                assert fetch(done,frozen,archive,'primary',1)[1]=='reused_verified'
                assert sha(proof)==read(marker)['proof_sha'] and len(read(proof)['attempts'])==1
        with closing(sqlite3.connect(done/'catalog.db')) as db, db:
            db.execute('CREATE TABLE IF NOT EXISTS obs(key TEXT PRIMARY KEY)')
            db.executemany('INSERT OR IGNORE INTO obs VALUES(?)',[(p,) for p in keys if p not in pending])
        durable_json(done/'phase_receipt.json',{'seed':previous.name})
        return {'catalog_run':done.name,'total':len(rows),'pending':pending,'settled':len(rows)-len(pending),'network_bytes':0}
    base_once=rolling.once
    def once(path,val):
        base_once(path,val)
        marker=root/'killed_receipt.json'
        if child_mode=='kill' and path==run/'receipts'/'b0001.json' and not marker.exists():
            # An actual child inherits the lease and outlives the killed controller.
            t=time.monotonic(); subprocess.Popen([sys.executable,'-c','import time; time.sleep(1.4)'],pass_fds=a.owner_fds)
            durable_json(marker,{'sha':sha(path),'head':read(run/'catalog_head.json')['catalog_run'],'killed_monotonic':t})
            os.kill(os.getpid(),signal.SIGKILL)
    if child_mode=='kill' and (root/'killed_receipt.json').is_file():
        durable_json(root/'replay_start.json',{'time':time.monotonic()})
    def save(path,value):
        durable_json(path,value)
        marker=root/'balance_killed.json'
        if child_mode=='balance-kill' and path==run/'retry_balance'/'b0002.json' and not marker.is_file():
            assert value['state']=='DRAINING' and value['next_partition']=='b0002'
            durable_json(marker,value)
            os.kill(os.getpid(),signal.SIGKILL)
    with patch.multiple(rolling,stage_batch=stage,retry_download=retry,once=once,write=save,
                        Archive=lambda *a,**kw: Archive(*a,**kw,local_test=True)):
        rolling.rolling(a)
    if child_mode is not None:
        return
    inside={p for p,in_stage,_ in wave if in_stage}
    assert inside&set(inputs[0]) and inside&set(inputs[1]) and len(inside)>=2, wave
    assert len(wave)==len({p for p,_,_ in wave})==4,wave
    receipts={p.stem:read(p) for p in (run/'receipts').glob('*.json')}
    assert sum(x.get('settled',0) for x in receipts.values())==5
    waves=[row['retry_selected'] for label,row in receipts.items() if '-r' in label]
    assert sorted(map(len,waves))==[2,2] and {p for row in waves for p in row}==set(inputs[0]+inputs[1])
    assert read(run/'retry_progress.json')['queued']==0
    return {'dispatched_while_one_main_stage_waited':len(inside),'cross_partition':True,'unique_retry_keys':len(wave),'ready_files_coalesced':True}


def real_supervised(root,url):
    def run_site(mode):
        where=root/mode; run=where/'runs'/'holo-m3-test-scheduler'; run.mkdir(parents=True)
        # The fixture owns its setup; the supervisor stores its state alongside it.
        args=[sys.executable,'-B',__file__,'--journal-child',str(where),url,mode]
        assert supervise(run,args,limit=2,backoff=.05,lease_wait=4)==0
        status=read(run/'supervisor_status.json'); assert status['state']=='NORMAL_EXIT'
        expected='PAUSED_RETRY_DRAIN_TIMEOUT' if mode=='balance-timeout' else 'PAUSED_OBSERVATION_LIMIT'
        assert read(run/'rolling_progress.json')['status']==expected
        return where,run,status
    where,run,status=run_site('kill'); evidence=read(where/'killed_receipt.json')
    assert status['restarts']==1 and evidence['head']=='holo-m3-test-seed'
    assert read(where/'replay_start.json')['time']-evidence['killed_monotonic']>=1.3
    assert sha(run/'receipts'/'b0001.json')==evidence['sha']
    with closing(sqlite3.connect(run.parent/(run.name+'-b0001')/'catalog.db')) as db:
        assert db.execute('SELECT COUNT(*) FROM obs').fetchone()[0]==2
    assert read(run/'catalog_head.json')['catalog_run']==run.name+'-b0003'
    _,stop,state=run_site('stop'); assert state['restarts']==0 and len(list((stop/'receipts').glob('b*.json')))==1
    download,trial,state=run_site('download')
    assert state['restarts']==1 and read(trial/'rolling_progress.json')['new_settled_entries']==5
    assert (download/'killed_download.json').is_file()
    balanced,balance_run,state=run_site('balance-kill')
    assert state['restarts']==1 and (balanced/'balance_resumed.json').is_file()
    assert read(balance_run/'rolling_progress.json')['new_settled_entries']==5
    assert read(balance_run/'retry_progress.json')['queued']==0
    _,held,state=run_site('balance-timeout')
    assert state['restarts']==0 and read(held/'rolling_progress.json')['partitions_completed']==1
    assert read(held/'retry_balance'/'b0002.json')['state']=='DRAIN_TIMEOUT'
    assert not (held/'receipts'/'b0002.json').is_file()
    return {'durable_balance_resumed_after_SIGKILL':True,'balance_timeout_stops_primary':True,
            'verified_download_proof_reused_after_real_exit':True,'signal':'SIGKILL','automatic_restarts':1,'waited_for_descendant_lease':True,
            'receipt_hash_preserved':True,'unique_records_after_restart':2,'intentional_stop_restarts':0}


def verify(output):
    with tempfile.TemporaryDirectory(prefix='aphcp-recovery-') as tmp, ThreadingHTTPServer(('127.0.0.1',0),Server) as srv:
        Thread(target=srv.serve_forever,daemon=True).start()
        root=Path(tmp); url=f'http://127.0.0.1:{srv.server_port}'
        evidence={'transport':transport(root,url),'journal':[journal(root,url,cut) for cut in ('sqlite-before-phase','phase-before-batch','batch-before-head')],
                  'continuous_retry':scheduler(root/'schedule',url),'supervised':real_supervised(root/'supervised',url)}
        srv.shutdown()
    durable_json(output,evidence)
    print(evidence,flush=True)


if __name__=='__main__':
    if len(sys.argv)==5 and sys.argv[1]=='--journal-child':
        scheduler(Path(sys.argv[2]),sys.argv[3],child_mode=sys.argv[4]); raise SystemExit(0)
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--output',required=True,type=Path)
    verify(p.parse_args().output)
