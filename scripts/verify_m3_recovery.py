"""Run isolated Linux controls for the rolling downloader and its commit journal."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
import fcntl
import gzip
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import multiprocessing as mp
import os
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
    return {'timeouts':deadlines,'max_parallel_across_processes':max_parallel,'max_ancillary':max_ancillary,
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


def verify(output):
    with tempfile.TemporaryDirectory(prefix='aphcp-recovery-') as tmp, ThreadingHTTPServer(('127.0.0.1',0),Server) as srv:
        Thread(target=srv.serve_forever,daemon=True).start()
        root=Path(tmp); url=f'http://127.0.0.1:{srv.server_port}'
        evidence={'transport':transport(root,url),'journal':[journal(root,url,cut) for cut in ('sqlite-before-phase','phase-before-batch','batch-before-head')]}
        srv.shutdown()
    durable_json(output,evidence)
    print(evidence,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--output',required=True,type=Path)
    verify(p.parse_args().output)
