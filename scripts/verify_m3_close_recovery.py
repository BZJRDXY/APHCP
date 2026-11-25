"""Isolated controls for M3 terminal accounting, interrupted backups and supervisor handoff."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
from unittest.mock import patch

from aphcp_control.registry import backup, connect, recover_empty_registry, RegistryDiskReserve
from holo_archive_runtime import durable_json as write
from holo_m3_register import read, sha
from holo_m3_supervise import supervise
from verify_m3_terminal import audit_keys


def wait_for(path, timeout=8):
    until=time.monotonic()+timeout
    while not path.is_file() and time.monotonic()<until:
        time.sleep(.03)
    assert path.is_file(), path


def fails(action, error, marker=None):
    try:
        action()
    except error as exc:
        assert marker is None or marker in str(exc), str(exc)
    else:
        raise AssertionError(f'expected {error.__name__}: {marker}')


def backup_control(root):
    root.mkdir(); source=root/'source.db'; target=root/'next.db'; proof=root/'sealed.json'
    with connect(source) as db:
        db.execute('PRAGMA journal_mode=WAL')
        db.execute('PRAGMA user_version=31')
        db.execute('CREATE TABLE observations(k TEXT PRIMARY KEY)')
        db.execute("INSERT INTO observations VALUES('previously-sealed')")
    old=sha(source)
    crashed=subprocess.run([sys.executable,'-B',__file__,'--cut-backup',str(source),str(target)])
    assert crashed.returncode==-signal.SIGKILL and not target.exists()
    abandoned=list(root.glob('.next.db.backup-*')); assert len(abandoned)==1 and not abandoned[0].stat().st_size
    target.touch(); recover_empty_registry(target,(proof,)); assert not target.exists()
    with patch.dict(os.environ,{'APHCP_ROLLING_RESERVE_BYTES':str(10**18)}):
        fails(lambda:backup(source,target),RegistryDiskReserve)
    backup(source,target)
    with connect(target) as db:
        assert list(map(tuple,db.execute('SELECT k FROM observations')))==[('previously-sealed',)]
        db.execute("INSERT INTO observations VALUES('new-observation')")
    write(proof,{'seed_sha256':old,'result_sha256':sha(target)})
    fails(lambda:backup(source,target),FileExistsError)
    recover_empty_registry(target,(proof,)); assert sha(source)==old and sha(target)==read(proof)['result_sha256']
    damaged=root/'damaged.db'; damaged.touch()
    fails(lambda:recover_empty_registry(damaged,(proof,)),ValueError,'commit receipt')
    assert damaged.is_file()
    return {'killed_before_sqlite_write':True,'destination_never_exposed':True,'legacy_empty_recovered':True,
            'seed_preserved':True,'no_clobber':True,'sealed_empty_refused':True,'disk_protection':True}


def source(pdb):
    low=pdb.lower(); base='pdb_0000'+low
    return {'pdb_id':pdb,'revision':'1.2','file_type':'pdbx_mmcif_gzip',
            'source_url':f'https://files-versioned.rcsb.org/pdb_versioned/data/entries/{low[1:3]}/{base}/{base}_xyz_v1-2.cif.gz'}


def terminal_control(root):
    run=root/'runs'/'holo-m3-terminal-fixture'; seed=run.with_name('holo-m3-seed')
    (run/'partitions').mkdir(parents=True); (run/'receipts').mkdir(); seed.mkdir()
    write(run/'excluded.json',[]); write(run/'rolling_config.json',{'seed_run':str(seed)})
    partitions=[]; previous=seed; all_ids=[]
    # Exercise the absolute ordinal, not a `max(batches)` proxy. No archive requests.
    for n in range(1,450):
        bid=f'b{n:04d}'; pdb=f'1{n:03X}'; all_ids.append(pdb); pin=source(pdb)
        inp=run/'partitions'/f'{bid}.json'; write(inp,[{'source':pin}]); partitions.append({'id':bid,'entries':1,'sha256':sha(inp)})
        phase=run.with_name(f'{run.name}-{bid}-phase-001'); phase.mkdir()
        m2=run.with_name(f'holo-m2-terminal-fixture-{bid}'); m2.mkdir(); write(m2/'frozen.json',{'files':[pin]})
        write(phase/'phase_receipt.json',{'registry_run':str(phase),'m2_run':str(m2),'previous_catalog_run':str(previous),
                                          'selected':[pdb],'m2_pending':[],'m3_pending':[]})
        write(run/'receipts'/f'{bid}.json',{'sequence':n,'previous_catalog':previous.name,'catalog_run':phase.name,
                                           'm3_run':phase.name,'total':1,'settled':1,'pending':[],'archive_deferred':[]})
        previous=phase
    write(run/'rolling_plan.json',{'denominator':449,'partitions':partitions,'exclusions':{},'excluded_sha256':sha(run/'excluded.json')})
    write(run/'catalog_head.json',{'catalog_run':previous.name})
    retry={'running':0,'retryable_not_ready':0,'verified_waiting_for_catalog':0,'exhausted':0,'queued':0}
    write(run/'retry_progress.json',retry)
    expected=sha(run/'rolling_plan.json'); check=lambda:audit_keys(run,expected,449,3)
    result=check(); assert result['rolling_entries']==result['settled']==449 and result['absolute_stop']==449
    assert check()==result
    last=run/'receipts'/'b0449.json'; original=read(last); last.unlink()
    try: fails(check,ValueError,'terminal_unattempted_partitions')
    finally: write(last,original)
    phase=previous/'phase_receipt.json'; correct=read(phase); write(phase,correct | {'selected':[]})
    try: fails(check,ValueError,'terminal_per_entry_disposition_differs')
    finally: write(phase,correct)
    write(run/'retry_progress.json',retry | {'verified_waiting_for_catalog':1})
    try: fails(check,ValueError,'terminal_actionable_retries_remaining')
    finally: write(run/'retry_progress.json',retry)
    pre=run/'preflight'; pre.mkdir(); stray=pre/'b0450.json'; write(stray,{})
    try: fails(check,ValueError,'terminal_dispatched_beyond_stop')
    finally: stray.unlink()
    assert check()==result
    return {'last_ordinal':449,'early_exit_rejected':True,'missing_entry_rejected':True,
            'unmerged_rejected':True,'dispatch_beyond_stop_rejected':True,'journal_replay_equivalent':True}


def child(run, mode):
    with (run/'rolling.lock').open('a') as lease:
        fcntl.flock(lease,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if mode=='outer' and not (run/'writer_entered').is_file():
            (run/'writer_entered').touch(); time.sleep(1.3)
            (run/'writer_released').touch()
            return 75
        if mode=='outer':
            assert (run/'writer_released').is_file()
            return 75
        if mode=='stall':
            n=0
            while True:
                write(run/'rolling_progress.json',{'status':'ROLLING','heartbeat':n,'partitions_completed':0})
                n+=1; time.sleep(.03)
        raise AssertionError(mode)


def outer_control(root):
    run=root/'holo-m3-outer-control'; run.mkdir(parents=True)
    args=[sys.executable,'-B',__file__,'--outer',str(run)]
    with subprocess.Popen(args) as outer:
        wait_for(run/'writer_entered'); outer.kill(); assert outer.wait(timeout=4)==-signal.SIGKILL
    assert read(run/'supervisor_recovery.json')['inflight']
    start=time.monotonic(); assert subprocess.run(args,timeout=10).returncode==75
    assert time.monotonic()-start>.8
    evidence=read(run/'supervisor_recovery.json'); assert evidence['used']==1 and evidence['exhausted']
    assert [e['reason'] for e in evidence['events']]==['OUTER_SUPERVISOR_LOST','RECOVERABLE_CONTROLLER_EXIT']
    assert subprocess.run(args,timeout=4).returncode==75 and read(run/'supervisor_recovery.json')==evidence
    stalled=root/'stall'; stalled.mkdir()
    assert supervise(stalled,[sys.executable,'-B',__file__,'--child',str(stalled),'stall'],limit=1,
                     backoff=.03,lease_wait=1,stall_s=.25,poll_s=.05)==1
    assert read(stalled/'supervisor_status.json')['state']=='STALLED_NO_BUSINESS_PROGRESS'
    pause=root/'pause'; pause.mkdir(); (pause/'pause.request.json').touch()
    assert supervise(pause,['should-not-start'])==0 and read(pause/'supervisor_status.json')['state']=='USER_PAUSE'
    terminal=root/'early'; terminal.mkdir()
    cmd=[sys.executable,'-B','-c','import sys; from pathlib import Path; Path(sys.argv[1]).write_text(\'{"status":"PAUSED_OBSERVATION_LIMIT"}\')',
         str(terminal/'rolling_progress.json'),'--terminal-plan-sha256']
    assert supervise(terminal,cmd)==1 and read(terminal/'supervisor_status.json')['state']=='UNACCEPTED_CONTROLLER_EXIT'
    return {'outer_killed':True,'relaunch_waited_for_writer':True,'persistent_restart_budget':True,
            'heartbeat_without_business_progress_stopped':True,'user_pause_respected':True,'exit_zero_without_acceptance_refused':True}


if __name__=='__main__':
    if sys.argv[1:2]==['--cut-backup']:
        import aphcp_control.registry as registry
        real=registry.connect; source_db,target_db=map(Path,sys.argv[2:])
        def cut(path,**kw):
            if Path(path)!=source_db: os.kill(os.getpid(),signal.SIGKILL)
            return real(path,**kw)
        with patch.object(registry,'connect',cut): backup(source_db,target_db)
    elif sys.argv[1:2]==['--child']:
        raise SystemExit(child(Path(sys.argv[2]),sys.argv[3]))
    elif sys.argv[1:2]==['--outer']:
        run=Path(sys.argv[2]); cmd=[sys.executable,'-B',__file__,'--child',str(run),'outer']
        raise SystemExit(supervise(run,cmd,limit=1,backoff=.06,lease_wait=4))
    else:
        p=argparse.ArgumentParser(description=__doc__); p.add_argument('--output',type=Path,required=True); a=p.parse_args()
        with tempfile.TemporaryDirectory(prefix='aphcp-m3-close-isolated-') as tmp:
            root=Path(tmp); result={'backup':backup_control(root/'backup'),'terminal':terminal_control(root/'terminal'),
                                     'supervision':outer_control(root/'supervision')}
        write(a.output,result); print(json.dumps(result),flush=True)
