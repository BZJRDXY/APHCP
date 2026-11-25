"""Supervise a rolling controller with a recovery budget that survives process loss."""
from __future__ import annotations

import argparse
import fcntl
import signal
import subprocess
import sys
import time
from pathlib import Path

from holo_archive_runtime import durable_json
from holo_m3_rolling import DISK_RESERVE_EXIT, TEMPORARY_INTERRUPTION
from holo_m3_register import read, require, utc_now

PLANNED = {'PAUSED_OBSERVATION_LIMIT', 'PAUSED_RETRY_ONLY', 'STAGED_PENDING_RETRY'}


def state_file(path):
    return read(path) if path.is_file() else {}


def business_progress(run):
    """Use actual receipts/checkpoints; log timestamps and retry heartbeats do not count."""
    progress = state_file(run/'rolling_progress.json')
    retry = state_file(run/'retry_progress.json')
    dispatch = state_file(run/'active_dispatch.json') or {}
    label = dispatch.get('label') or progress.get('active_partition')
    checkpoints = []
    if label:
        m2 = run.parent/f'holo-m2-{run.name.removeprefix("holo-m3-")}-{label}'
        m3 = run.parent/f'{run.name}-{label}'
        for folder in (m2/'source_proofs', m2/'entry_results', m2/'system_results'):
            checkpoints.append(sum(1 for _ in folder.glob('*.json')))
        phase = m3.with_name(m3.name+'-phase-001')
        checkpoints.extend((phase/name).is_file() for name in ('candidate_evidence.json','review_results.json',
                                                               'catalog_receipt.json','phase_receipt.json'))
    return (progress.get('partitions_completed'), progress.get('new_settled_entries'), progress.get('catalog_run'),
            retry.get('verified_waiting_for_catalog'), retry.get('settled_by_retry'), label, *checkpoints)


def supervise(run, command, limit=3, backoff=5, lease_wait=600, stall_s=3600, poll_s=10):
    snapshot=run/'supervisor_status.json'; budget=run/'supervisor_recovery.json'
    terminal=any(x=='--terminal-plan-sha256' or x.startswith('--terminal-plan-sha256=') for x in command)
    pause=run/'pause.request.json'; protective=run/'protective_stop.json'
    with (run/'supervisor.lock').open('a') as supervisor:
        fcntl.flock(supervisor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        recovery=state_file(budget) or {'used':0,'limit':limit,'inflight':False,'events':[]}
        require(recovery['limit']==limit, 'supervisor_recovery_budget_changed')

        def save(**changes):
            recovery.update(changes)
            durable_json(budget,recovery)

        def record(state, **evidence):
            durable_json(snapshot, {'state':state,'updated_utc':utc_now(),'restarts':recovery['used'],**evidence})

        def reserve(reason, **evidence):
            used=recovery['used']; event={'reason':reason,'utc':utc_now(),**evidence}
            if used>=limit:
                save(inflight=False,not_before=0,exhausted=True,
                     events=[*recovery['events'],event | {'accepted':False}])
                record('RESTART_LIMIT',**event)
                return False
            delay=backoff*2**used
            save(used=used+1,inflight=False,not_before=time.time()+delay,
                 events=[*recovery['events'],event | {'accepted':True,'delay_s':delay}])
            record('RECOVERY_BACKOFF',**event,delay_s=delay)
            return True

        while True:
            if pause.is_file():
                record('USER_PAUSE'); return 0
            if recovery.get('stall') or recovery.get('protective'):
                record(recovery.get('protective') or 'STALLED_NO_BUSINESS_PROGRESS'); return 1
            if recovery.get('exhausted'):
                record('RESTART_LIMIT'); return TEMPORARY_INTERRUPTION
            # A lost supervisor can leave a real writer behind. Wait for its inherited lease.
            until=time.monotonic()+lease_wait
            with (run/'rolling.lock').open('a') as lease:
                while True:
                    if pause.is_file():
                        record('USER_PAUSE'); return 0
                    try:
                        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        break
                    except BlockingIOError:
                        if time.monotonic()>=until:
                            record('WRITER_LEASE_STILL_HELD'); return 1
                        time.sleep(min(.25,max(.01,until-time.monotonic())))
            if recovery['inflight']:
                previous_status=state_file(run/'rolling_progress.json').get('status')
                if (reason:=state_file(protective).get('status')):
                    save(inflight=False,protective=reason); record(reason); return DISK_RESERVE_EXIT
                if (terminal and previous_status=='STAGED_COMPLETE' and (run/'terminal_receipt.json').is_file() and
                    state_file(run/'terminal_receipt.json').get('accepted')):
                    save(inflight=False,not_before=0)
                    record('COMPLETED',controller_status=previous_status); return 0
                if not reserve('OUTER_SUPERVISOR_LOST', controller_status=previous_status):
                    return TEMPORARY_INTERRUPTION
            while (remaining:=recovery.get('not_before',0)-time.time())>0:
                if pause.is_file():
                    record('USER_PAUSE'); return 0
                time.sleep(min(remaining,.25))
            if protective.is_file():
                record(state_file(protective)['status']); return DISK_RESERVE_EXIT
            save(inflight=True,not_before=0)
            child=subprocess.Popen(command)
            record('RUNNING',controller_pid=child.pid)
            last=business_progress(run); since=time.monotonic()
            while True:
                try:
                    status=child.wait(timeout=poll_s)
                    break
                except subprocess.TimeoutExpired:
                    actual=business_progress(run)
                    if actual!=last:
                        last=actual; since=time.monotonic()
                    if time.monotonic()-since>=stall_s:
                        save(inflight=False,stall=True)
                        # Park this Run. Do not kill descendants or ever take over a held writer lease.
                        child.terminate()
                        record('STALLED_NO_BUSINESS_PROGRESS',controller_pid=child.pid,
                               unchanged_s=round(time.monotonic()-since),business_progress=last)
                        return 1
            reason=state_file(run/'rolling_progress.json').get('status')
            if pause.is_file():
                save(inflight=False)
                record('USER_PAUSE',controller_exit=status); return 0
            if status in (TEMPORARY_INTERRUPTION,-signal.SIGKILL):
                if not reserve('RECOVERABLE_CONTROLLER_EXIT',controller_exit=status,controller_status=reason):
                    return TEMPORARY_INTERRUPTION
                continue
            save(inflight=False,not_before=0)
            if status==0 and reason=='STAGED_COMPLETE':
                if not terminal or state_file(run/'terminal_receipt.json').get('accepted'):
                    record('COMPLETED',controller_status=reason); return 0
            if status==0 and reason in PLANNED and not terminal:
                record('PLANNED_STOP',controller_status=reason); return 0
            if status==DISK_RESERVE_EXIT:
                record('PAUSED_DISK_RESERVE',controller_exit=status); return status
            if reason=='PAUSED_UNSEALED_PARTITION' or reason=='VERIFYING_TERMINAL':
                save(protective='CHECKPOINT_OR_INTEGRITY_STOP')
                record('CHECKPOINT_OR_INTEGRITY_STOP',controller_exit=status,controller_status=reason)
            else:
                save(protective='UNACCEPTED_CONTROLLER_EXIT')
                record('UNACCEPTED_CONTROLLER_EXIT',controller_exit=status,controller_status=reason)
            return status if status>0 else (128-status if status<0 else 1)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--restart-limit',type=int,default=3)
    p.add_argument('--backoff-s',type=float,default=5)
    p.add_argument('--lease-wait-s',type=float,default=600)
    p.add_argument('--stall-s',type=float,default=3600)
    options,args=p.parse_known_args()
    if args[:1]==['--']: args=args[1:]
    if ('--run' not in args or '--max-main-batches' in args or '--resume-paused' in args or
            not any(x in args for x in ('--retry-only','--stop-after-partition')) or
            not 0<=options.restart_limit<=5 or not 0<options.backoff_s<=60 or
            not 0<options.lease_wait_s<=3600 or not 900<=options.stall_s<=86400):
        p.error('supply --run and a persistent --retry-only or --stop-after-partition; unsafe restart options refused')
    run=Path(args[args.index('--run')+1]).resolve(strict=True)
    cmd=[sys.executable,'-B',str(Path(__file__).with_name('holo_m3_rolling.py')),*args]
    raise SystemExit(supervise(run,cmd,options.restart_limit,options.backoff_s,options.lease_wait_s,options.stall_s))


if __name__=='__main__':
    main()
