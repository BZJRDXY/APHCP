"""Restart a rolling controller after bounded, recognized interruption."""
from __future__ import annotations

import argparse
import fcntl
import signal
import subprocess
import sys
import time
from pathlib import Path

from holo_archive_runtime import durable_json
from holo_m3_rolling import TEMPORARY_INTERRUPTION
from holo_m3_register import utc_now


def supervise(run, command, limit=3, backoff=5, lease_wait=600):
    snapshot=run/'supervisor_status.json'
    def record(state, **evidence):
        durable_json(snapshot, {'state':state,'updated_utc':utc_now(),'restarts':restarts,**evidence})
    restarts=0
    with (run/'supervisor.lock').open('a') as supervisor:
        fcntl.flock(supervisor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        while True:
            if (run/'pause.request.json').is_file():
                record('USER_PAUSE'); return 0
            # Descendants inherit the writer lease. Never wait without a bound or kill them here.
            until=time.monotonic()+lease_wait
            with (run/'rolling.lock').open('a') as lease:
                while True:
                    if (run/'pause.request.json').is_file():
                        record('USER_PAUSE'); return 0
                    try:
                        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        break
                    except BlockingIOError:
                        if time.monotonic()>=until:
                            record('WRITER_LEASE_STILL_HELD'); return 1
                        time.sleep(min(.25,max(.01,until-time.monotonic())))
            with subprocess.Popen(command) as child:
                record('RUNNING',controller_pid=child.pid)
                status=child.wait()
            if status==0:
                record('NORMAL_EXIT'); return 0
            if (run/'pause.request.json').is_file():
                record('USER_PAUSE',controller_exit=status); return 0
            if status not in (TEMPORARY_INTERRUPTION,-signal.SIGKILL):
                record('CONTROLLER_STOP',controller_exit=status)
                return status if status>0 else 128-status
            if restarts>=limit:
                record('RESTART_LIMIT',controller_exit=status); return TEMPORARY_INTERRUPTION
            delay=backoff*2**restarts
            record('RECOVERY_BACKOFF',controller_exit=status,delay_s=delay)
            until=time.monotonic()+delay
            while (remaining:=until-time.monotonic())>0:
                if (run/'pause.request.json').is_file():
                    record('USER_PAUSE',controller_exit=status); return 0
                time.sleep(min(remaining,.25))
            restarts+=1


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--restart-limit',type=int,default=3)
    p.add_argument('--backoff-s',type=float,default=5)
    p.add_argument('--lease-wait-s',type=float,default=600)
    options,args=p.parse_known_args()
    if args[:1]==['--']: args=args[1:]
    if ('--run' not in args or '--max-main-batches' in args or '--resume-paused' in args or
            not any(x in args for x in ('--retry-only','--stop-after-partition')) or
            not 0<=options.restart_limit<=5 or not 0<options.backoff_s<=60 or
            not 0<options.lease_wait_s<=3600):
        p.error('supply --run and a persistent --retry-only or --stop-after-partition; unsafe restart options refused')
    run=Path(args[args.index('--run')+1]).resolve(strict=True)
    cmd=[sys.executable,'-B',str(Path(__file__).with_name('holo_m3_rolling.py')),*args]
    raise SystemExit(supervise(run,cmd,options.restart_limit,options.backoff_s,options.lease_wait_s))


if __name__=='__main__':
    main()
