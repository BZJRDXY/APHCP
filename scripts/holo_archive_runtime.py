"""Archive request limits and crash-safe publication for the rolling M3 runner."""
from contextlib import contextmanager
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
from urllib.parse import urlsplit


def sync_path(path):
    flags = os.O_RDONLY | (getattr(os, 'O_DIRECTORY', 0) if path.is_dir() else 0)
    fd = os.open(path, flags)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def durable_json(path, value):
    """Publish fully written JSON; on Linux persist the directory rename too."""
    with tempfile.NamedTemporaryFile(mode='w', encoding='utf8', newline='\n',
                                     prefix=f'.{path.name}.', dir=path.parent, delete=False) as f:
        temp = Path(f.name)
        try:
            json.dump(value, f, ensure_ascii=False, indent=2, sort_keys=True)
            f.write('\n'); f.flush(); os.fsync(f.fileno())
        except BaseException:
            f.close(); temp.unlink(missing_ok=True)
            raise
    try:
        temp.replace(path)
        if os.name == 'posix':
            sync_path(path.parent)
    finally:
        temp.unlink(missing_ok=True)


class Archive:
    """A shared process-safe pool. HEAD/retry/range cannot occupy the main reserve."""
    def __init__(self, folder, limit, ancillary, *, local_test=False, lease_fds=()):
        if not 1 <= ancillary < limit:
            raise ValueError('archive slots must leave a main GET reserve')
        self.folder = Path(folder); self.folder.mkdir(parents=True, exist_ok=True)
        self.limit = limit; self.ancillary = ancillary; self.local_test = local_test
        self.lease_fds = tuple(lease_fds)

    @contextmanager
    def slot(self, kind, deadline=None):
        import fcntl
        if kind not in {'get', 'head', 'retry', 'range'}:
            raise ValueError('unknown archive request category')
        slots = range(self.limit-1, -1, -1) if kind == 'get' else range(self.ancillary)
        while deadline is None or time.monotonic() < deadline:
            for number in slots:
                with (self.folder/f'{number:03d}.lock').open('ab') as f:
                    try:
                        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    except BlockingIOError:
                        continue
                    try:
                        yield number, f.fileno()
                    finally:
                        fcntl.flock(f, fcntl.LOCK_UN)
                    return
            time.sleep(.012)
        raise TimeoutError('Archive slot wait exhausted file budget')

    def request(self, url, kind, output, seconds, trace, *, headers=(), deadline=None):
        parsed = urlsplit(url)
        protocol = '=https'
        if self.local_test and parsed.scheme == 'http' and parsed.hostname == '127.0.0.1':
            protocol = '=http'
        elif (parsed.scheme != 'https' or parsed.netloc != 'files-versioned.rcsb.org'
              or not parsed.path.startswith('/pdb_versioned/data/entries/') or parsed.query or parsed.fragment):
            raise ValueError('archive request outside the fixed RCSB archive')
        output = Path(output)
        queued = time.monotonic()
        with self.slot(kind, deadline) as (number, slot_fd):
            started = time.monotonic()
            remaining = min(seconds, deadline-started) if deadline is not None else seconds
            if remaining < .005:
                raise TimeoutError('Archive file budget expired before request')
            trace.update(kind=kind, slot=number, queued_s=round(started-queued, 4), started_monotonic=started)
            fmt = '%{http_code}\t%{size_download}\t%{time_starttransfer}\t%{time_total}'
            with tempfile.TemporaryFile() as header:
                cmd = ['curl', '-q', '--silent', '--show-error', '--globoff', '--proto', protocol,
                       '--retry', '0', '--max-time', f'{remaining:.3f}',
                       '--connect-timeout', f'{min(10, remaining):.3f}', '--output', str(output),
                       '--dump-header', '/dev/fd/'+str(header.fileno()), '--write-out', fmt,
                       '--header', 'Accept-Encoding: identity', '--user-agent', 'APHCP-fixed-archive/2']
                if kind == 'head':
                    cmd += ['--head']
                for h in headers:
                    cmd += ['--header', h]
                try:
                    result = subprocess.run(cmd+['--url',url], capture_output=True,
                                            pass_fds=(header.fileno(), slot_fd, *self.lease_fds), timeout=remaining+.20)
                except subprocess.TimeoutExpired as exc:
                    trace.update(exit_code=28, body_bytes=output.stat().st_size if output.is_file() and kind!='head' else 0)
                    raise TimeoutError('Archive transport watchdog expired') from exc
                finally:
                    trace.update(finished_monotonic=time.monotonic(), elapsed_s=round(time.monotonic()-started,4))
                trace['exit_code'] = result.returncode
                fields = result.stdout.decode('ascii', errors='replace').strip().split('\t')
                if len(fields) == 4:
                    code,size,ttfb,wall = fields
                    trace.update(http_status=int(code), body_bytes=int(float(size)), ttfb_s=float(ttfb), transfer_s=float(wall))
                if result.returncode == 28:
                    raise TimeoutError('Archive transfer exhausted its wall-clock budget')
                if result.returncode:
                    raise OSError(f'Archive transfer failed (curl {result.returncode})')
                header.seek(0)
                raw = header.read().decode('iso-8859-1').replace('\r\n','\n').strip()
                lines = raw.split('\n\n')[-1].splitlines()[1:]
                return {k.lower():v.strip() for line in lines if ':' in line for k,v in [line.split(':',1)]}

    def head(self, pin, seconds=12):
        from datetime import datetime, timezone
        trace = {}; out = pin | {'checked_utc':datetime.now(timezone.utc).isoformat()}
        try:
            fields = self.request(pin['source_url'], 'head', '/dev/null', seconds, trace)
            n = int(fields.get('content-length', 0)); status = trace.get('http_status')
            if status==200 and n>0:
                return out | {'status':status,'content_length':n,'final_url':pin['source_url'],
                              'availability':'VERIFIED','request':trace}
        except (OSError, ValueError) as e:
            out['error']=f'{type(e).__name__}: {e}'
        return out | {'status':trace.get('http_status'),'availability':'DEFERRED_BEFORE_FREEZE','request':trace}
