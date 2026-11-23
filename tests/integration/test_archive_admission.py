"""Exercise HEAD admission with real curl, POSIX file leases and a local archive."""
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import os
from pathlib import Path
import shutil
import sys
import tempfile
from threading import Event, Lock, Thread
import time
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from holo_archive_runtime import Archive


@unittest.skipUnless(os.name == 'posix' and shutil.which('curl'), 'requires POSIX leases and curl')
class ArchiveAdmission(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.archive = Archive(self.root / 'locks', 8, 4, local_test=True)
        self.release_get = Event()
        self.changed = Event()
        self.lock = Lock()
        self.active = 0
        self.peak = 0
        self.hits = {}
        test = self

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *_args):
                pass

            def reply(self, has_body):
                with test.lock:
                    test.active += 1
                    test.peak = max(test.active, test.peak)
                    test.hits[self.path] = test.hits.get(self.path, 0) + 1
                    test.changed.set()
                try:
                    if has_body and self.path.startswith('/hold'):
                        test.release_get.wait(10)
                    elif self.path == '/real-timeout':
                        time.sleep(.45)
                    elif not has_body:
                        time.sleep(.12)
                    self.send_response(206 if self.headers.get('Range') else 200)
                    self.send_header('Content-Length', '5')
                    self.end_headers()
                    if has_body:
                        self.wfile.write(b'hello')
                finally:
                    with test.lock:
                        test.active -= 1
                        test.changed.set()

            def do_HEAD(self):
                self.reply(False)

            def do_GET(self):
                self.reply(True)

        class Server(ThreadingHTTPServer):
            request_queue_size = 32

        self.server = Server(('127.0.0.1', 0), Handler)
        self.thread = Thread(target=self.server.serve_forever, kwargs={'poll_interval': .02}, daemon=True)
        self.thread.start()
        self.addCleanup(self.close_server)

    def close_server(self):
        self.release_get.set()
        self.server.shutdown()
        self.thread.join(2)
        self.server.server_close()

    def url(self, path):
        return f'http://127.0.0.1:{self.server.server_port}/{path}'

    def get(self, index, kind):
        trace = {}
        headers = ('Range: bytes=0-4',) if kind == 'range' else ()
        path = self.root / str(index)
        began = time.monotonic()
        self.archive.request(self.url(f'hold-{index}'), kind, path, 8, trace,
                             headers=headers, deadline=began+8)
        self.assertEqual(path.read_bytes(), b'hello')
        return trace

    def head(self, path, request=.32, admission=3):
        return self.archive.head({'source_url':self.url(path)}, seconds=request, queue_seconds=admission)

    def wait_holders(self, count):
        until = time.monotonic() + 4
        while time.monotonic() < until:
            with self.lock:
                if len(self.hits) >= count:
                    return
                self.changed.clear()
            self.changed.wait(.1)
        self.fail(f'Only {len(self.hits)} requests in HTTP; expected {count}')

    def test_head_gets_full_http_budget_under_mixed_load(self):
        with ThreadPoolExecutor(max_workers=16) as pool:
            main = [pool.submit(self.get, k, 'get') for k in range(6)]
            self.wait_holders(6)
            retry = pool.submit(self.get, 6, 'range')
            self.wait_holders(7)
            try:
                # Under this saturation only the protected HEAD lane is free.
                # The old combined deadline waits out some requests locally.
                def legacy(i):
                    trace = {}
                    try:
                        self.archive.request(self.url(f'legacy-{i}'), 'head', '/dev/null', .32,
                                             trace, deadline=time.monotonic()+.32)
                    except TimeoutError:
                        pass
                    return trace
                old = list(pool.map(legacy, range(6)))
                self.assertTrue(any(x.get('phase') == 'QUEUE' for x in old), old)
                # With an expired client, the HTTP handler may continue to
                # serve for a moment; measure new requests after that drains.
                until = time.monotonic() + 2
                while True:
                    with self.lock:
                        if self.active == 7:
                            self.peak = self.active
                            break
                        self.changed.clear()
                    if time.monotonic() >= until:
                        self.fail('Legacy HEAD did not drain')
                    self.changed.wait(.05)
                new = list(pool.map(lambda i:self.head(f'new-{i}'), range(6)))
                self.assertTrue(all(x['availability']=='VERIFIED' for x in new), new)
                self.assertTrue(any(x['request']['queued_s']>.32 for x in new), new)
                self.assertTrue(all(x['request']['slot']==0 and x['request']['http_started'] for x in new))
                self.assertEqual(self.peak, 8)
            finally:
                self.release_get.set()
            self.assertEqual({f.result()['slot'] for f in main}, set(range(2,8)))
            self.assertEqual(retry.result()['slot'], 1)

    def test_local_admission_is_not_remote_timeout(self):
        with self.archive.slot('head'), self.archive.slot('head'), self.archive.slot('head'):
            result = self.head('never-sent', admission=.06)
        self.assertEqual(result.get('failure_domain'), 'LOCAL_ADMISSION')
        self.assertFalse(result['request']['http_started'])
        self.assertNotIn('/never-sent', self.hits)
        result = self.head('real-timeout', request=.16)
        self.assertEqual(result.get('failure_domain'), 'ARCHIVE_RESPONSE')
        self.assertTrue(result['request']['http_started'])
        self.assertEqual(result['request']['exit_code'], 28)
        self.assertEqual(self.hits['/real-timeout'], 1)


if __name__ == '__main__':
    unittest.main()
