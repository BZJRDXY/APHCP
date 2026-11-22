import hashlib
import json
import os
import re
import socket
import time
from contextlib import contextmanager
from pathlib import Path
from uuid import uuid4

from .registry import connect, database, migrate, utc_now


def digest(data):
    return hashlib.sha256(data).hexdigest()


def atomic_write(path, data):
    temp = path.with_suffix(path.suffix + ".partial")
    temp.write_bytes(data)
    temp.replace(path)


@contextmanager
def coordinator(root):
    lock_path = root / "registry" / "coordinator.lock"
    with lock_path.open("ab") as stream:
        if os.name == "nt":
            import msvcrt
            stream.seek(0)
            msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def prepare(db, items, task_key, commit, batch, task):
    db.execute("INSERT OR IGNORE INTO batch(id,module_id,status,candidate_count) VALUES "
               "('M0-EMPTY','holo','PLANNED',0)")
    db.execute("""
        INSERT OR IGNORE INTO batch(id,module_id,status,source_release,input_manifest,config_version,
                                    git_commit,candidate_count,next_step)
        VALUES (?,'holo','PROCESSING','synthetic-m0','tests/fixtures/m0_work_items.json',
                'm0-fixture-v1',?,?,'完成三条合成测试工作项')
    """, (batch, commit, len(items)))
    db.execute("UPDATE batch SET status='PROCESSING' WHERE id=?", (batch,))
    db.execute("INSERT OR IGNORE INTO task VALUES (?,?,'fixture-check',?,'RUNNING')", (task, batch, task_key))
    db.execute("UPDATE task SET status='RUNNING' WHERE id=?", (task,))
    for item in items:
        pdb = item["pdb_id"]
        version = f"{pdb}-fixture-v1"
        assembly = f"{version}-assembly-1"
        db.execute("INSERT OR IGNORE INTO pdb_entry(id,title) VALUES (?,'合成测试记录')", (pdb,))
        db.execute("INSERT OR IGNORE INTO entry_version VALUES (?,?, 'test','synthetic-m0',?,'1')",
                   (version, pdb, pdb))
        db.execute("INSERT OR IGNORE INTO work_item(id,entry_version_id,kind,label) "
                   "VALUES (?,?,'assembly','测试组装体 1')", (assembly, version))
        db.execute("INSERT OR IGNORE INTO work_item(id,entry_version_id,parent_id,kind,label) "
                   "VALUES (?,?,?,'system',?)", (item["id"], version, assembly, item["id"]))
        chain = f"{assembly}-chain-A"
        ligand = f"{item['id']}-ligand"
        db.execute("INSERT OR IGNORE INTO work_item(id,entry_version_id,parent_id,kind,label) "
                   "VALUES (?,?,?,'chain','测试蛋白链 A')", (chain, version, assembly))
        db.execute("INSERT OR IGNORE INTO work_item(id,entry_version_id,parent_id,kind,label) "
                   "VALUES (?,?,?,'ligand','测试配体')", (ligand, version, assembly))
        db.execute("INSERT OR IGNORE INTO system_member VALUES (?,?, 'protein')", (item["id"], chain))
        db.execute("INSERT OR IGNORE INTO system_member VALUES (?,?, 'ligand')", (item["id"], ligand))


def recover(db, task):
    abandoned = db.execute("SELECT * FROM run WHERE task_id=? AND status='RUNNING'", (task,)).fetchall()
    for record in abandoned:
        run_id = record["id"]
        items = [dict(row) for row in db.execute("SELECT * FROM run_item WHERE run_id=?", (run_id,))]
        recovered_at = utc_now()
        summary = {"fixture_only": True, "run_id": run_id, "status": "PARTIAL",
                   "recovered_at": recovered_at, "known_items": items,
                   "reason": "协调器中断，恢复时收口"}
        atomic_write(Path(record["report_path"]), json.dumps(summary, ensure_ascii=False, indent=2).encode())
        db.execute("""
            UPDATE run SET status='PARTIAL',end_time_utc=?,error_summary='协调器中断，恢复时收口',
                success_records=(SELECT COUNT(*) FROM run_item WHERE run_id=? AND status IN ('PASS','REUSED'))
            WHERE id=?
        """, (recovered_at, run_id, run_id))


def trial(root, fixture, commit, env_lock, interrupt_after=None):
    root = Path(root).resolve()
    if root.name != "m0-test":
        raise ValueError("合成测试仅允许写入名为 m0-test 的独立数据根")
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("需要明确的 40 位 Git Commit")
    source = fixture.read_bytes()
    items = json.loads(source)
    if not isinstance(items, list) or len(items) != 3 or any(
        not re.fullmatch(r"[A-Z0-9]{4}", item["pdb_id"]) or
        not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,50}", item["id"]) for item in items
    ) or len({item["id"] for item in items}) != len(items):
        raise ValueError("夹具需要三条不同的合成测试工作项")

    input_hash = digest(source)
    config_hash = digest(b"m0-fixture-v1")
    env_hash = digest(env_lock.read_bytes())
    task_key = digest("|".join(("holo", input_hash, "synthetic-m0", commit, config_hash, env_hash)).encode())
    batch = f"M0-RECOVERY-{task_key[:12]}"
    task = f"{batch}-verify"
    run = f"M0-RUN-{uuid4().hex[:12]}"
    db_path = database(root)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    with coordinator(root):
        migrate(db_path)
        out_dir = root / "artifacts" / task_key
        out_dir.mkdir(parents=True, exist_ok=True)
        report = root / "runs" / run / "report.json"
        report.parent.mkdir(parents=True)
        start = time.monotonic()
        cpu_start = time.process_time()
        with connect(db_path) as db:
            prepare(db, items, task_key, commit, batch, task)
            recover(db, task)
            db.execute("""
                INSERT INTO run(id,batch_id,task_id,status,pid,git_commit,code_version,config_version,
                  config_hash,source_release,input_manifest_hash,environment_lock_hash,host,start_time_utc,
                  input_records,report_path)
                VALUES (?,?,?,'RUNNING',?,?,?,'m0-fixture-v1',?,'synthetic-m0',?,?,?,?,?,?)
            """, (run, batch, task, os.getpid(), commit, commit, config_hash,
                  input_hash, env_hash, socket.gethostname(), utc_now(), len(items), str(report)))
        atomic_write(report, json.dumps({"fixture_only": True, "run_id": run, "status": "RUNNING"}).encode())
        results = []
        bytes_read = len(source) + env_lock.stat().st_size
        bytes_written = 0
        try:
            for index, item in enumerate(items, start=1):
                output = out_dir / f"{item['id']}.txt"
                expected = f"m0-fixture-v1:{item['pdb_id']}:{item['id']}\n".encode()
                prior = output.read_bytes() if output.is_file() else None
                bytes_read += len(prior or b"")
                if prior == expected:
                    status = "REUSED"
                else:
                    atomic_write(output, expected)
                    bytes_written += len(expected)
                    status = "PASS"
                result = {"id": item["id"], "pdb_id": item["pdb_id"], "status": status,
                          "sha256": digest(expected), "output_uri": str(output)}
                results.append(result)
                with connect(db_path) as db:
                    db.execute("""
                        INSERT INTO run_item(run_id,object_type,object_id,pdb_id,stage,status,output_uri,updated_at)
                        VALUES (?,'system',?,?,'FIXTURE',?,?,?)
                    """, (run, item["id"], item["pdb_id"], status, str(output), utc_now()))
                if index == interrupt_after:
                    os._exit(75)
        except Exception as exc:
            with connect(db_path) as db:
                db.execute("UPDATE run SET status='FAILED',end_time_utc=?,error_summary=? WHERE id=?",
                           (utc_now(), str(exc), run))
                db.execute("UPDATE task SET status='FAILED' WHERE id=?", (task,))
                db.execute("UPDATE batch SET status='FAILED' WHERE id=?", (batch,))
            raise

        manifest = json.dumps(results, ensure_ascii=False, sort_keys=True).encode()
        output_hash = digest(manifest)
        final = {"run_id": run, "batch_id": batch, "fixture_only": True, "status": "SUCCEEDED", "git_commit": commit,
                 "config_hash": config_hash, "environment_lock_hash": env_hash,
                 "input_manifest_hash": input_hash, "output_manifest_hash": output_hash, "items": results}
        manifest_path = report.parent / "output_manifest.json"
        atomic_write(manifest_path, manifest)
        atomic_write(report, json.dumps(final, ensure_ascii=False, indent=2).encode())
        bytes_written += manifest_path.stat().st_size + report.stat().st_size
        with connect(db_path) as db:
            db.execute("""
                UPDATE run SET status='SUCCEEDED',end_time_utc=?,elapsed_seconds=?, success_records=?,
                  failed_records=0,uncertain_records=0,bytes_read=?,bytes_written=?,
                  downloaded_bytes=0,cpu_seconds=?,output_manifest_hash=? WHERE id=?
            """, (utc_now(), time.monotonic() - start, len(results), bytes_read, bytes_written,
                  time.process_time() - cpu_start, output_hash, run))
            db.execute("UPDATE task SET status='SUCCEEDED' WHERE id=?", (task,))
            db.execute("UPDATE batch SET status='ACCEPTED',next_step='完成合成测试' WHERE id=?", (batch,))
        return final
