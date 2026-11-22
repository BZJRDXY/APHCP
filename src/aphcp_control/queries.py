import json
import os
import re
import shutil
import time
from collections import deque
from datetime import datetime, timezone
from pathlib import Path

from .registry import utc_now


LIST_LIMIT = 25
QUERY_IDS = {"batch", "run", "pdb"}
QUERIES = {"project", "modules", "batches", "batch", "runs", "run", "pdb", "storage", "deployment",
           "governance"}


def query_args(view, value=None):
    if view not in QUERIES or bool(value) != (view in QUERY_IDS):
        raise ValueError("无效的查询类型或标识符参数")
    if value:
        pattern = r"(?:[A-Za-z0-9]{4}|pdb_[A-Za-z0-9]{8})" if view == "pdb" else r"[A-Za-z0-9][A-Za-z0-9_.-]{0,95}"
        if not re.fullmatch(pattern, value, flags=re.ASCII | re.IGNORECASE):
            raise ValueError("标识符格式不正确")
        return view, value.upper() if view == "pdb" else value
    return (view,)


def rows(db, sql, *parameters):
    return [dict(row) for row in db.execute(sql, parameters)]


def one(db, sql, *parameters):
    row = db.execute(sql, parameters).fetchone()
    return dict(row) if row else None


def limited(db, sql, *parameters):
    return rows(db, sql + " LIMIT ?", *parameters, LIST_LIMIT)


def deployment(db, code_root):
    marker = code_root / "current" / "deployment.json"
    deployed = json.loads(marker.read_text(encoding="utf-8")) if marker.is_file() else None
    history = []
    event_file = code_root / "deploy" / "history.jsonl"
    if event_file.is_file():
        with event_file.open(encoding="utf-8") as records:
            history = [json.loads(line) for line in deque(records, maxlen=LIST_LIMIT)]
    if history and deployed and history[-1]["commit"] == deployed["commit"]:
        deployed = history[-1]
    return {"current": deployed, "history": history,
            "run_commits": limited(db, "SELECT id,git_commit,status FROM run ORDER BY start_time_utc DESC")}


def project(db, data_root, code_root):
    status = one(db, "SELECT * FROM project")
    return {
        "project": status,
        "modules": rows(db, "SELECT * FROM module ORDER BY id"),
        "active_batches": limited(db, """
            SELECT * FROM batch WHERE status IN ('INVENTORY_RUNNING','DOWNLOADING','PROCESSING')
            ORDER BY id
        """),
        "active_tasks": limited(db, """
            SELECT t.*,r.id AS run_id,r.start_time_utc FROM task t
            LEFT JOIN run r ON r.task_id=t.id AND r.status='RUNNING'
            WHERE t.status='RUNNING' ORDER BY t.id
        """),
        "last_success": one(db, "SELECT * FROM run WHERE status='SUCCEEDED' ORDER BY end_time_utc DESC LIMIT 1"),
        "last_failure": one(db, "SELECT * FROM run WHERE status='FAILED' ORDER BY end_time_utc DESC LIMIT 1"),
        "blockers": limited(db, "SELECT * FROM issue WHERE status IN ('OPEN','BLOCKED') ORDER BY id"),
        "resources": resources(db),
        "disks": disks(data_root, code_root),
        "deployment": deployment(db, code_root)["current"],
    }


def batches(db, batch_id=None):
    sql = """
        SELECT b.*,
            (SELECT COUNT(*) FROM run WHERE batch_id=b.id) AS run_count
        FROM batch b
    """
    if batch_id is None:
        return {"batches": limited(db, sql + " ORDER BY b.id DESC")}
    batch = one(db, sql + " WHERE b.id=?", batch_id)
    tasks = limited(db, "SELECT * FROM task WHERE batch_id=? ORDER BY id", batch_id)
    runs = limited(db, "SELECT * FROM run WHERE batch_id=? ORDER BY start_time_utc DESC", batch_id)
    outcomes = rows(db, """
        WITH latest AS (
            SELECT ri.object_type,ri.stage,ri.status, row_number() OVER (
                PARTITION BY ri.object_type,ri.object_id,ri.stage
                ORDER BY ri.updated_at DESC,ri.run_id DESC
            ) AS rank
            FROM run_item ri JOIN run r ON r.id=ri.run_id WHERE r.batch_id=?
        )
        SELECT object_type,stage,status, COUNT(*) AS amount FROM latest WHERE rank=1
        GROUP BY object_type,stage,status ORDER BY stage,object_type,status
    """, batch_id)
    usage = one(db, """
        SELECT SUM(downloaded_bytes) AS downloaded_bytes, MAX(peak_disk_usage) AS peak_disk_usage,
            (SELECT COUNT(*) FROM asset_file a JOIN run r ON r.id=a.fetch_run_id
               WHERE r.batch_id=? AND a.fetched_at IS NOT NULL) AS recorded_downloaded_files
        FROM run WHERE batch_id=?
    """, batch_id, batch_id)
    return {"batch": batch, "tasks": tasks, "runs": runs, "latest_stage_outcomes": outcomes, "usage": usage}


def runs(db, run_id=None, data_root=None):
    if run_id is None:
        return {"runs": limited(db, "SELECT * FROM run ORDER BY start_time_utc DESC")}
    record = one(db, "SELECT * FROM run WHERE id=?", run_id)
    preview = None
    if record and record["report_path"]:
        path = Path(record["report_path"]).resolve()
        if path.is_relative_to(data_root / "runs") and path.is_file():
            with path.open(encoding="utf-8", errors="replace") as report:
                preview = report.read(2800)
    return {
        "run": record,
        "report_preview": preview,
        "items": limited(db, "SELECT * FROM run_item WHERE run_id=? ORDER BY updated_at DESC", run_id),
        "item_count": db.execute("SELECT COUNT(*) FROM run_item WHERE run_id=?", (run_id,)).fetchone()[0],
    }


def pdb(db, pdb_id):
    entry = one(db, "SELECT * FROM pdb_entry WHERE id=?", pdb_id)
    versions = limited(db, "SELECT * FROM entry_version WHERE pdb_id=? ORDER BY id", pdb_id)
    files = limited(db, """
        SELECT asset_file.* FROM asset_file JOIN entry_version ON entry_version.id=entry_version_id
        WHERE pdb_id=? ORDER BY asset_file.id
    """, pdb_id)
    for asset in files:
        asset["exists_now"] = Path(asset["uri"]).is_file()
    objects_sql = """
        SELECT wi.*, ri.stage AS latest_stage,ri.status AS latest_status,
            ri.error_code,ri.error_message,ri.run_id,ri.output_uri,
            (SELECT status FROM run_item WHERE object_id=wi.id AND object_type=wi.kind AND stage='QC'
             ORDER BY updated_at DESC,run_id DESC LIMIT 1) AS qc_status,
            (SELECT status FROM run_item WHERE object_id=wi.id AND object_type=wi.kind AND stage='SIMILARITY'
             ORDER BY updated_at DESC,run_id DESC LIMIT 1) AS similarity_status
        FROM work_item wi JOIN entry_version ev ON ev.id=wi.entry_version_id
        LEFT JOIN run_item ri ON ri.rowid=(
            SELECT rowid FROM run_item x WHERE x.object_id=wi.id AND x.object_type=wi.kind
            ORDER BY updated_at DESC,run_id DESC LIMIT 1
        ) WHERE ev.pdb_id=? AND wi.kind=? ORDER BY wi.id
    """
    by_kind = {kind: limited(db, objects_sql, pdb_id, kind)
               for kind in ("assembly", "chain", "ligand", "system")}
    for system in by_kind["system"]:
        for role in ("protein", "ligand"):
            members = limited(db, """
                SELECT wi.id,wi.label FROM system_member sm JOIN work_item wi ON wi.id=sm.member_id
                WHERE sm.system_id=? AND sm.role=? ORDER BY wi.id
            """, system["id"], role)
            system[f"{role}_instances"] = ", ".join(item["id"] for item in members)
    history = limited(db, """
        SELECT ri.*,r.batch_id,r.git_commit,r.config_version FROM run_item ri
        JOIN run r ON r.id=ri.run_id WHERE ri.pdb_id=? ORDER BY ri.updated_at DESC
    """, pdb_id)
    return {"entry": entry, "versions": versions, "files": files, "objects": by_kind, "history": history}


def disks(data_root, code_root):
    return {key: shutil.disk_usage(path)._asdict() for key, path in
            (("data", data_root), ("code", code_root)) if path.is_dir()}


def resource_counts(db):
    counts = one(db, """
        SELECT (SELECT COUNT(*) FROM (
              SELECT DISTINCT source_name, source_release, source_record_id FROM entry_version
            )) AS source_records,
            (SELECT COUNT(DISTINCT sha256) FROM asset_file) AS unique_files,
            (SELECT COUNT(*) FROM asset_file) AS registered_files,
            (SELECT COUNT(*) FROM pdb_entry) AS entries,
            (SELECT COUNT(DISTINCT pdb_id) FROM run_item WHERE stage IN ('PARSE','SYSTEM','QC')
                 AND status IN ('PARSED','IDENTIFIED','ACCEPTED','REJECTED','UNCERTAIN')) AS processed_entries,
            (SELECT COUNT(*) FROM entry_version) AS entry_versions,
            (SELECT COUNT(*) FROM release WHERE kind='resource') AS resource_releases,
            (SELECT COUNT(DISTINCT component_id) FROM work_item WHERE kind='ligand') AS ligand_components
    """)
    for kind in ("assembly", "chain", "ligand", "system"):
        counts[kind] = db.execute("SELECT COUNT(*) FROM work_item WHERE kind=?", (kind,)).fetchone()[0]
    for status in ("ACCEPTED", "UNCERTAIN", "REJECTED"):
        counts[f"QC_{status}"] = db.execute("""
            SELECT COUNT(*) FROM work_item wi
            WHERE wi.kind='system' AND ?=(
                SELECT status FROM run_item
                WHERE object_type='system' AND object_id=wi.id AND stage='QC'
                ORDER BY updated_at DESC,run_id DESC LIMIT 1
            )
        """, (status,)).fetchone()[0]
    if counts["source_records"] == 0:
        return {"status": "未盘点", "counts": None, "basis": "等待登记真实来源快照"}
    return {"status": "已登记", "counts": counts,
            "basis": "登记库范围。来源按来源名、发布版和来源记录号去重。文件按 SHA-256 去重，QC 取各 System 最新 QC 记录。"}


def resources(db):
    snapshot = one(db, "SELECT * FROM storage_snapshot ORDER BY sampled_at DESC LIMIT 1")
    if snapshot:
        summary = json.loads(snapshot["report_json"]).get("resources")
        if summary:
            return dict(summary, sampled_at=snapshot["sampled_at"])
    has_entries = db.execute("SELECT EXISTS(SELECT 1 FROM entry_version)").fetchone()[0]
    return {"status": "待汇总" if has_entries else "未盘点", "counts": None,
            "basis": "等待执行端发布资源统计快照" if has_entries else "等待登记真实来源快照"}


def host_sample():
    stat = Path("/proc/stat")
    mem = Path("/proc/meminfo")
    result = {"cpu_count": os.cpu_count()}
    if stat.is_file() and mem.is_file():
        with stat.open() as stream:
            before = [int(value) for value in stream.readline().split()[1:9]]
        time.sleep(0.10)
        with stat.open() as stream:
            after = [int(value) for value in stream.readline().split()[1:9]]
        elapsed = sum(after) - sum(before)
        result["cpu_used_percent"] = round(100 * (1 - ((after[3] + after[4]) -
                                                     (before[3] + before[4])) / elapsed), 1) if elapsed else None
        with mem.open() as stream:
            for line in stream:
                key, value = line.split(":", 1)
                if key in ("MemTotal", "MemAvailable"):
                    result[key + "_bytes"] = int(value.split()[0]) * 1024
    return result


def storage(db, data_root, code_root):
    snapshot = one(db, "SELECT * FROM storage_snapshot ORDER BY sampled_at DESC LIMIT 1")
    active = limited(db, """SELECT id,batch_id,task_id,pid,start_time_utc,git_commit FROM run
                           WHERE status='RUNNING' ORDER BY start_time_utc DESC""")
    for run in active:
        run["pid_visible"] = Path(f"/proc/{run['pid']}").is_dir() if Path("/proc").is_dir() else None
        run["elapsed_seconds"] = round((datetime.now(timezone.utc) -
                                       datetime.fromisoformat(run["start_time_utc"])).total_seconds(), 1)
    data = json.loads(snapshot["report_json"]) if snapshot else {}
    layer_snapshot = {"sampled_at": snapshot["sampled_at"],
                      **data.get("occupied_bytes", {})} if snapshot else None
    return {"sampled_at": utc_now(), "disks": disks(data_root, code_root), "host": host_sample(),
            "active_runs": active, "layer_snapshot": layer_snapshot, "resources": resources(db),
            "transfers": data.get("transfers", {})}


def query(db, data_root, code_root, view, value=None):
    query_args(view, value)
    if view == "project":
        return project(db, data_root, code_root)
    if view == "modules":
        return {"modules": rows(db, "SELECT * FROM module ORDER BY id")}
    if view in ("batches", "batch"):
        return batches(db, value)
    if view in ("runs", "run"):
        return runs(db, value, data_root)
    if view == "pdb":
        return pdb(db, value.upper())
    if view == "storage":
        return storage(db, data_root, code_root)
    if view == "deployment":
        return deployment(db, code_root)
    return {table: limited(db, f"SELECT * FROM {table} ORDER BY id DESC")
            for table in ("decision", "threshold", "issue", "config_version")}
