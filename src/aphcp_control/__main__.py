import argparse
import json
import sqlite3
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from .queries import QUERIES, LIST_LIMIT, one, query, query_args, resource_counts
from .registry import backup, connect, database, migrate, utc_now
from .trial import trial


def cli():
    parser = argparse.ArgumentParser(description="APHCP 项目管理")
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--code-root", type=Path)
    commands = parser.add_subparsers(dest="command", required=True)
    status = commands.add_parser("status")
    status.add_argument("view", choices=sorted(QUERIES))
    status.add_argument("identifier", nargs="?")
    status.add_argument("--json", action="store_true", help="JSON 是默认输出格式")
    commands.add_parser("init")
    commands.add_parser("backup")
    restore = commands.add_parser("restore")
    restore.add_argument("from_file", type=Path)
    restore.add_argument("to_file", type=Path)
    commands.add_parser("measure-storage")
    smoke = commands.add_parser("trial")
    smoke.add_argument("--fixture", required=True, type=Path)
    smoke.add_argument("--git-commit", required=True)
    smoke.add_argument("--env-lock", required=True, type=Path)
    smoke.add_argument("--interrupt-after", type=int, choices=(1, 2))
    return parser.parse_args()


def measure_storage(db_path, data_root, code_root):
    paths = {name: data_root / name for name in (
        "raw", "normalized", "annotations", "relations", "manifests", "runs", "releases"
    )}
    paths.update({"data_root": data_root, "code_root": code_root, "scratch": code_root / "work" / "scratch"})
    values = {}
    for name, path in paths.items():
        if path.is_dir():
            output = subprocess.check_output(["du", "-sB1", "--", str(path)], encoding="utf-8", timeout=120)
            values[name] = int(output.split(maxsplit=1)[0])
    sampled_at = utc_now()
    with connect(db_path) as db:
        summary = {"occupied_bytes": values, "resources": resource_counts(db),
                   "transfers": one(db, """SELECT SUM(downloaded_bytes) AS downloaded_bytes,
                       (SELECT COUNT(*) FROM run_item WHERE object_type='asset' AND status='REUSED') AS cache_hits
                       FROM run""")}
        db.execute("INSERT INTO storage_snapshot VALUES (?,?)", (sampled_at, json.dumps(summary)))
    return {"sampled_at": sampled_at, "occupied_bytes": values}


def run(args):
    data_root = args.data_root.resolve()
    code_root = (args.code_root or data_root).resolve()
    db_path = database(data_root)
    if args.command == "init":
        return {"database": str(db_path), "schema_version": migrate(db_path)}
    if args.command == "backup":
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        path = backup(db_path, db_path.parent / "backups" / f"registry-{stamp}.sqlite3")
        return {"backup": str(path)}
    if args.command == "restore":
        source, target = args.from_file.resolve(), args.to_file.resolve()
        if any(not path.is_relative_to(db_path.parent) for path in (source, target)):
            raise ValueError("恢复文件须位于项目 registry 内")
        return {"restored_file": str(backup(source, target))}
    if args.command == "trial":
        return trial(data_root, args.fixture, args.git_commit, args.env_lock, args.interrupt_after)
    if args.command == "measure-storage":
        return measure_storage(db_path, data_root, code_root)
    arguments = query_args(args.view, args.identifier)
    with connect(db_path, readonly=True) as db:
        result = query(db, data_root, code_root, *arguments)
    return {"queried_at": utc_now(), "detail_limit": LIST_LIMIT, "data": result}


def main():
    try:
        print(json.dumps(run(cli())))
    except (OSError, ValueError, sqlite3.Error, subprocess.SubprocessError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
