#!/usr/bin/env bash
set -euo pipefail
umask 027
code=$1
data=$2
base_python=$3

for root in "$code" "$data"; do
    if [[ -e "$root" || -L "$root" ]]; then
        [[ ! -L "$root" && -f "$root/.aphcp-m0" ]] || {
            printf 'Existing path is outside this initialization: %s\n' "$root" >&2
            exit 1
        }
    fi
done
for root in "$code" "$data"; do
    if [[ ! -d "$root" ]]; then
        mkdir -- "$root"
        printf 'APHCP M0\n' > "$root/.aphcp-m0"
    fi
done
mkdir -p "$code/releases" "$code/deploy" "$code/work/scratch"
mkdir -p "$data/registry/backups" "$data/raw" "$data/normalized" "$data/annotations" \
    "$data/relations" "$data/manifests" "$data/batches" "$data/runs" "$data/reports" "$data/releases/resource"
if [[ ! -f "$code/repo.git/HEAD" ]]; then
    git init --bare "$code/repo.git"
fi
if [[ ! -x "$code/venv/bin/python" ]]; then
    "$base_python" -m venv "$code/venv"
fi
"$code/venv/bin/python" -B - "$data/registry" <<'PY'
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

with tempfile.TemporaryDirectory(prefix="sqlite-preflight-", dir=sys.argv[1]) as scratch:
    path = Path(scratch) / "probe.sqlite3"
    holder = sqlite3.connect(path)
    try:
        holder.execute("CREATE TABLE probe (id INTEGER)")
        holder.execute("BEGIN IMMEDIATE")
        contending_writer = """
import sqlite3, sys
with sqlite3.connect(sys.argv[1], timeout=0.25) as db:
    db.execute('BEGIN IMMEDIATE')
    db.execute('INSERT INTO probe VALUES (1)')
"""
        command = [sys.executable, "-B", "-c", contending_writer, str(path)]
        attempt = subprocess.run(command, capture_output=True, text=True, timeout=10)
        if attempt.returncode == 0 or "database is locked" not in attempt.stderr:
            raise SystemExit("SQLite: cross-process write lock failed")
        holder.rollback()
        subprocess.run(command, check=True, timeout=10)
        if holder.execute("SELECT COUNT(*) FROM probe").fetchone()[0] != 1:
            raise SystemExit("SQLite: post-lock write or durable read failed")
    finally:
        holder.close()
print("SQLite on-disk cross-process lock and recovery passed")
PY
{
    date -Is
    printf 'code=%s\ndata=%s\n' "$code" "$data"
    "$code/venv/bin/python" -B --version
    git --git-dir="$code/repo.git" rev-parse --is-bare-repository
    printf 'sqlite_on_disk_lock=passed\n'
} >> "$data/reports/M0_server_initialization.txt"
