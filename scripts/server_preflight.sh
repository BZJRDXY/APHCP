#!/usr/bin/env bash
set -euo pipefail
code=$1
data=$2
base_python=$3
upstream=$4
id -un
hostname
date -Is
findmnt -T "$(dirname "$code")"
findmnt -T "$(dirname "$data")"
df -hT "$(dirname "$code")" "$(dirname "$data")"
"$base_python" -B -c 'import sys,sqlite3; print(sys.version.split()[0],sqlite3.sqlite_version); sqlite3.connect(":memory:").execute("select 1")'
test -w "$(dirname "$code")" && test -w "$(dirname "$data")"
test -r "$upstream"
ls -ld "$upstream"
for path in "$code" "$data"; do
    if test -e "$path"; then ls -ld -- "$path"; else printf 'absent %s\n' "$path"; fi
done
