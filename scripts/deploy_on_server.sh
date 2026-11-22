#!/usr/bin/env bash
set -euo pipefail
umask 027
code=$1
data=$2
commit=$3
[[ "$commit" =~ ^[a-f0-9]{40}$ ]]
test -f "$code/.aphcp-m0" && test -f "$data/.aphcp-m0"
python="$code/venv/bin/python"
repo="$code/repo.git"
release="$code/releases/$commit"
test "$(git --git-dir="$repo" rev-parse "${commit}^{commit}")" = "$commit"
if [[ -e "$code/current" && ! -L "$code/current" ]]; then
    printf 'current must be a symbolic link\n' >&2
    exit 1
fi

candidate="$release"
if [[ ! -d "$release" ]]; then
    candidate=$(mktemp -d "$code/releases/.staging-${commit:0:12}.XXXXXX")
    git --git-dir="$repo" archive "$commit" | tar -xf - -C "$candidate"
fi

dependency_receipt="$code/deploy/requirements-server.applied.txt"
if [[ -f "$dependency_receipt" ]]; then
    cmp "$candidate/requirements-server.txt" "$dependency_receipt" || {
        printf 'This commit needs a separately prepared server environment\n' >&2
        exit 1
    }
else
    PIP_CACHE_DIR="$code/work/scratch/pip-cache" TMPDIR="$code/work/scratch" \
        "$python" -m pip --disable-pip-version-check install -r "$candidate/requirements-server.txt"
    "$python" -m pip check
    {
        "$python" --version
        "$python" -m pip freeze
    } > "$code/deploy/environment.lock"
    cp -- "$candidate/requirements-server.txt" "$dependency_receipt"
fi

export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$candidate/src"
"$code/venv/bin/snakemake" --version
cli=("$python" -m aphcp_control --data-root "$data" --code-root "$code")
db="$data/registry/aphcp_registry.sqlite3"
if [[ ! -f "$db" || ! -e "$data/registry/.initialized" ]]; then
    "${cli[@]}" init
    "${cli[@]}" backup
    touch "$data/registry/.initialized"
fi
"${cli[@]}" status project --json

stamp=$(date -u +%Y-%m-%dT%H:%M:%SZ)
deployment=$(printf '{"commit":"%s","deployed_at_utc":"%s"}' "$commit" "$stamp")
if [[ "$candidate" != "$release" ]]; then
    printf '%s\n' "$deployment" > "$candidate/deployment.json"
    chmod -R a-w "$candidate"
    mv -T -- "$candidate" "$release"
fi
"$python" -B -c 'import json,sys; expected=sys.argv[2]; actual=json.load(open(sys.argv[1]))["commit"]; sys.exit(0 if actual == expected else 1)' "$release/deployment.json" "$commit"
link="$code/.current-$$"
ln -s "releases/$commit" "$link"
mv -Tf -- "$link" "$code/current"
printf '%s\n' "$deployment" >> "$code/deploy/history.jsonl"
printf '%s\n' "$deployment"
