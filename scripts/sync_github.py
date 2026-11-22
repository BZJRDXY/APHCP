"""Sync APHCP's public source snapshot without staging or publishing the private local tree."""
import argparse
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import tempfile

from public_review_status import SOURCES, STATUS, render, runtime_fingerprint

ROOT = Path(__file__).resolve().parents[1]
REMOTE = 'https://github.com/BZJRDXY/APHCP.git'
BRANCH = 'refs/heads/main'
SYNCED = 'refs/aphcp/github-public'
RULE = 'docs/governance/代码同步规则.md'
REVIEW_GUIDE = 'docs/public/review.md'
TOP = {'.gitattributes', '.gitignore', 'pyproject.toml', 'requirements-local.txt', 'requirements-server.txt', 'LICENSE', '.streamlit/config.toml'}
CONFIG = {'deployment.example.yaml', 'project.yaml', 'environment.yaml', 'temporal_split.yaml'}
SECRET = re.compile(rb'-----BEGIN [A-Z ]*PRIVATE KEY-----|gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}|'
                    rb'AKIA[A-Z0-9]{16}|https?://[^/\s:]+:[^/\s@]+@|/data/(?:exp|fast2)/users/|'
                    rb'(?i:[A-Z]:[/\\]Users[/\\][^/\\\s]+)|codex:' rb'//threads/')
PLACEHOLDERS = {
    'config/deployment.example.yaml': {
        'host': 'aphcp.example.org', 'user': 'aphcp', 'base_python': '/opt/aphcp/env/bin/python',
        'code_root': '/srv/aphcp/code', 'data_root': '/srv/aphcp/data', 'upstream_root': '/srv/aphcp/upstream',
    },
    'config/project.yaml': {'physical_data_root': '/srv/aphcp/data', 'environment_lock': '/srv/aphcp/code/deploy/environment.lock'},
}


def git(*args, env=None, data=None, allow_missing=False):
    p = subprocess.run(['git', *args], cwd=ROOT, env=env, input=data, stdout=subprocess.PIPE)
    if allow_missing and p.returncode == 1:
        return b''
    p.check_returncode()
    return p.stdout


def included(name):
    p = PurePosixPath(name)
    if name in TOP or name in {RULE, REVIEW_GUIDE, 'workflow/Snakefile'}:
        return True
    if len(p.parts) < 2 or any(part.startswith('.') for part in p.parts):
        return False
    group = p.parts[0]
    return (group == 'scripts' and p.suffix in {'.py', '.sh', '.ps1'} or
            group == 'src' and p.suffix in {'.py', '.html', '.css', '.js'} or
            group == 'tests' and (p.suffix == '.py' or name == 'tests/fixtures/m0_work_items.json') or
            group == 'migrations' and p.suffix == '.sql' or
            group == 'config' and len(p.parts) == 2 and p.name in CONFIG)


def snapshot(require_receipts=False):
    source = set(git('ls-files', '--cached', '--others', '--exclude-standard', '-z').decode('utf8').split('\0'))
    readme = 'README.github.md' if (ROOT/'README.github.md').is_file() else 'README.md'
    source.add(readme)
    payload = {}
    for name in sorted(source):
        if name != readme and not included(name):
            continue
        path = ROOT/name
        if not path.exists():
            continue
        if not path.is_file() or path.is_symlink() or path.resolve() != path.absolute():
            raise ValueError(f'Not an ordinary workspace file: {name}')
        data = path.read_bytes()
        if name in PLACEHOLDERS:
            text = data.decode('utf8')
            for key, val in PLACEHOLDERS[name].items():
                text, count = re.subn(rf'(?m)^(\s*{key}:)[ \t]*[^\r\n]*$', lambda m: f'{m[1]} {val}', text)
                if count != 1:
                    raise ValueError(f'Expected one deployment template field: {name}: {key}')
            data = text.encode('utf8')
        if name == 'scripts/m1_holo_join.py':
            data = re.sub(rb"(/data/(?:exp|fast2)/users/)[^'\r\n]+/PLINDER(?=')", b'/srv/aphcp/upstream/PLINDER', data)
        if SECRET.search(data) or len(data) > 2_000_000:
            raise ValueError(f'Review this file before publishing: {name}')
        payload['README.md' if name == readme else name] = data
    required = {'README.md', RULE, REVIEW_GUIDE, 'pyproject.toml', 'scripts/sync_github.py',
                'src/aphcp_control/__main__.py'}
    if not required <= payload.keys():
        raise ValueError(f'Missing public source: {sorted(required - payload.keys())}')
    fingerprint = runtime_fingerprint(payload)
    if (ROOT / SOURCES).is_file():
        payload[STATUS] = render(ROOT, fingerprint).encode('utf-8')
    elif require_receipts:
        raise ValueError('Publish from the project workspace with its accepted private receipt selection.')
    else:
        payload[STATUS] = (ROOT / STATUS).read_bytes()
        if f'Public runtime fingerprint (SHA-256): `{fingerprint}`'.encode() not in payload[STATUS]:
            raise ValueError('Published status does not cover this checkout; regenerate from accepted receipts.')
    if SECRET.search(payload[STATUS]):
        raise ValueError('Private value in generated public status')
    return payload


def tree(payload):
    with tempfile.TemporaryDirectory(prefix='aphcp-public-index-') as temp:
        env = dict(os.environ, GIT_INDEX_FILE=str(Path(temp)/'index'))
        git('read-tree', '--empty', env=env)
        index = b''.join(b'100644 '+git('hash-object', '-w', '--stdin', data=data).strip()+b'\t'+name.encode('utf8')+b'\0'
                         for name, data in sorted(payload.items()))
        git('update-index', '-z', '--index-info', env=env, data=index)
        return git('write-tree', env=env).decode().strip()


def tip():
    data = git('ls-remote', '--refs', REMOTE, BRANCH).split()
    return data[0].decode() if data else ''


def sync(new_tree, message):
    parent = tip()
    last = git('for-each-ref', '--format=%(objectname)', SYNCED).decode().strip()
    if last and parent != last:
        raise ValueError('GitHub main changed independently; reconcile those changes before publishing again.')
    if parent:
        git('fetch', '--no-tags', REMOTE, BRANCH)
        if git('rev-parse', 'FETCH_HEAD').decode().strip() != parent:
            raise ValueError('GitHub main changed during fetch; rerun the sync.')
        if git('rev-parse', f'{parent}^{{tree}}').decode().strip() == new_tree:
            git('update-ref', SYNCED, parent)
            print(f'GitHub main already matches: {parent}')
            return
    env = dict(os.environ, GIT_AUTHOR_NAME='BZJRDXY', GIT_COMMITTER_NAME='BZJRDXY',
               GIT_AUTHOR_EMAIL='201316341+BZJRDXY@users.noreply.github.com',
               GIT_COMMITTER_EMAIL='201316341+BZJRDXY@users.noreply.github.com')
    commit = git('commit-tree', new_tree, *(['-p', parent] if parent else []), '-m', message, env=env).decode().strip()
    git('push', REMOTE, f'{commit}:{BRANCH}')
    if tip() != commit:
        raise ValueError('GitHub main did not return the commit just pushed.')
    git('update-ref', SYNCED, commit)
    print(f'GitHub main verified: {commit}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--push', action='store_true', help='publish the verified source snapshot')
    mode.add_argument('--show-status', action='store_true', help='preview the exact generated public summary')
    parser.add_argument('--expect-tree', help='abort if the current snapshot differs from a reviewed tree')
    parser.add_argument('--message', default='Sync APHCP project code')
    args = parser.parse_args()
    payload = snapshot(require_receipts=args.push)
    if args.show_status:
        print(payload[STATUS].decode('utf-8'), end='')
        return
    new_tree = tree(payload)
    if args.expect_tree and args.expect_tree != new_tree:
        parser.error('working snapshot changed since preview')
    print(f'Public tree: {new_tree} ({len(payload)} files)', flush=True)
    if args.push:
        sync(new_tree, args.message)
        (ROOT / STATUS).write_bytes(payload[STATUS])
    else:
        print('\n'.join(sorted(payload)))


if __name__ == '__main__':
    main()
