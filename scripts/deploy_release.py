import argparse
import json
import os
import shlex
import subprocess
from pathlib import Path

from aphcp_control.transport import fetch, git_state, in_sync, load_settings, ssh_command


ROOT = Path(__file__).resolve().parents[1]


def git(*args):
    return ["git", "-C", str(ROOT), *args]


def main():
    parser = argparse.ArgumentParser(description="APHCP 部署预审与按 Commit 部署")
    parser.add_argument("action", choices=("preflight", "initialize", "deploy", "verify"))
    parser.add_argument("--commit", default="HEAD")
    parser.add_argument("--approved", action="store_true", help="确认用户已明确授权本次 Pro5000 写操作")
    args = parser.parse_args()
    if args.action in ("initialize", "deploy") and not args.approved:
        parser.error("服务器写入须经用户明确批准")
    if args.action == "deploy":
        try:
            pending = subprocess.check_output(git("status", "--porcelain", "--untracked-files=all"))
        except (OSError, subprocess.CalledProcessError):
            parser.error("无法检查本地 Git 工作区，已停止部署")
        if pending:
            parser.error("本地 Git 工作区存在未提交的变更（含未跟踪文件），已停止部署")
        commit = subprocess.check_output(
            git("rev-parse", "--verify", "--end-of-options", f"{args.commit}^{{commit}}"), text=True
        ).strip()
    settings = load_settings(ROOT)
    server = settings["server"]
    code, data, upstream = (server[key] for key in ("code_root", "data_root", "upstream_root"))
    base_python = server["base_python"]
    ssh = ssh_command(settings)
    if args.action == "verify":
        local = git_state(ROOT)
        deployed = fetch(settings, "deployment")["data"]
        print(json.dumps({"local": local, "server": deployed, "in_sync": in_sync(local, deployed["current"])}))
        return
    if args.action in ("preflight", "initialize"):
        name = {"preflight": "server_preflight.sh", "initialize": "server_initialize.sh"}[args.action]
        script = (ROOT / "scripts" / name).read_bytes()
        cmd = shlex.join(["bash", "-s", "--", code, data, base_python, upstream])
        subprocess.run([*ssh, cmd], input=script, check=True)
        return
    print(f"部署 Commit SHA: {commit}", flush=True)
    host = ssh[-1]
    environment = dict(os.environ, GIT_SSH_COMMAND=shlex.join(ssh[:-1]))
    repo = f"{code}/repo.git"
    subprocess.run(git("push", f"ssh://{host}{repo}", f"{commit}:refs/heads/deploy/{commit}"),
                   env=environment, check=True)
    deployer = shlex.join(["git", f"--git-dir={repo}", "show", f"{commit}:scripts/deploy_on_server.sh"])
    parameters = shlex.join(["bash", "-s", "--", code, data, commit])
    subprocess.run([*ssh, f"set -o pipefail; {deployer} | {parameters}"], check=True)


if __name__ == "__main__":
    main()
