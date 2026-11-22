import json
import re
import shlex
import subprocess
from pathlib import Path, PurePosixPath

from .queries import query_args


class QueryError(RuntimeError):
    def __init__(self, state, message):
        super().__init__(message)
        self.state = state


def load_settings(root):
    import yaml

    local = root / "config" / "deployment.local.yaml"
    settings_file = local if local.is_file() else local.with_name("deployment.example.yaml")
    settings = yaml.safe_load(settings_file.read_text(encoding="utf-8"))
    ssh = settings["ssh"]
    if not re.fullmatch(r"[A-Za-z0-9.-]+", ssh["host"]) or not re.fullmatch(r"[A-Za-z0-9_-]+", ssh["user"]):
        raise ValueError("SSH 主机或用户格式不正确")
    for key in ("code_root", "data_root"):
        path = PurePosixPath(settings["server"][key])
        if not path.is_absolute() or ".." in path.parts:
            raise ValueError(f"服务器路径必须为明确的绝对路径：{key}")
    return settings


def ssh_command(settings):
    ssh = settings["ssh"]
    identity = Path(ssh["identity"]).expanduser()
    return ["ssh", "-T", "-i", str(identity), "-o", "IdentitiesOnly=yes", "-o", "BatchMode=yes",
            "-o", "StrictHostKeyChecking=yes", "-o", "UpdateHostKeys=no", "-o", "ClearAllForwardings=yes",
            "-o", f"ConnectTimeout={int(ssh['connect_timeout'])}", f"{ssh['user']}@{ssh['host']}"]


def fetch(settings, view, identifier=None):
    arguments = query_args(view, identifier)
    server = settings["server"]
    code = PurePosixPath(server["code_root"])
    python = str(code / "venv" / "bin" / "python")
    current = str(code / "current")
    cli = ["env", f"PYTHONPATH={current}/src", python, "-B", "-m", "aphcp_control",
           "--data-root", server["data_root"], "--code-root", str(code), "status", *arguments, "--json"]
    command = f"test -x {shlex.quote(python)} && test -d {shlex.quote(current)} || exit 78; exec {shlex.join(cli)}"
    return request_json(settings, command)


def fetch_inventory(settings, identifier):
    """Query an approved frozen snapshot without a formal server deployment."""
    from .inventory_lookup import canonical_id

    _, identifier = query_args("pdb", identifier)
    command = shlex.join(["python3", "-B", "-", settings["server"]["data_root"], canonical_id(identifier)])
    program = Path(__file__).with_name("inventory_lookup.py").read_text(encoding="utf-8-sig")
    return request_json(settings, command, stdin=program, limit=32_000)


def request_json(settings, command, *, stdin=None, limit=512_000):
    try:
        result = subprocess.run([*ssh_command(settings), command], input=stdin, capture_output=True, text=True,
                                encoding="utf-8", errors="replace", timeout=9, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise QueryError("离线", str(exc)) from exc
    if result.returncode == 78:
        raise QueryError("SSH 已连接 · 待部署", "Pro5000 尚未提供 APHCP 查询程序")
    if result.returncode == 255:
        raise QueryError("离线", result.stderr.strip()[:240] or "SSH 连接失败")
    if result.returncode:
        raise QueryError("查询失败", result.stderr.strip()[:240] or f"查询退出码 {result.returncode}")
    if len(result.stdout.encode("utf-8")) > limit:
        raise QueryError("查询失败", "摘要超过大小上限，请缩小查询范围")
    try:
        payload = json.loads(result.stdout)
        if not isinstance(payload["data"], dict) or not payload["queried_at"]:
            raise ValueError("摘要格式错误")
        return payload
    except (ValueError, KeyError, TypeError) as exc:
        raise QueryError("查询失败", "服务器没有返回有效的 JSON 摘要") from exc


def git_state(root):
    def read(*args):
        process = subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True,
                                 encoding="utf-8", errors="replace", timeout=4, check=False)
        return process.stdout.strip() if process.returncode == 0 else None

    mods = read("status", "--porcelain")
    return {"branch": read("branch", "--show-current"), "commit": read("rev-parse", "HEAD"),
            "dirty": bool(mods) if mods is not None else None}


def in_sync(local, remote):
    return bool(local["commit"] and local["dirty"] is False and remote and local["commit"] == remote.get("commit"))
