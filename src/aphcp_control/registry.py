import sqlite3
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path


MIGRATIONS = Path(__file__).resolve().parents[2] / "migrations"


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def database(root):
    return Path(root) / "registry" / "aphcp_registry.sqlite3"


@contextmanager
def connect(path, *, readonly=False):
    path = Path(path).resolve()
    target = path.as_uri() + "?mode=ro" if readonly else str(path)
    connection = sqlite3.connect(target, uri=readonly, timeout=3)
    connection.row_factory = sqlite3.Row
    try:
        connection.execute("PRAGMA foreign_keys=ON")
        if readonly:
            connection.execute("PRAGMA query_only=ON")
        with connection:
            yield connection
    finally:
        connection.close()


def backup(source, target):
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("xb"):
        pass
    try:
        with connect(source, readonly=True) as live, connect(target) as saved:
            live.backup(saved)
            if saved.execute("PRAGMA quick_check").fetchone()[0] != "ok":
                raise ValueError("SQLite 备份完整性检查失败")
    except Exception:
        target.unlink()
        raise
    return target


def migrate(path, directory=MIGRATIONS):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with connect(path) as connection:
        version = connection.execute("PRAGMA user_version").fetchone()[0]
        migrations = sorted(directory.glob("[0-9][0-9][0-9]_*.sql"))
        if not migrations or version > int(migrations[-1].name[:3]):
            raise ValueError("迁移文件缺失，或数据库版本比当前代码新")
        pending = [file for file in migrations if int(file.name[:3]) > version]
        if version and pending:
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
            backup(path, path.parent / "backups" / f"pre-migration-{stamp}.sqlite3")
        for file in pending:
            next_version = int(file.name[:3])
            if next_version != version + 1:
                raise ValueError(f"缺少迁移版本 {version + 1}")
            script = file.read_text(encoding="utf-8")
            connection.executescript("BEGIN IMMEDIATE;\n" + script)
            connection.execute(
                "INSERT INTO schema_migrations VALUES (?, ?, ?)", (next_version, file.name, utc_now())
            )
            connection.execute(f"PRAGMA user_version={next_version}")
            if connection.execute("PRAGMA foreign_key_check").fetchone():
                raise ValueError(f"迁移 {file.name} 产生外键错误")
            connection.commit()
            version = next_version
    return version
