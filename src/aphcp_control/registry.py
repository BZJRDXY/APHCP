import os
import shutil
import sqlite3
import tempfile
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


def sqlite_identity(db):
    return (db.execute("PRAGMA application_id").fetchone()[0],
            db.execute("PRAGMA user_version").fetchone()[0],
            db.execute("PRAGMA page_count").fetchone()[0],
            [tuple(r) for r in db.execute("SELECT type,name,tbl_name,sql FROM sqlite_schema ORDER BY type,name")])


def publish_registry(staged, target):
    """Durably publish a complete SQLite file without replacing another writer's file."""
    if not staged.stat().st_size:
        raise ValueError("Empty SQLite backup")
    with staged.open("r+b" if os.name == "nt" else "rb") as stream:
        os.fsync(stream.fileno())
    target.hardlink_to(staged)
    if os.name == "posix":
        fd = os.open(target.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)


class RegistryDiskReserve(RuntimeError):
    pass


def require_backup_space(source, target):
    reserve = int(os.environ.get('APHCP_ROLLING_RESERVE_BYTES','0'))
    if reserve:
        needed = Path(source).stat().st_size
        wal = Path(f'{source}-wal')
        if wal.exists():
            needed += wal.stat().st_size
        if shutil.disk_usage(target.parent).free < reserve + needed:
            raise RegistryDiskReserve('registry_copy_would_reach_disk_reserve')


def backup(source, target):
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() or target.is_symlink():
        raise FileExistsError(target)
    require_backup_space(source,target)
    fd, name = tempfile.mkstemp(dir=target.parent, prefix=f".{target.name}.backup-")
    os.close(fd)
    staged = Path(name)
    try:
        with connect(source, readonly=True) as live, connect(staged) as saved:
            live.execute("BEGIN")  # Pin the source snapshot across inspection and backup.
            original = sqlite_identity(live)
            if not original[2]:
                raise ValueError("Empty SQLite backup source")
            live.backup(saved)
            if sqlite_identity(saved) != original or saved.execute("PRAGMA quick_check").fetchone()[0] != "ok":
                raise ValueError("SQLite backup source or integrity differs")
            if saved.execute("PRAGMA journal_mode=DELETE").fetchone()[0].lower() != "delete":
                raise ValueError("SQLite backup must be a self-contained file")
        publish_registry(staged, target)
    finally:
        staged.unlink(missing_ok=True)
    return target


def recover_empty_registry(target, sealed):
    """Called only under the writer lease for a distinct staging destination."""
    if target.is_file() and not target.stat().st_size:
        if any(p.exists() for p in sealed) or any(Path(f"{target}{s}").exists() for s in ("-wal", "-shm", "-journal")):
            raise ValueError("Empty registry has a commit receipt or active journal")
        target.unlink()


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
