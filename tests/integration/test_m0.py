import json
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from streamlit.testing.v1 import AppTest

from aphcp_control import observation_lookup
from aphcp_control import transport
from aphcp_control.queries import query, query_args
from aphcp_control.registry import backup, connect, database, migrate


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests" / "fixtures" / "m0_work_items.json"


class M0Tests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / "m0-test"
        self.db_path = database(self.root)

    def run_cli(self, *args):
        cmd = [sys.executable, "-B", "-m", "aphcp_control", "--data-root", str(self.root), *args]
        return subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", timeout=30)

    def json_cli(self, *args):
        proc = self.run_cli(*args)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        return json.loads(proc.stdout)

    def test_registry_migration_backup_and_rollback(self):
        self.assertEqual(migrate(self.db_path), 1)
        with connect(self.db_path) as db:
            db.execute("INSERT INTO issue VALUES ('test','M0','OPEN','fixture',NULL)")
        with connect(self.db_path, readonly=True) as db:
            with self.assertRaises(sqlite3.OperationalError):
                db.execute("DELETE FROM issue")
        source = self.db_path.parent / "backups" / "test.sqlite3"
        backup(self.db_path, source)
        restored = self.db_path.parent / "restored.sqlite3"
        self.json_cli("restore", str(source), str(restored))
        with connect(restored, readonly=True) as db:
            self.assertEqual(db.execute("SELECT id FROM issue").fetchone()[0], "test")
        migrations = self.root / "migration-fixture"
        migrations.mkdir()
        for file in (ROOT / "migrations").glob("*.sql"):
            (migrations / file.name).write_bytes(file.read_bytes())
        (migrations / "002_broken.sql").write_text("CREATE TABLE rolled_back(x);\nINVALID SQL;", encoding="utf-8")
        with self.assertRaises(sqlite3.Error):
            migrate(self.db_path, migrations)
        with connect(self.db_path, readonly=True) as db:
            self.assertEqual(db.execute("PRAGMA user_version").fetchone()[0], 1)
            self.assertIsNone(db.execute("SELECT 1 FROM sqlite_master WHERE name='rolled_back'").fetchone())
        (migrations / "002_broken.sql").unlink()
        (migrations / "002_valid.sql").write_text("CREATE TABLE migration_passed(x);", encoding="utf-8")
        self.assertEqual(migrate(self.db_path, migrations), 2)
        self.assertEqual(migrate(self.db_path, migrations), 2)
        with connect(self.db_path, readonly=True) as db:
            self.assertEqual(db.execute("SELECT COUNT(*) FROM schema_migrations").fetchone()[0], 2)
        self.assertGreaterEqual(len(list((self.db_path.parent / "backups").glob("pre-migration-*"))), 1)

    def test_abrupt_exit_resume_and_two_systems(self):
        opts = ("trial", "--fixture", str(FIXTURE), "--git-commit", "a" * 40,
                "--env-lock", str(ROOT / "requirements-server.txt"))
        interrupted = self.run_cli(*opts, "--interrupt-after", "1")
        self.assertEqual(interrupted.returncode, 75, interrupted.stderr)
        runs = self.json_cli("status", "runs")["data"]["runs"]
        first = runs[0]["id"]
        self.assertEqual(runs[0]["status"], "RUNNING")
        detail = self.json_cli("status", "run", first)["data"]
        self.assertEqual(detail["item_count"], 1)
        a_path = Path(detail["items"][0]["output_uri"])
        first_mtime = a_path.stat().st_mtime_ns
        result = self.json_cli(*opts)
        self.assertEqual([item["status"] for item in result["items"]], ["REUSED", "PASS", "PASS"])
        self.assertEqual(a_path.stat().st_mtime_ns, first_mtime)
        prior = self.json_cli("status", "run", first)["data"]
        self.assertEqual(prior["run"]["status"], "PARTIAL")
        self.assertIn("PARTIAL", prior["report_preview"])
        detail = self.json_cli("status", "pdb", "0tst")["data"]
        self.assertEqual({system["id"] for system in detail["objects"]["system"]},
                         {"TEST-SYSTEM-A", "TEST-SYSTEM-B"})
        self.assertEqual(len({item["ligand_instances"] for item in detail["objects"]["system"]}), 2)
        self.assertEqual(len({item["protein_instances"] for item in detail["objects"]["system"]}), 1)
        self.assertEqual(len(detail["history"]), 3)
        batch = self.json_cli("status", "batch", result["batch_id"])["data"]
        self.assertEqual(batch["batch"]["status"], "ACCEPTED")
        self.assertEqual(sum(row["amount"] for row in batch["latest_stage_outcomes"]), 3)
        self.assertEqual(self.json_cli("status", "batch", "M0-EMPTY")["data"]["batch"]["candidate_count"], 0)
        corrupted = Path(result["items"][1]["output_uri"])
        corrupted.write_text("corrupt", encoding="utf-8")
        result = self.json_cli(*opts)
        self.assertEqual([item["status"] for item in result["items"]], ["REUSED", "PASS", "REUSED"])
        self.assertIn("m0-fixture-v1", corrupted.read_text(encoding="utf-8"))

    def test_ssh_allowlist_and_dirty_git(self):
        for invalid in ("1M17;uname", "1M17\nls", "../", "--json", "中文"):
            with self.assertRaises(ValueError):
                query_args("pdb", invalid)
        settings = transport.load_settings(ROOT)
        with patch("aphcp_control.transport.subprocess.run") as shell:
            shell.return_value = subprocess.CompletedProcess([], 255, "", "unreachable")
            with self.assertRaisesRegex(transport.QueryError, "unreachable") as error:
                transport.fetch(settings, "pdb", "1m17")
            self.assertEqual(error.exception.state, "离线")
            argv = shell.call_args.args[0]
            self.assertIn("StrictHostKeyChecking=yes", argv)
            self.assertIn("status pdb 1M17 --json", argv[-1])
        self.assertFalse(transport.in_sync({"commit": "a", "dirty": True}, {"commit": "a"}))
        self.assertTrue(transport.in_sync({"commit": "a", "dirty": False}, {"commit": "a"}))
        self.assertFalse(transport.in_sync({"commit": "a", "dirty": False}, {"commit": "b"}))

    def test_dashboard_pages_refresh_and_disconnection(self):
        migrate(self.db_path)

        def fixture_query(settings, view, identifier=None):
            with connect(self.db_path, readonly=True) as db:
                data = query(db, self.root, self.root, view, identifier)
            return {"data": data, "queried_at": "测试 UTC 时间", "detail_limit": 25}

        path = ROOT / "src" / "aphcp_dashboard" / "app.py"
        with patch.object(transport, "fetch", side_effect=fixture_query) as read, patch.object(
            transport, "fetch_inventory", return_value={"data": {"state": "absent"}, "queried_at": "测试 UTC 时间"}
        ) as inventory, patch.object(observation_lookup, "fetch_observations", return_value={
                "data": {"state": "absent"}, "queried_at": "测试 UTC 时间"}) as observations:
            app = AppTest.from_file(str(path), default_timeout=12).run()
            self.assertFalse(list(app.exception))
            self.assertEqual(read.call_count, 0, "Receipt-first homepage should not prefetch the registry")
            self.assertEqual(len(app.radio(key="nav").options), 8)
            for page in ("home", "m1", "m2", "pdb", "batches", "runs", "server", "govern"):
                app.radio(key="nav").set_value(page).run()
                self.assertFalse(list(app.exception), page)
            app.radio(key="nav").set_value("pdb").run()
            app.text_input[0].set_value("1M17")
            next(button for button in app.button if button.label == "查询 PDB").click().run()
            self.assertTrue(any("已查询的来源均未找到" in msg.value for msg in app.info))
            self.assertTrue(inventory.called)
            self.assertTrue(observations.called)
            observations.return_value = {"queried_at": "测试 UTC 时间", "data": {
                "state": "found", "catalog_run": "holo-m3-test-phase-001", "matching": 1, "omitted": 0,
                "records": [{"status": "REVIEW", "system_id": "test-system", "chemical_type": "UNKNOWN",
                             "files": {"system_cif": "system.cif"}}]}}
            next(button for button in app.button if button.label == "刷新本页").click().run()
            self.assertFalse(list(app.exception))
            self.assertTrue(any("holo-m3-test-phase-001" in line.value for line in app.caption))
            self.assertTrue(any("待核" in str(table.value) and "UNKNOWN" in str(table.value) for table in app.dataframe))
            observations.side_effect = transport.QueryError("离线", "目录断开")
            next(button for button in app.button if button.label == "刷新本页").click().run()
            self.assertTrue(any("上次缓存" in line.value and "目录断开" in line.value for line in app.caption))
            app.radio(key="nav").set_value("batches").run()
            before = read.call_count
            read.side_effect = transport.QueryError("离线", "测试断开")
            next(button for button in app.button if button.label == "刷新本页").click().run()
            self.assertGreater(read.call_count, before)
            self.assertTrue(any("上次缓存" in line.value for line in app.caption))
            self.assertFalse(list(app.exception))
        with patch.object(transport, "load_settings", side_effect=ValueError("测试配置异常")):
            app = AppTest.from_file(str(path), default_timeout=12).run()
            self.assertTrue(any("本地验收回执仍可用" in msg.value for msg in app.warning))
            app.radio(key="nav").set_value("m1").run()
            self.assertFalse(list(app.exception))


if __name__ == "__main__":
    unittest.main()
