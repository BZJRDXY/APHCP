"""Run or resume one APHCP M1 metadata inventory."""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone

from m1_holo_config import (metadata_complete, read_config, read_json, reusable_window,
                           search_complete, sha256, window, write_json)


def save_input(source, destination):
    if destination.exists():
        if source.read_bytes() != destination.read_bytes():
            raise ValueError(f"Frozen input changed: {destination}. Use a new run directory.")
        return
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    shutil.copyfile(source, temporary)
    temporary.replace(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--config", type=Path, help="Inventory query JSON for a new run")
    parser.add_argument("--selection", type=Path, help="Reviewed M2 selection JSON")
    parser.add_argument("--register", action="store_true", help="Record the finished batch and runs in SQLite")
    args = parser.parse_args()
    run = args.run_directory.resolve()
    inputs = run / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    if args.config:
        save_input(args.config, inputs / "query.json")
    cfg = read_config(run)
    start, end, _ = window(cfg)
    receipt_path = run / "run_receipt.json"
    progress_path = run / "progress.json"
    search_path = run / "search" / "summary.json"
    receipt = read_json(receipt_path) or {}
    if args.selection:
        if reusable_window(run) and not (inputs / "m2_selection.json").is_file():
            parser.error("Inventory is finalized; use a separate nomination run.")
        save_input(args.selection, inputs / "m2_selection.json")

    def execute(stage, *options):
        print(f"M1 {stage} {' '.join(options)}".rstrip(), flush=True)
        subprocess.run([sys.executable, "-B", str(code / f"m1_holo_{stage}.py"), str(run), *options], check=True)

    search = read_json(search_path) or {}
    # Completed historical runs can be reused without requiring identical current code.
    def complete_empty():
        if not receipt:
            write_json(receipt_path, {"status":"INVENTORY_COMPLETE", "batch_id":cfg["batch_id"],
                "run_directory":str(run), "window":{"start_inclusive":start.isoformat(),"end_exclusive":end.isoformat()},
                "empty_window":True,"search_summary_sha256":sha256(search_path),"coordinate_downloads":0})
        progress_path.unlink(missing_ok=True)
        print(f"No entries in this window: {search_path}", flush=True)
        return 0

    if search_complete(search) and search["unique_ids"] == 0:
        return complete_empty()
    if reusable_window(run):
        if args.register and not receipt.get("registry"):
            code = inputs / "resume-v2" if (inputs / "resume-v2" / "m1_holo_receipt.py").is_file() else inputs
            execute("receipt")
        print(f"Inventory complete: {receipt_path} (reused)", flush=True)
        return 0

    # Keep original stage code and any incorrectly finalized legacy receipt as evidence.
    code = inputs / "resume-v2"
    code.mkdir(exist_ok=True)
    source_dir = Path(__file__).resolve().parent
    for stage in ("config", "inventory", "fetch", "join", "budget", "receipt"):
        name = f"m1_holo_{stage}.py"
        save_input(source_dir / name, code / name)
    legacy = run / "superseded_run_receipt.json"
    if receipt and not legacy.exists():
        receipt_path.replace(legacy)

    def partial(stage, **issues):
        write_json(progress_path, {"status":"PARTIAL", "stage":stage,"batch_id":cfg["batch_id"],
            "updated_utc":datetime.now(timezone.utc).isoformat(), **issues})
        print(f"M1 partial ({stage}): {progress_path}", flush=True)
        return 2

    if not search_complete(search):
        try:
            execute("fetch", "search")
        except subprocess.CalledProcessError:
            search = read_json(search_path) or {}
            return partial("search", refresh_offsets=search.get("refresh_offsets", []))
        search = read_json(search_path)
    if not search["entry_ids"]:
        return complete_empty()

    data_path = run / "data" / "summary.json"
    entries_path = run / "data" / "entries.json"
    data = read_json(data_path) or {}
    if not metadata_complete(search, data) or not entries_path.is_file() or (
            data.get("entries_sha256") and data["entries_sha256"] != sha256(entries_path)):
        try:
            execute("fetch", "metadata")
        except subprocess.CalledProcessError:
            return partial("metadata", error="Metadata stage did not finish; see request receipts")
        data = read_json(data_path)
    if not metadata_complete(search, data):
        return partial("metadata", missing_ids=data.get("missing_ids", []))

    # A new metadata snapshot invalidates only the dependent stages. Last written
    # markers keep a stopped run resumable even if it died during regeneration.
    join_path = run / "inventory" / "inventory_summary.json"
    join = read_json(join_path) or {}
    entries_hash = sha256(entries_path)
    period_path = run / "inventory" / "plinder_period_only_reconcile.json"
    budget_path = run / "budget" / "budget_summary.json"
    if join.get("entry_metadata_sha256") != entries_hash:
        budget_path.unlink(missing_ok=True)
        period_path.unlink(missing_ok=True)
        execute("join")
    period = read_json(period_path)
    if period is None or any(r.get("_technical_missing") for r in period):
        execute("fetch", "period")
        period = read_json(period_path)
    unresolved = [r["entry_pdb_id"].upper() for r in period if r.get("_technical_missing")]
    if unresolved:
        return partial("period", period_missing_ids=unresolved)

    selection_supplied = (inputs / "m2_selection.json").is_file()
    budget = read_json(budget_path)
    if budget is None or (selection_supplied and not budget["m2_nomination_count"]):
        execute("budget")
    execute("receipt", *(() if args.register else ("--no-register",)))
    progress_path.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
