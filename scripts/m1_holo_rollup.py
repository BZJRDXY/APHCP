"""Publish a single set of M1 tables from reconciled window runs."""

from collections import Counter
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from m1_holo_config import read_json, write_json


ENTRY_ROWS = "candidate_entries.jsonl"
LIGAND_ROWS = "plinder_ligand_rows.jsonl"
FILE_ROWS = "monthly_file_plan.csv"
BANDS = ("p10", "p50", "p90")
COUNT_FIELDS = (
    "plinder_annotation_entries", "plinder_no_annotation_entries",
    "plinder_annotation_rows", "plinder_candidate_systems", "rcsb_assembly_pairs",
)
COVERAGE_FIELDS = (
    "methods", "discovery_margins", "discovery_combinations", "family_entry_evidence",
    "protein_entity_evidence", "indexed_ligand_row_classes", "pli_binding_evidence",
    "all_source_revision_comparisons", "metadata_review",
)


def append_jsonl(source, target, expected):
    count = 0
    with source.open("rb") as stream:
        for line in stream:
            if not line.strip():
                continue
            target.write(line if line.endswith(b"\n") else line + b"\n")
            count += 1
    if count != expected:
        raise ValueError(f"{source}: expected {expected} rows, found {count}")


def load(path):
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def write_inventory(destination, report):
    if report["status"] != "INVENTORY_COMPLETE":
        raise ValueError("Inventory export requires a complete cross-window reconciliation")
    summary_path = destination / "summary.json"
    prior = read_json(summary_path)
    if prior and prior["plan_sha256"] == report["plan_sha256"]:
        if all((destination / name).is_file() for name in prior["tables"]):
            return str(summary_path)
    if destination.exists():
        raise ValueError(f"Choose a new inventory directory: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    coverage = {key: Counter() for key in COVERAGE_FIELDS}
    counts = Counter()
    budget = {tag: Counter() for tag in BANDS}
    file_count = 0
    release_ids = set()

    with TemporaryDirectory(prefix=".m1-rollup-", dir=destination.parent) as temporary:
        stage = Path(temporary) / "inventory"
        stage.mkdir()
        with (
            (stage / ENTRY_ROWS).open("wb") as entries,
            (stage / LIGAND_ROWS).open("wb") as ligands,
            (stage / FILE_ROWS).open("w", encoding="utf-8", newline="") as file_table,
        ):
            writer = None
            for window in report["windows"]:
                expected = window["search_total_count"]
                if not expected:
                    continue
                run = Path(window["run_directory"])
                cfg = load(run / "inputs/query.json")
                release_ids.add(cfg["plinder_release"])
                stats = load(run / "inventory/inventory_summary.json")["counts"]
                costs = load(run / "budget/budget_summary.json")
                data = load(run / "data/summary.json")
                counts["entries"] += expected
                counts["metadata_received"] += data["received"]
                counts["confirmed_unavailable"] += len(data.get("confirmed_unavailable", {}))
                for key in COUNT_FIELDS:
                    counts[key] += stats[key]
                for key in COVERAGE_FIELDS:
                    coverage[key].update(stats[key])
                append_jsonl(run / "inventory" / ENTRY_ROWS, entries, expected)
                append_jsonl(run / "inventory" / LIGAND_ROWS, ligands, stats["plinder_annotation_rows"])

                file_ids = set()
                scoped_ids = set(load(run / "search/summary.json")["entry_ids"])
                with (run / "budget" / FILE_ROWS).open(encoding="utf-8", newline="") as stream:
                    rows = csv.DictReader(stream)
                    for row in rows:
                        if writer is None:
                            writer = csv.DictWriter(file_table, fieldnames=rows.fieldnames)
                            writer.writeheader()
                        pid = row["pdb_id"]
                        if pid in file_ids or pid not in scoped_ids:
                            raise ValueError(f"Unexpected or duplicate planned PDB file for {pid}: {run}")
                        file_ids.add(pid)
                        writer.writerow(row)
                        for tag in BANDS:
                            budget[tag]["compressed_original_bytes"] += int(row[f"compressed_bytes_{tag}"])
                            budget[tag]["download_bytes"] += int(row[f"download_bytes_{tag}"])
                missing = costs["files_without_current_revision"]
                if len(file_ids) + missing != expected:
                    raise ValueError(f"Coordinate file plan does not reconcile: {run}")
                file_count += len(file_ids)
                counts["files_without_current_revision"] += missing
                for tag in BANDS:
                    scenario = costs["month_all_indexed_systems_scenario"]["bands"][tag]
                    budget[tag]["derived_cif_scenario_bytes"] += scenario["derived_cif"]
                    budget[tag]["largest_window_serial_scratch_bytes"] = max(
                        budget[tag]["largest_window_serial_scratch_bytes"], scenario["serial_scratch_peak"]
                    )

        if counts["entries"] != report["unique_ids"] or len(release_ids) > 1:
            raise ValueError("Merged inventory totals or PLINDER releases do not reconcile")
        summary = {
            "plan_sha256": report["plan_sha256"],
            "scope": report["scope"], "denominator": report["denominator"],
            "window_count": len(report["windows"]), "plinder_releases": sorted(release_ids),
            "counts": counts, "coverage": coverage, "unique_coordinate_files": file_count,
            "coordinate_budget": budget,
            "budget_note": "Sum of window planning scenarios. Scratch is the largest single serial window estimate.",
            "tables": [ENTRY_ROWS, LIGAND_ROWS, FILE_ROWS],
        }
        write_json(stage / "summary.json", summary)
        stage.rename(destination)
    return str(summary_path)
