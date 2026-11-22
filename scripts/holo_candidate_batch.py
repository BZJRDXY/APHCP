"""Run an isolated, resumable PLINDER/source validation over a candidate CSV."""

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
import time
import traceback
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import gemmi
from plinder.core import PlinderSystem, query_table
from plinder.core.release import PlinderRelease

import holo_pilot
from holo_native_candidates import assemble_candidate_ledger
from holo_source_mapping import MappingError, verify_system


def now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")
    temp.replace(path)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def revision_from_cif(path, entry_id):
    block = gemmi.cif.read(str(path)).sole_block()
    if block.find_value("_entry.id").lower() != entry_id.lower():
        raise ValueError("source entry ID differs from manifest")
    rows = block.find("_pdbx_audit_revision_history.",
                      ["data_content_type", "major_revision", "minor_revision"])
    revisions = [(int(row[1]), int(row[2])) for row in rows
                 if gemmi.cif.as_string(row[0]).casefold() == "structure model"]
    if not revisions:
        raise ValueError("source has no Structure model revision history")
    return block, ".".join(map(str, max(revisions)))


def source_inventory(block, structure):
    model_count = len(structure)
    assembly_ids = [assembly.name for assembly in structure.assemblies]
    nonpolymer = Counter()
    polymer_residues = Counter()
    chains = set()
    covalent_connections = 0
    for model in structure:
        for chain in model:
            chains.add(chain.name)
            for residue in chain:
                if residue.is_water():
                    continue
                if residue.entity_type == gemmi.EntityType.Polymer:
                    polymer_residues[residue.name] += 1
                else:
                    nonpolymer[residue.name] += 1
    conn = block.get_mmcif_category("_struct_conn.")
    covalent_connections = sum(kind == "covale" for kind in conn.get("conn_type_id", []))
    return {
        "model_count": model_count,
        "assembly_ids": assembly_ids,
        "asymmetric_unit_chain_count": len(chains),
        "nonwater_nonpolymer_residue_counts": dict(sorted(nonpolymer.items())),
        "polymer_residue_name_counts": dict(sorted(polymer_residues.items())),
        "source_covalent_connection_count": covalent_connections,
    }


def parse_system_id(system_id):
    parts = system_id.split("__")
    if len(parts) != 4:
        raise ValueError(f"unrecognized PLINDER system ID: {system_id}")
    assembly = parts[1]
    instance_pattern = r"[0-9]+\.[A-Za-z0-9]+"
    proteins = re.findall(instance_pattern, parts[2])
    ligands = re.findall(instance_pattern, parts[3])
    if not proteins or not ligands:
        raise ValueError(f"missing chain instances in PLINDER system ID: {system_id}")
    return assembly, tuple(proteins), tuple(ligands)


def verify_system_source(structure, source_block, derivative_block, system_id, annotation_rows):
    return verify_system(structure, source_block, derivative_block, *parse_system_id(system_id), annotation_rows)


def download_source(row, raw_dir):
    pdb = row["pdb_id"].lower()
    revision = row["revision"]
    path = raw_dir / f"{pdb}-v{revision}.cif.gz"
    if path.exists():
        if row.get("source_sha256") and sha256(path) == row["source_sha256"]:
            return path, "reused_verified", sha256(path), path.stat().st_size
        path.unlink()
    temp = path.with_suffix(path.suffix + ".part")
    error = None
    for attempt in range(3):
        try:
            request = urllib.request.Request(row["source_url"], headers={"User-Agent": "APHCP-HoloCandidateValidation/1.0"})
            with urllib.request.urlopen(request, timeout=20) as response, temp.open("wb") as output:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    output.write(chunk)
            with gzip.open(temp, "rb") as compressed:
                while compressed.read(1024 * 1024):
                    pass
            if row.get("expected_bytes") and temp.stat().st_size != int(row["expected_bytes"]):
                raise ValueError(f"download size mismatch: {temp.stat().st_size} != {row['expected_bytes']}")
            temp.replace(path)
            return path, "downloaded", sha256(path), path.stat().st_size
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"
            temp.unlink(missing_ok=True)
            time.sleep(1 + attempt)
    raise RuntimeError(error or "download failed")


def process_entry(candidate, source_meta, root, release, code_commit, *, role_rules=None):
    started = now()
    started_clock = time.monotonic()
    pdb = candidate["pdb_id"].lower()
    output_path = root / "entries" / f"{pdb}.json"
    existing = None
    if output_path.exists():
        try:
            existing = json.loads(output_path.read_text(encoding="utf-8"))
            if existing.get("run_status") == "completed":
                return existing
        except Exception:
            pass
    result = {
        "pdb_id": pdb.upper(),
        "candidate_no": candidate["candidate_no"],
        "primary_stratum": candidate["primary_stratum"],
        "selection_rationale": candidate["selection_rationale"],
        "control_role": candidate.get("control_role", ""),
        "control_role_rationale": candidate.get("control_role_rationale", ""),
        "run_started_at_utc": started,
        "run_status": "running",
        "source_revision": source_meta["revision"],
        "source_url": source_meta["source_url"],
        "source_file_status": "not_attempted",
        "stages": [],
        "system_results": [],
    }
    atomic_json(output_path, result)
    source = None
    block = None
    structure = None
    annotations = []
    try:
        source, download_state, digest, byte_count = download_source(source_meta, root / "raw")
        result.update({"source_file_status": download_state, "source_sha256": digest, "source_bytes": byte_count,
                       "source_relative_path": source.relative_to(root).as_posix()})
        result["stages"].append({"stage": "download_and_gzip_check", "status": "pass", "completed_at_utc": now()})
    except Exception as exc:
        result["source_file_status"] = "failed"
        result["candidate_outcome"] = "source_download_failure"
        result["stages"].append({"stage": "download_and_gzip_check", "status": "technical_failure",
                                 "error_type": type(exc).__name__, "error": str(exc)})
        result["run_status"] = "completed"
        result["run_finished_at_utc"] = now()
        result["elapsed_seconds"] = round(time.monotonic() - started_clock, 3)
        atomic_json(output_path, result)
        return result

    try:
        block, parsed_revision = revision_from_cif(source, pdb)
        structure = gemmi.read_structure(str(source))
        inventory = source_inventory(block, structure)
        if not inventory["model_count"]:
            raise ValueError("fixed mmCIF has no coordinate model")
        result["source_inventory"] = inventory
        result["parsed_revision"] = parsed_revision
        result["source_revision_match"] = parsed_revision == source_meta["revision"]
        if not result["source_revision_match"]:
            raise ValueError(f"fixed revision mismatch: expected {source_meta['revision']} got {parsed_revision}")
        result["source_parse_status"] = "pass"
        result["stages"].append({"stage": "gemmi_source_parse_and_revision", "status": "pass",
                                 "completed_at_utc": now()})
    except Exception as exc:
        result["source_parse_status"] = "technical_failure"
        result["candidate_outcome"] = "source_parse_or_revision_failure"
        result["stages"].append({"stage": "gemmi_source_parse_and_revision", "status": "technical_failure",
                                 "error_type": type(exc).__name__, "error": str(exc)})
        result["run_status"] = "completed"
        result["run_finished_at_utc"] = now()
        result["elapsed_seconds"] = round(time.monotonic() - started_clock, 3)
        atomic_json(output_path, result)
        return result

    try:
        annotations = query_table("annotation", release=release,
                                  filters=[("entry_pdb_id", "==", pdb)],
                                  columns=["system_id", "ligand_id", "ligand_ccd_code", "ligand_is_proper"])
        annotation_rows = annotations.to_dict(orient="records")
        entry_sources = query_table("entry_sources", release=release,
                                    filters=[("entry_pdb_id", "==", pdb)])
        revisions = {f"{row.source_mmcif_major_revision}.{row.source_mmcif_minor_revision}"
                     for row in entry_sources.itertuples()}
        plinder_revision = next(iter(revisions)) if len(revisions) == 1 else None
        result["plinder_annotation_row_count"] = len(annotation_rows)
        result["plinder_system_count"] = len({row["system_id"] for row in annotation_rows})
        result["plinder_proper_row_count"] = sum(bool(row["ligand_is_proper"]) for row in annotation_rows)
        result["plinder_source_revision"] = plinder_revision
        result["plinder_source_revisions"] = sorted(revisions)
        result["plinder_revision_match"] = plinder_revision == source_meta["revision"]
        result["plinder_candidate_ligands"] = annotation_rows
        result["plinder_annotation_status"] = "present" if annotation_rows else "absent"
        result["stages"].append({"stage": "plinder_annotation_lookup", "status": "pass",
                                 "completed_at_utc": now()})
    except Exception as exc:
        result["plinder_annotation_status"] = "technical_failure"
        result["candidate_outcome"] = "plinder_annotation_lookup_failure"
        result["candidate_ledger"] = assemble_candidate_ledger(
            block, structure, pdb.upper(), source_meta["revision"], [], role_rules or {}, direct=True)
        result["stages"].append({"stage": "plinder_annotation_lookup", "status": "technical_failure",
                                 "error_type": type(exc).__name__, "error": str(exc)})
        result["run_status"] = "completed"
        result["run_finished_at_utc"] = now()
        result["elapsed_seconds"] = round(time.monotonic() - started_clock, 3)
        atomic_json(output_path, result)
        return result

    proper_rows = [row for row in annotation_rows if bool(row["ligand_is_proper"])]
    if not annotation_rows:
        result["plinder_case"] = "no_annotation_rows"
    elif not proper_rows:
        result["plinder_case"] = "annotation_present_no_proper_ligand"
    else:
        result["plinder_case"] = "proper_ligand_annotations_present"

    rows_by_system = {}
    for row in annotation_rows if result["plinder_revision_match"] else []:
        rows_by_system.setdefault(row["system_id"], []).append(row)
    for system_id, system_rows in rows_by_system.items():
        key = hashlib.sha1(system_id.encode("utf-8")).hexdigest()[:12]
        system_result = {"system_id": system_id, "annotation_rows": system_rows,
                         "reconstruction_status": "not_attempted", "readback_status": "not_attempted",
                         "coordinate_validation_status": "not_attempted"}
        output_dir = root / "systems" / pdb / key
        stage = "reconstruction_status"
        try:
            plinder_system = PlinderSystem(system_id=system_id, release=release, source_mmcif=source,
                                           reconstruction_dir=output_dir)
            derivative = Path(plinder_system.system_cif)
            if not derivative.is_file():
                raise FileNotFoundError(f"system CIF missing: {derivative}")
            system_result["derivative_relative_path"] = derivative.relative_to(root).as_posix()
            system_result["derivative_bytes"] = derivative.stat().st_size
            system_result["derivative_sha256"] = sha256(derivative)
            system_result["reconstruction_status"] = "generated"
            stage = "readback_status"
            derivative_block = gemmi.cif.read(str(derivative)).sole_block()
            reread = gemmi.make_structure_from_block(derivative_block)
            if not len(reread) or not reread[0].count_atom_sites():
                raise ValueError("reconstructed CIF has no coordinate atoms")
            system_result["readback_status"] = "pass"
            stage = "coordinate_validation_status"
            system_result["coordinate_validation"] = verify_system_source(
                structure, block, derivative_block, system_id, system_rows)
            system_result["coordinate_validation_status"] = system_result["coordinate_validation"]["status"]
        except Exception as exc:
            system_result[stage] = "failed"
            system_result["failure_stage"] = stage
            system_result["failure_type"] = type(exc).__name__
            system_result["failure"] = str(exc)
            if isinstance(exc, MappingError):
                system_result["mapping_error"] = {"kind": exc.kind, "evidence": exc.evidence}
            system_result["traceback_tail"] = traceback.format_exc().splitlines()[-8:]
        result["system_results"].append(system_result)
        atomic_json(output_path, result)

    result["run_status"] = "completed"
    result["run_finished_at_utc"] = now()
    result["elapsed_seconds"] = round(time.monotonic() - started_clock, 3)
    statuses = Counter(item["coordinate_validation_status"] for item in result["system_results"])
    result["system_status_counts"] = dict(statuses)
    result["system_attempt_count"] = len(result["system_results"])
    candidates_for_ledger = [(s["system_id"], parse_system_id(s["system_id"])[0], s["annotation_rows"],
                              s.get("coordinate_validation")) for s in result["system_results"]]
    result["candidate_ledger"] = assemble_candidate_ledger(
        block, structure, pdb.upper(), source_meta["revision"], candidates_for_ledger, role_rules or {},
        direct=not annotation_rows or not result["plinder_revision_match"])
    if result["plinder_case"] == "no_annotation_rows":
        result["candidate_outcome"] = "source_parsed_plinder_no_annotation"
    elif not result["plinder_revision_match"]:
        result["candidate_outcome"] = "source_parsed_plinder_revision_blocked"
    elif statuses.get("coordinate_identity_pass") == len(result["system_results"]):
        result["candidate_outcome"] = "all_indexed_system_coordinates_validated"
    elif statuses.get("coordinate_identity_pass"):
        result["candidate_outcome"] = "partial_indexed_system_coordinates_validated"
    elif result["system_results"]:
        result["candidate_outcome"] = "system_reconstruction_or_coordinate_failures"
    else:
        result["candidate_outcome"] = "source_parsed_without_reconstructable_system"
    atomic_json(output_path, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--plinder-release", type=Path, required=True)
    parser.add_argument("--git-commit", required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--role-evidence", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    if root.exists() and not args.resume and any(root.iterdir()):
        raise SystemExit(f"run root is not empty (pass --resume): {root}")
    (root / "raw").mkdir(parents=True, exist_ok=True)
    (root / "entries").mkdir(parents=True, exist_ok=True)
    with args.candidates.open(encoding="utf-8-sig", newline="") as stream:
        candidates = list(csv.DictReader(stream))
    with args.source_manifest.open(encoding="utf-8-sig", newline="") as stream:
        sources = {row["pdb_id"]: row for row in csv.DictReader(stream)}
    roles = json.loads(args.role_evidence.read_text(encoding="utf-8"))["entries"] if args.role_evidence else {}
    if len(candidates) != 100 or len({row["pdb_id"] for row in candidates}) != 100:
        raise SystemExit("candidate manifest must contain exactly 100 unique entries")
    if set(sources) != {row["pdb_id"] for row in candidates}:
        raise SystemExit("source manifest IDs differ from candidate manifest")
    metadata = {
        "run_id": root.name,
        "scope": "100 real candidate-entry validations; four fixed regressions tracked separately",
        "started_at_utc": now(),
        "git_commit": args.git_commit,
        "python": sys.version.split()[0],
        "gemmi": gemmi.__version__,
        "plinder": __import__("importlib.metadata", fromlist=["version"]).version("plinder"),
        "plinder_release": str(args.plinder_release),
        "candidate_entries": len(candidates),
        "resume_enabled": args.resume,
        "input_sha256": {"candidates": sha256(args.candidates), "source_manifest": sha256(args.source_manifest),
                         "batch_script": sha256(Path(__file__)), "pilot_script": sha256(Path(holo_pilot.__file__)),
                         "mapping": sha256(Path(__file__).with_name("holo_source_mapping.py")),
                         "native": sha256(Path(__file__).with_name("holo_native_candidates.py")),
                         "roles": sha256(args.role_evidence) if args.role_evidence else None},
    }
    if args.resume:
        previous = json.loads((root / "run.json").read_text(encoding="utf-8"))
        for key in ("git_commit", "python", "gemmi", "plinder", "plinder_release", "input_sha256"):
            if previous.get(key) != metadata[key]:
                raise SystemExit(f"resume identity changed: {key}")
        metadata["started_at_utc"] = previous["started_at_utc"]
    atomic_json(root / "run.json", metadata | {"status": "running"})
    release = PlinderRelease(data_dir=args.plinder_release)
    selected = candidates[:args.limit] if args.limit else candidates
    results = []
    for index, candidate in enumerate(selected, 1):
        print(json.dumps({"event": "entry_started", "index": index, "total_selected": len(selected),
                          "pdb_id": candidate["pdb_id"]}), flush=True)
        try:
            result = process_entry(candidate, sources[candidate["pdb_id"]], root, release, args.git_commit,
                                   role_rules=roles.get(candidate["pdb_id"], {}))
        except Exception as exc:
            result = {"pdb_id": candidate["pdb_id"], "candidate_no": candidate["candidate_no"],
                      "run_status": "completed", "candidate_outcome": "unhandled_entry_failure",
                      "failure_type": type(exc).__name__, "failure": str(exc),
                      "traceback_tail": traceback.format_exc().splitlines()[-8:], "run_finished_at_utc": now()}
            atomic_json(root / "entries" / f"{candidate['pdb_id'].lower()}.json", result)
        results.append(result)
        print(json.dumps({"event": "entry_completed", "index": index, "pdb_id": candidate["pdb_id"],
                          "outcome": result.get("candidate_outcome"), "systems": result.get("system_attempt_count", 0),
                          "elapsed_seconds": result.get("elapsed_seconds")}), flush=True)
        atomic_json(root / "run_progress.json", {
            "updated_at_utc": now(),
            "processed_in_this_call": index,
            "selected_in_this_call": len(selected),
            "completed_entry_files": len(list((root / "entries").glob("*.json"))),
            "system_attempts_in_completed_files": sum(
                len(row.get("system_results", [])) for row in results
            ),
        })
    all_rows = [json.loads(path.read_text(encoding="utf-8")) for path in sorted((root / "entries").glob("*.json"))]
    summary = Counter(row.get("candidate_outcome", "unknown") for row in all_rows)
    system_summary = Counter()
    for row in all_rows:
        system_summary.update(row.get("system_status_counts", {}))
    finished = metadata | {"status": "completed" if len(all_rows) == len(candidates) else "partial",
                          "finished_at_utc": now(),
                          "selected_in_this_call": len(selected),
                          "entry_files_present": len(all_rows),
                          "candidate_outcome_counts": dict(summary),
                          "system_coordinate_status_counts": dict(system_summary)}
    atomic_json(root / "run.json", finished)
    atomic_json(root / "summary.json", finished)
    print(json.dumps({"event": "run_completed", "entry_files": len(all_rows),
                      "outcomes": dict(summary), "system_statuses": dict(system_summary)}), flush=True)


if __name__ == "__main__":
    main()
