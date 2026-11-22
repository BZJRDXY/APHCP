"""Recheck concrete assembly/residue mapping failures without replacing baseline results."""

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import gemmi

import holo_candidate_batch as batch
from holo_source_mapping import assembly_instances




def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--entries", nargs="*")
    args = parser.parse_args()
    rows = []
    for path in sorted((args.root / "entries").glob("*.json")):
        if args.entries and path.stem.upper() not in args.entries:
            continue
        entry = json.loads(path.read_text(encoding="utf-8"))
        if entry.get("run_status") != "completed":
            continue
        for s in entry.get("system_results", []):
            before = s.get("failure", "")
            if not before.startswith(("source assembly lacks instance", "Ambiguous atom correspondence",
                                      "Ambiguous residue correspondence")) and not (
                entry["pdb_id"] == "3HS4" and s.get("coordinate_validation_status") == "coordinate_identity_pass"
            ):
                continue
            if s.get("readback_status") != "pass" or entry.get("source_inventory", {}).get("model_count") != 1:
                continue
            source = args.root / entry["source_relative_path"]
            derivative = args.root / s["derivative_relative_path"]
            asm, proteins, ligands = batch.parse_system_id(s["system_id"])
            block = gemmi.cif.read(str(source)).sole_block()
            structure = gemmi.make_structure_from_block(block)
            instances = assembly_instances(structure, asm, structure[0])
            mapping = {name: instances[name][0].name for name in sorted(set(proteins + ligands)) if name in instances}

            row = {"pdb_id": entry["pdb_id"], "system_id": s["system_id"],
                   "baseline_status": s["coordinate_validation_status"], "baseline_error": before,
                   "assembly": asm, "copy_instance_to_mmcif_operation": mapping,
                   "derivative_sha256": s["derivative_sha256"]}
            try:
                if batch.sha256(derivative) != s["derivative_sha256"]:
                    raise ValueError("derivative checksum differs from baseline")
                result = batch.verify_system_source(structure, block, gemmi.cif.read(str(derivative)).sole_block(),
                                                    s["system_id"], s["annotation_rows"])
                row.update(status=result["status"], matched_heavy_atoms=result["matched_heavy_atoms"],
                           differences=result["differences"])
            except Exception as exc:
                row.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            rows.append(row)
    result = {"checked_at_utc": datetime.now(timezone.utc).isoformat(), "systems": rows,
              "counts": dict(Counter(r["status"] for r in rows)),
              "new_passes": sum(r["status"] == "coordinate_identity_pass" and r["baseline_error"] != "" for r in rows)}
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2)
    print(json.dumps({k: v for k, v in result.items() if k != "systems"}))


if __name__ == "__main__":
    main()
