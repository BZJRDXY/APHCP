"""Export bounded Entry/System evidence from a completed Holo pilot."""

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import gemmi


COMPONENT_CHECKS = {
    "2Q9S": ("EIC",), "1LN1": ("DLP",), "1CCG": ("IMD", "HEM"),
    "4CVL": ("ACP",), "5I6Z": ("NAG",), "7MFF": ("215",),
    "4PYP": ("BNG",), "2CFQ": ("HG",), "2H8I": ("PEG",),
    "1FE3": ("OLA",), "1URN": ("GOL",), "1JQA": ("GOL",), "8RH4": ("EST",),
    "2OAL": ("FAD",), "8UIV": ("FAD",), "1XUC": ("PB3",),
}
NEGATIVE_COMPONENT_SETS = {
    "1BRS": set(), "3NCY": set(), "6V3G": set(), "2CFQ": {"HG"}, "4XRK": set(),
}


def write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def source_evidence(path, pdb, inspect_instances):
    block = gemmi.cif.read(str(path)).sole_block()
    structure = gemmi.make_structure_from_block(block)
    entities = block.get_mmcif_category("_entity.")
    types = block.get_mmcif_category("_entity_poly.")
    entity_types = dict(zip(types.get("entity_id", []), types.get("type", [])))
    polymers = " / ".join(f"{entity_id}: {entity_types[entity_id]} ({name})"
                          for entity_id, name in zip(entities.get("id", []), entities.get("pdbx_description", []))
                          if entity_id in entity_types)
    seen = Counter()
    selected = []
    focus = COMPONENT_CHECKS.get(pdb, ())
    for chain in structure[0]:
        for residue in chain:
            if residue.is_water() or residue.entity_type == gemmi.EntityType.Polymer:
                continue
            seen[residue.name] += 1
            if (residue.name in focus or inspect_instances) and len(selected) < 60:
                heavy = sum(not atom.element.is_hydrogen for atom in residue)
                selected.append(f"{residue.name}@{residue.subchain}({chain.name}:{residue.seqid};{heavy} heavy)")
    assemblies = ";".join(a.name for a in structure.assemblies)
    return seen, polymers, assemblies, selected


def counts_text(counts):
    return "; ".join(f"{key}={value}" for key, value in sorted(counts.items()))


def failure_group(system):
    text = system.get("failure", "")
    for prefix, group in (
        ("Ambiguous residue correspondence", "residue_identity"),
        ("Residue count differs", "residue_count"),
        ("Ambiguous atom correspondence", "atom_identity_or_xyz"),
        ("Missing atom without", "missing_atom"),
        ("current verifier requires one", "multiple_models"),
        ("system chain instances differ", "chain_membership"),
        ("source assembly lacks instance", "assembly_instance"),
    ):
        if text.startswith(prefix):
            return group
    if system["coordinate_validation_status"] == "component_identity_mismatch":
        return "component_identity"
    if system.get("failure_stage"):
        return f"{system['failure_stage']}:{system.get('failure_type', 'unknown')}"
    return ""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--display-candidates", type=Path, required=True)
    parser.add_argument("--revision-recheck", type=Path)
    args = parser.parse_args()
    run = json.loads((args.root / "summary.json").read_text(encoding="utf-8"))
    if run["status"] != "completed" or run["entry_files_present"] != 100:
        raise SystemExit("run is not yet complete")
    entries = sorted((json.loads(path.read_text(encoding="utf-8"))
                      for path in (args.root / "entries").glob("*.json")),
                     key=lambda row: int(row["candidate_no"]))
    if len(entries) != 100 or len({row["pdb_id"] for row in entries}) != 100:
        raise SystemExit("expected 100 unique Entry records")
    recheck = None
    if args.revision_recheck:
        recheck = json.loads(args.revision_recheck.read_text(encoding="utf-8"))
        row_no = next(i for i, row in enumerate(entries) if row["pdb_id"] == recheck["pdb_id"])
        prior = entries[row_no]
        if (prior["candidate_outcome"] != "source_parse_or_revision_failure"
            or recheck["source_sha256"] != prior["source_sha256"]
            or recheck["source_parse_status"] != "pass"):
            raise SystemExit("revision recheck does not resolve the original source receipt")
        entries[row_no] = recheck
    with args.display_candidates.open(encoding="utf-8-sig", newline="") as handle:
        display = {row["pdb_id"]: row for row in csv.DictReader(handle)}
    for entry in entries:
        current = display[entry["pdb_id"]]
        for key in ("candidate_no", "primary_stratum", "control_role", "control_role_rationale"):
            if key in entry and entry[key] != current[key]:
                raise SystemExit(f"report roster changed {entry['pdb_id']} {key}")

    entry_rows, system_rows = [], []
    for entry in entries:
        pdb = entry["pdb_id"]
        corrected = recheck is not None and pdb == recheck["pdb_id"]
        evidence_root = args.revision_recheck.parent.parent if corrected else args.root
        receipt_prefix = f"../{evidence_root.name}/" if corrected else ""
        stages = {row["stage"]: row["status"] for row in entry.get("stages", [])}
        inventory = entry.get("source_inventory", {})
        systems = entry.get("system_results", [])
        raw = entry.get("source_relative_path", "")
        role = display[pdb]["control_role"]
        no_annotation = entry.get("plinder_annotation_status") == "absent"
        source_recheck, seen, polymers, assemblies, selected = "not_attempted", {}, "", "", []
        role_result, reason = "", ""
        if no_annotation or role or pdb in COMPONENT_CHECKS:
            if stages.get("gemmi_source_parse_and_revision") == "pass":
                try:
                    seen, polymers, assemblies, selected = source_evidence(evidence_root / raw, pdb, no_annotation)
                    source_recheck = "pass"
                except Exception as exc:
                    source_recheck = f"failed: {type(exc).__name__}: {exc}"
        if pdb in NEGATIVE_COMPONENT_SETS and source_recheck == "pass":
            if set(seen) == NEGATIVE_COMPONENT_SETS[pdb]:
                role_result = "expected_exclusion_for_declared_target"
                reason = "First model has only the expected non-target components, or none. Control is limited to the declared small-molecule/substrate scope."
        if not role_result and role == "\u6b63\u5411\u8fb9\u754c" and source_recheck == "pass":
            target = {"1FE3": "OLA", "8RH4": "EST", "2OAL": "FAD", "8UIV": "FAD", "1XUC": "PB3", "1JQA": "GOL"}[pdb]
            if seen.get(target):
                role_result = "expected_positive_component_present"
                reason = f"Fixed source model 1 contains {seen[target]} {target} instance(s); biological role derives from the declared control."
        if role and not role_result:
            role_result = "role_pending"
            reason = "Component coordinates or index presence alone do not establish the declared biological role."

        reasons = []
        for s in systems:
            verification = s.get("coordinate_validation", {})
            mismatches = [f"{a['ligand_id']}: CCD={a['ligand_ccd_code']}, source={'-'.join(a['source_component_sequence'])}"
                          for a in verification.get("target_annotations", []) if a.get("component_match") is False]
            mismatches.extend(d["kind"] for d in verification.get("differences", []))
            problem = s.get("failure") or " | ".join(mismatches)
            group = failure_group(s)
            if group:
                reasons.append(group)
            system_rows.append({
                "pdb_id": pdb, "system_id": s["system_id"],
                "annotation_rows": len(s["annotation_rows"]),
                "nonproper_rows": sum(not a["ligand_is_proper"] for a in s["annotation_rows"]),
                "ligand_annotations": "; ".join(f"{a['ligand_id']}={a['ligand_ccd_code']}"
                                                for a in s["annotation_rows"]),
                "reconstruction": s["reconstruction_status"], "readback": s["readback_status"],
                "coordinate_comparison": s["coordinate_validation_status"], "problem_group": group,
                "problem": problem, "matched_heavy_atoms": verification.get("matched_heavy_atoms", verification.get("matched_heavy_atom_count", "")),
                "omitted_positive_altloc_atoms": verification.get("omitted_source_altloc_atoms", verification.get("omitted_source_altloc_atom_count", "")),
                "changed_occupancies": verification.get("changed_occupancies", verification.get("occupancy_change_count", "")),
                "derivative_path": receipt_prefix + s["derivative_relative_path"] if s.get("derivative_relative_path") else "",
                "derivative_sha256": s.get("derivative_sha256", ""),
            })
        entry_rows.append({
            "candidate_no": entry["candidate_no"], "pdb_id": pdb, "stratum": display[pdb]["primary_stratum"],
            "selection_rationale": display[pdb]["selection_rationale"], "control_role": role,
            "control_role_rationale": display[pdb]["control_role_rationale"], "control_result": role_result,
            "control_reason": reason, "fixed_revision": entry.get("source_revision", ""),
            "download_integrity": stages.get("download_and_gzip_check", "not_attempted"),
            "source_parse_revision": stages.get("gemmi_source_parse_and_revision", "not_attempted"),
            "original_bytes": entry.get("source_bytes", ""), "original_sha256": entry.get("source_sha256", ""),
            "source_models": inventory.get("model_count", ""),
            "source_assemblies": assemblies or ";".join(inventory.get("assembly_ids", [])),
            "source_nonwater_nonpolymer_counts_all_models": counts_text(inventory.get("nonwater_nonpolymer_residue_counts", {})),
            "source_native_recheck": source_recheck, "source_polymers": polymers,
            "source_native_instances_first_model": "; ".join(selected),
            "plinder_annotation": entry.get("plinder_annotation_status", "not_attempted"),
            "plinder_annotation_rows": entry.get("plinder_annotation_row_count", ""),
            "plinder_system_count": entry.get("plinder_system_count", ""),
            "plinder_revision_match": entry.get("plinder_revision_match", ""),
            "systems_attempted": len(systems),
            "systems_generated": sum(s["reconstruction_status"] == "generated" for s in systems),
            "systems_readback": sum(s["readback_status"] == "pass" for s in systems),
            "systems_coordinate_pass": sum(s["coordinate_validation_status"] == "coordinate_identity_pass" for s in systems),
            "problem_groups": counts_text(Counter(reasons)),
            "entry_outcome": entry["candidate_outcome"], "stage_errors": entry.get("failure") or " | ".join(
                s["error"] for s in entry.get("stages", []) if s.get("error")),
            "entry_receipt": f"{receipt_prefix}entries/{pdb.lower()}.json",
        })

    args.output.mkdir(exist_ok=False)
    write_csv(args.output / "entry_results.csv", entry_rows)
    write_csv(args.output / "system_results.csv", system_rows)
    run.update({
        "baseline_candidate_outcome_counts": run["candidate_outcome_counts"],
        "baseline_system_coordinate_status_counts": run["system_coordinate_status_counts"],
        "candidate_outcome_counts": dict(Counter(row["entry_outcome"] for row in entry_rows)),
        "system_coordinate_status_counts": dict(Counter(row["coordinate_comparison"] for row in system_rows)),
        "revision_rechecked_entry": recheck["pdb_id"] if recheck else None,
        "display_candidate_sha256": hashlib.sha256(args.display_candidates.read_bytes()).hexdigest(),
        "strata": dict(Counter(row["stratum"] for row in entry_rows)),
        "download_integrity": dict(Counter(row["download_integrity"] for row in entry_rows)),
        "source_parse_revision": dict(Counter(row["source_parse_revision"] for row in entry_rows)),
        "control_results": dict(Counter(row["control_result"] for row in entry_rows if row["control_result"])),
        "source_native_rechecks": dict(Counter(row["source_native_recheck"] for row in entry_rows)),
        "compressed_original_bytes": sum(int(row["original_bytes"] or 0) for row in entry_rows),
        "systems_attempted": len(system_rows),
        "systems_generated": sum(row["systems_generated"] for row in entry_rows),
        "systems_readback": sum(row["systems_readback"] for row in entry_rows),
        "nonproper_rows_checked": sum(row["nonproper_rows"] for row in system_rows),
        "system_problem_groups": dict(Counter(row["problem_group"] for row in system_rows if row["problem_group"])),
    })
    (args.output / "receipt.json").write_text(json.dumps(run, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(run, indent=2))


if __name__ == "__main__":
    main()
