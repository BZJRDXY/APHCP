"""Validate selected fixed-revision Holo samples and their PLINDER representative views."""

import argparse
import hashlib
import json
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import gemmi
import numpy as np
from plinder.core import PlinderSystem, query_table
from plinder.core.release import PlinderRelease


CASES = (
    ("1m17", "1.3", "m1", "AQ4", "1m17__1__1.A__1.B", ("1.A",), "1.B", ("AQ4",)),
    ("1b86", "2.3", "b8", "DG2", "1b86__1__1.B_1.D__1.K", ("1.B", "1.D"), "1.K", ("DG2",)),
    ("1hew", "2.1", "he", "NAG-NAG-NAG", "1hew__1__1.A__1.B", ("1.A",), "1.B", ("NAG",) * 3),
    ("3tct", "1.2", "tc", "3MI", "3tct__1__1.A_1.B_2.A_2.B__1.C", ("1.A", "1.B", "2.A", "2.B"), "1.C", ("3MI",)),
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def point(position):
    return tuple(round(value, 3) for value in (position.x, position.y, position.z))


def subchains(structure):
    if len(structure) != 1:
        raise ValueError("Pilot requires a single model")
    chains = defaultdict(list)
    for chain in structure[0]:
        for residue in chain:
            if not residue.is_water():
                chains[residue.subchain].append((chain.name, residue))
    return chains


def source_instances(structure, assembly_id):
    chains = subchains(structure)
    assembly = next(item for item in structure.assemblies if item.name == assembly_id)
    instances = {}
    for generator in assembly.generators:
        for operator in generator.operators:
            for label in generator.subchains:
                if label in chains:
                    instances[f"{operator.name}.{label}"] = (operator, label, chains[label])
    return instances


def match_residues(source, derived, *, single_ligand):
    """Match within a known instance before comparing atoms."""
    if len(source) != len(derived):
        raise ValueError("Residue count differs from source instance")
    matched = set()
    for _, target in derived:
        if single_ligand:
            choices = [0] if source[0][1].name == target.name else []
        else:
            choices = [
                i for i, (_, original) in enumerate(source) if original.name == target.name and (
                    original.label_seq == target.label_seq if target.label_seq is not None
                    else str(target.seqid) in (str(original.seqid), str(original.label_seq))
                )
            ]
        if len(choices) != 1 or choices[0] in matched:
            raise ValueError(f"Ambiguous residue correspondence for {target.subchain} {target.seqid}")
        matched.add(choices[0])
        yield source[choices[0]], target


def altloc(atom):
    return None if atom.altloc in ("\0", " ", ".", "?") else atom.altloc


def compare_residue(original, target, operator):
    """Resolve within the mapped residue and atom name; reject ambiguous coordinates."""
    source_atoms = {}
    for atom in original:
        if atom.element.is_hydrogen:
            continue
        key = (atom.name, altloc(atom))
        if key in source_atoms:
            raise ValueError(f"Duplicate source atom identity: {original.seqid} {key}")
        source_atoms[key] = (point(operator.transform.apply(atom.pos)), round(atom.occ, 3))

    matches = {}
    target_keys = set()
    for atom in target:
        if atom.element.is_hydrogen:
            continue
        target_key = (atom.name, altloc(atom))
        xyz = point(atom.pos)
        candidates = [key for key in source_atoms if key[0] == atom.name
                      and (target_key[1] is None or key[1] == target_key[1])
                      and source_atoms[key][0] == xyz]
        if len(candidates) != 1 or candidates[0] in matches or target_key in target_keys:
            raise ValueError(f"Ambiguous atom correspondence: {target.subchain} {target.seqid} {target_key}")
        target_keys.add(target_key)
        matches[candidates[0]] = (round(atom.occ, 3), target_key[1])

    kept_alts = {key for key in matches if source_atoms[key][1] > 0 and matches[key][0] > 0}
    for key, (_, occupancy) in source_atoms.items():
        if occupancy > 0 and key not in kept_alts:
            if key[1] is None or not any(other[0] == key[0] and other[1] not in (None, key[1])
                                         for other in kept_alts):
                raise ValueError(f"Missing atom without a retained alternate: {target.subchain} {target.seqid} {key}")
    return source_atoms, matches


def ligand_links(block, label, members):
    category = block.get_mmcif_category("_struct_conn.")
    links = []
    for i, kind in enumerate(category.get("conn_type_id", [])):
        if kind != "covale" or any(category[f"ptnr{side}_label_asym_id"][i] != label for side in (1, 2)):
            continue
        endpoints = []
        for side in (1, 2):
            prefix = f"ptnr{side}_"
            atom = category[f"{prefix}label_atom_id"][i]
            code = category[f"{prefix}label_comp_id"][i]
            seq = category[f"{prefix}auth_seq_id"][i]
            insertion = category.get(f"pdbx_{prefix}PDB_ins_code", [None] * (i + 1))[i]
            seqid = str(seq) + (str(insertion) if insertion else "")
            endpoint = [member for member in members
                        if member["component"] == code and member["source_author_seq_id"] == seqid
                        and any(item["atom"] == atom for item in member["atoms"])]
            if len(endpoint) != 1:
                raise ValueError(f"Unresolved ligand connection: {seqid} {code} {atom}")
            match = endpoint[0]
            endpoints.append({"source_author_seq_id": seqid, "component": code, "atom": atom,
                              "output_author_seq_id": match["output_author_seq_id"]})
        links.append({"source_struct_conn_id": category["id"][i], "endpoints": endpoints})
    if len(members) > 1 and not links:
        raise ValueError("Multi-residue pilot ligand has no mapped source connection")
    return links


def verify_coordinates(source, assembly, proteins, ligand, components, derivative):
    structure = gemmi.read_structure(str(source))
    instances = source_instances(structure, assembly)
    observed = subchains(gemmi.read_structure(str(derivative)))
    if set(observed) != {*proteins, ligand} or not set(observed) <= instances.keys():
        raise ValueError(f"Unexpected non-water instances in reconstruction: {sorted(observed)}")
    if tuple(res.name for _, res in instances[ligand][2]) != tuple(components):
        raise ValueError("Unexpected source ligand residue composition")

    protein_members = {name: (operator, label, residues) for name, (operator, label, residues) in instances.items()
                       if any(res.entity_type == gemmi.EntityType.Polymer for _, res in residues)}
    if not set(proteins) <= protein_members.keys():
        raise ValueError("Retained protein chains do not belong to the source assembly")
    members, ligand_xyz, per_instance = [], [], {}
    for name, target in sorted(observed.items()):
        operator, label, original = instances[name]
        evidence = {"label_asym_id": label, "operator": operator.name, "mapped_heavy_atoms": 0,
                    "omitted_alternate_positions": [], "occupancy_changes": []}
        single_ligand = name == ligand and len(components) == 1
        for (author, res), out_res in match_residues(original, target, single_ligand=single_ligand):
            atoms, matched = compare_residue(res, out_res, operator)
            evidence["mapped_heavy_atoms"] += len(matched)
            identities = []
            for (atom, alt), (xyz, occupancy) in atoms.items():
                normalized, output_alt = matched.get((atom, alt), (None, None))
                identity = {"component": res.name, "source_author_chain": author, "source_author_seq_id": str(res.seqid),
                            "source_label_seq_id": res.label_seq, "output_author_seq_id": str(out_res.seqid),
                            "atom": atom, "altloc": alt, "source_occupancy": occupancy,
                            "output_altloc": output_alt, "output_occupancy": normalized}
                if normalized is None and occupancy > 0:
                    evidence["omitted_alternate_positions"].append(identity)
                if normalized is not None and normalized != occupancy:
                    evidence["occupancy_changes"].append(identity)
                if name == ligand:
                    identities.append(identity)
                    if normalized and occupancy > 0:
                        ligand_xyz.append(xyz)
            if name == ligand:
                members.append({"component": res.name, "source_author_seq_id": str(res.seqid),
                                "output_author_seq_id": str(out_res.seqid), "atoms": identities})
        per_instance[name] = evidence

    for member in members:
        retained_alts = {atom["altloc"] for atom in member["atoms"]
                         if atom["altloc"] is not None and atom["output_occupancy"] is not None}
        if len(retained_alts) > 1:
            raise ValueError("Representative ligand residue mixes source alternate conformers")
    if not ligand_xyz:
        raise ValueError("Representative ligand has no matched observed heavy atoms")
    distances = {}
    for name, (operator, _, residues) in sorted(protein_members.items()):
        xyz = np.array([point(operator.transform.apply(atom.pos)) for _, res in residues for atom in res
                        if not atom.element.is_hydrogen and atom.occ > 0])
        distance = float(np.linalg.norm(xyz[:, None, :] - np.array(ligand_xyz)[None, :, :], axis=2).min())
        distances[name] = round(distance, 3)
    return {"model": structure[0].num, "assembly_member_chains": sorted(protein_members),
            "retained_protein_chains": sorted(proteins),
            "contact_rule": "pilot: source observed heavy atoms within 4 Å of selected representative ligand pose",
            "ligand_contact_chains": [name for name, value in distances.items() if value <= 4],
            "minimum_distance_A": distances, "instances": per_instance, "ligand_members": members,
            "ligand_links_from_source": ligand_links(gemmi.cif.read(str(source)).sole_block(),
                                                     instances[ligand][1], members)}


def sample(root, release, case):
    pdb, revision, shard, component, system_id, proteins, ligand, components = case
    source = root / "raw" / f"{pdb}-v{revision}.cif.gz"
    raw = gemmi.cif.read(str(source)).sole_block()
    revisions = raw.find(["_pdbx_audit_revision_history.major_revision",
                          "_pdbx_audit_revision_history.minor_revision"])
    seen = ".".join(revisions[len(revisions) - 1])
    if raw.find_value("_entry.id").lower() != pdb or seen != revision:
        raise ValueError(f"Source identity or revision mismatch: {source}")

    filters = [("entry_pdb_id", "==", pdb)]
    annotations = query_table("annotation", release=release, filters=filters,
                              columns=["system_id", "ligand_ccd_code", "ligand_id"])
    revision_row = query_table("entry_sources", release=release, filters=filters).iloc[0]
    indexed_revision = f"{revision_row.source_mmcif_major_revision}.{revision_row.source_mmcif_minor_revision}"
    if indexed_revision != revision:
        raise ValueError(f"PLINDER revision differs from source for {pdb}")
    ligand_id = f"{pdb}__1__{ligand}"
    matches = annotations[(annotations.system_id == system_id) & (annotations.ligand_id == ligand_id)
                          & (annotations.ligand_ccd_code == component)]
    if len(matches) != 1:
        raise ValueError(f"PLINDER target is missing or ambiguous: {system_id}, {ligand_id}")
    hem_labels = sorted({str(row[0]) for row in raw.find(["_atom_site.label_asym_id", "_atom_site.label_comp_id"])
                         if row[1] == "HEM"})
    if pdb == "1b86":
        indexed_hem = annotations.loc[annotations.ligand_ccd_code == "HEM", "ligand_id"]
        if len(hem_labels) != 4 or indexed_hem.nunique() != 4:
            raise ValueError("Repeated HEM instances did not match source counts")

    system = PlinderSystem(system_id=system_id, release=release, source_mmcif=source,
                           reconstruction_dir=root / "reconstructed" / pdb)
    was_cached = (root / "reconstructed" / pdb / "system.cif").is_file()
    derivative = Path(system.system_cif)
    checks = verify_coordinates(source, "1", proteins, ligand, components, derivative)
    source_url = ("https://files-versioned.rcsb.org/pdb_versioned/data/entries/"
                  f"{shard}/pdb_0000{pdb}/pdb_0000{pdb}_xyz_v{revision.replace('.', '-')}.cif.gz")
    return {
        "pdb_id": pdb.upper(), "revision": revision, "assembly": "1", "plinder_system_id": system_id,
        "ligand_component": component, "ligand_components": components, "ligand_instance": ligand,
        "plinder_ligand_id": ligand_id, "protein_instances": proteins,
        "source": source.relative_to(root).as_posix(), "source_url": source_url, "source_sha256": digest(source),
        "source_bytes": source.stat().st_size, "derivative": derivative.relative_to(root).as_posix(),
        "derivative_bytes": derivative.stat().st_size,
        "plinder_candidate_rows_for_entry": len(annotations),
        "plinder_candidate_systems_for_entry": int(annotations.system_id.nunique()),
        "source_hem_instances": hem_labels, "reused_derivative": was_cached,
        "coordinate_checks": checks, "pilot_qc": "PASS_WITH_SOURCE_PROVENANCE",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("plinder_release", type=Path)
    parser.add_argument("--git-commit", required=True)
    parser.add_argument("--case", action="append", choices=[case[0] for case in CASES])
    args = parser.parse_args()
    root = args.root.resolve()
    release = PlinderRelease(data_dir=args.plinder_release)
    started_at = datetime.now(timezone.utc)
    clock = time.monotonic()
    selected = args.case or ["1m17", "1b86"]
    systems = [sample(root, release, case) for case in CASES if case[0] in selected]
    manifest = root / "holo_systems.json"
    with manifest.open("x", encoding="utf-8") as stream:
        json.dump(systems, stream, ensure_ascii=False, indent=2)
    reread = json.loads(manifest.read_text(encoding="utf-8"))
    for system in reread:
        source = root / system["source"]
        if digest(source) != system["source_sha256"]:
            raise ValueError("Source checksum differs on readback")
        checks = verify_coordinates(source, system["assembly"], system["protein_instances"],
                                    system["ligand_instance"], system["ligand_components"],
                                    root / system["derivative"])
        if checks != system["coordinate_checks"]:
            raise ValueError("Coordinate evidence differs on readback")
    result = {
        "scope": "isolated Holo representative-view verification", "status": "PASS", "git_commit": args.git_commit,
        "started_at_utc": started_at.isoformat(), "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(time.monotonic() - clock, 3), "python": sys.version.split()[0],
        "gemmi": version("gemmi"), "plinder": version("plinder"), "plinder_release": "2026-09",
        "source_files": len(reread), "source_bytes": sum(item["source_bytes"] for item in reread),
        "pilot_systems": len(reread), "readback_passed": True, "manifest": manifest.name,
    }
    with (root / "report.json").open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
