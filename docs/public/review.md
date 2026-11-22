# APHCP independent review guide

Start with the [latest published verification snapshot](status.md). It separates public code provenance from the date and scope of accepted scientific evidence.

## Review one exact version

1. Resolve `main` to a full Git commit SHA when starting a review. Open the files and this page at that commit, or use the two commit SHAs supplied for a change review. If you can only see an older version, say so before reviewing.
2. For each finding, return `commit → file and line → trigger or counterexample → impact → smallest verifying test`. Separate a code finding from something that needs an actual structure, server run or scientific decision.
3. Report the denominator and evidence scope for numeric claims. The maintainer verifies findings against protected receipts or an isolated run and synchronizes any tested code fix. A reply in a separate GPT conversation is not automatically delivered by GitHub.

## What the evidence means

- A source entry is a fixed PDB revision. A system identity also distinguishes model, Assembly operation, bound protein chain instances and target partner instances. Identical CCD codes do not imply identical instances. Keep explicitly supported partner members and source mappings; a connection or contact alone does not establish experimental role.
- Coordinate identity, experimental fidelity, topology, partner scope, experimental role and chemical type are separate evidence. `UNKNOWN` chemistry can be annotated later without discarding an otherwise eligible structure. `REVIEW` keeps the evidence available. Neither a successful readback nor an index hit is a scientific confirmation by itself.
- The inventory counts metadata candidates, not structures. The observation catalog includes retrievable but unconfirmed systems. `REGISTERED` in current snapshots means confirmed in unpublished staging, not a published Holo Master resource.
- A receipt may predate the current source snapshot. The runtime fingerprint changes with public executable inputs and excludes this review tooling, docs and tests. A stage only claims coverage of the current public runtime if accepted evidence is explicitly pinned to that exact fingerprint.

## Where to begin

| Question | Public code |
| --- | --- |
| Inventory, association and cross-window reconciliation | [Window runner](../../scripts/m1_holo_windows.py), [rollup](../../scripts/m1_holo_rollup.py) |
| Discovery, structure export and provenance | [Native candidates](../../scripts/holo_native_candidates.py), [source mapping](../../scripts/holo_source_mapping.py), [V1 export](../../scripts/holo_v1_export.py) |
| Chemistry, registration and readback | [Chemistry](../../scripts/holo_m2_chemtype.py), [new-batch driver](../../scripts/holo_m3_increment.py), [registration and lookup](../../scripts/holo_m3_register.py), [decision history](../../scripts/holo_m3_decisions.py), [observation catalog](../../scripts/holo_m3_observations.py) |
| Small reproducible checks | [Synthetic operations input](../../tests/fixtures/m0_work_items.json), [integration test](../../tests/integration/test_m0.py), [synthetic receipt example](../../tests/integration/test_public_review_status.py) |

Run the synthetic receipt-to-summary check with Python 3.10 or newer:

```sh
python -m unittest discover -s tests/integration -p test_public_review_status.py -v
```

Production coordinates, entry-level ledgers, complete receipts, and private deployment settings stay outside this repository. A request that needs them should identify the specific evidence needed so the maintainer can check it or supply a bounded extract.
