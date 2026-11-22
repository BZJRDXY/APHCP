# APHCP public verification snapshot

Generated during GitHub synchronization from accepted local receipts. The dates below are **evidence dates**, not the time this page was generated. This snapshot describes inventory and unpublished staging; it is not a dataset release or a live service status.

Public runtime fingerprint (SHA-256): `62fab046341ac930ea5099109f641b8ee1ad38325f72e0256172b5ad3e90eb1f`

## M1 · Metadata inventory

- Independently verified: **2026-11-21T07:43:47+00:00**.
- **257,840 unique experimental protein PDB entries** across **43 date windows**.
- PLINDER ligand annotations: **140,236** entries with annotations; **117,604** without. **609,131** indexed candidate systems, before structure QC.
- Verification and current public code: Unknown: the historical Run is not pinned to a public runtime fingerprint.

## M3 · Staging registration and structure observation catalog

- Catalog independently verified: **2026-11-22 (Asia/Shanghai)**.
- Registration review: **527** proposals; **472** confirmed in staging, **51** held for review with files retained, **4** out of scope.
- Separate observation catalog: **43,013** records, **43,007** retrievable, including the **472** confirmed staging systems. The catalog has **6** out-of-scope records without supplied assets. Its denominator differs from the registration review.
- **2,344** retrievable observations have an unknown primary chemical type. Unknown type does not by itself block eligible structures.
- Verification and current public code: Unknown: the historical Run is not pinned to a public runtime fingerprint.

## M3 · New-batch pilot (unpublished)

- Independently verified: **2026-11-22T17:22:50+00:00**. **9** new entries: **3** same-revision, **3** different-revision, and **3** without PLINDER ligand annotations.
- **29** candidate systems, **21** structures read back. **4** newly confirmed in this isolated staging batch and **17** held with files retained. These numbers are separate from the historical totals above.
- A test-only `REVIEW` → `REGISTERED` transition preserved the same system and files, retained both decisions, and was independent of run order. Replaying the real batch reused its outputs.

## Open review

- Measure sustained throughput on a separately approved production increment; the bounded pilot is not a scale estimate.
- Keep partner-connection and role exceptions under review; the bounded pilot does not broaden scientific acceptance.

See the [review guide](review.md) for scope, evidence definitions and review instructions.
