# APHCP public verification snapshot

Generated during GitHub synchronization from accepted local receipts. The dates below are **evidence dates**, not the time this page was generated. This snapshot describes inventory and unpublished staging; it is not a dataset release or a live service status.

Public runtime fingerprint (SHA-256): `84883fe9ae26f0f092857269ae0e196e9344d0ac11c1dbd3c056af07f3dfcf23`

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

## M3 · First continuous production batch (unpublished)

- Independently verified: **2026-11-22T18:28:55+00:00**. First release dates from **2018-02-07** (included) to **2018-02-15** (excluded).
- **359** exact-revision originals verified; **3,844** candidates, **3,497** read back and **3,497** new retrievable systems. **69** confirmed in staging; **3,428** held for review. Unified catalog: **46,525** retrievable systems.
- **145,717,825 bytes** of archival response bodies; **1456.0 s** for the production pipeline in this window. Replaying the same input reused committed phases. In a separate control, the resume command repaired one failed Entry without modifying prior results.

## M3 · Scale gate: second continuous production batch (unpublished)

- Independently verified: **2026-11-23T06:18:02+00:00**. First release dates from **2018-02-22** (included) to **2018-03-01** (excluded).
- **325** exact-revision originals verified; **2,174** candidates, **2,035** read back and **2,035** new retrievable systems. **72** confirmed in staging; **1,963** held for review. Unified catalog: **48,560** retrievable systems.
- **115,424,513 bytes** of archival response bodies; **456.5 s** for the production pipeline in this window. Replaying the same input reused committed phases. In a separate control, the resume command repaired one failed Entry without modifying prior results.

## Open review

- Measure sustained whole-catalog throughput and allocated storage; the bounded scale gate is not a completion date.
- Keep partner-connection and role exceptions under review; the bounded pilot does not broaden scientific acceptance.

See the [review guide](review.md) for scope, evidence definitions and review instructions.
