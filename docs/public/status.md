# APHCP public verification snapshot

Generated during GitHub synchronization from accepted local receipts. The dates below are **evidence dates**, not the time this page was generated. This snapshot describes inventory and unpublished staging; it is not a dataset release or a live service status.

Public runtime fingerprint (SHA-256): `5f8b54f65d5bd47226e74a3bc2f96884dc13d960cb75bcdf626303468c733534`

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

## M3 · Rolling recovery control (isolated, unpublished)

- Verified: **2026-11-23T08:51:36+00:00**. The interrupted batch reconciled **550** entries in an isolated catalog: **547** reused, **3** recovered. The live catalog pointer was left in place.
- A separate **34**-entry control sealed its next primary batch before the slow-file retry merged. A bounded Range retry recovered the previously partial archive bodies. The retry wave merged once, and replay did not change committed receipts.
- This is engineering verification on already scoped entries. Sustained production throughput and the live controller handoff are outside this control.

## M3 · Bounded production recovery gate (unpublished)

- Independently verified: **2026-11-23T09:21:00+00:00**. The previous controller was closed, its sealed history was preserved, and the interrupted batch was completed by reusing **547** settled entries and recovering **3** missing ones.
- Exactly two new production batches attempted **1,100** entries; **57** were still pending in those batches at the observation limit. Successive batch seals took **398.2 s / 441.1 s**, including handoff/preflight where applicable. These observations do not establish a full-catalog ETA.
- The staging catalog contained **59,036** retrievable structures, including **809** confirmed in staging. Overall, **3,243 / 233,639** frozen rolling entries were settled; **57** attempted entries remained pending. Production paused at the gate; no dataset release was created.

## M3 · Engineering verification after ten main batches (isolated)

- Evidence snapshot **2026-11-23T13:29:10+00:00**: **160** of **160** existing technical gaps recovered; **0** remained at the snapshot. The main catalog pointer was unchanged and no new main batch started.
- Verified that retries continue during a main-stage wait. In an isolated process interruption, recovery waited for the inherited writer lease and reused the committed receipt and source proof. Intentional exit did not restart.
- Under the same 30-second transfer budget, whole GET verified **3 of 24** and Range verified **3 of 24**, with **24** and **66** requests respectively. The small-file Range threshold was retained.
- On the same existing heavy entry, structure build took **183.4 → 74.9 s**; the separate independent review control took **149.4 → 49.9 s**, with equivalent candidate decisions and structure evidence. These are isolated timings, not whole-pipeline throughput.
- With the same 34-system staging input, registration took **23.1 → 5.8 s**. The database integrity and foreign-key checks still ran, and the review result hash and registry counts matched.
- Same-input replay preserved **174** checked evidence files, including the catalog, receipts, and source proofs.

## M3 · Three consecutive production batches (unpublished)

- Verified **2026-11-23T14:08:48+00:00** against execution commit `e13dbfebe264886e8e37d167e709b66d668a3eac`. Current public runtime is not covered by this evidence (verified runtime: `790d6af2f4dc1a26969f054c0bdadfea7feccecbf6f9a6a94b1573e577c81ce6`).
- Batches 11–13 attempted **1,369** entries. **237** technical gaps appeared; **31** were recovered and merged by the last main seal and **19** more at shutdown, leaving **187**. Batch 14 was neither preflighted nor dispatched. Three batches do not establish sustained recovery capacity or a full-catalog ETA.
- **6,682 / 233,639** rolling entries were settled and **226,770** had not started. The staging catalog contained **73,031** retrievable systems, including **1,357** confirmed in staging. Database integrity, foreign keys and file retrieval were checked. No dataset release was created.

## M3 · Archive admission and retry capacity (unpublished)

- Verified **2026-11-23T16:10:58+00:00** against execution commit `8000770af198e23a9afae98a6b7a4ff09d0d0127`. From **187** existing gaps, **185** were recovered and merged before the next main batch.
- Batches 14–16 attempted **1,650** entries. **1,650** HEAD requests reached HTTP; **0** timed out waiting locally. Main GET had **170** partial-body timeouts. The shared archive slot limit was unchanged.
- **180** new technical gaps appeared. **15** were merged during production and **60** more at shutdown, leaving **107** pending overall. Expansion was held because retries did not keep up.
- **8,412 / 233,639** rolling entries were settled. The staging catalog contained **80,058** retrievable systems; no dataset release was created.

## M3 · Bounded rolling-production launch (unpublished)

- Batches 17–19 passed the gate at **2026-11-23T17:14:34+00:00** using `49ef3c74d7ee72fcababdd9f4764389df24d4c58`. Automatic pause and resumption occurred under production load. File retrieval, same-input replay, and the closing backlog check passed.
- **1,650** new entries were attempted; **110** retry entries were recovered. At the gate, **10,166 / 233,639** rolling entries were settled, **3** technical entries had exhausted retries, and no actionable retries remained.
- The staging catalog then contained **88,623** retrievable systems, including **1,933** confirmed in staging. No dataset release was created.
- The subsequent run was independently observed after launch at **2026-11-23T17:15:56+00:00** and is capped at absolute batch **100**. This is launch evidence, not a live status or completion claim.

## M3 · Pipeline and source acquisition comparison

- At **2026-11-24T06:11:14+00:00**, two isolated runs from the same catalog and **170** fixed entries compared serial and overlapped scheduling. Time decreased from **51.5 s** to **42.8 s** with **1,079** identical reviewed system identities, QC outcomes and CIF hashes. This is a cached-input test, not a production throughput forecast.
- On a separate **24**-file sample enriched for historical slow downloads, the alternate route verified **24** files and the server route verified **5** on the first attempt with the same 30-second budget. The alternate-source cache verified input identities and reused all transferred files on replay. No dataset release was made.

## M3 · Batch 100 closure and subsequent bounded launch

- At the batch-100 stop, **53,173/233,639** rolling entries were settled, **5** had exhausted retries, and **180,461** were not yet attempted. The unpublished catalog had **326,582** retrievable systems, **10,655** confirmed in staging.
- A targeted recheck reused frozen files and resolved **1** technical entry exception. **2** additional systems are retrievable with structural evidence; their scientific disposition is still REVIEW. No new coordinates were downloaded for this recheck.
- Batch **101** was observed in progress at **2026-11-24T11:53:30+00:00**, with an internal gate at **105** and an absolute stop at **200** for **51,003** approved entries. Execution code: `29691cdfadd5576639f52652692682fe77a38cdb`. Neither the gate nor completion is attested by this launch snapshot. No dataset release was created.

## M3 · Batches 101–200 accepted and stopped

- Final acceptance: **2026-11-25T03:22:55+00:00**. Batch **200** sealed and production stopped as planned. No batch 201 dispatch or actionable retry backlog. File retrieval, provenance, foreign keys and same-input replay were verified.
- Batches 101–200: **50,995/51,003** frozen PDB entries settled; **8** exhausted technical retries. Rolling total: **104,169/233,639** settled, **12** exhausted, **129,458** not yet attempted. The larger M1 inventory also includes excluded and already processed same-revision entries.
- The unpublished catalog has **681,283** retrievable systems, including **22,143** confirmed in staging. PDB entries and systems use different denominators. Formal releases: **0**.
- Net throughput including retries and catalog merges: **2,588.29 entries/hour** in batches 101–105 with the backup feed; **3,340.21 entries/hour** in batches 106–200 after switching to the server-only feed. Execution code: `29691cdfadd5576639f52652692682fe77a38cdb`.

## Open review

- Measure sustained whole-catalog throughput and allocated storage; the bounded scale gate is not a completion date.
- Keep partner-connection and role exceptions under review; the bounded pilot does not broaden scientific acceptance.

See the [review guide](review.md) for scope, evidence definitions and review instructions.
