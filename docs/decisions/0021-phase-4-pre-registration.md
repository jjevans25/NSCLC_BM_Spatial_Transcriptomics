# ADR 0021 — Phase 4 pre-registration: the pairing rule, the detection filter, and the two nulls

**Date:** 2026-09-14
**Status:** Accepted
**Task:** P4-T1 (before any Phase 4 result was computed)
**Discharges:** PROJECT_PLAN §6 P4-T1 ("write it in an ADR before running
anything") and the P4-T2 detection filter
**Constrained by:** ADR 0014 (A5 is exploratory), ADR 0007 (the one detection
rule), ADR 0012 (primary + sensitivity, both reported), ADR 0016 §3 (a threshold
chosen after seeing the result is not a threshold), ADR 0018 (a restricted
primary with a declared exploratory secondary)

## Context

Phase 4 nominates inferred ligand–receptor crosstalk between adjacent
compartments, within site, across two adjacencies: lung `L ↔ TIME-L` and brain
`LB ↔ TIME-B`. Every quantity it produces is a ranked correlation, and a ranked
correlation is the single easiest thing in this project to manufacture by
choosing a threshold late. PROJECT_PLAN §6 says so in terms for P4-T1, and
ADR 0016 §3 exists because the same failure was live in Phase 3.

Eight decisions fix what Phase 4 is allowed to say. All eight are recorded here,
in the same commit that opens `config/config.yaml → phases.crosstalk`, and
before `06_crosstalk.smk` contained a single rule.

**A5 is exploratory (ADR 0014) and nothing below changes that.** The demotion is
about what may be claimed, not about what is built. No Phase 4 result may be a
headline claim; every nomination states the per-site detection status of both
partners; and the anticipated null is an assay-sensitivity limit, never evidence
that the crosstalk is absent.

## Decision

### 1. Duplicate AOIs are averaged on the log2 scale; the higher-QC pick is a sensitivity

Verified against `results/tables/design_matrix.tsv` rather than taken from the
plan's estimate — the two agree, and now the agreement is checked:

| adjacency | patients | duplicates |
|---|---|---|
| lung, `L ↔ TIME-L` | **13** — P5, P12, P15, P18, P19, P24, P26, P29, P30, P32, P35, P40, P43 | P12 and P24, **two `TIME-L` each** |
| brain, `LB ↔ TIME-B` | **8** — P5, P12, P14, P15, P19, P20, P31, P35 | **none** |

Both duplicates sit on the same side of the same adjacency. P15's duplicate is
`TBME`, which Phase 4 does not pair. So this ADR decides exactly one thing.

```yaml
crosstalk:
  pairing:
    duplicate_rule: mean_log2
    duplicate_sensitivity: highest_detection_rate
```

**Primary: the arithmetic mean of `log2(q3 + 1)` across a patient's duplicate
AOIs**, pre-registered before any rho is computed.

**Rejected: taking the higher-QC AOI as the primary.** It needs a tiebreak rule
that does not exist anywhere in this project, and inventing one after seeing the
correlations is precisely the failure ADR 0016 §3 was written to prevent. It is
kept as a **declared sensitivity** — highest `gene_detection_rate` from
`results/tables/qc_metrics.tsv` — run in full and reported beside the primary
with a `delta_vs_primary` column, the shape ADR 0012 established for `dsp_run`.

**The choice cannot be load-bearing and the ADR says so up front: it moves 2 of
13 lung patients and 0 of 8 brain ones.** If the sensitivity ever disagrees
materially with the primary, that is a finding about two AOIs from one patient,
not about crosstalk.

This is flag-don't-drop applied to a collapse rather than an exclusion: no AOI
is discarded, both contribute, and `results/tables/qc_excluded.tsv` gains
nothing.

### 2. The database is CellChatDB **v2** (jinworks), pinned at a commit, acquired by a sanctioned writer

PROJECT_PLAN §6 offers CellChatDB or CellPhoneDB. **CellChatDB, because of one
column.** Its `annotation` field classifies every interaction as Secreted
Signaling, ECM-Receptor or Cell-Cell Contact, and ADR 0014's case for demoting
A5 was a claim about exactly that axis — "secreted ligands and chemokines are
systematically undetected". With the class column, P4-T2's before-and-after
filtering counts can be reported **by interaction class**, which turns ADR 0014's
sentence into a measured number in the one place it matters most. CellPhoneDB
has no equivalent.

**The version is v2**, `jinworks/CellChat` at commit
`75253cd0c9e68410e6e721a6d3a0419a1d7e358f`, `data/CellChatDB.human.rda`,
1,519,587 bytes, sha256 `582e99db…3ad8`.

*Rejected: CellChatDB v1*, `sqjin/CellChat` at its final archived commit
`e4f68625b074247d619c2e488d33970cc531e17c`. The argument for v1 was that the
repository is archived, that every one of its 1,939 interactions is
protein-coding so no exclusion rule is needed, and that its smaller family
matches PROJECT_PLAN §6's own "~1000 pairs at n = 13" sizing. **All three
arguments fail on inspection and the decision was changed accordingly.** A
commit SHA is already immutable, so archiving buys nothing a pin does not; §6's
empirical FDR is computed by permutation, which self-calibrates to whatever the
family size is, so a larger family is not the penalty it would be under a
Bonferroni-style correction; and §6's "~1000" was the plan's estimate, not a
constraint.

The one argument that could not be settled by reasoning was that v2's larger
complex vocabulary (338 entries of up to 5 subunits, against v1's 157 of up to
4) might push *more* interactions out of the 1:1 primary under §4, leaving a
smaller primary despite the larger database. **It was measured rather than
argued, against the 18,694 measured symbols, and it is false:**

| | v1 | v2 |
|---|---|---|
| interactions | 1,939 | 3,233 |
| non-protein, excluded | 0 | 994 |
| usable (protein-coding) | 1,939 | 2,239 |
| 1:1 candidates | 1,011 | 1,256 |
| **1:1 with both symbols measured — the primary** | **963** | **1,149** |
| complex, all subunits measured — the exploratory | 913 | 967 |

**And the gain sits where this assay has resolution.** Of v2's 186 additional
measurable 1:1 interactions, 138 are **Cell-Cell Contact** (257 → 395) — the
membrane-bound class `docs/NEXT_STEPS.md`'s reconnaissance finds detected in
nearly every AOI — against +3 ECM-Receptor and +45 Secreted Signaling, the class
ADR 0014 established is largely below background at both sites. v2 adds coverage
almost entirely to the interactions Phase 4 can actually measure, which is the
opposite of the concern that motivated v1.

**v2's 994 `Non-protein Signaling` interactions are excluded before anything
else, as a pre-registered rule.** Their ligands are metabolites and
neurotransmitters with no transcript to measure, so leaving them in would let
them be counted as "filtered out for low detection" and silently inflate the
before/after gap this task exists to report honestly. The exclusion is reported
as its own count, separately from the detection filter.

v2 also carries per-partner annotation v1 lacks — `ligand.secreted_type`,
`ligand.transmembrane`, `receptor.transmembrane`, `receptor.surfaceome_*` — which
lets Phase 4's likely central finding be stated in the database's own molecular
terms rather than inferred from the three-way class column alone. Those columns
are **descriptive and never a filter**: nothing enters or leaves an analysis
table on the strength of them.

**Acquisition.** The CellChat paper's Springer supplementary files were probed
first, on the hope of reusing P3-T2b's dependency-free xlsx reader.
`41467_2021_21246_MOESM4_ESM.xlsx` is a **method comparison** (sheets: CellChat,
CellPhoneDB, iTALK, SingleCellSignalR) and MOESM5 is developmental trajectory
data; **the paper does not deposit the database**. Recorded here so the route is
not re-probed.

So the pin is the `.rda`, read with base R `load()` — a base function, no
CellChat installation — in the **existing, unmodified `r-geomx` env**, which
already reads `.RData` in `deconvolve.R`. An R rule declares
`results/interim/env_repair/r-geomx.ok` as an input like every other R rule
(ADR 0013), exports the interaction, complex and annotation tables as TSV across
the ADR 0001 boundary, and a Python rule does the resolution and the detection
filter. **No conda env YAML is edited and no dependency is added anywhere** —
`py-analysis.yaml` is `p0t2_fetch_geo`'s env and its `protected()` outputs abort
the DAG on a software-env trigger, and an edit to `r-stats.yaml` or
`r-geomx.yaml` would needlessly re-run Phase 2 and Phase 3's fits.

`resources/` is read-only to rules (hard constraint 1, ADR 0005), so the release
is acquired the way P0-T2, P2-T5 and P3-T2b acquired theirs — **the fourth
sanctioned writer**: a `protected()` fetch rule reusing `acquire_geo.py`
unchanged, a `record_provenance.py` rule writing
`resources/ligand_receptor_provenance.tsv` and
`resources/ligand_receptor_checksums.sha256`, and a pinned `sha256` with
`stability: fixed`. `base_url` pins a **commit**, never a branch, for the reason
`config.yaml`'s `reference.base_url` does: a branch ref would silently change the
database under a rerun and the digest would become a failure rather than a
control.

The pin lives in `config/ligand_receptor.yaml`, not in `config/config.yaml`,
for the reason `config/external_validation.yaml` does: `config.yaml`'s SHA-256
is recorded in the `.h5ad` `uns`, so a re-pin must not cost a Phase 1–3 rebuild.

### 3. The detection filter is 0.5 — Phase 2's and Phase 3's value, not a new one

```yaml
crosstalk:
  detection:
    detected_in_aoi_fraction: 0.5
```

The per-value rule underneath is unchanged and is **not** re-derived: a gene is
detected in an AOI when `layers['q3']` exceeds `qc.detection_background_multiple`
(2.0) × that AOI's `NegProbe-WTX` level. That is ADR 0007's definition,
implemented in `score_signatures.py` and lifted by `checkpoint_detection.py`.
Phase 4 lifts it a third time rather than writing a third copy of it.

**A pair enters only if both members clear the floor, each in the compartment
that member is measured in** — the ligand in the source compartment, the
receptor in the target — and, because CellChatDB is directed and both directions
are tested (§9), admission is evaluated **per direction**, not per pair.

Detection is computed over **all** AOIs of an `aoi_code` (`L` 30, `TIME-L` 15,
`LB` 27, `TIME-B` 8), not over the paired subset. That keeps Phase 4's coverage
numbers directly comparable with Phase 2's and Phase 3's, which are the only
coverage numbers a reader has to calibrate against.

0.5 is chosen because it is `contexture.scoring.detected_in_aoi_fraction` and
`checkpoints.detection.detected_in_aoi_fraction`, which are in turn ADR 0008
point 4 verbatim. **Phase 4 introduces no threshold the project did not already
operate under.**

**Rejected: a stricter 0.75.** More conservative against both the background
gradient and the high-expression bias §6 addresses, but it is a number with no
precedent anywhere in this project, and a *different* floor here than in Phase 3
would be a decision presented as a detail.

**Rejected: a looser 0.3.** It would admit more of the ligand side, including
some brain-limited chemokines, at the cost of re-opening the ADR 0014 artefact
that Phases 2 and 3 each spent a task closing.

### 4. Multi-subunit interactions: a 1:1 primary, complexes in a declared exploratory table

CellChatDB contains heteromeric entries (`TGFB1 → TGFBR1_TGFBR2` and similar).
A detection floor means something different for a complex than for a single
gene, and putting both in one table would make one column mean two things.

**Primary — one-to-one only.** An interaction enters the primary if its ligand
and its receptor each resolve to exactly one measured gene symbol.

**Exploratory — complexes, expanded to subunits** (CellChatDB v2 carries 338
complex entries of up to five subunits), admitted only when **every** subunit
clears the floor in its compartment, carrying **no FDR of any kind** and every
row labelled. This is ADR 0018's two-table shape applied to a second
question, and it carries the same prohibition: no Phase 4 sentence may rest on
it.

`n_complex_dropped` is a reported number in `crosstalk_lr_summary.json`, not
bookkeeping — it is part of the same before/after accounting §2 exists to
produce.

**Rejected: expanding complexes into the primary**, valuing each as the mean of
its subunits' log2 values (CellChat's geometric-mean convention on the linear
scale). It admits more real biology — TGF-β, the integrins, the heteromeric
interleukin receptors — but the detection floor then binds one subunit in a
complex row and one whole gene in a 1:1 row, inside one table with one q column.

**Rejected: dropping complexes entirely.** Simplest and most conservative, but
the loss would be recorded only as a count, and the exploratory table costs one
extra output to record it as data.

### 5. The primary correlation keeps every paired patient; a detected-both refit is the sensitivity

A gene can clear the compartment-level floor of §3 and still be below background
in one particular patient's AOI.

**Primary: rho is computed over all paired patients — lung n = 13, brain n = 8 —
with `n_detected_both` reported on every row** beside `n`, so a reader sees
immediately how much of a correlation rests on values at background.

**Sensitivity: the detected-both-only refit**, bound by
`crosstalk.correlation.min_n_patients: 6`, with a `delta_vs_primary` column.

Two reasons the restriction is the sensitivity and not the primary. Dropping on
detection is **expression-dependent selection** — it removes the low values of
both members preferentially, which is a mechanism for inducing correlation, not
for removing one. And it makes n vary per pair, so §6's permutation null would
have to be recomputed at every distinct n rather than once per site; a null
whose calibration depends on the quantity being calibrated is worth avoiding
when the alternative is a column.

**Rejected: dropping undetected patients from the primary.** Closer to "the
measurement exists", which is the instinct ADR 0018 is built on — but ADR 0018's
restriction operates on a *pre-specified group*, not on individual observations
chosen by their own value, and that difference is the whole of why one is
defensible and the other is not.

### 6. Two nulls, both pre-registered, and only one of them is the gate

**Null A — the gate.** Permute **the pairing**, not the expression: shuffle
patient labels on the immune side, **within site**, so the null preserves both
n's exactly. `n_permutations: 10000`. Pool the null |rho| and report an
**empirical FDR** column — expected null pairs at or above a threshold divided by
observed pairs at or above it. `p1t3_variance_partition`'s permutation null is
the in-repo precedent for how this project structures one.

**Null B — the abundance-matched control.** `docs/NEXT_STEPS.md`'s
reconnaissance warns that a nomination list dominated by VEGFA, TGFB1 and CD47
may be measuring "high-expression genes correlate with high-expression genes"
rather than crosstalk — these are broadly expressed, and their correlation across
13 patients may track shared AOI quality. Null B answers that directly: for each
candidate, draw `n_matched_draws: 10000` random gene pairs whose members sit in
the same mean-expression decile (`n_expression_deciles: 10`) as the real ligand
in the source compartment and the real receptor in the target compartment, with
the **true** pairing intact, and report the fraction reaching the observed |rho|.

**Gate 4 turns on Null A alone.** Null B is a control that qualifies a
nomination; it is not a second gate, it introduces no second FDR family, and a
pair that clears A but not B is reported as such rather than removed.

Both nulls seed from `config["seed"]`, passed explicitly through
`params: seed=SEED` and echoed into the summary JSON. No implicit RNG anywhere.

**Rejected: Null A only**, which is all PROJECT_PLAN §6 P4-T4 mandates. It
satisfies Gate 4 exactly and leaves the abundance confound stated in prose. The
reconnaissance makes that confound the single most likely way this phase
produces a wrong answer, and prose is not a measurement.

**Rejected: adding Null B only if the ranking turns out to be
abundance-dominated.** That is ADR 0016 §3's failure with a different subject.

**Multiplicity: one family per site, and no other family anywhere.** The
empirical FDR is computed across **all primary direction-rows within a site**,
separately for lung and brain. Two reasons, and the first is not a convention:
the permutation null *is* the correction, and it is calibrated on a permutation
that preserves that site's n — so a null built at n = 13 cannot be applied to
the brain's n = 8. The second is ADR 0016 §4's rule applied unchanged — a
family is the set of tests actually run to answer one question, and "which
interactions correlate in lung" is one question.

Both directions of an interaction sit in the *same* family. They are distinct
hypotheses, but they are hypotheses asked in one sweep of one database, and
splitting them would create two families where one was asked.

**No second family is created anywhere.** Null B produces an abundance-matched
p per nomination and no q — it qualifies a row, it does not test a new
hypothesis. The complex table (§4) carries no q of any kind. The detected-both
refit (§5) is a sensitivity on an existing hypothesis, and ADR 0012 already
established that re-adjusting a sensitivity would treat the same hypothesis
under a different adjustment as a new one.

### 7. The nomination cap is the plan's number, carried as a rule param

PROJECT_PLAN §6 P4-T5 says "Nominate 5–10 pairs", so the cap of 10 per site is
**not a threshold this phase invented** and does not belong in the
pre-registered config block beside the ones that are. It is carried as a rule
param on `p4t5_nominations`, exactly as `p3t4_checkpoint_paired_check` carries
`expected_immune_paired=5` from PROJECT_PLAN §2.3.

It truncates a long list and **never licenses a longer one**. Gate 4 turns on
P4-T4's FDR, not on this table, and the cap is one of the mechanisms keeping
P4-T5 out of the gate's subject.

There is a second, practical reason it is not a config key: `config.yaml` is a
declared input of `p0t7_assemble_h5ad`, so adding one would rebuild the `.h5ad`
and re-run Phases 1–3 for a display cap.

### 8. `LGALS9` is not measurable, and no paralog is substituted

P4-T6 cross-references three pairs the Nat Commun 2024 Visium study reports as
adjacent. **`LGALS9` is not among the 18,694 measured genes** — the panel carries
`LGALS1`, `LGALS3`, `LGALS8` and `LGALS9C`, but not `LGALS9` itself — so
`LGALS9–HAVCR2` is reported `not_measurable`.

**`LGALS9C` is not substituted for it.** Adding a gene to an analysis panel is a
stop-and-ask (CLAUDE.md), ADR 0011 and ADR 0015 govern membership, and ADR 0015
exists *because* citation checking caught four errors a schema could not. A
paralog swapped in to keep a comparator alive would be a membership change made
to protect a result.

### 9. Both directions of every interaction are tested and labelled

CellChatDB is directed. "Ligand on the tumour AOI, receptor on the paired immune
AOI" and its reverse are different biological claims, and Phase 4 reports both
with an explicit `direction` column. Neither is a control for the other.

## Consequences

- `config/config.yaml` gains a `crosstalk:` block and `phases.crosstalk` goes
  true **in a single commit, with every Phase 4 key present**. That file is a
  declared input of `p0t7_assemble_h5ad`, so the edit rebuilds the `.h5ad` and
  re-runs Phases 1–3. Phase 4 pays that cost exactly once, here, the way P3-T1
  paid it at `f1f55b8` — and verifies it the same way, by confirming the existing
  tables return byte-identical except `h5ad_summary.json`'s `git_sha` and
  `config_sha256`.
- `workflow/schemas/config.schema.yaml` has `additionalProperties: false` at the
  root, so `crosstalk` joins both `properties` and `required` in the same edit or
  the workflow fails at parse.
- **Only the 1:1 primary table carries an empirical FDR.** The complex table and
  both sensitivities carry none, by §4 and ADR 0012.
- **The expected Gate 4 outcome is not a bare null.** The reconnaissance
  indicates a structurally biased subset will survive — the membrane-bound and
  matrix-associated axis, with the soluble-immune axis absent at both sites — and
  **that bias is the finding**, the same shape Gate 3's answer took. §2's
  by-class counts are what let it be stated as a measurement.
- Changing any threshold in this ADR remains a stop-and-ask, and a result that
  makes one look wrong is a conversation, not an edit.
- ADR 0014 §4 recorded that "any obligation an ADR defers to a later gate must
  land in that phase's task table". This ADR is owned by P4-T1 and is the first
  item in Phase 4's commit sequence, so the obligation it creates cannot go
  unowned the way the A5 re-decision did.
