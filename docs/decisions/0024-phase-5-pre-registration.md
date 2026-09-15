# ADR 0024 — Phase 5 pre-registration: the cohort, the endpoints, the split, and what a null may say

**Date:** 2026-09-15
**Status:** Accepted
**Task:** P5-T1 (before any Phase 5 result was computed)
**Discharges:** PROJECT_PLAN §6 P5-T1 ("hard gate: if it doesn't join cleanly on
patient ID, stop the phase"), §6 P5-T5 ("you tested 5 signatures × 2 cohorts —
say so, with correction")
**Constrained by:** ADR 0005 (acquisition provenance, `resources/` read-only),
ADR 0007 (the one detection rule), ADR 0008 (A4 exploratory; a near-background
gene is "not assessable", never "lower"), ADR 0009 (patient random intercept is
mandatory on non-independence), ADR 0012 (primary + sensitivity, both reported),
ADR 0016 §3 (a threshold chosen after seeing the result is not a threshold),
ADR 0023 (environments are pinned; an env edit is a transition, not a casual act)

## Context

Phase 5 asks whether an immune-contexture signature measured in a `TIME` AOI
carries prognostic information. PROJECT_PLAN calls it a stretch goal on a
one-week timebox, Gate 5 is *"timebox respected"*, and P5-T1's own acceptance
criterion permits **"a documented abandonment"**.

Its central quantity is a median-split Kaplan–Meier at small n. That is the
easiest quantity in this project to manufacture after the fact: the split point,
the endpoint, the censoring convention and the covariate set are each a degree of
freedom, and every one of them can be chosen to make a curve separate. ADR 0016
§3 exists because that failure was live in Phase 3, and ADR 0021 was written for
the same reason in Phase 4. Phase 5 has strictly more exposure than either,
because its n is smaller than both.

Everything below is fixed **before any association between a signature score and
an outcome has been computed**. What had been measured when this ADR was written
is the design only — the cohort table's shape, the join, the event counts — and
those are pre-registration inputs in the sense that n always is. No
score-versus-outcome quantity of any kind existed.

## Decision

### 1. The clinical table is a fifth sanctioned writer, with its own manifest

`resources/` gains `clinical/`, alongside `raw/` (P0-T2, GEO), `reference/`
(P2-T5, SpatialDecon), `supplementary/` (P3-T2b, the source publication's
supplementary) and `ligand_receptor/` (P4-T2, CellChatDB). ADR 0005's rules apply
unchanged: `protected()` outputs, a pinned digest per artifact, a
`*_checksums.sha256` / `*_provenance.tsv` pair, and verification in a second rule
because the pins cannot be params of a rule whose outputs are write-protected.

The manifest is **`config/clinical.yaml`**, not `config/external_validation.yaml`.
Two reasons, and the second is the load-bearing one:

- That file is named for external *validation*. Supplementary Data 1 is analysis
  *input* — the thing Phase 5 models, not a comparator this project is checked
  against. Filing it there would make the manifest's name a lie about what the
  bytes are for.
- `EXTERNAL_VALIDATION` is a `params:` of `p3t2b_external_validation`
  (`05_checkpoints.smk:255`), so adding an artifact to it changes that rule's
  params hash and re-parses the 45 MB Source Data workbook for no reason
  connected to Phase 5.

The pin lives outside `config/config.yaml` for the reason
`config/ligand_receptor.yaml` does (ADR 0021 §2): `config.yaml` is a declared
input of the `.h5ad` build and its SHA-256 is recorded in `uns`, so re-pinning an
artifact must never cost a Phases 1–4 rebuild. `config.yaml` names the file; the
digest sits in the file.

TCGA (P5-T4) becomes a **sixth** sanctioned writer on identical terms —
`config/tcga.yaml` → `resources/tcga/`. Its path is named in `config.yaml` now,
in the same single edit, so P5-T4 costs no second rebuild. The source is
**cBioPortal**, because PMID 36216799's own TCGA LUAD analyses came from
cBioPortal and a comparator should share provenance with the claim it
contextualises. UCSC Xena is the fallback if the datahub tarball proves awkward;
that substitution changes `config/tcga.yaml` only and is free.

### 2. The cohort is the 16 distinct `TIME` patients, and the arms are within site

Phase 2 scored signatures on `TIME-L` and `TIME-B` only. `signature_scores.tsv`
therefore covers **13 lung + 8 brain patients, 5 of them shared — 16 distinct**.

PROJECT_PLAN §6 P5-T3 says "at ~35 patients this is descriptive". That 35 is the
number of patients with an `L` or `LB` AOI, and **none of them has a signature
score**. The plan's estimate and Phase 2's actual scope disagree, and the plan is
wrong here rather than this file.

**Phase 5 does not extend signature scoring to `L`/`LB` to reach 35.** Doing so
would be a scope change of the kind CLAUDE.md lists as a stop-and-ask, and it
would re-open the detection floor in compartments where the panel clears far less
than it does in `TIME-L` (P3-T2: `L` 2/9 genes, `LB` 3/9, against `TIME-L` 8/9).
A larger n bought with a less measurable compartment is not a larger n.

Arms are **within site**, never pooled: compartment is 41% of variance (Gate 1)
and `site` is the contrast Phases 2–4 all estimate within compartment. Pooling 13
lung and 8 brain patients into one 21-AOI survival model would mix the two.

### 3. Endpoints are fixed here, and the censoring convention is measured, not assumed

Supplementary Data 1 carries two time-to-event columns. They are assigned by
site, once, before any curve:

| arm | endpoint column (normalised) |
|---|---|
| lung (`TIME-L`) | `Primary lung cancer diagnosis to death (Months)` |
| brain (`TIME-B`) | `Brain metastasis diagnosis to death (Months)` |

**Censoring is the literal string `Alive` in the lung column, and an EMPTY CELL
in the brain column.** This is not a guess and it is not what
`docs/data-provenance.md` §Q5 implies. Measured from the deposited file: patients
6, 11 and 35 carry `Alive` in the lung column, and **exactly those three** carry
an empty brain cell. The sets are identical.

That identity is what licenses reading a blank brain value as **censored rather
than missing**, and it is asserted by the rule, not trusted — if the two sets ever
differ, the phase stops. Treating those blanks as missing instead would silently
drop the only censored observations in the study.

Event counts follow, and they are small enough to state in the pre-registration:
**lung 13 patients, 12 events, 1 censored; brain 8 patients, 7 events, 1
censored.** The censored patient is P35 in both arms.

### 4. Age is an ordered factor; free text and missingness are normalised, never trusted

`docs/data-provenance.md` §Q5's three constraints, with the exact tokens observed
in the deposited file rather than the recollection of them:

- **Age is banded** — exactly five levels, `40s`, `50s`, `60s`, `70s`, `90s`, no
  missing. It enters any model as an **ordered factor**, never as a continuous
  covariate. There is no `80s` and the gap is real, so the levels are declared,
  not inferred from what appears.
- **Missing tokens are `N/A`, `NA`, `n/a`, `Unspecified` and the empty cell.** All
  five occur. The set is declared in the manifest and applied uniformly; a token
  outside it is a hard failure, because a sixth spelling that silently reads as
  data is exactly the failure this rule exists to prevent.
- **Free text is inconsistently cased and padded, including the headers.** The
  header cells carry **embedded newlines** — `Age at \nNSCLC diagnosis`,
  `Primary lung cancer \ndiagnosis to death (Months)` — and
  ` Location of BrM ` has a leading *and* trailing space.
  `docs/data-provenance.md` §Q5 transcribed these already-normalised, so its
  column list is a normalised list and must not be matched against verbatim.
  **Normalisation is: collapse every whitespace run to one space, then strip.**
  Applied to headers and to string values alike.
  `Gender` additionally carries a lowercase `'m'` beside `'M'`/`'F'` — a fourth
  inconsistency §Q5 does not mention — so categorical values are casefolded
  before comparison.

### 5. The join is a hard gate whose failure action is to stop the phase

Two assertions, from `docs/data-provenance.md` §Q5, both of which must pass
before the join counts as established:

1. Every GEO patient number appears in Supplementary Data 1.
2. The AOI-code numeric suffix equals the patient number.

**No fuzzy matching, at all.** PROJECT_PLAN §6 is explicit — "do not spend three
days on fuzzy ID matching for a stretch goal" — and P5-T1's acceptance criterion
is "a joined table **or** a documented abandonment". If either assertion fails,
Phase 5 stops and the abandonment is the deliverable. Phase 6, which the plan
calls the primary deliverable, is already unblocked.

The assertions are re-derived by the rule on every run. They are not satisfied by
having been checked once by hand.

### 6. The split is the median, and n per arm is reported inline

Median split per signature, within site. Not a tertile, not an optimal cutpoint,
not a maximally-selected rank statistic — those are the standard routes to a
separation that does not replicate, and at this n they are not defensible.

The median is computed **within the arm being split**, so the arms are as close
to balanced as an odd n allows: 6/7 in lung, 4/4 in brain.

Every reported curve states its n per arm inline, and every brain statement
states **`TIME-B` n = 8** inline (CLAUDE.md hard constraint 8). Every reported
p-value carries n, effect size and a confidence interval (hard constraint 7) —
for a log-rank that means the hazard ratio and its CI beside the test, never the
p-value alone.

### 7. A null is uninformative, not negative — declared before the first curve

With 12 events in lung and 7 in brain, this analysis detects only very large
hazard ratios. **A null result in Phase 5 is uninformative, not evidence of no
association**, in exactly the sense the 1.1–1.3 SD power floor makes a Phase 2
null uninformative (P0-T8), and in exactly the sense ADR 0022 makes Phase 4's
empty nomination table an assay-sensitivity limit rather than a claim about
biology.

This sentence is fixed here so it cannot be softened later if a curve happens to
separate. **A separation at 6 versus 7 patients is a description of this cohort,
never an inferential claim**, and no Phase 5 result may be a headline claim —
the same status ADR 0008 gave A4 and ADR 0014 gave A5.

### 8. The multiplicity family is declared before the first p-value

**6 signatures × 2 cohorts = 12 tests**, corrected with Benjamini–Hochberg at
α = 0.05, the method and level every other phase uses.

Six, not the plan's five: `config/signatures.yaml` resolves to six sets
(`antigen_presentation`, `cytotoxicity`, `exhaustion`, `myeloid_m1`,
`myeloid_m2`, `tls`). PROJECT_PLAN §6 P5-T5 says "5 signatures"; the repository
says six, and the repository is what will be tested.

**Three of those six failed their Phase 2 detection-coverage floor** —
`exhaustion` 1/6 and `tls` 1/5 in brain, `myeloid_m1` 2/10 at *both* sites. They
are scored, corrected within the family, and reported as **"not assessable"**
wherever they failed that floor, never as a prognostic null. ADR 0008's rule is
not suspended by a change of outcome variable: a score built from genes at
background measures background, whatever it is regressed against.

The two cohorts are this study and TCGA LUAD (P5-T4). They are **one family**,
because the plan's own framing — "does the same signature stratify a properly
powered independent cohort" — is a single question asked twice, and correcting
each cohort separately would understate the search.

## Consequences

- `config/config.yaml` is edited **once**, opening `phases.survival` and the
  `survival:` block with every Phase 5 key present including P5-T4's. That edit
  changes the file's SHA-256, which is recorded in the `.h5ad` `uns`, so it
  rebuilds Phases 1–4 (~40 min). Verified the way P3-T1 and P4-T1 verified
  theirs: every pre-existing table byte-identical except `h5ad_summary.json`'s
  `git_sha` and `config_sha256`.
- `workflow/schemas/config.schema.yaml` is `additionalProperties: false` at the
  root, so `survival` is added to `properties` in the same commit or the config
  does not load.
- **No conda environment is edited.** `scikit-survival` is already in
  `workflow/envs/py-analysis.yaml`, so Phase 5 avoids ADR 0023's
  pin-regeneration transition entirely — and must keep avoiding it. Adding a
  survival package would invalidate all three `*.pin.txt` files and fire the
  software-env trigger on `p0t2_fetch_geo`'s `protected()` outputs.
- The stdlib `.xlsx` reader in `workflow/scripts/external_validation.py` is
  **lifted, not imported**, following the idiom `checkpoint_detection.py:145`
  uses for the detection rule. This project has no cross-script imports,
  Snakemake script-mode sibling imports are not reliable across versions, and
  editing `external_validation.py` would fire its code rerun-trigger. `openpyxl`
  stays absent from `py-analysis.yaml` for the reason that script's docstring
  gives.
- Gate 5 is *"timebox respected"*. This ADR does not commit the project to a
  result; it commits it to which result would count.
