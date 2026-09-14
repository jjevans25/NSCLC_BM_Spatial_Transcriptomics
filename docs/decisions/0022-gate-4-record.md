# ADR 0022 — Gate 4 record: the null was computed and honoured, and nothing survived it

**Date:** 2026-09-14
**Status:** Accepted
**Task:** P4-T7 (records Phase 4)
**Gate:** this is the Gate 4 record
**Follows:** ADR 0014 (A5 demoted at Gate 2), ADR 0021 (Phase 4 pre-registration)

## 1. Gate 4's question, and the answer

> **GATE 4** — Empirical FDR computed and nominations survive it. If nothing
> survives, report that: "no LR pair exceeded chance expectation at n=13" is an
> honest, useful, publishable-to-blog result, and it's a better outcome than a
> ranked list you can't defend.

**The FDR was computed and honoured. Nothing survived it. The gate passes on
its own second clause.**

| | lung | brain |
|---|---|---|
| paired patients | 13 | **8** |
| primary direction-rows tested | 268 | 260 |
| **clearing the empirical FDR at 0.05** | **0** | **0** |
| smallest empirical FDR | 0.497 | 0.228 |
| largest observed \|rho\| | 0.742 | 0.952 |

The observed distribution sits **on top of** the null rather than beside it —
median \|rho\| 0.231 observed against 0.203 null in lung, 0.262 against 0.286 in
brain (`crosstalk_null.tsv`, `crosstalk_null.png`). The strongest single result
in each arm illustrates why the null matters: brain's `SEMA4C → PLXNB2` at
rho −0.952 on n = 8 is the kind of number a ranked list would lead with, and
10,000 permutations of the pairing produce values that large routinely at that n.

**The gate turns on P4-T4, not P4-T5**, and that distinction did the work it was
meant to: `p4t5_nominations` was written with no fallback path and no second
threshold, so an empty table was reachable rather than something to be escaped.

**This is an assay-sensitivity limit, never evidence that the crosstalk is
absent.** A correlation that cannot be resolved at n = 13 and n = 8, on an assay
whose ligand side is largely at background, is not a correlation of zero. That
sentence is the one a reader gets wrong from silence, and it is the same
distinction Phases 2 and 3 both turned on.

## 2. The finding the gate question did not ask about, and it is the larger one

**The detection gap by interaction class is ADR 0014's claim made numeric.** Of
the interactions CellChatDB v2 contains, the fraction admitted — both partners
above background, each in the compartment it is measured in:

| class | brain | lung |
|---|---|---|
| **ECM-Receptor** | **45.3%** | **39.6%** |
| Cell-Cell Contact | 13.7% | 17.3% |
| **Secreted Signaling** | **6.4%** | **5.3%** |

ADR 0014 demoted A5 on the claim that "secreted ligands and chemokines are
systematically undetected", measured on three cytokine signatures. **It
generalises, and the gradient is eightfold**: the matrix-associated axis is
measurable in roughly two of five interactions, the soluble axis in one in
twenty. `docs/NEXT_STEPS.md`'s reconnaissance predicted the direction; this is
the measurement, on 2,239 interactions rather than a handful of marker genes.

**That bias is the phase's result, and it has the same shape as Gate 3's** —
where the audit was larger than the one interpretable estimate. What Phase 4
can say is not "there is no crosstalk" but "this assay resolves matrix and
membrane interactions and does not resolve secreted ones, and no surviving
signal was detectable in the part it does resolve."

## 3. Background travels with the patient, which is why Null B existed

ADR 0018 established that detection and expression move opposite ways under one
background gradient. Phases 2 and 3 met that as a *group-contrast* problem.
**Phase 4 correlates across patients, so the relevant question was different and
had to be measured rather than inherited**: does a patient's tumour AOI share
its background with that patient's immune AOI?

| | negprobe_log2 | gene_detection_rate |
|---|---|---|
| brain (n = 8) | ρ **+0.548** | ρ **+0.667** |
| lung (n = 13) | ρ **+0.352** | ρ +0.231 |

**Three of four exceed the 0.3 threshold.** AOI quality is a patient-level
property, so two genes that merely track background would correlate across
patients for no biological reason at all. P3-T4 found the *site* gradient
negligible; that is a different gradient and inheriting its verdict would have
been assuming exactly what needed checking.

The pre-registered response (ADR 0021 §6) was the abundance-matched Null B, and
it was decided **before** this number was seen — which is the point. 14 rows per
site beat their abundance-matched null; **zero rows beat both nulls**, because
zero beat the first.

## 4. What was found by a check rather than by eye

Recorded because each is a defect a plausible-looking result would have hidden.

- **The FDR step-down was wrong under ties.** Spearman on 8 patients takes 42
  distinct values across 260 rows. Ranking rows 1..N and dividing by the running
  index gives tied rows different denominators, which is simply the wrong count:
  at a threshold, every row reaching it is observed. Now computed once per
  distinct \|rho\| and mapped back by value. **The monotonicity assertion caught
  it**; without that assertion the phase would have shipped a plausible ranked
  table.
- **`CD96–NECTIN1` was entered into the comparator backwards** and reported as
  absent from a database that contains it. The comparator measures co-expression
  within a spot, which is **undirected**, so the lookup now matches either
  orientation — the orientation was a hand-transcription step with nothing
  checking it.
- **Both figures clipped on first render.** The network's right-hand labels
  showed "CEACAM" for both CEACAM1 and CEACAM5, and the colourbar label rendered
  off-canvas. Fixed with the `main | colourbar` gridspec `checkpoint_dotplot.py`
  already uses. **Caught only by opening the PNGs** — P3-T5's lesson, which the
  plan warned about and which still recurred.

## 5. Two defects in PROJECT_PLAN §6, recorded because the plan is gitignored

Same reason ADR 0008 and ADR 0014 §4 exist.

- **P4-T6 names four comparator pairs; the source reports five.** De Zuani et
  al. list `NRP1-VEGFA`, `NECTIN2-TIGIT`, `LGALS9-HAVCR2` **and `CD96-NECTIN1`**
  as enriched, with `PD1-PDL1` absent from that list. `CD96-NECTIN1` is included
  here: a comparator list that drops one of the source's own results is a
  comparison chosen after the fact.
- **P4-T7's "grep the entire repo and fix every hit" is unrunnable as written.**
  It fails on `CLAUDE.md`'s own hard constraint 6, on
  `config/compartment_map.yaml`'s header and on `docs/limitations.md` §3 —
  because **the files that state the rule have to name what they forbid.** A
  check that cannot tell a prohibition from a claim would force the rule's own
  statement to be deleted to make the check pass. `p4t7_language_audit`
  classifies the two and fails the build only on the latter: **0 violations, 22
  prohibitions, 4 quotations.** The criterion should read "no committed file
  *asserts* a spatial relationship", not "no file contains the string".

  **A third category emerged while the audit was being written, and it is the
  one that matters most.** `config/ligand_receptor.yaml` quotes the P4-T6
  comparator's own sentence, which says "colocalized" — because that is what
  De Zuani et al. wrote, and ADR 0021 §2 quotes rather than paraphrases
  precisely so the verdict is not written by us. **Editing a quotation to
  satisfy this project's style rule would falsify the source.** So a match
  inside a declared quotation field is counted separately and does not fail the
  build. The cue is deliberately narrow — the field must be named for quoting,
  or the surrounding lines must say the text is quoted — so it cannot be reached
  by writing a claim and calling it a quote.

  The audit carries its own rerun trigger, a digest of every file it scans
  computed at parse time (`_tracked_text_digest` in `common.smk`). Without one
  it would run once and never again, and a check that goes stale the moment
  someone edits a doc is worse than no check: the green output keeps asserting
  something nobody re-tested.

## 6. Two things that did not go as the handoff expected

- **`VEGFA–NRP1` is absent from CellChatDB v2 entirely.**
  `docs/NEXT_STEPS.md` called it "the one comparator where both partners are
  detected nearly everywhere and the concordance check is genuinely
  informative". Both partners *are* detected — but CellChatDB carries VEGFA to
  FLT1 and to KDR, with NRP1 appearing only as a SEMA3x co-receptor. The
  comparator used CellPhoneDB, which pairs them. **This is a database
  difference, not a detection one, and not a biological one** — exactly the
  category ADR 0021 §2 said to keep separate, arriving in the one place nobody
  predicted it. Of the five comparators, four are unevaluable for four
  different reasons: not in our database (2), not on the panel (1, `LGALS9`),
  below the floor (1 arm of `NECTIN2–TIGIT`).
- **Phase 4 was not starved of genes, only of the right ones.** 1,062 of the
  1,148 symbols CellChatDB names are on the panel, and 528 primary
  direction-rows were admitted. The constraint was never coverage of the
  database; it was which *class* survives background.

## Consequences

- **Phase 4 is complete and A5 remains exploratory.** Nothing here may be a
  headline claim (ADR 0014), and the honest one-sentence result is: *no inferred
  ligand–receptor pair exceeded chance expectation in either adjacency —
  lung n = 13, brain `TIME-B` n = 8 — and the assay resolves matrix-associated
  interactions roughly eight times as often as secreted ones.*
- **`docs/limitations.md` §10 gains Phase 4's version of the detection limit**,
  which is the assay-wide, class-resolved form of what Phase 2 measured on six
  signatures and Phase 3 on nine genes.
- **The empirical-FDR machinery is reusable and now tested.** The identity
  permutation reproduces P4-T3b's rho to 1.1e-16 and 2.2e-16.
- **A5's demotion is vindicated rather than merely honoured.** ADR 0014 demoted
  it on the prediction that the ligand side would not support the analysis;
  Phase 4 ran the analysis properly and measured exactly that, which is the
  methods value the demotion decision said was real.
- Phase 5 (survival, Q5 = yes) and Phase 6 (clean-room repro, DOI) remain. The
  slip rule in PROJECT_PLAN §timeline says cut Phase 5 before Phase 6.
