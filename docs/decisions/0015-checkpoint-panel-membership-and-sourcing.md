# ADR 0015 — the checkpoint panel is the plan's nine, and the source is cited for the clinical warrant

**Date:** 2026-09-02
**Status:** Accepted
**Task:** P3-T1
**Follows:** ADR 0011, which set the same contract for `config/signatures.yaml`
**Constrained by:** ADR 0008, which demoted A4 and forbade choosing this panel
from the reconnaissance data

## Context

`config/checkpoints.yaml` was left empty at Phase 0 on purpose. ADR 0008 ran a
reconnaissance detection table over ten canonical checkpoint genes and then said,
in terms:

> **`config/checkpoints.yaml` is still a Phase 3 decision.** Whichever genes go
> in it, the detection reporting in point 3 applies. The table above is
> reconnaissance and confers no licence to populate the panel without asking.

That constraint is the whole reason P3-T1 exists as a separate task. The panel
has to be fixed **before** the detection audit, or the audit stops being a
measurement and becomes a selection.

Two questions had to be answered: which genes, and what counts as a source.

## Decision

### 1. Membership is PROJECT_PLAN §6's nine, unchanged

`CD274`, `PDCD1`, `CTLA4`, `LAG3`, `HAVCR2`, `TIGIT`, `IDO1`, `VSIR`, `CD276`.

All nine resolve against the 18,694 genes in
`results/interim/aoi_normalised.h5ad`; `p3t1_resolve_checkpoints` re-checks on
every run and stops if one ever fails to. Detection over all 120 AOIs at the
project's 2× background rule, for scale: CD276 105, CD274 72, VSIR 64, PDCD1 50,
HAVCR2 31, IDO1 28, CTLA4 23, LAG3 17, TIGIT 15.

**Rejected: adding `BTLA` and `PDCD1LG2` as in-panel negative controls.** ADR
0008 measured BTLA at 0/15 and 0/8 and PDCD1LG2 at 3/15 and 1/8, so both would
serve as a visible assay floor inside the deliverable table. That is a real
benefit and it was considered. It was declined because the nine are what the
plan pre-registered, because ADR 0008's table is exactly the evidence that must
not shape the panel, and because the detection floor is already made visible by
the P3-T2 audit and the P3-T5 figure rather than needing genes chosen to
demonstrate it.

**Rejected: restricting the panel to the three genes ADR 0008 found reliably
detected** (VSIR, CD274, HAVCR2). This is the failure mode P3-T1 exists to
prevent. A panel selected on measured detection would make the P3-T2 audit
circular and would delete the phase's most useful finding — that most of the
clinically important checkpoint panel is not measurable in 8 brain immune AOIs.

Changing membership from here is a stop-and-ask (`CLAUDE.md`), **including
after seeing a result**. Especially then.

### 2. Each gene carries `alias`, `trial_stage`, `rationale`, `source`

§6 P3-T1 asks for "one line of clinical rationale per gene (approved agent /
trial stage), sourced", and its Accept condition is "YAML validates; every gene
annotated". ADR 0014 §4 recorded that P2-T7's criterion failed because
"publication-grade" had no operational definition. So `trial_stage` is a
**controlled enum** — `approved`, `phase_3`, `phase_2`, `investigational` —
rather than prose, which makes "every gene annotated" something the schema can
actually check.

`trial_stage` records the most advanced stage **the cited source establishes**,
not the most advanced trial in existence. Where those differ, `source_note` says
so. TIGIT is `phase_2` because CITYSCAPE is a phase 2 study, even though phase 3
tiragolumab work exists; VSIR is `investigational` because the cited paper is
mechanistic, even though first-in-human agents exist. Understating what a
citation supports is the safe direction.

### 3. `source` admits `PMID:` and `NCT`, and deliberately not `MSigDB:`

`workflow/schemas/signatures.schema.yaml` admits `PMID:` or `MSigDB:`. The
checkpoint schema admits `PMID:` or `NCT[0-9]{8}` instead. A checkpoint's
warrant is clinical rather than a gene-set database — §6 names
ClinicalTrials.gov as a P3-T1 source — so MSigDB would be the wrong currency
here and a gene-set identifier could not support a `trial_stage` claim.

### 4. Three of the nine citations were wrong on the first pass

This is recorded because it is the second time and it will be the third.

ADR 0011 noted that two of four PMIDs proposed for `signatures.yaml` pointed at
unrelated papers. At P3-T1 the same check was run against NCBI eutils before
anything was written, and **three of the nine proposed PMIDs pointed at entirely
unrelated papers**: a digital-health regulatory guide (proposed for TIGIT), a
review of autoimmune encephalitides (VSIR), and a paper on the
hydroaminomethylation of α-olefins (CD276). A fourth error was caught one level
down — the CD276 replacement's first author is Malapelle, not the name first
written.

Every citation in the file has now been verified against PubMed by title,
journal, year **and first author**. The eight surviving attributions were
confirmed the same way.

**A schema can enforce that a source is present. It can never enforce that a
source is correct.** The panel is nine genes and the check took one query.

## Consequences

- The panel is fixed for Phase 3. P3-T2 measures it; it does not choose it.
- Five panel genes (`PDCD1`, `CTLA4`, `LAG3`, `HAVCR2`, `TIGIT`) are also
  members of the `exhaustion` signature. `common.smk` computes and records that
  overlap and `p3t1_resolve_checkpoints` writes it into the summary — recorded,
  never rejected, exactly as `SIGNATURE_OVERLAP` handles overlap within
  `signatures.yaml`. Phase 2 measured `exhaustion` at 1 of 6 genes above
  background in `TIME-B`, so the overlap is also the honest prior for what P3-T2
  will find, and no Phase 2 and Phase 3 result about those genes may be
  presented as independent evidence.
- `config/checkpoints.yaml` is **not** a declared input of `p0t7_assemble_h5ad`,
  for the reason ADR 0011 gives for `signatures.yaml`: editing `config.yaml`
  rebuilds the `.h5ad` and re-runs Phases 1 and 2, and a panel edit must not
  cost that.
- Any future addition must be verified against PubMed before it is written, not
  after. The schema is a floor, not a check.
