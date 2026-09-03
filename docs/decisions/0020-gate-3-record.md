# ADR 0020 — Gate 3 record: one interpretable pattern, and an audit that is the larger result

**Date:** 2026-09-03
**Status:** Accepted
**Task:** P3-T7 (records Phase 3)
**Gate:** this is the Gate 3 record
**Follows:** ADR 0008 (A4 demoted at Gate 0), ADR 0016 (Phase 3 pre-registration),
ADR 0018 (Option C), ADR 0019 (the P3-T6 slider)

## 1. Gate 3's question, and the answer

> At least one checkpoint gene shows an interpretable compartment-resolved
> pattern with a CI you're willing to show a stranger. If everything is
> not-assessable, that's still a legitimate finding.

**Both clauses are satisfied, and they were satisfied by different tasks.**

**First clause — `CD274` (PD-L1) in the lung immune compartment.**
**+0.612 SD [+0.361, +0.863], q = 0.0001**, immune minus tumour within lung,
`expression ~ compartment + (1|patient_id)`, 45 AOIs across 30 patients,
detected 19/30 in `L` and 13/15 in `TIME-L` (`checkpoint_carrier_models.tsv`,
P3-T3, primary). It survives both pre-registered sensitivities — batch +0.599,
background +0.635 — the fit is not singular, and the `statsmodels.MixedLM`
refit agrees with `lme4` to 1.1e-6.

Two properties make it interpretable rather than merely significant. The
**background adjustment raises it** (+0.023), and the immune compartments are
precisely where background is higher, so the one gradient that could have
manufactured an immune-side enrichment is measured and is not doing so
(ADR 0018). And **`CD276`, the only other gene clearing the primary
restriction, shows no compartment preference at either site** (+0.147 lung,
+0.159 brain, both null), so this is not the panel-wide immune-side drift an
artefact would produce.

**Second clause — the audit, and it is the larger result.** **40 of 63 gene ×
compartment cells are below the pre-registered floor, and the panel clears it in
`TIME-L` alone** (`checkpoint_detection.tsv`, P3-T2). That is a legitimate
finding *because the detection table licenses it*, and the table is now
cross-validated three ways: against `export_tsv.py`'s independent
implementation, against ADR 0008's Gate 0 reconnaissance, and — since P3-T2b —
against the source publication's own deposited Source Data, where **all
2,243,280 values of `layers['q3']` match exactly** (ADR 0017).

## 2. What the gate deliberately does not rest on

**The lung-vs-brain shift is null in every assessable gene, and those nulls are
uninformative rather than negative.** The largest is `CD274` at
−0.260 [−0.595, +0.075], q = 0.462, brain minus lung within the immune
compartment, **`TIME-B` n = 8**. The design detects roughly **1.1–1.3 SD at 80%
power** (P0-T8); every P3-T4 estimate sits far below that. Reporting them as
evidence of no difference would fail this gate more surely than having no result
at all.

**Seven of nine genes never entered a primary fit**, and the exploratory table
that carries them has no q-value and licenses nothing (ADR 0018). Its largest
apparent effects — `CTLA4` +1.425, `VSIR` +1.336, `TIGIT` +1.036 in lung — are
every one of them a gene detected in 0–9 of 30 tumour AOIs.

**`TBME`, `mLN` and `BC` are audited and plotted but never modelled**, because
each sits wholly within one DSP run and no covariate recovers that
(ADR 0009 §3, ADR 0016 §2). §6's "tumour vs. immune vs. glial" is answered as
tumour vs. immune, and the narrowing is recorded rather than silent.

**A4 remains exploratory** (ADR 0008). Nothing above is a headline claim, and
the clinical reading in `docs/analysis-notes.md` is framed as reconnaissance
that would justify a targeted study.

## 3. What Phase 3 established that Gate 3 did not ask for

1. **This project's inputs are demonstrably the published inputs** (ADR 0017).
   Every check before P3-T2b was internal.
2. **The 119-vs-120 gap is closed — the paper dropped `TBME15b`**, recorded for
   two phases as not recoverable from GEO, which was true of GEO.
3. **Detection and expression move opposite ways under the same background
   gradient** (ADR 0018), so the argument for relaxing the pre-registered
   restriction was inverted. Measured, not reasoned.
4. **The background artefact is a compartment-contrast problem, not a
   site-contrast one** — the opposite of what ADR 0008 anticipated. Every P3-T4
   gradient is `negligible` (site background deltas +0.03 and −0.06 against
   0.31–0.40 for the compartment contrasts).

## Consequences

- **Gate 3 is PASSED.** Phase 3 is complete: P3-T1 … P3-T7.
- **The reportable Phase 3 result is one gene in one compartment in one site**,
  plus an audit saying why it is one. Any future writing that promotes a second
  is reaching into the exploratory table, which ADR 0018 forbids.
- **`CD274`'s compartment localisation is a hypothesis about source, not about
  cell type.** A `TIME` AOI is the PanCK-negative segment of an ROI sited in a
  CD45-rich region, not a CD45-sorted population (Q2), and safeTME puts `TIME-L`
  at a mean 0.198 fibroblasts beside 0.200 macrophages. "Not tumour-cell-derived"
  is supported; "myeloid-derived" is not.
- **Phase 4 (A5) inherits the reason it was demoted.** The ligand side is
  largely undetected (ADR 0014), and Phase 3 has now confirmed the same pattern
  gene by gene for a panel overlapping `exhaustion` in five of nine genes.
