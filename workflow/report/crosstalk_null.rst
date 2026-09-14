P4-T4 empirical null calibration for inferred ligand–receptor crosstalk between
adjacent compartments (Aim A5). One panel per adjacency: lung ``L`` ↔
``TIME-L`` across 13 patients, brain ``LB`` ↔ ``TIME-B`` across **8**. The
filled histogram is the permuted-pairing null, the outlined one the observed
\|Spearman ρ\| across all admitted one-to-one direction-rows. The dashed line, if
present, is the \|ρ\| at which the empirical FDR reaches 0.05.

**This is the figure Gate 4 turns on, and the gate is about the null, not the
list.** PROJECT_PLAN §6 P4-T4 is explicit: with roughly a thousand interactions
at n = 13, some \|ρ\| > 0.7 arises by chance, and "without this, the ranking is
uninterpretable". The gate passes if the empirical FDR is computed and
**honoured** — a surviving list and an empty list both pass; a ranked list
without a null does not.

**What is permuted is the pairing, never the expression.** Each permutation
shuffles which patient's immune AOI is matched to which patient's tumour AOI,
**within site**, so both n's are preserved exactly. Every gene keeps its own
values, its distribution, its abundance and its detection status. What is
destroyed is precisely the claim the analysis makes — that the ligand and the
receptor were measured in the same person — which is what makes the overlap of
the two histograms readable as chance rather than as noise.

**The observed distribution sitting on top of the null is the expected
outcome, not a failure.** Aim A5 is **exploratory** (ADR 0014), demoted at
Gate 2 on measured evidence: a ligand–receptor analysis needs the ligand above
background, and Phase 2 measured the ligand side directly at both sites and
found it largely absent. "No LR pair exceeded chance expectation at n = 13" is
an honest and useful result and it passes this gate.

**A null here is an assay-sensitivity limit, never evidence that the crosstalk
is absent.** That distinction is the one a reader gets wrong from silence, and
it is the same one Phases 2 and 3 both turned on. A correlation that cannot be
computed is not a correlation of zero.

**The second null is not on this figure.** ``crosstalk_fdr.tsv`` also carries an
abundance-matched p per row, drawn from random gene pairs in the same
mean-expression decile as the real partners with the true pairing intact. It
asks whether a nomination beats "high-expression genes correlate with
high-expression genes" — the likeliest way this phase produces a wrong answer,
given that the ligand class that survives detection here is the broadly
expressed membrane and matrix one. It **qualifies** a nomination; it is not a
second gate and creates no second FDR family, so putting it on the same axes
would imply a symmetry the design does not have.

**The family is one per site.** The FDR is computed across all primary
direction-rows within a site, lung and brain separately, because a null built
at n = 13 cannot be applied to n = 8 (ADR 0021 §6). Both directions of an
interaction sit in the same family.

**Every quantity behind this figure is restricted to interactions whose two
partners both clear the detection floor**, each in the compartment that partner
is measured in, at the project's one detection rule — ``q3`` above 2.0 × that
AOI's NegProbe-WTX level (ADR 0007). Complex interactions are **not** here: they
are in a declared exploratory table that carries no FDR of any kind
(ADR 0021 §4).

**The binding constraint applies to the whole right-hand panel: ``TIME-B``
n = 8.** Nothing about this design supports a small effect at that n, and the
brain arm is reported alongside lung, never alone.

The word for what this figure is about is **inferred crosstalk between adjacent
compartments**. It is not colocalisation and not spatial proximity: this assay
carries no coordinates, and a ``TIME`` AOI is the PanCK-negative segment of an
ROI sited in a CD45-rich region — not a CD45-sorted population (Q2's caveat).
