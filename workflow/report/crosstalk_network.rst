P4-T7 inferred ligand–receptor crosstalk between adjacent compartments
(Aim A5). One panel per adjacency: lung ``L`` ↔ ``TIME-L`` across 13 patients,
brain ``LB`` ↔ ``TIME-B`` across **8**. Nodes are genes, placed on the left if
measured in the tumour compartment and on the right if measured in the immune
one; arrows run ligand → receptor. The layout is bipartite because that is what
the design is — a spring layout would invent a topology the data does not have.

**This is not colocalisation and it is not spatial proximity.** This assay
carries **no coordinates**: an AOI is a region a pathologist placed, and there
is no spot grid. An edge here means *the two compartments' expression correlated
across patients*, which is an inference from a pinned interaction database plus
a correlation — not a measured spatial relationship. The phrase is **"inferred
crosstalk between adjacent compartments"**, and ``p4t7_language_audit`` enforces
it across every committed file rather than trusting anyone to remember.

**Edge colour is ``RdBu_r``, and that is a deliberate reuse rather than a
borrow.** This project's colour language is already fixed: ``RdBu_r`` means a
signed effect (P2-T7) and ``PRGn`` centred on the detection multiple means a
background ratio (P3-T5, ADR 0016 §5). A Spearman ρ **is** a signed effect on
[−1, 1] — the same quantity class ``RdBu_r`` was established for — so speaking
that language is consistent, and inventing a third one here would imply ρ is a
third kind of quantity. Line width also tracks \|ρ\|.

**Edge style, not edge presence, carries the FDR.** Solid means the edge clears
the pre-registered empirical FDR from P4-T4; dashed means it does not. Drawing
only survivors would produce a blank panel for an outcome Gate 4 explicitly
licenses — "no LR pair exceeded chance expectation at n = 13" — and a blank
panel reads as a broken rule rather than as a finding. **A dashed edge is not a
weak result; it is a ρ that chance explains.**

**The strongest edges are drawn whether or not any of them is real.** The panel
shows the top 25 by \|ρ\| per site. That is a display cap and decides nothing
about which edges are real — the style does, from the pre-registered FDR.

**A node is a (gene, compartment) pair, not a gene.** The same gene measured in
the tumour and immune compartments is two different measurements, and merging
them would draw an edge across a compartment it never touched.

**Everything here is restricted to interactions whose two partners both clear
the detection floor**, each in the compartment that partner is measured in, at
the project's one detection rule — ``q3`` above 2.0 × that AOI's NegProbe-WTX
level (ADR 0007). **Complex (multi-subunit) interactions are absent by design**:
they carry no FDR of any kind (ADR 0021 §4), so they have no style to draw.

**Aim A5 is exploratory (ADR 0014).** No edge on this figure may be a headline
claim. The ligand side of this assay is largely at background — Phase 2 measured
it directly at both sites — so an absent or dashed edge is an
**assay-sensitivity limit, never evidence that the crosstalk is absent**.

**The binding constraint applies to the whole right-hand panel: ``TIME-B``
n = 8.** The brain arm is reported alongside lung, never alone, and the
lung-vs-brain contrast detects roughly 1.1–1.3 SD at 80% power (P0-T8).

A ``TIME`` AOI is the PanCK-negative segment of an ROI sited in a CD45-rich
region — **not** a CD45-sorted population (Q2's caveat). The compartment label
is not a cell-type label, so a node on the right is a gene measured in that
segment, not a gene expressed by a lymphocyte.
