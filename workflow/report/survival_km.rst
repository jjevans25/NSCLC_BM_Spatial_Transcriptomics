P5-T3 Kaplan–Meier by **median signature split, within arm** (Phase 5,
**stretch**; Gate 5 is *"timebox respected"*, not "a result was found"). Rows are
the six signatures of ``config/signatures.yaml``; columns are the two arms —
lung (``TIME-L``) and brain (``TIME-B``), never pooled (ADR 0024 §2). Curves are
the ``ssgsea`` primary; ``zscore`` is a declared sensitivity and lives in the
tables only.

**The split is** ``score > median``, **computed within the arm being split, and
the tie side is declared.** ``high`` is strictly above that arm's median,
``low`` at or below it. Stated rather than defaulted, because at n = 7 the median
*is* an observation and which side it falls on changes the arm sizes. This gives
**6/6 in lung and 3/4 in brain**, which ADR 0025 §3 recorded before any curve
existed. It is not a tertile, not an optimal cutpoint and not a
maximally-selected rank statistic — those are the standard routes to a separation
that does not replicate, and at these counts they are not defensible.

**Every p-value on this figure is RAW AND UNCORRECTED, and that is a decision,
not an omission.** ADR 0026: the multiplicity family spans this study *and* the
TCGA LUAD cohort (P5-T4), because "does the same signature stratify an
independent cohort" is one question asked twice. It therefore cannot be corrected
until both have been fitted, so **P5-T5 applies Benjamini–Hochberg once over the
union** and only that q is reportable. A q-value here would be a corrected value
that nothing corrected. ADR 0024 §8 originally declared a family of 12; that
arithmetic was wrong — it counted cohorts and forgot that this cohort has two
arms — and ADR 0026 corrects it to the **assessable** cells, **8 here** (lung 5,
brain 3) plus TCGA's.

**Four of twelve panels are hatched, and they are as much the finding as the
curves.** ``exhaustion`` and ``tls`` in brain, and ``myeloid_m1`` at **both**
sites, sit below Phase 2's coverage floor — the fraction of a set's genes
detected in at least half that site's AOIs, ``detected_in_aoi_fraction: 0.5``.
They are **not fitted**, and the finding is **"not assessable in that arm"**,
never "no prognostic association". ADR 0008's rule is not suspended by a change
of outcome variable: **a score built from genes at background measures
background, whatever it is regressed against.** They are hatched rather than
omitted so a reader can tell *not tested* from *tested, null* — the same
reasoning that hatches 40 of 63 cells in ``checkpoint_dotplot.png`` (P3-T5).

**The binding constraint is on every panel:** ``TIME-B`` **n = 8** (hard
constraint 8), of which **7 enter a fit**. The eighth is patient 35, censored
with **no follow-up time in either endpoint column** — ``Alive`` is a status, not
a duration, and the brain cell is blank, and Supplementary Data 1 carries no
last-contact date. It is excluded from fits and recorded in
``survival_excluded.tsv`` rather than imputed, because administrative censoring
at the cohort maximum would invent an observation *and* place it exactly where a
median split is most sensitive (ADR 0025). Both numbers belong in any reported
sentence: 8 is the cohort, 7 is what was fitted.

**Every patient entering a fit has an event.** With the one censored case
removed, the arms carry **no censoring at all** — 12 events in 12 lung patients
and 7 in 7 brain ones — so these curves fall to zero and are a comparison of
ordered survival times rather than a censored-data estimate.

**Read the confidence intervals, not the curves.** At n = 3 versus 4 a hazard
ratio's 95% interval spans roughly **thirty-fold**. A separation that looks
dramatic and an interval that spans 1 are the same observation drawn two ways,
and the interval is the honest one. **A null here is UNINFORMATIVE, NOT
NEGATIVE** (ADR 0024 §7, fixed before the first curve so it could not be softened
afterwards) — the same status the 1.1–1.3 SD power floor gives a Phase 2 null
(P0-T8) and ADR 0022 gives Phase 4's empty nomination table.

**A separation at 3 versus 4 patients describes this cohort. It is never an
inferential claim, however small its p-value, and no Phase 5 result may be a
headline claim.**

**What a compartment label is not.** A ``TIME`` AOI is the PanCK-negative segment
of an ROI sited in a CD45-rich region — **not** a CD45-sorted population (Q2).
The signature scores behind this split are compartment-level, and the compartment
label is not a cell-type label.
