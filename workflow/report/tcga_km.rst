P5-T4 **TCGA LUAD external check** on Phase 5's null (Phase 5, **stretch**).
Six panels, one per signature of ``config/signatures.yaml``, Kaplan–Meier by
**median split of the ``ssgsea`` score within this cohort** — the same split rule
and the same tie side (``score > median``) that P5-T3 applied, because
**comparability is the entire point**: a different estimator would confound
"different cohort" with "different method", which is the inference this task
exists to license. ``zscore`` is a declared sensitivity and lives in the tables.

**Why this cohort exists in the plan.** P5-T3 found no immune signature
stratified survival here — eight assessable cells, every hazard-ratio interval
spanning 1, at lung n = 12 and brain n = 7. PROJECT_PLAN §6 P5-T4: *"a negative
here mostly tells you about your n, not about the biology."* At ~500 patients
TCGA can separate a power limit from an absent effect. **A null in TCGA means
something a null in this study does not** — and that asymmetry, not any single
p-value, is the deliverable.

**Every p-value on this figure is RAW AND UNCORRECTED.** The Phase 5
multiplicity family is **8 (P5-T3) + 6 (here) = 14** and **P5-T5 applies
Benjamini–Hochberg once over all of it** (ADR 0026 §2, ADR 0027 §2); only that q
is reportable. A q-value here would be a corrected value that nothing corrected.
Where the log-rank and the hazard-ratio interval disagree, **the interval
governs** — the log-rank is a score test and the interval a Wald interval, and a
p below α beside an interval spanning 1 is not a detected effect. Any such panel
says so on its face.

**All six signatures are tested, and that is a decision.** Phase 2's coverage
floor is ``q3 > 2 ×`` that AOI's ``NegProbe-WTX`` (ADR 0007) — a property of the
**GeoMx assay**. TCGA is bulk RNA-seq with **no negative-probe channel**, so
nothing in it corresponds to that threshold, and restricting to the five
signatures assessable in the GeoMx lung arm would import a limitation the
external cohort does not have. ``myeloid_m1`` is therefore tested here and its
result is reportable, but it is labelled **TCGA-only: it has no GeoMx
comparator**, because it was not assessable at either GeoMx site. What this
figure computes instead is a **gene-presence** table — whether a signature's
genes are in the Xena matrix at all. **Presence is not detection**: "not in this
matrix" and "at background in this assay" are unrelated limitations and must
never be read as the same measurement.

**THIS IS NOT THE SAME MEASUREMENT, and the comparison is bounded by that.** A
GeoMx ``TIME-L`` score comes from the **PanCK-negative segment** of an ROI sited
in a CD45-rich region — and per Q2's live caveat that is **not** a CD45-sorted
population either. A TCGA score comes from **whole bulk tumour**: tumour cells,
stroma and immune cells together, in one homogenate. "The same signature" is
computed on two different things. **A disagreement between the cohorts is
therefore not necessarily a disagreement about biology, and an agreement is not
necessarily a replication.** This is the same class of caveat as Phase 2's
reference-matrix confound, and it is a property of the design, not a shortcoming
of either cohort.

**There is no brain comparator, at all.** TCGA LUAD is primary lung. The brain
arm — the one carrying ``TIME-B`` **n = 8**, the binding constraint on the whole
project — receives **no external check from this task**. This figure bounds the
**lung** null only, and nothing on it speaks to the brain metastasis compartment.

**The cohort is constructed, not taken as given.** TCGA barcode positions 14–15
are the sample type: only ``01`` primary solid tumour is kept, and this matrix
contains ``11`` solid-tissue normals and ``02`` recurrent tumours that would
otherwise be scored as though they were primary tumours. Cases with a missing or
non-positive follow-up time are **excluded and counted, never imputed** — the
same rule ADR 0025 applied to this study's own censored patient. ``OS.time``
arrives in **days** and is converted at ``365.25 / 12``, a factor declared in
``config/tcga.yaml`` rather than inlined. Endpoints are the Liu et al. 2018
TCGA-CDR curated set, not raw TCGA follow-up.

**The source is UCSC Xena, not cBioPortal, and that was forced.** ADR 0024 §1
preferred cBioPortal because the source publication's own TCGA analyses came from
there. Measured 2026-09-15, its bulk files return **403** (datahub) and **404**
(LFS mirror); only the live REST API answers, and a per-study JSON endpoint
cannot carry a pinned SHA-256 under ADR 0005. ADR 0024 §1 pre-authorised this
fallback in terms, and ADR 0027 records the substitution with the evidence.

**A continuous Cox on the score is reported beside the split**, in the tables,
as a declared sensitivity (ADR 0012's shape). At this n a median split discards
real information, and the continuous fit is where the cohort's power actually
shows — so **a null split beside a significant continuous fit would be a
statement about the split, not about the biology**, and must be readable as such.
Neither is the tiebreak for the other.
