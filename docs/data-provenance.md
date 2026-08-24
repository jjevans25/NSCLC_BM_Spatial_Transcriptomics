# Data provenance and the answers to Q1–Q5

Answers to the five open questions in `Markdowns/PROJECT_PLAN.md` §3. Owner task
is **P0-T4**; Q2, Q4 and Q5 are answered here ahead of it because they need only
the paper and the GEO metadata, and Q2 gates the project's vocabulary (§3: "Q2
gets a written answer in `docs/data-provenance.md` before Phase 2 starts").

Each question carries a **Status**. `Resolved` means answered with a citation or
an explicit "not available in GEO" statement — P0-T4's acceptance criterion.
`Provisional` means the evidence is strong but comes from metadata rather than
from the data itself, and P0-T4 must confirm it against the downloaded file
before anything is hardened against it.

## Sources

| | |
|---|---|
| GEO series | [GSE200563](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200563) — public 2022-04-25, last updated 2023-02-26 |
| Publication | Zhang Q, Abdo R, Iosef C, Kaneko T, Cecchini M, Han VK, Li SS. *The spatial transcriptomic landscape of non-small cell lung cancer brain metastasis.* Nat Commun 13:5983 (2022) |
| PMID | 36216799 |
| PMC | [PMC9551067](https://pmc.ncbi.nlm.nih.gov/articles/PMC9551067/) (open access, CC-BY) |
| DOI | 10.1038/s41467-022-33365-y |
| Platform | GPL21697 (NanoString GeoMx DSP, WTA; 18,694 genes) |
| Accessed | 2026-08-24 |

Quotations below are verbatim from the PMC full text or from the GEO SOFT
records, as marked. Nothing here is paraphrase.

---

## Q1 — Is the GEO matrix raw counts, Q3-normalised, or already log-transformed?

**Status: Provisional — Q3-normalised. Confirm empirically at P0-T4.**

Every one of the 120 GSM records carries the same two `!Sample_data_processing`
lines:

> GeoMx® DSP counts from each AOI were scaled to the 75th percentile of expression

> Supplementary files format and content: tab-delimited text files include 75th
> percentile of expression values for each Sample

This matches the paper:

> The sequencing data were normalized using the third quartile expression (Q3)
> and validated to ensure quality; and the 0.75 quantile-scaled data were used
> for all subsequent analysis

So `GSE200563_processed_data.txt.gz` is **Q3-normalised, not raw counts, and not
log-transformed**. Raw counts are separately available as per-AOI `.dcc` files
inside `GSE200563_RAW.tar` (see Q3).

Confirm at P0-T4 with the checks §3 names — values non-integer, no negatives,
distribution right-skewed rather than symmetric, and the column-wise third
quartile near-constant across AOIs (that last one is the direct test: Q3 scaling
is exactly the operation that makes it so).

**Consequence, once confirmed:** `config.yaml → normalisation.method` moves
`q3` → `verify`; the schema already permits that value. Do not make the change
on the metadata alone — the point of Q1 is to look at the file.

## Q2 — Antibody-segmented compartments, or geometric ROIs with descriptive names?

**Status: Resolved — antibody-segmented. `compartment` is retained.**

This was the highest-impact unknown in the project. It resolves in favour of the
existing vocabulary: **no `region` rename, and Aim A4 is not reframed.**

Collection was by antibody-driven segmentation. From the GeoMx DSP methods, and
repeated verbatim in all 120 GSM records under `!Sample_extract_protocol_ch1` —
independent corroboration, since the sample records were deposited separately
from the manuscript:

> Oligoes from PanCK+ and PanCK- regions were collected separately by UV-cleavage.

ROI *placement* was a separate, earlier step: pathologist review of H&E plus
immunofluorescence.

> Regions-of-interest (ROI) for DSP were annotated based on histology by a
> pathologist and immunofluorescence staining with the morphological markers
> PanCK (for epithelial cells), CD45 (for hematopoietic cells), and GFAP (for
> brain cells).

> sections of the TMA were stained simultaneously with antibodies against the
> leukocyte marker CD45 to demarcate the tumor-immune microenvironment (TIME),
> the epithelial cell marker PanCK to mark the tumor cores (L and LB), GFAP
> (glial fibrillary acidic protein) to identify the tumor brain microenvironment
> (TBME), and SYTO83 to mark the cell nuclei

The answer is therefore: **two-channel antibody segmentation (PanCK⁺ / PanCK⁻)
within marker-guided ROIs.**

### The caveat, which constrains what may be claimed

**PanCK was the only collection mask.** CD45 and GFAP were guides for where a
human placed the ROI — they were not segmentation channels, and no AOI was
collected on a CD45 or GFAP mask.

Concretely, a `TIME-L` AOI is *the PanCK-negative segment of an ROI sited by a
pathologist in a CD45-rich region*. It is **not** a CD45-sorted population. The
same holds for `TBME` and GFAP.

This is a materially weaker warrant than "CD45⁺-segmented immune compartment",
and the difference matters wherever a compartment label is treated as a cell-type
label. Two consequences that are not optional:

1. **P0-T6's marker sanity check is doing real work, not ceremony.** The
   compartment labels rest on a human's placement judgement plus a PanCK⁻ mask;
   `PTPRC` enrichment in `TIME` AOIs and `GFAP`/`AQP4` enrichment in `TBME` AOIs
   are the evidence that the placement worked. If it fails, the labels are
   wrong, not the marker panel.
2. **Deconvolution (`SpatialDecon`) is not confirmatory here.** A PanCK⁻ segment
   contains whatever non-epithelial cells were in the ROI — immune, stromal,
   glial, vascular. Cell-type fractions describe that mixture; they do not
   validate the label.

Same reasoning as hard constraint 6, which forbids "colocalisation": state the
claim the assay supports, not the one the label suggests.

### The caveat was tested, and the labels hold (P0-T6)

Because `TIME` and `TBME` rest on human ROI placement rather than a sorting
mask, the transcriptome is the only independent evidence that the placement
worked. P0-T6 supplies it — all four gating criteria pass, on all 120 AOIs,
against a `NegProbe-WTX` background floor of log2 5.23:

| Criterion | median (high) | median (rest) | Δ log2 | 95% CI | n |
|---|---|---|---|---|---|
| tumour epithelial-high (EPCAM) | 7.92 | 5.50 | +2.42 | [+1.96, +2.66] | 70 / 50 |
| tumour epithelial-high (KRT19) | 9.02 | 6.06 | +2.97 | [+2.51, +3.53] | 70 / 50 |
| `TIME` `PTPRC`-high | 8.19 | 6.05 | +2.14 | [+1.51, +2.38] | 23 / 97 |
| `TBME` `GFAP`-high | 11.89 | 5.84 | +6.05 | [+1.67, +6.64] | 20 / 100 |

Two details worth carrying forward:

- **Epithelial signal in non-tumour compartments is at background, not merely
  lower.** `BC` sits on the floor for both EPCAM and KRT19. The PanCK⁻ mask
  did what it claims.
- **`TBME` GFAP (11.89) exceeds `BC` GFAP (10.18).** Tumour-adjacent brain is
  *more* astrocytic than normal brain, consistent with the reactive astrocytes
  the source paper reports — the label is not merely correct, it is capturing
  the biology it was drawn for.

Evidence: `results/figures/marker_sanity.png`,
`results/tables/marker_sanity_verdict.tsv`.

## Q3 — Is slide / TMA / batch identifiable per AOI?

**Status: Resolved — yes, from the DCC filenames. Extracted by P0-T3 into
`samples.tsv` (`dsp_run`, `dsp_well`) and `results/tables/batch_crosstab.tsv`.
But it is partially confounded with compartment — see below.**

**Not from the GEO sample metadata.** The SOFT family file carries exactly one
characteristics field per sample:

```
!Sample_characteristics_ch1 = cell type: <compartment description>
```

There is no slide field, no TMA block field, no batch field, and no
annotation/LabWorksheet file in the series.

**But it is recoverable from the raw files.** The series `filelist.txt` names one
`.dcc.gz` per AOI, and the DSP run identifier is embedded in each filename
(`GSM6573697_DSP-1012300141221-A-A02.dcc.gz`). Across the 120 DCCs:

| DSP run ID | AOIs |
|---|---|
| `DSP-1012300141221` | 91 |
| `DSP-1012310141221` | 29 |

Plus a plate well position (`A02`, `A03`, …) per AOI. **This is why P0-T2
acquires `GSE200563_RAW.tar` and `filelist.txt`, which the original task sketch
did not scope** — without them Q3 is unanswerable and the QC model has no batch
term.

**This is not the paper's TMA blocks.** The paper describes four:

> NSCLC patients with metastases to the brain (n = 44) were represented in four
> tissue microarray (TMA) blocks (LB-D1 to D4)

Two DSP runs ≠ four TMA blocks, and there is no mapping between them in GEO. So
the available batch variable is the **sequencing/DSP run**, not the TMA block.

The plate letter (`A` / `B`) is perfectly collinear with the run id, so this is
one two-level variable, not two.

### The confounding, which is the part that matters

P0-T3 joined the filenames to the design (`results/tables/batch_crosstab.tsv`):

| Code | `DSP-1012300141221` | `DSP-1012310141221` | Total |
|---|---|---|---|
| L | 22 | 8 | 30 |
| LB | 21 | 6 | 27 |
| mLN | **13** | **0** | 13 |
| TBME | **20** | **0** | 20 |
| TIME-L | 9 | 6 | 15 |
| TIME-B | 6 | 2 | 8 |
| BC | **0** | **7** | 7 |
| **Total** | **91** | **29** | **120** |

**`mLN`, `TBME` and `BC` each sit entirely inside one run.** For those
compartments batch and biology are inseparable: any `BC`-vs-anything difference
is also a run-B-vs-run-A difference, and no model can tell them apart. This is a
limitation to state (P0-T8), not a defect to correct.

**The core comparison survives.** `TIME-L` (9/6) and `TIME-B` (6/2) both span
both runs, so a batch term is estimable for the primary contrast. Note it rests
on only **2 `TIME-B` AOIs** in the smaller run, so that estimate is thin — worth
a sensitivity check at P0-T5 rather than blind inclusion.

Aim A4 and anything else resting on `TBME` inherits the `TBME` confounding in
full.

## Q4 — Does per-AOI nuclei count / surface area survive into GEO metadata?

**Status: Resolved — no.**

Not available in GEO. `cell type` is the only per-sample characteristics field
(see Q3), and the series contains no DSP annotation worksheet. The paper reports
only a cohort-level average:

> A total of 119 ROIs (average 0.2 mm² each) were analyzed.

**Consequences:**

- Models **cannot** be weighted by AOI area or nuclei count. State this in
  `docs/limitations.md` at P0-T8.
- It also removes the usual normalisation cross-check — with no nuclei count
  there is no way to ask whether library size tracks cellularity, so the P0-T5
  QC thresholds rest on detection rate and library size alone.

**One open thread:** the `.dcc` headers themselves carry per-AOI sequencing QC
fields (raw / trimmed / stitched / aligned / deduplicated reads, sequencing
saturation). Those are not cellularity, but they may serve as QC covariates.
Check once `RAW.tar` lands at P0-T2.

## Q5 — Is patient-level clinical/survival metadata extractable and joinable?

**Status: Resolved — yes. Phase 5 is viable.**

> The clinical–histological characteristics of brain metastasis patients are
> described in Supplementary Data 1.

Supplementary Data 1 is `41467_2022_33365_MOESM4_ESM.xlsx`. Note the retrieval
quirk: `pmc.ncbi.nlm.nih.gov/articles/instance/9551067/bin/…` returns an HTML
interstitial rather than the file; `static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-33365-y/MediaObjects/…`
serves the xlsx itself.

Columns:

`Patient ID`, `Age at NSCLC diagnosis`, `Gender`, `Histological type of NSCLC`,
`Grade`, `Largest dimension (cm)`, `Visc Pleural Invasion`, `Margins`,
`Treatment of NSCLC`, `Metastasis Intervals to the brain (Months)`,
`Location of BrM`, `Histologic type of metastatic NSCLC/subtype`,
`Treatment of BrMs`, **`Primary lung cancer diagnosis to death (Months)`**,
**`Brain metastasis diagnosis to death (Months)`**.

Two time-to-event variables, with `Alive` appearing as a value — i.e. censoring
is representable. That the cohort has usable outcome data is corroborated by the
paper's own use of it:

> Kaplan–Meier survival analysis and Cox proportional hazards of the current
> cohort (n = 30) and the TCGA LAUD cohort (n = 501)

(The paper's *primary* survival analyses were on TCGA LGG and LUAD from
cBioPortal; the `n = 30` and `n = 23` panels of Fig. 8 are this cohort.)

### The join

GEO sample titles have the form `Patient <n> [<AOI code><n>]`, e.g.
`Patient 1 [LB01]`, `Patient 12 [TIME-L12a]`. Patient numbers run up to 43, so
**GEO preserved the paper's original 1–44 patient numbering** rather than
renumbering its 35-case subset. The join to Supplementary Data 1's `Patient ID`
should therefore be direct.

*Expected clean, not verified.* **Assert it at P0-T3** — every GEO patient number
must appear in Supplementary Data 1, and the AOI-code numeric suffix must equal
the patient number (§6 already requires that second assertion). Do not treat the
join as established until those pass.

### Three constraints on how the clinical data may be used

1. **Age is banded, not exact** (`40s`, `50s`, `60s`, `70s`, `90s`) —
   de-identified. It enters a model as an **ordered factor**, never as a
   continuous covariate.
2. **Missingness is real and inconsistently coded**: `N/A`, `NA`, `n/a`,
   `Unspecified` all appear, as does inconsistent capitalisation in the
   histology strings (`Adenocarcinoma/solid`, `adenocarcinoma/solid `, with
   trailing spaces). The P5 parser must normalise these, not trust them.
3. **Coverage is 35 of 44.** The table covers the full 44-patient cohort; GEO
   carries 35. Nine supplementary rows have no expression data.

**`config.yaml → phases.survival` stays `false`** until Phase 5 has rules. Q5 no
longer gates it; the remaining risk is n, not availability.

---

## Acquisition first-look (P0-T2, 2026-08-24)

The four artifacts are acquired, pinned and verified (`resources/provenance.tsv`).
What follows is **reconnaissance, not a result** — it comes from ad-hoc commands
over the downloaded files, so per hard constraint 5 nothing here is reportable
until P0-T3/P0-T4 produce it from a rule. It is recorded because it changes what
those tasks must handle.

- **Q3 confirmed in hand.** The 120 DCC entries in `GSE200563_filelist.txt` split
  91 / 29 across `DSP-1012300141221` and `DSP-1012310141221`, as the pre-download
  survey predicted. The batch variable is real and available.
- **Q1 consistent with Q3-normalised.** The matrix is 121 columns (`Gene#` + 120
  AOIs). Of 2,243,400 values, **0 are negative and 99.3% are non-integer** — so
  not raw counts. Magnitudes are in the tens, not the units, so not log₂. Still
  **Provisional**: the discriminating test is whether the column-wise third
  quartile is near-constant across AOIs, and that belongs in P0-T4's rule.
- **The matrix has 18,695 rows, not 18,694.** The extra row is **`NegProbe-WTX`**.
  Two consequences:
  1. **P0-T5 gets its negative-probe control.** §6 lists "negative-probe geomean
     *if available*" — it is available, in the processed matrix.
  2. **It must be separated from the gene matrix.** A naive row count will report
     18,695 genes and every gene-wise operation will carry a non-gene row.
     P0-T7's `var.tsv` must exclude it or flag it.
- **Columns are AOI codes, not GSM IDs.** The header reads `L01`, `L02`,
  `TIME-L12a`, … — i.e. the bracketed portion of the SOFT sample titles
  (`Patient 12 [TIME-L12a]`). **The join key between the expression matrix and
  `samples.tsv` is the AOI code**, and GSM ID reaches the matrix only via the
  SOFT titles. P0-T3 must carry both.

---

## Background is not constant across compartments (P0-T6)

`NegProbe-WTX` — the whole-transcriptome negative probe, present in the
processed matrix — is **not flat** across the design:

| Code | n | median log2 |
|---|---|---|
| L | 30 | 4.96 |
| LB | 27 | 4.99 |
| mLN | 13 | 5.03 |
| TIME-B | 8 | 5.30 |
| TIME-L | 15 | 5.36 |
| TBME | 20 | 5.62 |
| BC | 7 | 5.67 |

A monotone ~0.71 log2 (≈1.6×) rise from lung tumour cores to brain
parenchyma. Q3 normalisation scales each AOI by its own third quartile, so it
does not remove a background difference that tracks tissue type.

**Consequence for P0-T5:** a gene expressed near background will look
differentially expressed between lung and brain compartments on background
alone. Either subtract or model the negative-probe level, or restrict claims to
genes comfortably above it. This is a threshold decision, so it needs an ADR
and is a "stop and ask" — flagged here, not decided.

It does not touch the P0-T6 verdicts: every effect there is 2–6 log2, far above
a 0.7 log2 background gradient.

---

## Discrepancies and reconciliations

Both of these will be rediscovered as apparent bugs at P0-T3 if not recorded.

### 120 AOIs (GEO) vs 119 ROIs (paper)

GEO deposits **120** samples and 120 `.dcc.gz` files. The paper analyses **119**.

The difference is entirely in `TBME`. The paper's Fig. 1c legend:

> (L) = 30 samples, (LB) = 27 samples, TBME = 19 samples, TIME-L = 15 samples,
> TIME-B = 8 samples, mLN = 13 samples, BC = 7 samples.

That sums to 119, with `TBME = 19`. Tallying the 120 GEO sample titles gives:

| Code | GEO | Paper Fig. 1c | `PROJECT_PLAN` §2.1 |
|---|---|---|---|
| L | 30 | 30 | 30 |
| LB | 27 | 27 | 27 |
| mLN | 13 | 13 | 13 |
| TBME | **20** | **19** | **20** |
| TIME-L | 15 | 15 | 15 |
| TIME-B | 8 | 8 | 8 |
| BC | 7 | 7 | 7 |
| **Total** | **120** | **119** | **120** |

**GEO matches the design table exactly.** One TBME AOI is present in GEO and
excluded from the paper's analysis, with no stated reason. So:

- `PROJECT_PLAN` §2.1 stays authoritative and `expected_n_aoi: 120` is correct.
- `common.smk`'s hard failure on a row-count mismatch stays as written — it is
  right to fail if P0-T3 produces 119.
- Which TBME AOI the authors dropped is not recoverable from GEO. If it matters
  later (it would only matter to a direct replication of their TBME figures),
  it is a limitation, not a bug. Note it at P0-T8.

### 44 patients (paper) vs 35 (GEO)

Not a discrepancy — a stated subset. `!Series_overall_design`:

> In the study presented here, 35 cases from a well-defined cohort of 44 NSCLC
> cases with paired brain metastases was used to acquire expression profiles of
> a total of 18,694 unique genes.

This is why the Q5 join covers 35 of 44 supplementary rows.

### Design-table facts confirmed from the SOFT titles

Checked while answering Q3/Q5; all three load-bearing claims in `CLAUDE.md` hold.

- **Replicates exist and are suffixed `a`/`b`**: `TIME-L12a`/`TIME-L12b` and
  `TIME-L24a`/`TIME-L24b`. P0-T3's regex must handle a trailing letter after the
  patient number.
- **Exactly five patients — 5, 12, 15, 19, 35 — have both `TIME-L` and
  `TIME-B`.** The paired-immune consistency check has n = 5, as stated.
- **`TIME-L` / `TIME-B` are distinguishable two ways**: from the title code
  itself, and independently from `cell type` (`tumor immune microenviroment in
  the primary lung cancer` vs `… in the metastatic brain` — note the source's
  spelling of "microenviroment", which any string match must reproduce). P0-T3
  should parse the title and **cross-check** against `cell type` rather than
  relying on either alone.

---

## Known gaps

- **The WTA PKC file is not in GEO.** `GeomxTools::readNanoStringGeoMxSet()`
  requires a probe-kit configuration file (`Hs_R00000001_WTA.pkc`) to map probe
  IDs to targets. The series carries DCCs but no PKC, so it must come from
  NanoString — a non-GEO source needing its own provenance record. Raise at
  **P0-T5**, when the R path first needs it. Note the DCCs may be usable for QC
  metrics without a PKC; it is the expression matrix that needs the mapping.
- **TMA block (LB-D1..D4) is not recoverable per AOI.** See Q3.
- **Which TBME AOI the paper dropped is not recoverable.** See above.
