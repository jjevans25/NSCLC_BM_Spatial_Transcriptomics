# Next steps — Phase 5 (prognostic association) or Phase 6 (FAIR packaging)

**Last session:** 2026-09-14
**Branch:** `P4`, **not yet merged**. Two commits on top of `main` at `3ccabe0`.
**Gate 0:** PASSED. **Gate 1:** PASSED. **Gate 2:** PASSED (ADR 0014).
**Gate 3:** PASSED (ADR 0020). **Gate 4:** **PASSED (ADR 0022).**
**Next free ADR number: 0023.**

**Phase 4 is complete and everything it produced is committed.** `snakemake
--lint` is clean and a re-run reports "nothing to be done". Squash-merge `P4`
into `main` with the Gate 4 verdict in the message (ADR 0004) before starting
anything new.

---

## READ THIS FIRST — the one sentence Phase 4 earned

> No inferred ligand–receptor pair exceeded chance expectation in either
> adjacency — lung n = 13, brain **`TIME-B` n = 8** — and the assay resolves
> matrix-associated interactions roughly **eight times** as often as secreted
> ones.

Both halves matter and the second is the larger one. **The null is an
assay-sensitivity limit, never evidence that the crosstalk is absent.** Anyone
writing this up (P6-T7) who states the first half without the second has
converted a measurement of the assay into a claim about biology.

---

## Which phase is next — decide this before writing code

PROJECT_PLAN's slip rule: *"if you're behind at the end of week 4, cut Phase 5
entirely and protect Phase 6. A finished, reproducible, well-documented
four-phase project beats a sprawling six-phase one that nobody can run."*

**Phase 5 is a stretch goal with a one-week timebox and Gate 5 is "timebox
respected".** Q5 resolved **yes** — Supplementary Data 1 carries two
time-to-event columns — so it is viable. Three constraints from
`docs/limitations.md` §8 bind it before P5-T1 starts:

- **Age is banded** (`40s`, `60s`), so it enters a model as an ordered factor,
  never as a continuous covariate.
- **Missingness is inconsistently coded** — `N/A`, `NA`, `n/a`, `Unspecified`.
- **The table covers all 44 cohort patients while GEO carries 35**, so nine rows
  have no expression data. P5-T1's hard gate is that it joins cleanly on patient
  ID; if it does not, stop the phase rather than spend days on fuzzy matching.

**Phase 6 is arguably the primary deliverable** (PROJECT_PLAN §6's own words),
and **P6-T1 has a known blocker that will not fix itself**: see below.

**Recommendation — do P6-T1 first even if you intend to do Phase 5.** It is the
only task whose failure invalidates everything else, and it is currently blocked
on a stale lock file. Finding that out in week 6 is worse than finding it out
now.

---

## The P6-T1 blocker is RESOLVED (ADR 0023), and the note about it was wrong

The long-standing note read "`py-analysis.conda-lock.yml` is stale and
`conda-lock` is not installed". **Both halves were false and the truth was
worse:** `conda-lock` was installed all along, and Snakemake never read
`*.conda-lock.yml` at all — it reads `<env>.<platform>.pin.txt`
(`deployment/conda.py:128`). **Nothing was pinned.** Two of the three locks were
also stale in content, missing `xarray` and `r-reformulas`, both of which were
added to fix real failures — so wiring them up would have reproduced two known
bugs in the very run meant to prove the pipeline works.

**Now fixed.** `workflow/envs/<env>.osx-arm64.pin.txt` exists for all three
envs, generated with `conda list --explicit --md5` from the environments that
actually produced every committed result. Verified by forcing one rule per env
to re-execute: all three outputs byte-identical, including the R mixed-model
fit. The three dead `*.conda-lock.yml` files are removed.

**What remains for P6-T1:**

- **The pins are `osx-arm64` only.** On `linux-64` Snakemake finds no pin and
  solves the YAML, exactly as it did before — nothing is lost, but nothing is
  guaranteed either. **P6-T6's limitations document must say so.** If
  cross-platform reproduction is wanted, render a `linux-64.pin.txt` from a
  `conda-lock` solve; `conda-lock` is still in `environment.yml` for that.
- **Regenerating a pin is now mandatory after any env YAML edit**, because a
  stale pin installs the wrong environment silently instead of failing.
- **P6-T1 must be run on a path containing a space**, or it will not exercise
  `p2t0_repair_r_env` and will report a false pass. This working directory has
  one ("Biomedical Data Science"); a clean clone elsewhere may not.
- **`r-geomx.yaml`'s header is stale**: it says "Phase 0" but
  `p4t2c_export_lr_database` now runs in it too. Editing it fires the
  software-env trigger on every rule in that env AND invalidates the pin, so it
  needs the ADR 0023 transition (rebuild, re-repair, `--touch`) — do it in a
  commit where that is the point, not casually.

---

## What Phase 4 left for Phase 6 to use

| need | where |
|---|---|
| the Gate 4 record and every Phase 4 decision | `docs/decisions/0022-gate-4-record.md`, `0021-phase-4-pre-registration.md` |
| the detection-by-class table (the headline) | `results/tables/crosstalk_lr_summary.json` → `by_site` |
| the honest null | `results/tables/crosstalk_fdr.tsv`, `crosstalk_null_summary.json` |
| the empty nomination table, which is a result | `results/tables/crosstalk_nominations.tsv` |
| the concordance prose for P6-T7 | `docs/analysis-notes.md`, P4-T6 section |
| figures, both `report()`-wrapped with captions | `crosstalk_null.png`, `crosstalk_network.png` |
| a build-time language check | `p4t7_language_audit` — **runs on every build**, so P6-T7's prose is checked automatically |
| a permutation-null implementation | `workflow/scripts/lr_null_calibration.py` — identity permutation reproduces the observed rho to 1e-16 |

---

## Things Phase 4 learned that the next phase should not relearn

1. **Open the PNG.** Both Phase 4 figures clipped on first render — the network
   showed "CEACAM" for both CEACAM1 and CEACAM5, and the colourbar label went
   off-canvas. P3-T5 had the same failure and the plan warned about it. Use the
   `main | colourbar` gridspec from `checkpoint_dotplot.py`; **do not use
   `tight_layout` with an attached colourbar** — it warns and then clips.
2. **Assertions catch what review does not.** The FDR step-down was wrong under
   ties (Spearman on 8 patients takes 42 distinct values across 260 rows) and
   would have shipped a plausible ranked table. The monotonicity assertion
   caught it. Write the assertion that would fail if you were wrong.
3. **A hand-transcribed orientation has nothing checking it.** `CD96–NECTIN1`
   was entered backwards into the comparator and reported as absent from a
   database that contains it. Where the underlying measure is undirected, match
   both orientations.
4. **A grep that cannot tell a prohibition from a claim is unrunnable.**
   PROJECT_PLAN §6 P4-T7's literal instruction fails on `CLAUDE.md`'s own hard
   constraint 6, because the files that state a rule must name what they forbid.
5. **`config/config.yaml` has been edited for the last time in the analysis
   phases.** Every further threshold belongs in a sibling YAML. The Phase 4 edit
   cost a full Phases 1–3 rebuild and was verified at 58 of 59 tables
   byte-identical; there is no reason to pay it again.

---

## State you'll have forgotten

- **`--conda-prefix "$HOME/nsclc-envs"` is MANDATORY on this machine**, alongside
  `--use-conda`. Recreate with `ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"`.
  **Do not** relocate the envs: that fires the software-env trigger on
  `p0t2_fetch_geo` and its `protected()` outputs abort the DAG. ADR 0013.
- **A bare `snakemake -n` is not a clean dry run.** Always run from the
  `nsclc_bm_spatial` env (Snakemake 8.30); base anaconda's 9.20 writes
  `.snakemake/metadata` in a format 8.30 reads as stale.
- **`resources/` now has four sanctioned writers** — `raw/` (P0-T2),
  `reference/` (P2-T5), `supplementary/` (P3-T2b) and `ligand_receptor/`
  (P4-T2). Each has a `*_checksums.sha256` and `*_provenance.tsv` pair, and
  P6-T2's RO-Crate must describe all four.
- **`results/reports/landscape_explorer` still renders on a black background** —
  P1-T5's notebook, left alone because editing it fires that rule's rerun
  trigger. The fix is the `plt.rcParams` reset in
  `notebooks/apps/checkpoint_explorer.py`. Worth doing when something else
  touches Phase 1, and P6-T4b is that occasion.
- **`Markdowns/PROJECT_PLAN.md` is gitignored** — ADRs are the durable record.
  ADR 0022 §5 records two further defects in §6's Phase 4 text.

## Useful commands

```bash
conda activate nsclc_bm_spatial
ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"   # once per machine
git switch main && git merge --squash P4             # ADR 0004, Gate 4 verdict in the message
snakemake -n --use-conda --conda-prefix "$HOME/nsclc-envs"   # must stay clean
snakemake --lint                                             # must stay clean
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" --cores 4
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" \
  --report results/reports/workflow-report.html
/gate 5
```

---

## Phase 4, for reference

**Gate 4 (ADR 0022) passed on its second clause.** The empirical FDR was
computed and honoured; **0 of 268 lung and 0 of 260 brain** primary
direction-rows clear FDR 0.05, smallest 0.499 and 0.224. Observed |rho| sits
**on** the permuted-pairing null. The gate turns on P4-T4, not P4-T5:
`lr_nominations.py` was written with no fallback and no second threshold, so an
empty table was reachable rather than escapable.

**The larger result is the detection gap by class**, ADR 0014's claim made
numeric on 2,239 interactions: ECM-Receptor **45.3%** brain / **39.6%** lung
admitted, Cell-Cell Contact 13.7% / 17.3%, Secreted Signaling **6.4%** / **5.3%**.

**Background is a patient-level property** — across a patient's two paired AOIs,
`negprobe_log2` ρ +0.548 brain / +0.352 lung, `gene_detection_rate` ρ +0.667 /
+0.231. Three of four above threshold, which is why ADR 0021 §6's
abundance-matched Null B was pre-registered. P3-T4's *site* gradient was
negligible; that is a different gradient.

**`VEGFA–NRP1`, which this file previously called the one informative
comparator, is absent from CellChatDB v2 entirely** — it carries VEGFA to FLT1
and KDR, with NRP1 only as a SEMA3x co-receptor. A database difference, not a
detection one. Four of the five comparators are unevaluable for four different
reasons.

## Phase 3, for reference

`CD274` (PD-L1) is enriched in the **lung immune** compartment, **+0.612 SD
[+0.361, +0.863], q = 0.0001**, 45 AOIs across 30 patients. **40 of 63 gene ×
compartment cells are below the pre-registered floor** and the panel clears it in
`TIME-L` alone. The lung-vs-brain shift is null in every assessable gene —
largest `CD274` at −0.260 [−0.595, +0.075], q = 0.462, **`TIME-B` n = 8** — and
against the 1.1–1.3 SD power floor those nulls are **uninformative, not
negative**.

## Phase 2, for reference

Brain minus lung, within `TIME`, `TIME-B` n = 8: antigen presentation
**−0.477 SD [−1.030, +0.075]** (q = 0.098, 13/13 detected at both sites),
cytotoxicity **−0.946 SD [−1.590, −0.302]** (q = 0.018), myeloid M2 null. Three
of six signatures failed their coverage floor; the genes that failed are the
cytokine genes, which Phase 4 has now measured across a whole LR database.
