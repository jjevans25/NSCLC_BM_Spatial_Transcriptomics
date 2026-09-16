# Next steps — there are none. The project is complete.

**Last session:** 2026-09-16
**Branch:** `P6` was the last phase branch, squash-merged into `main` with the
Gate 6 verdict in the message (ADR 0004). There is no next phase branch.
*(No SHA is quoted here on purpose: a commit hash in this file is stale the
moment the next commit lands, and a stale one reads as authoritative. Run
`git log --oneline -1`.)*
**Gates 0–6: all PASSED.** Gate 6 is ADR 0031.
**Next free ADR number: 0033.**

---

## START HERE — nothing is outstanding

**All six phases are done and every gate has passed.** The last open item,
P6-T5's Zenodo DOI, was **abandoned deliberately** (ADR 0032) rather than left
pending. There is no backlog, no half-finished task and no TODO waiting in a
script.

**If you come back to this repository, the useful entry points are:**

| Want | Read |
|---|---|
| what it found, with the constraint on each finding | the results table at the bottom of this file |
| what it cannot support | `docs/limitations.md` — §9 first, then §1 |
| why any threshold is what it is | `docs/decisions/` — 32 ADRs, each committed before the result it governs |
| how to run it | `README.md` § How to run |
| what a reader gets without running it | `ro-crate-metadata.json`, `results/tables/`, `results/reports/workflow-report.html` |
| the write-up | `docs/substack-draft.md` (draft, unpublished) |

### If you ever do want the DOI

The mechanism is still wired and enforced; only the artifact was declined. About
an hour, and the **ordering matters**:

1. Turn on the **GitHub–Zenodo integration** for the repository **first**. If the
   release goes out before the toggle, Zenodo never sees it and you cut a second
   one.
2. Tag `v1.0.0` and publish a GitHub **release** — a pushed tag alone is not
   enough; Zenodo hooks the release event.
3. Set `fair.identifiers.doi` in `config/config.yaml`, add `doi:` to
   `CITATION.cff` and `identifier` to `codemeta.json`, mention it in
   `README.md` § Citing, then re-run.

The crate picks it up on its own — `build_ro_crate.py` writes `sameAs` only when
the field is populated — and `p6t3_citation_audit` **enforces the policy in both
directions**, so a half-finished backfill cannot ship. Note that editing
`config/config.yaml` costs the Phases 1–5 rebuild, as every phase switch has.

**Do not invent a DOI to make something look finished.** A persistent identifier
that does not resolve is worse than none: it looks like provenance, survives
copy-paste into someone else's bibliography, and fails silently years later.

---

## What Phase 6 produced

| Artifact | Checked by | Result |
|---|---|---|
| `ro-crate-metadata.json` — RO-Crate 1.1, **281 entries, 327 entities, 126 MB** | `p6t2b_validate_ro_crate` **and** an out-of-band `rocrate-validator` run | **24/24 internal, 38/38 external REQUIRED** |
| `metadata/ontology-terms.tsv` — **15 terms** | `p6t2a_resolve_ontology_terms`, live against EBI OLS4 every build | all 15 resolve at the pinned CURIE |
| `CITATION.cff` + `codemeta.json` | `p6t3_citation_audit` | **28/28** |
| `results/reports/workflow-report.html` — **14 MB, zero external references** | `p6t4_report_audit` | **9/9**, 17 figures |
| Clean-room reproduction | P6-T1, see below | **green, and it found four defects** |

---

## The four defects P6-T1 found, because this is the point of the exercise

A clean-room run is not a formality and this one earned its three hours. Every
one of these passed in the working directory and failed on a fresh clone:

1. **The workflow did not parse.** `config/ro_crate.yaml` globbed
   `resources/raw/`, `reference/` and `supplementary/` — all three gitignored —
   so the `inputs` group matched 19 files against an `expect_min` of 25. The
   threshold did its job; the group was wrong in *kind*. Those artifacts are
   declared in six manifests, so they are read from the manifests now, exactly as
   the results group is read from the phase TARGETS. **A glob cannot see a file
   that does not exist yet, and a clean room is the state where nothing exists
   yet.**
2. **The repository-URL check was a property of the checkout.**
   `p6t3_citation_audit` compared the declared URL to `git remote get-url
   origin`, and a clone taken from a local path fails that — as would a fork, a
   mirror, or a tarball with no git at all. **A rule that only passes in one
   working directory makes the workflow unrunnable for exactly the stranger
   P6-T1 exists to serve.** The three declarations are now checked against each
   other unconditionally and against git only where origin is a web URL.
3. **The README's own setup instruction fails on a fresh clone.** `ln -sfn
   "$PWD/.snakemake/conda" …` makes a *dangling* symlink when `.snakemake/` does
   not exist yet, and Snakemake then dies with a bare `FileExistsError` naming
   the symlink rather than the missing target. Fixed with an `mkdir -p` and an
   explanation of why the symlink exists at all.
4. **Before any of that: the pin did not carry `rocrate`** (ADR 0029). It was
   declared under a `pip:` section, and `conda list --explicit` — which is what
   generates the pins — records conda packages only. It worked here because the
   environment already existed. **A dependency no rule exercises is not pinned;
   it is declared.**

**The shape they share is worth more than any of them individually.** All four
are checks or instructions that were correct *about this machine* and silently
wrong about anywhere else — and three of the four were written in Phase 6 itself,
by someone who had just finished writing the argument for why that class of error
matters.

---

## State you'll have forgotten

- **`--conda-prefix "$HOME/nsclc-envs"` is MANDATORY on this machine**, alongside
  `--use-conda`. Recreate with
  `mkdir -p .snakemake/conda && ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"`.
  **Do not** relocate the envs: that fires the software-env trigger on
  `p0t2_fetch_geo` and its `protected()` outputs abort the DAG. ADR 0013.
- **A bare `snakemake -n` is not a clean dry run.** Always run from the
  `nsclc_bm_spatial` env (Snakemake 8.30); base anaconda's 9.20 writes
  `.snakemake/metadata` in a format 8.30 reads as stale.
- **No env YAML may contain a `pip:` section** (ADR 0029), and a pin must be
  regenerated after ANY env YAML edit. `p6t3_citation_audit` now enforces both
  on every build, so this is checked rather than remembered.
- **A module-level `_pattern` in any `.smk` kills every shell-using rule in the
  workflow.** An include's globals are the namespace Snakemake formats `shell:`
  strings against, and `snakemake.utils.format` has its own argument by that
  name. The error — "format() got multiple values for argument '_pattern'" —
  names a rule nowhere near the file that caused it.
- **The fetch rules declare no `retries:`.** A transient network error fails the
  whole run; re-running resumes it, and the digests are pinned and re-verified so
  nothing is at risk. Seen once during P6-T1 on `p5t4a_fetch_tcga`. Left as-is
  deliberately: adding `retries:` changes the code of rules whose outputs are
  `protected()`, which is a transition not worth paying for a condition that
  resolves by typing the command again.
- **A first-ever run needs a second `snakemake` invocation** before the dry run
  reports "Nothing to be done". `p4t7_language_audit`'s rerun trigger is a digest
  computed at PARSE time, and a first run is the only one that creates most of
  the repository's tracked outputs. The second pass re-runs four jobs and then it
  is stable. Not a defect to chase further — measured afterwards, zero digest
  inputs change — but know it before assuming something is wrong.
- **The crate is rebuilt whenever anything it describes changes**, which
  includes every doc. Edit prose, re-run, *then* commit — or the committed crate
  describes the previous version of the repository.
- **The crate carries `CITATION.cff`'s version, not a commit SHA**, and
  deliberately: a file cannot contain the hash of the commit that contains it,
  and a rerun trigger on the SHA would put the repository on a treadmill.
- **`Markdowns/PROJECT_PLAN.md` is gitignored** — ADRs are the durable record.
  ADR 0022 §5 and ADR 0030 §7 each record a defect in §6's text.
- **`results/reports/landscape_explorer`'s black background is FIXED** (P6-T4b),
  verified by decoding the figure out of the built `index.html`.

## Useful commands

```bash
conda activate nsclc_bm_spatial
mkdir -p .snakemake/conda && ln -sfn "$PWD/.snakemake/conda" "$HOME/nsclc-envs"
snakemake -n --use-conda --conda-prefix "$HOME/nsclc-envs"   # must stay clean
snakemake --lint                                             # must stay clean
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" --cores 4
snakemake --use-conda --conda-prefix "$HOME/nsclc-envs" \
  --report results/reports/workflow-report.html
uvx --from roc-validator rocrate-validator -y validate \
  --profile-identifier ro-crate-1.1 .          # out-of-band, ADR 0030 §4
uvx marimo check notebooks/**/*.py
/gate 6 && /fairscan
```

---

## The results, for reference — each with the constraint that bounds it

**Nothing below is a headline claim without its qualifier.** The qualifiers are
the result.

| Phase | Finding | Constraint |
|---|---|---|
| 1 | **41% of variance is compartment-attributable** (PVCA over 14 PCs spanning 60.4%; per-gene median 21%; permutation null 0.3%). Clustering ARI: compartment **0.741**, patient **0.003** | the random intercept for patient is justified by **non-independence**, never by variance share (ADR 0009) |
| 2 | Brain − lung within `TIME`: cytotoxicity **−0.946 SD** [−1.590, −0.302], q = 0.018; antigen presentation **−0.477 SD** [−1.030, +0.075], q = 0.098 | **`TIME-B` n = 8**; the 1.1–1.3 SD power floor makes the second uninformative, not negative |
| 3 | `CD274` (PD-L1) enriched in the lung immune compartment, **+0.612 SD** [+0.361, +0.863], q = 0.0001 | **40 of 63 gene × compartment cells are below the detection floor**; the panel clears it in `TIME-L` alone |
| 4 | Detection by interaction class: ECM-Receptor **45.3%** brain / 39.6% lung, Cell-Cell Contact 13.7% / 17.3%, Secreted Signaling **6.4%** / 5.3% | **0 of 268 lung and 0 of 260 brain** pairs exceed chance — an **assay-sensitivity limit, never evidence the crosstalk is absent** |
| 5 | **0 of 14** tests survive BH at FDR 0.05; smallest q **0.256** (TCGA antigen presentation, HR 0.70 [0.52, 0.94], n = 502, 182 events) | **the two nulls are not the same null** — 12 and 7 patients here against 502 there. That asymmetry is the deliverable |

**The sentence Phase 4 earned, both halves of it:**

> No inferred ligand–receptor pair exceeded chance expectation in either
> adjacency — lung n = 13, brain **`TIME-B` n = 8** — and the assay resolves
> matrix-associated interactions roughly **eight times** as often as secreted
> ones.

Anyone writing this up who states the first half without the second has
converted a measurement of the assay into a claim about biology.
