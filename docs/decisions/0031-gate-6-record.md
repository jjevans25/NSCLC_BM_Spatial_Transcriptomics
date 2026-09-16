# ADR 0031 — Gate 6 record: the clean room was green, and it found four defects

**Date:** 2026-09-16
**Status:** Accepted
**Task:** Phase 6 (P6-T1, T2, T2a, T2b, T3, T4, T4b, T6, T7)
**Related:** ADR 0023 (pins), ADR 0029 (a pin must carry every dependency),
ADR 0030 (the crate is manifest-driven)

## Context

**PROJECT_PLAN §6 defines no GATE line for Phase 6.** Every other phase has one.
The `/gate` command's own rule applies — "if an acceptance criterion cannot be
checked, the plan says it is not an acceptance criterion; flag it as a defect in
the plan rather than quietly passing it" — so the gate is defined here rather
than assumed, and the omission is recorded as the third defect this project has
found in §6's text (after ADR 0022 §5 and ADR 0030 §7).

## Decision — what Gate 6 is

Phase 6 passes when all five hold:

1. **A clean-room run is green from a fresh checkout**, on a path containing a
   space, with environments built from the pins.
2. **The RO-Crate validates against RO-Crate 1.1**, internally *and* against an
   external validator.
3. **The report opens standalone with no broken links.**
4. **Every notebook runs from one `uvx` command**, and `marimo check` is clean.
5. **`/fairscan` returns no failure**, with every partial named and justified.

## Verdict: PASSED

### 1. Clean-room reproduction — green, 72 of 72

Fresh `git clone` into `~/Documents/Biomedical Data Science/cleanroom-P6T1` —
**the path contains a space**, without which `p2t0_repair_r_env` is never
exercised and the run reports a false pass (ADR 0013) — with every derived output
deleted so nothing could be mistaken for up to date.

| | |
|---|---|
| jobs | **72, exit 0** |
| three environments, built from the pins | **62 s**, 3.4 GB (warm package cache; a cold machine downloads it) |
| pipeline, `--cores 4` | **19 min 24 s** |
| total, clone to finished | **20 min 26 s** |
| compute across all rules (benchmarks) | **21.4 min**, of which `p1t3_variance_partition` is 15.8 |

**Reproduction, measured rather than asserted:**

- **All 15 figures byte-identical** to the working tree, including the PCA
  ordination, the variance partition, the permutation null and both KM plots.
- **100 of 102 committed tables byte-identical.**
- **11 of 109 re-derived tracked files differ, and every one is a provenance
  field**: six `*_provenance.tsv` access timestamps (the artifacts were
  re-fetched and their digests re-verified), the git SHA in `h5ad_summary.json`,
  the tracked-prose digest in `crosstalk_language_audit.json`, the git-remote
  string and the built-from commit in the two audit summaries, and the crate's
  own byte count. **No estimate, interval, q-value or count moved.**

**One transient failure, recorded because it is not a defect and the difference
matters.** A confirmation run on the final commit lost `p5t4a_fetch_tcga` to a
network error; the same URL returned 200 by hand seconds later and the rule
succeeded on the next attempt. **The fetch rules declare no `retries:`, so a
flaky network fails the run** — recoverable by re-running, since the digests are
pinned and re-verified, but worth knowing before blaming the pipeline. It is
left as-is rather than fixed: adding `retries:` changes the code of rules whose
outputs are `protected()`, and that transition is not worth paying for a
condition that resolves by typing the command again.

### 2. RO-Crate — valid, internally and externally

**24/24** internal checks (`p6t2b_validate_ro_crate`) and **38/38 REQUIRED**
from `rocrate-validator` 0.11.4: *"RO-Crate is a valid ro-crate-1.1"*. 281
entries, 327 graph entities, 126 MB of described payload. ADR 0030 §4 records
the two REQUIRED failures the external validator found on the first attempt,
both real, one of which this project's own validator had encoded backwards.

### 3. The report — standalone

14 MB, **zero external `src`/`href` references**, all seven categories present
(Phase 0–5 plus Interactive), 17 figures each with a caption file that exists.
`p6t4_report_audit` passes **9/9**. Opened and read.

### 4. Notebooks — five of five

`marimo check` clean on all five, and **one notebook per tier was actually
launched** with `uvx marimo edit --sandbox` — review, app and explore — because
"a stranger can run any committed notebook with one command" is the stated
acceptance criterion and a lint pass is not that test. All three resolved their
PEP 723 dependencies and served.

`notebooks/explore/00_first_look.py` has PEP 723 metadata and an explicit SCRATCH header recording that Q1 was answered by
`p0t4_normalisation_check`, not by it. `landscape_explorer.py`'s black-background
WASM render is **fixed and verified by decoding the figure out of the built
`index.html`** — the only way that bug has ever been visible.

### 5. FAIR self-assessment — 12 pass, 3 partial, 0 fail

Recorded in `CHANGELOG.md` per PROJECT_PLAN §7.3. The three partials are named
and justified there; none is a gap that Phase 6 could close without either a DOI
or a second platform. **P6-T5 was subsequently abandoned by decision**
(ADR 0032), so two of the three are now permanently partial rather than pending
— which changes the wording, not the score.

## The part worth keeping: four defects, one shape

**Every one of these passed in the working directory and failed on a fresh
clone.** They are listed in full in `docs/NEXT_STEPS.md`; what matters here is
what they have in common.

1. The crate's `inputs` group globbed three gitignored directories, so the
   workflow **did not parse** on a clean checkout.
2. `p6t3_citation_audit` compared the declared repository URL against `git remote
   get-url origin` — **a property of the checkout, not of the project** — so any
   fork, mirror, local clone or tarball failed it.
3. The README's own setup instruction produced a **dangling symlink** on a fresh
   clone, and Snakemake died with a `FileExistsError` naming the symlink rather
   than the missing target.
4. `rocrate` was declared under a `pip:` section and **absent from the pin**
   (ADR 0029), so it existed here and would not have existed there.

**Three of the four were written during Phase 6 itself**, by an author who had
just finished writing the argument for why this class of error matters. That is
the finding. A check can only be trusted where it has run, and until P6-T1 every
check in this repository had only ever run in one directory, on one machine, with
one history. **"It passes here" is a statement about here.**

A fifth, smaller one belongs with them: `p4t7_language_audit` scans every tracked
file, and `ro-crate-metadata.json` is both tracked and generated by a rule
*downstream* of the audit. So the audit read it in one run and not in the other,
and reported different scan counts — **a check whose scanned set depends on
scheduling gives different answers on different runs.** The crate is now excluded
on the principle the audit already stated for `results/`: its prose comes from
`config/ro_crate.yaml`, which is scanned.

**And a sixth, which is the one that would have hurt most.** The language
audit's rerun trigger is a digest of every tracked text file, and
`ro-crate-metadata.json` is tracked, generated, and produced by a rule
*downstream* of the audit. So the two chased each other forever: rebuild the
crate, the digest changes, the audit re-runs, its output is a crate input,
rebuild the crate. **A completed run left four jobs pending, every time** —
`snakemake -n` could never come back clean, which is the standard every gate in
this project has been checked against.

**It took three attempts, and that is the point.** Excluding the crate alone was
not enough — `metadata/ro_crate_summary.json` and
`metadata/ro_crate_validation.json` are rewritten by the same rules and are
tracked too. Excluding `metadata/` as well settled the working directory, and
**the clean room still needed a second `snakemake` invocation**: the fetch rules
rewrite `resources/*_provenance.tsv` with a fresh access timestamp *during* the
run, after the digest was computed at parse time. Invisible here, where nothing
is ever re-fetched; unavoidable there, where everything is fetched for the first
time.

The digest now excludes **every committed rule output** — `results/`,
`metadata/`, the provenance TSVs and the crate — which is the category the
`results/` exclusion was always one instance of. Nothing is lost: the provenance
files are URLs and digests, and the rest are generated copies of prose that lives
in `config/ro_crate.yaml` and in the scripts, all of which are covered.

**The exclusion has to name the category, not the file that revealed it.** Three
fixes, each correct about the case in front of it and wrong about the next one,
is most of the lesson.

**What remains, stated exactly, because it is a residual rather than a fix.** In
the working directory the DAG now settles in one pass: run, then
`snakemake -n` reports "Nothing to be done", stably. **On a first-ever run in a
clean room it takes two invocations.** The second one re-runs four jobs — the
audit and the crate pair — and then it is stable, verified by running again and
getting "Nothing to be done" twice.

The mechanism is structural rather than mysterious: the digest is computed **at
parse time**, before any rule has run, and a first-ever run is the only run that
creates most of the repository's tracked outputs. Every file identified as
changing across that boundary is now excluded, and the residual pass was not
traced to a specific remaining file — measuring it afterwards shows **zero**
changed digest inputs, which is why the second invocation is the last one.

**It is recorded rather than claimed away.** A bounded, deterministic, one-pass
residual on a first-ever run is a fair thing for a stranger to be told; "the dry
run is clean after one run" would have been a claim about this working directory,
which is the exact error the four defects above are all instances of.

## Consequences

- **Phase 6 is complete.** P6-T5 (the Zenodo DOI) was the one remaining item and
  it was **abandoned by decision** immediately after this gate — ADR 0032.
  `fair.identifiers.doi` stays `null` and `p6t3_citation_audit` still enforces
  the policy in both directions, so the mechanism remains wired and a
  half-finished backfill still cannot ship.
- **The project is reproducible on `osx-arm64` and merely runnable on
  `linux-64`**, and `docs/limitations.md` §11 says so. The pins describe
  environments that exist, and these three exist for one platform.
- **A clean-room run is now the acceptance test of record.** It is cheap — 20
  minutes — and it is the only test in this project that runs somewhere its
  author is not.
