# Notebooks

**Every notebook here is a [marimo](https://marimo.io) notebook (`.py`). There are no
`.ipynb` files in this project — not created, not committed, not accepted as
deliverables** (PROJECT_PLAN §4.5). A `PreToolUse` hook enforces this
(`.claude/hooks/no-ipynb.sh`, ADR 0002); it is not an honour system.

Jupyter material arriving from a paper or tutorial gets converted before
anything is built on it:

```bash
marimo convert old.ipynb -o notebooks/explore/old.py
```

## The hard rule

> **No number, table, or figure that appears in the final report may originate
> in a notebook.**

Notebooks *read* pipeline outputs and make them explorable; `workflow/scripts/`
*produces* them. If an exploration turns out to matter, it gets promoted into a
script and a Snakemake rule, and the notebook then reads that rule's output.
This is the single constraint that keeps Snakemake the source of truth — without
it, notebooks quietly become a second, undocumented pipeline.

## The three tiers

| Tier | Directory | In the DAG? | Reads from | May write to | Committed? |
|---|---|---|---|---|---|
| **Explore** | `explore/` | No | anything | `results/interim/` only | Yes — but never cited as the source of a result |
| **Review** | `review/` | No — a human runs it at a gate | `results/` | nothing | Yes |
| **App** | `apps/` | **Yes** — a rule exports it headlessly | `results/`, declared as rule inputs | its own exported HTML | Yes |

## Conventions

- **Parameters, not hardcoded paths.** Every app-tier notebook opens with the
  §A.5 block: `mo.cli_args()` with sensible defaults, so it opens standalone
  with a bare `marimo edit` and the Snakemake rule can pass explicit paths.
- **Declared inputs.** Anything an app-tier notebook reads must also appear in
  its rule's `input:`. An undeclared read is a silent provenance hole.
- **Seeds come from `config/config.yaml`**, same as the scripts. No
  notebook-local RNG.
- **One notebook, one purpose.** marimo forbids redefining a variable across
  cells; treat that as a design aid, not an obstacle.
- **Prose carries the reasoning.** A review notebook with no `mo.md` cells is a
  script wearing a costume.
- **`uvx marimo check` must pass before commit** — enforced by the
  `marimo-check` `PostToolUse` hook (§5.6).

## Running them

```bash
# Edit, with the browser reloading on every agent-side write:
marimo edit --watch notebooks/review/qc_review.py

# Read-only:
marimo run notebooks/apps/checkpoint_explorer.py

# Lint (what the hook runs):
uvx marimo check notebooks/**/*.py

# A marimo notebook is also just a Python script:
python notebooks/explore/00_first_look.py
```

App-tier notebooks are exported to self-contained interactive HTML by a rule, so
they land in `snakemake --report` alongside the static figures:

```bash
# What the rules run -- `p1t5_landscape_explorer` (03_landscape.smk) and
# `p3t6_checkpoint_explorer` (05_checkpoints.smk). The output is a DIRECTORY
# containing index.html, not a single file, and the notebook's declared inputs
# are staged into a sibling `public/` so the export can fetch them over HTTP.
# That data path is ADR 0010 -- read it before adding a third app.
marimo export html-wasm --execute --mode run -f <staged>/<nb>.py -o results/reports/<nb>
```

The export **must be served over HTTP**, never opened as a `file://` URL:

```bash
python -m http.server --directory results/reports/<nb>
```

## Current inventory

| Notebook | Tier | Phase | Status |
|---|---|---|---|
| `explore/00_first_look.py` | Explore | P0 | Written — Q1 reconnaissance. **Q1 is resolved** (Q3-normalised; `normalisation_check.tsv`), so this is scratch. No PEP 723 block yet — a P6-T4b item |
| `review/design_review.py` | Review | P0 | Written (P0-T3) |
| `review/qc_review.py` | Review | P0 | Written (P0-T5) |
| `apps/landscape_explorer.py` | App | P1 | Written (P1-T5) — exported by `p1t5_landscape_explorer` to `results/reports/landscape_explorer`. Known: renders on a black background; left alone so the rule's rerun trigger does not fire |
| `apps/checkpoint_explorer.py` | App | P3 | Written (P3-T6, ADR 0019) — exported by `p3t6_checkpoint_explorer` to `results/reports/checkpoint_explorer`. **Browser-verified 2026-09-11** |

Statuses are checked against the tree, not remembered: `find notebooks -name "*.py"`
must return exactly the rows above. This table sat three phases out of date once.
