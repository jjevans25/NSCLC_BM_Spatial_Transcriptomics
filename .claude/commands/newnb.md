---
description: Scaffold a marimo notebook in the right tier directory
argument-hint: <tier: explore|review|apps> <name>
allowed-tools: Bash(uvx marimo check:*), Bash(ls:*), Read, Write
---

Scaffold a marimo notebook: **$ARGUMENTS**

First re-read `notebooks/README.md` for the tier rules, then create
`notebooks/<tier>/<name>.py` with:

1. A header `mo.md` cell naming the notebook, its **tier**, its **phase**, and
   its purpose in one sentence.
2. The §A.5 parameter cell — `mo.cli_args()` with sensible defaults, so the
   notebook opens standalone under a bare `marimo edit` and a Snakemake rule can
   pass explicit paths.
3. A load cell, reactively downstream of the parameters.
4. An `mo.stop()` guard when the inputs may not exist yet, carrying a message
   that names the task which produces them.

Tier rules — enforce these, do not just mention them:

| Tier | In the DAG? | Reads | Writes |
|---|---|---|---|
| `explore` | no | anything | `results/interim/` only |
| `review` | no — a human runs it at a gate | `results/` | **nothing** |
| `apps` | **yes** — a rule exports it | `results/`, declared as rule inputs | its own exported HTML |

For an **apps**-tier notebook, also add the `export_notebook` rule wiring to the
right `.smk`, listing every file the notebook reads as a rule input. An
undeclared read is a silent provenance hole (§4.5).

Never put a reported number, table, or figure's *computation* in a notebook —
notebooks read pipeline outputs; `workflow/scripts/` produces them.

Finish by running `uvx marimo check --strict` on the new file. Do not report the
notebook as scaffolded until it passes.
