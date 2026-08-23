---
description: Run marimo check on a notebook and fix only what it reports
argument-hint: <path to notebook.py>
allowed-tools: Bash(uvx marimo check:*), Bash(marimo check:*), Read, Edit
---

Run `uvx marimo check --fix $ARGUMENTS`, then read the remaining diagnostics and
fix **only** what the checker reports.

Constraints:

- **Fix what is reported, nothing else.** Do not refactor, rename, reorder
  cells, or "improve" code the checker did not flag. Unrequested edits to a
  notebook are how a reviewed notebook silently becomes an unreviewed one.
- `multiple-definitions` is the common one. marimo forbids redefining a
  variable across cells. Prefer renaming to reflect what the value actually is;
  fall back to an underscore prefix (`_tmp`) only for genuinely cell-local
  throwaways.
- Re-run the check after fixing and report the final state.
- If a diagnostic implies a real design change (splitting a cell that does two
  unrelated things, say), **say so and stop** rather than guessing.

Finish with `uvx marimo check --strict $ARGUMENTS` — warnings are fatal there,
and that is the bar the pre-commit hook holds notebooks to.
