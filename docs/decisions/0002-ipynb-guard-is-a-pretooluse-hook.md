# ADR 0002 — The `.ipynb` guard is a PreToolUse hook, and both hooks match Bash

**Date:** 2026-08-23
**Status:** Accepted
**Task:** P0-T1
**Supersedes:** the hook registration sketched in PROJECT_PLAN §5.6

## Context

§5.6 specifies two `PostToolUse` hooks: `marimo-check` and `no-ipynb`. The
acceptance criterion for P0-T1 is that "a throwaway `.ipynb` write is **blocked**
by the hook."

Two problems with the sketch as written:

1. **`PostToolUse` cannot block.** It fires *after* the tool has run. A
   `PostToolUse` hook returning exit 2 tells the agent it did something wrong,
   but the `.ipynb` is already on disk. That is a complaint, not a guard, and
   it fails the acceptance criterion on the plain meaning of "blocked".

2. **Matching only `Edit|Write` leaves a hole.** Files in this project are
   frequently written through `Bash` heredocs rather than the `Write` tool. A
   policy enforced only on `Write|Edit` is silently unenforced along the path
   actually being used — which is precisely the failure mode §5.6 invokes when
   it says "a constraint that's only in CLAUDE.md is a suggestion".

## Decision

- `no-ipynb` is registered as **`PreToolUse`** on `Write|Edit|NotebookEdit|Bash`.
  Exit 2 refuses the call, so the file is never created.
- `marimo-check` stays **`PostToolUse`** (it must run against the saved file),
  but matches `Write|Edit|Bash`. When the tool is `Bash` and no `file_path` is
  present, it extracts `.py` paths from the command line and checks any that
  turn out to be marimo notebooks.
- The `Bash` branch of `no-ipynb` allows `marimo convert`, which is the
  sanctioned route for Jupyter material arriving from a paper or tutorial
  (§4.5), and otherwise blocks redirections and `touch`/`cp`/`mv`/`tee` that
  would put an `.ipynb` into the tree.

## Consequences

- The §4.5 policy is enforced rather than advertised, along both write paths.
- A `Bash` command that merely *mentions* `.ipynb` without creating one (say,
  `ls *.ipynb`) is allowed through; the guard targets creation, not mention.
  A determined bypass is still possible (e.g. a python one-liner writing the
  file); the hook is a guard rail, not a sandbox, and `CLAUDE.md` carries the
  intent.
- Deviating from the plan text here is deliberate and recorded rather than
  silent, per the project's own convention that real decisions get an ADR.
