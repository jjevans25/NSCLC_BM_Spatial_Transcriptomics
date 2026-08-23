---
description: Dry-run the workflow and summarise what would change
allowed-tools: Bash(snakemake -n:*), Bash(snakemake --lint:*), Bash(snakemake --summary:*), Read
---

Run the workflow dry and report what *would* happen — never execute it.

```
snakemake -n --rerun-triggers mtime
snakemake --lint
```

Then summarise for the user:

1. **Job count by rule.** If zero jobs, say so plainly and say why (empty DAG
   because a phase switch is off, or everything is up to date — these are
   different situations and the user needs to know which).
2. **What would be created vs. re-created.** Re-creation of an existing output
   is the interesting case: name the trigger (changed input, changed params,
   changed code, mtime).
3. **Anything alarming.** A rule that would re-run the entire pipeline because
   a config value moved, a missing input that is not yet produced by any rule,
   or a lint warning.

Report `--lint` output verbatim if it is non-empty. A workflow that does not
dry-run and lint cleanly does not reach `main` (§5.5).

$ARGUMENTS
