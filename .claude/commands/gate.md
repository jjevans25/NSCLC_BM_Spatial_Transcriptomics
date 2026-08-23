---
description: Print a phase's acceptance criteria and check which outputs exist
argument-hint: <phase number, e.g. 0>
allowed-tools: Bash(ls:*), Bash(test:*), Bash(wc:*), Bash(snakemake -n:*), Read, Grep
---

Check gate readiness for **Phase $ARGUMENTS**.

1. Read the phase's task table and its **GATE** line from
   `Markdowns/PROJECT_PLAN.md` §6.
2. For every task in that phase, list its **Accept** condition and check it
   against reality:
   - Does each declared output path exist?
   - For table outputs, does the row count match what the plan states?
   - For decisions, is there an ADR in `docs/decisions/`?
3. Produce a table: task ID | acceptance condition | **met / not met / cannot
   check** | evidence.

Then state the gate verdict in one sentence.

**Be strict about this:** "cannot check" is a real answer and must not be
rounded up to "met". If an acceptance criterion cannot be checked, the plan
says it is not an acceptance criterion (§0) — flag it as a defect in the plan
rather than quietly passing it.

Do not run the pipeline. Read state; do not create it.
