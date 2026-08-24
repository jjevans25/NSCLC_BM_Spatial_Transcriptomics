# ADR 0004 — One branch per phase, not one branch per task

**Date:** 2026-08-24
**Status:** Accepted
**Task:** P0-T1 (decided during), P0-T4 (recorded)
**Supersedes:** the branch-name clause of `Markdowns/PROJECT_PLAN.md` §0

## Context

The plan specifies git branch granularity twice, and the two specifications
contradict each other.

§0 ("How to read this plan"), in the table defining the shape of a task:

> | **ID** | `P<phase>-T<task>`, used as the Snakemake rule prefix and the git branch name |

That makes the branch name `P0-T1` — one branch per **task**.

§5.5 ("Working style"):

> - **One phase per branch**, squash-merge with the phase gate result in the commit message.

That makes the branch name `P0` — one branch per **phase**.

Both cannot hold, and the conflict is not cosmetic. The §5.5 form carries a
requirement the §0 form cannot satisfy: **a gate result cannot go in a
squash-merge message if the branch is one task long**, because gates are
phase-level. §6 defines a single GATE 0 for all eight of P0-T1..T8; there is no
such thing as a P0-T3 gate result to put in a P0-T3 merge message. Following §0
literally would mean squash-merging eight branches into `main` with nothing to
say in seven of them, and the mechanism §5.5 describes — "the gate is the
decision point that keeps a project this size from sprawling" — would have no
place to leave its record.

Per-task branching has one real advantage, which is why §0's version is
tempting: a failed gate discards one task's work rather than a phase's. That
advantage is what per-phase branching gives up, and it is priced below in
Consequences.

The decision was made in practice during P0-T1 but written down nowhere.
`CLAUDE.md` §Conventions asserts the outcome ("One phase per branch") without
the reasoning, so the next reader hits the same §0/§5.5 ambiguity with no
record that it was ever resolved — which is what this ADR exists to stop.

## Decision

**One branch per phase.** The branch is named for the phase (`P0`, `P1`, …) and
squash-merges into `main` with the phase gate result in the commit message.

`P<phase>-T<task>` survives in its two other roles, which §0 also assigns it and
which are not in dispute: the **Snakemake rule prefix** (`p0t2_fetch_geo`) and
the **commit message prefix** (`P0-T2: acquire GSE200563 with pinned
provenance`). Only the "git branch name" clause of §0's ID row is superseded.

**Exception — infrastructure tasks may branch and merge early.** A task that
produces no scientific result cannot be invalidated by a gate, so holding it
behind one buys nothing and costs review latency. This exception has already
been exercised: P0-T1 (scaffold, Snakemake skeleton, notebook infra) was
branched as `P0-T1` and merged via PR #1 before Gate 0 (`57b7531`). A task
qualifies only if a failing gate would leave its output untouched — scaffolding,
hooks, env specs, CI. Anything that produces a number, a table, a figure or a
threshold does not qualify, however mechanical it looks.

## Consequences

What this costs, first, because these are the ones that will actually be felt:

- **A phase branch is long-lived.** `P0` spans eight tasks and 12–16 hours of
  work. It diverges from `main` for that whole span, so conflicts accumulate and
  rebases are more painful than they would be against an eight-branch sequence.
- **A failing gate discards more.** If Gate 0's marker sanity check (P0-T6)
  fails and the compartment mapping has to be re-derived, the branch carrying
  that failure also carries T2–T5 and T7–T8. Per-task branching would have
  landed those already. This is the real price, and it is highest for the phases
  with the most tasks.
- **Review arrives late and large.** There is one review surface per phase, not
  per task, so a wrong turn at T3 can go unreviewed until T8.

What it buys:

- **The gate result has somewhere to live.** This is the whole point. One
  squash-merge per phase, its message stating whether the gate passed and on
  what evidence, gives `main`'s history one entry per gated decision — which is
  the granularity at which this project's decisions are actually made.
- **The phase is reviewable whole.** Phase 0's output is "a validated, annotated,
  normalised `.h5ad` you trust". That is a property of the phase, not of any
  task in it, and it can only be reviewed as a unit.
- **`main` stays gate-clean.** Nothing reaches `main` that has not passed a gate,
  except the infrastructure exception above — which by construction cannot carry
  a scientific claim.

The mitigation for the first two costs is not more branches, it is **merging
`main` into the phase branch regularly** rather than letting it diverge, and
using the per-task commit prefix so that `git log --oneline` on the phase branch
still reads as a task sequence and individual tasks stay revertable within it.
