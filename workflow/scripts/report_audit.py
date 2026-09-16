"""P6-T4 — check that `snakemake --report` will be worth opening.

Owner task: P6-T4. Driven by rule p6t4_report_audit.

PROJECT_PLAN §6 P6-T4: "`snakemake --report`. Every figure captioned,
categorised by phase, provenance-linked; app notebooks under an 'Interactive'
category. **Accept:** report opens standalone with no broken links."

**The report itself is a COMMAND, not a rule, and that is the one sanctioned
exception to hard constraint 5 (ADR 0030 §3).** `snakemake --report` consumes
the completed DAG — it renders what the workflow produced, after it produced it
— so it cannot be a node inside that DAG. Making it one would mean a rule whose
input is the result of running every rule including itself.

**What CAN be a rule is the condition that makes the report worth generating**,
and that is this. A figure without a caption still appears in the report; it
just appears uselessly, and nobody notices until a reader asks what they are
looking at. So:

  * every figure output in every `.smk` is wrapped in `report(...)`;
  * each carries a caption file that EXISTS and is not empty;
  * each carries a category, and the categories are phase-shaped;
  * no caption file is orphaned — a `.rst` nobody references is a caption that
    was written for a figure that has since been renamed, and it will not
    appear in the report at all;
  * every reported figure is a declared TARGET, because a figure no target
    asks for is a figure the report has a hole where.

**Keyed on a digest of what it scans.** `p4t7_language_audit` shipped with no
rerun trigger and would have run once and then reported a green result forever
(NEXT_STEPS lesson 6). A check that cannot go stale is the only kind worth
having, so this rule declares every file it reads as an input.
"""

import json
import re
from pathlib import Path

# A category is phase-shaped or "Interactive". A pattern rather than a list, so
# adding Phase 7 does not require editing this file — and so a typo
# ("Phase 3 - checkpoints", with a hyphen) still fails.
CATEGORY = re.compile(r"^(Phase \d+ — .+|Interactive)$")

# Outputs the report renders as figures. The app notebooks export to a
# directory of WASM assets whose entry point is index.html, and those are the
# "Interactive" category §6 asks for.
FIGURE_SUFFIX = (".png", ".pdf", ".svg", ".html")

paths = snakemake.params.paths
targets = set(snakemake.params.targets)

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

failures = []
checks = []


def check(name, condition, detail=""):
    checks.append({"check": name, "passed": bool(condition), "detail": detail})
    if not condition:
        failures.append(f"{name}: {detail}")
    return bool(condition)


def resolve(expression):
    """Turn an f-string output expression into the path it will produce.

    The .smk files write paths as f"{PATHS['figures']}/x.png". Substituting the
    config's own values is what makes this audit agree with the DAG rather than
    with a second guess at where things land.
    """
    resolved = expression
    for key, value in paths.items():
        resolved = resolved.replace(f"{{PATHS['{key}']}}", value)
        resolved = resolved.replace(f'{{PATHS["{key}"]}}', value)
    return resolved


def report_blocks(text):
    """Every report(...) CALL in a .smk, as (body, start_line).

    Anchored on `<name>=report(`, the form every output entry takes, so the
    audit does not match its own rule's docstring where it explains what it
    looks for. `p4t7_language_audit` assembles its patterns from fragments to
    solve the same problem; here an anchor is enough, because a prose mention
    is never an assignment.
    """
    for match in re.finditer(r"\w+\s*=\s*report\(", text):
        depth = 1
        index = match.end()
        while depth and index < len(text):
            if text[index] == "(":
                depth += 1
            elif text[index] == ")":
                depth -= 1
            index += 1
        yield text[match.end():index - 1], text[: match.start()].count("\n") + 1


with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    emit("P6-T4 — report readiness audit")
    emit("=" * 70)
    emit()
    emit("  `snakemake --report` is a COMMAND, not a rule (ADR 0030 §3): it")
    emit("  consumes the completed DAG, so it cannot be a node inside it. This")
    emit("  is the part that CAN be checked automatically.")
    emit()

    reported = {}
    unwrapped = []

    for smk in sorted(snakemake.input.rules):
        text = Path(smk).read_text(encoding="utf-8")

        wrapped_spans = []
        for body, line in report_blocks(text):
            wrapped_spans.append(body)

            path_match = re.search(r'f?"([^"]+)"', body)
            caption = re.search(r'caption\s*=\s*"([^"]+)"', body)
            category = re.search(r'category\s*=\s*"([^"]+)"', body)
            labels = re.search(r"labels\s*=", body)

            if not path_match:
                failures.append(f"{smk}:{line}: report() with no output path")
                continue

            resolved = resolve(path_match.group(1))
            reported[resolved] = {
                "file": smk,
                "line": line,
                "caption": caption.group(1) if caption else None,
                "category": category.group(1) if category else None,
                "has_labels": bool(labels),
            }

        # A figure output NOT inside any report() block. Found by removing the
        # blocks and looking at what is left: a figure that renders into the
        # report uncaptioned is the failure this catches.
        remainder = text
        for body in wrapped_spans:
            remainder = remainder.replace(body, "")
        for match in re.finditer(r'f"(\{PATHS\[[^\]]+\]\}/[^"]+)"', remainder):
            candidate = resolve(match.group(1))
            if not candidate.endswith(FIGURE_SUFFIX):
                continue
            if candidate.startswith(paths["figures"]) and candidate not in reported:
                unwrapped.append(f"{smk}: {candidate}")

    # Anything listed as a target but also unwrapped is a genuine miss; the
    # same path appearing in a TARGETS list AND in a report() is normal.
    unwrapped = sorted({u for u in unwrapped if u.split(": ")[1] not in reported})

    emit(f"1. Coverage: {len(reported)} reported figures")
    check(
        "every_figure_is_reported",
        not unwrapped,
        f"{len(unwrapped)} figure output(s) not wrapped in report(): {unwrapped}",
    )
    emit()

    emit("2. Captions")
    missing_caption = []
    missing_file = []
    empty_caption = []
    referenced = set()
    for path, entry in sorted(reported.items()):
        if not entry["caption"]:
            missing_caption.append(path)
            continue
        # caption= is relative to the .smk's directory, as Snakemake resolves it.
        caption_path = (Path(entry["file"]).parent / entry["caption"]).resolve()
        referenced.add(caption_path)
        if not caption_path.exists():
            missing_file.append(f"{path} -> {entry['caption']}")
        elif not caption_path.read_text(encoding="utf-8").strip():
            empty_caption.append(f"{path} -> {entry['caption']}")

    check("every_figure_has_a_caption", not missing_caption, f"{missing_caption}")
    check("every_caption_file_exists", not missing_file, f"{missing_file}")
    check("no_caption_is_empty", not empty_caption, f"{empty_caption}")

    orphans = sorted(
        str(p.relative_to(Path.cwd()))
        for p in (Path.cwd() / "workflow" / "report").glob("*.rst")
        if p.resolve() not in referenced
    )
    check(
        "no_orphaned_caption_files",
        not orphans,
        f"{orphans} are referenced by no report(): a caption written for a "
        "figure that has since been renamed never appears in the report.",
    )
    emit(f"   {len(referenced)} caption files referenced, {len(orphans)} orphaned")
    emit()

    emit("3. Categories")
    bad_category = [
        f"{path}: {entry['category']!r}"
        for path, entry in sorted(reported.items())
        if not entry["category"] or not CATEGORY.match(entry["category"])
    ]
    check("every_category_is_phase_shaped", not bad_category, f"{bad_category}")

    categories = sorted({e["category"] for e in reported.values() if e["category"]})
    for category in categories:
        n = sum(1 for e in reported.values() if e["category"] == category)
        emit(f"   {n:2d}  {category}")
    check(
        "interactive_category_is_used",
        any(c == "Interactive" for c in categories),
        "no figure is categorised Interactive — §6 asks for the app notebooks "
        "to appear there",
    )
    emit()

    emit("4. Every reported figure is a declared target")
    not_targeted = sorted(path for path in reported if path not in targets)
    check(
        "every_reported_figure_is_a_target",
        not not_targeted,
        f"{not_targeted} are wrapped in report() but no TARGETS list asks for "
        "them, so the report would have a hole where each one should be.",
    )
    emit(f"   {len(reported) - len(not_targeted)}/{len(reported)} are targets")
    emit()

    check(
        "every_figure_has_labels",
        all(e["has_labels"] for e in reported.values()),
        "a report() call carries no labels=, so the report's table view cannot "
        "group it by task",
    )

    passed = sum(1 for c in checks if c["passed"])
    summary = {
        "task": "P6-T4",
        "n_checks": len(checks),
        "n_passed": passed,
        "n_failed": len(checks) - passed,
        "n_reported_figures": len(reported),
        "categories": {
            category: sum(1 for e in reported.values() if e["category"] == category)
            for category in categories
        },
        "n_caption_files": len(referenced),
        "figures": {
            path: {"category": e["category"], "caption": e["caption"]}
            for path, e in sorted(reported.items())
        },
        "checks": checks,
        "valid": not failures,
        "note": (
            "`snakemake --report` is a command, not a rule (ADR 0030 §3): it "
            "consumes the completed DAG. This audit is the part of P6-T4 that "
            "can be a rule, and it fails the build."
        ),
    }
    Path(snakemake.output.summary).write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    emit(f"wrote {snakemake.output.summary}")
    emit()

    if failures:
        for failure in failures:
            emit(f"  FAIL  {failure}")
        raise RuntimeError(
            f"{len(failures)} of {len(checks)} report checks failed — see "
            f"{snakemake.log[0]}. An uncaptioned figure still renders; it just "
            "renders uselessly, and nobody notices until a reader asks what "
            "they are looking at."
        )

    emit(f"CLEAN. {passed}/{len(checks)} checks passed: {len(reported)} figures,")
    emit("each captioned from a file that exists, categorised by phase, labelled,")
    emit("and asked for by a target.")
