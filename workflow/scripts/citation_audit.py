"""P6-T3 — audit CITATION.cff and codemeta.json against reality.

Owner task: P6-T3. Driven by rule p6t3_citation_audit.

PROJECT_PLAN §6 P6-T3 asks for `CITATION.cff` + `codemeta.json` citing the data,
the workflow, every major package and each skill used. Those two files are
hand-maintained — they are project metadata, not results — so this rule does not
generate them. It checks them, and it fails the build when they are wrong.

**Why a rule rather than care.** Both files carried
`https://github.com/jjevans25/nsclc-brainmet-time` for five phases. That is the
project *code*, not the repository name, and the URL 404s. The two files whose
entire job is identification pointed at nothing, through five gates, because
nothing ever compared them to `git remote`. Now something does.

**Versions are checked against the pins, never trusted.** This is ADR 0026's
lesson in a different costume: a number that is quoted rather than counted from
its source is a number that goes stale silently, and a reader has no way to tell
a current version string from a stale one. Every cited package version is
matched against `workflow/envs/<env>.<platform>.pin.txt` — the file that decides
what actually gets installed.

**The DOI policy is enforced, not remembered.** `fair.identifiers.doi` is null —
P6-T5 was abandoned by decision rather than deferred (ADR 0032) — so neither
citation file may carry a DOI. A placeholder DOI would be a persistent
identifier that does not persist, and it is exactly the kind of
plausible-looking string that survives review.

The check runs **in both directions**, which is what keeps the decision
reversible: with a DOI in config, a citation file missing it fails too. So the
mechanism stays honest whether or not one is ever minted.

**ADR 0029 is made durable here.** No environment YAML may contain a `pip:`
section, and every environment must have a pin for the platform whose pins the
project maintains. `conda list --explicit` cannot express a pip dependency, so a
`pip:` section produces a pin that silently omits it — which is how `rocrate`
reached Phase 6 declared, installed here, and absent from the clean room. A
decision that lives only in an ADR is a decision that gets re-made.
"""

import json
import re
import subprocess
from pathlib import Path

import yaml

# Cited title -> the token the conda pin uses. A translation table, not a
# threshold: conda names differ from project names, and `matplotlib` resolves to
# `matplotlib-base` in a pin while `ro-crate-py` is packaged as `rocrate`.
PIN_NAME = {
    "R": "r-base",
    "GeomxTools": "bioconductor-geomxtools",
    "standR": "bioconductor-standr",
    "SpatialDecon": "bioconductor-spatialdecon",
    "limma": "bioconductor-limma",
    "edgeR": "bioconductor-edger",
    "lme4": "r-lme4",
    "lmerTest": "r-lmertest",
    "emmeans": "r-emmeans",
    "matplotlib": "matplotlib-base",
    "ro-crate-py": "rocrate",
}

CFF_REQUIRED = [
    "cff-version",
    "title",
    "message",
    "type",
    "authors",
    "license",
    "repository-code",
    "version",
    "date-released",
]
CODEMETA_REQUIRED = [
    "@context",
    "@type",
    "name",
    "version",
    "license",
    "codeRepository",
    "author",
]

fair = snakemake.params.fair
log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

failures = []
checks = []


def check(name, condition, detail=""):
    checks.append({"check": name, "passed": bool(condition), "detail": detail})
    if not condition:
        failures.append(f"{name}: {detail}")
    return bool(condition)


def normalise_repo(url):
    return url.rstrip("/").removesuffix(".git").lower()


with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    emit("P6-T3 — citation and software metadata audit")
    emit("=" * 70)
    emit()

    cff = yaml.safe_load(Path(snakemake.input.cff).read_text(encoding="utf-8"))
    codemeta = json.loads(
        Path(snakemake.input.codemeta).read_text(encoding="utf-8")
    )

    # --- 1. both files are complete ---------------------------------------
    emit("1. Required fields")
    for field in CFF_REQUIRED:
        check(f"cff_has_{field.replace('-', '_')}", bool(cff.get(field)),
              f"CITATION.cff has no {field}")
    check(
        "cff_version_is_1_2_0",
        str(cff.get("cff-version")) == "1.2.0",
        f"cff-version is {cff.get('cff-version')!r}",
    )
    for field in CODEMETA_REQUIRED:
        check(
            f"codemeta_has_{field.replace('@', '').replace('-', '_')}",
            bool(codemeta.get(field)),
            f"codemeta.json has no {field}",
        )
    check(
        "codemeta_context_is_codemeta_2_0",
        "codemeta-2.0" in str(codemeta.get("@context", "")),
        f"@context is {codemeta.get('@context')!r}",
    )
    emit(f"   {sum(1 for c in checks if c['passed'])}/{len(checks)} present")
    emit()

    # --- 2. the two files agree with each other and with git --------------
    emit("2. Identity")
    check(
        "version_agrees_between_files",
        str(cff.get("version")) == str(codemeta.get("version")),
        f"CITATION.cff {cff.get('version')!r} vs codemeta.json "
        f"{codemeta.get('version')!r}",
    )

    declared = {
        "cff": cff.get("repository-code", ""),
        "codemeta": codemeta.get("codeRepository", ""),
        "config": fair["identifiers"]["repository"],
    }

    # The three declarations must agree with EACH OTHER unconditionally. This
    # half is a property of the project and is checkable anywhere.
    check(
        "repository_agrees_across_all_three_declarations",
        len({normalise_repo(url) for url in declared.values()}) == 1,
        f"CITATION.cff, codemeta.json and config.yaml disagree: {declared}",
    )

    # Comparing against `git remote` is a property of the CHECKOUT, not of the
    # project, and P6-T1's clean room proved the difference: a clone taken from
    # a local path, a fork, a mirror, or a tarball export with no git at all
    # fails a comparison that says nothing about whether the metadata is right.
    # A rule that only passes in one working directory makes the workflow
    # unrunnable for exactly the stranger it is meant to serve.
    #
    # So the comparison runs only where it MEANS something — when origin is a
    # web URL — and otherwise reports "not checkable here", which is a real
    # answer (the /gate command says so in as many words) rather than a pass
    # rounded up. It still catches what it was written for: both citation files
    # carried a URL that 404s for five phases, and this is what found it.
    try:
        remote = subprocess.run(
            ["git", "remote", "get-url", "origin"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        remote = ""

    remote_is_a_url = remote.startswith(("http://", "https://", "git@"))
    if remote_is_a_url:
        for label, url in declared.items():
            check(
                f"{label}_repository_matches_git_remote",
                normalise_repo(url) == normalise_repo(remote),
                f"{label} says {url!r}, git remote is {remote!r}",
            )
        emit(f"   remote   {remote}")
    else:
        emit(f"   remote   {remote or '(none)'} — NOT a web URL, so the")
        emit("            declared repository is not checkable in this")
        emit("            checkout. Reported, not passed. The three")
        emit("            declarations were still checked against each other.")
    emit(f"   version  {cff.get('version')}")
    emit()

    # --- 3. cited versions match the pins ---------------------------------
    emit("3. Cited versions against the environment pins")

    pinned = {}
    for pin in snakemake.input.pins:
        for line in Path(pin).read_text(encoding="utf-8").splitlines():
            if not line.startswith("http"):
                continue
            filename = line.split("/")[-1].split("#")[0]
            filename = re.sub(r"\.(conda|tar\.bz2)$", "", filename)
            # <name>-<version>-<build>; names contain hyphens, builds do not
            # contain them after the last two splits.
            parts = filename.rsplit("-", 2)
            if len(parts) == 3:
                pinned.setdefault(parts[0], set()).add(parts[1])

    unpinned_but_cited = []
    mismatched = []
    for reference in cff.get("references", []):
        if reference.get("type") != "software":
            continue
        title = reference.get("title", "")
        version = reference.get("version")
        if version is None:
            continue
        name = PIN_NAME.get(title, title)
        if name not in pinned:
            unpinned_but_cited.append(f"{title} ({name})")
            continue
        if str(version) not in pinned[name]:
            mismatched.append(
                f"{title}: cited {version}, pinned {sorted(pinned[name])}"
            )

    check(
        "cited_versions_match_pins",
        not mismatched,
        f"{len(mismatched)}: {mismatched}",
    )

    # Snakemake is the one exception and it is a real one: it drives the
    # workflow from the bootstrap environment and is deliberately not in any
    # per-rule env, so no pin can carry it. Checked against environment.yml's
    # constraint instead, which is the only declaration there is.
    snakemake_cited = next(
        (
            r.get("version")
            for r in cff.get("references", [])
            if r.get("title") == "Snakemake"
        ),
        None,
    )
    driver = yaml.safe_load(
        Path(snakemake.input.environment).read_text(encoding="utf-8")
    )
    constraint = next(
        (
            d
            for d in driver["dependencies"]
            if isinstance(d, str) and d.startswith("snakemake")
        ),
        "",
    )
    check(
        "snakemake_version_satisfies_environment_yml",
        bool(snakemake_cited)
        and constraint.split("=")[-1].rstrip("*").strip(".")
        == str(snakemake_cited).split(".")[0],
        f"cited {snakemake_cited!r} against environment.yml {constraint!r}",
    )

    # A cited package that no pin carries is not automatically wrong — it may
    # be a driver-environment tool — but it is unverifiable, so it is REPORTED
    # rather than silently accepted.
    emit(f"   {len(pinned)} packages across {len(snakemake.input.pins)} pins")
    emit(f"   {len(mismatched)} mismatched, {len(unpinned_but_cited)} cited but "
         "not in any pin (reported, not failed):")
    for item in unpinned_but_cited:
        emit(f"     - {item}")
    emit()

    # --- 4. the DOI policy --------------------------------------------------
    emit("4. DOI policy")
    doi = fair["identifiers"].get("doi")
    cff_doi = cff.get("doi") or any(
        i.get("type") == "doi" for i in cff.get("identifiers", [])
    )
    codemeta_doi = bool(codemeta.get("identifier", "").startswith("10.")) or bool(
        codemeta.get("doi")
    )
    if doi:
        check("cff_carries_the_doi", bool(cff_doi), "config has a DOI, CITATION.cff does not")
        check("codemeta_carries_the_doi", bool(codemeta_doi),
              "config has a DOI, codemeta.json does not")
        emit(f"   {doi}")
    else:
        check(
            "no_placeholder_doi",
            not cff_doi and not codemeta_doi,
            "config.yaml has no DOI but a citation file carries one — a "
            "persistent identifier that does not resolve is worse than none",
        )
        emit("   none — abandoned by decision (ADR 0032), not pending — and")
        emit("   neither citation file claims one")
    emit()

    # --- 5. ADR 0029, made durable -----------------------------------------
    emit("5. Every environment is pinnable (ADR 0029)")
    pip_sections = []
    missing_pins = []
    for env_yaml in snakemake.input.envs:
        spec = yaml.safe_load(Path(env_yaml).read_text(encoding="utf-8"))
        for dependency in spec.get("dependencies", []):
            if isinstance(dependency, dict) and "pip" in dependency:
                pip_sections.append(env_yaml)
        stem = Path(env_yaml).stem
        if not any(Path(p).name.startswith(f"{stem}.") for p in snakemake.input.pins):
            missing_pins.append(env_yaml)

    check(
        "no_pip_section_in_any_env",
        not pip_sections,
        f"{pip_sections} declare a `pip:` section. `conda list --explicit` "
        "records conda packages only, so a pip dependency is invisible to the "
        "pin and absent from a clean-room environment (ADR 0029).",
    )
    check(
        "every_env_has_a_pin",
        not missing_pins,
        f"{missing_pins} have no .pin.txt — the clean room would solve them "
        "fresh against whatever the channel serves that day.",
    )
    emit(f"   {len(snakemake.input.envs)} environments, "
         f"{len(snakemake.input.pins)} pins, {len(pip_sections)} pip sections")
    emit()

    passed = sum(1 for c in checks if c["passed"])
    summary = {
        "task": "P6-T3",
        "n_checks": len(checks),
        "n_passed": passed,
        "n_failed": len(checks) - passed,
        "version": cff.get("version"),
        "repository": declared["config"],
        "git_remote": remote or None,
        "git_remote_checked": remote_is_a_url,
        "doi": doi,
        "n_references": len(cff.get("references", [])),
        "n_software_cited": sum(
            1 for r in cff.get("references", []) if r.get("type") == "software"
        ),
        "n_versions_checked_against_pins": sum(
            1
            for r in cff.get("references", [])
            if r.get("type") == "software"
            and r.get("version")
            and PIN_NAME.get(r["title"], r["title"]) in pinned
        ),
        "cited_but_not_pinned": sorted(unpinned_but_cited),
        "checks": checks,
        "valid": not failures,
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
            f"{len(failures)} of {len(checks)} citation checks failed — see "
            f"{snakemake.log[0]}. Both citation files pointed at a 404 URL for "
            "five phases because nothing compared them to anything."
        )

    emit(f"CLEAN. {passed}/{len(checks)} checks passed: both files are complete,")
    emit(
        "they agree with each other"
        + (" and with git's origin," if remote_is_a_url else
           " (git's origin is not a web URL here, so it was not compared),")
    )
    emit("every cited version matches the pin that installs it, and no file")
    emit("claims a DOI that does not exist.")
