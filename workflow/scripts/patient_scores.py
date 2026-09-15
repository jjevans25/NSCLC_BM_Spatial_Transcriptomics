"""P5-T2 — collapse Phase 2's per-AOI signature scores to one value per patient.

Owner task: P5-T2. Driven by rule p5t2_patient_scores.

Phase 2 scored six signatures on 23 `TIME` AOIs. A survival model needs one
value per patient, so something has to collapse the two `TIME-L` AOIs that P12
and P24 each contribute. PROJECT_PLAN §6 P5-T2 says "specify the aggregation
rule"; this script does not choose it — `config.yaml → survival.aggregation`
does, and ADR 0024 pre-registered it before any Phase 5 number existed:

  METHOD          mean of the per-AOI scores, matching ADR 0021 §1's
                  duplicate_rule so Phase 5 invents no second convention
  SCORING METHOD  ssgsea is PRIMARY; zscore is a declared sensitivity
                  (ADR 0012's shape). Phase 2 treated the AGREEMENT between
                  them as the robustness evidence at this n, and neither is the
                  tiebreaker for the other.

The collapse touches 2 of 21 patient-sites — P12 and P24, both lung, nothing in
the brain set — so it cannot be load-bearing, and the assertions below are what
make that checkable rather than asserted.

Output is the ANALYSIS-READY table P5-T3 fits directly: score, coverage and
outcome in one row, so no joining happens inside a modelling script.

THE CENSORED PATIENT WITH NO FOLLOW-UP TIME (ADR 0025)
------------------------------------------------------
A censored patient carries an event status and NO TIME. `Alive` is a status, not
a duration, and the brain cell is blank; Supplementary Data 1 has no
last-contact date, so this project holds no time-axis position for those three
patients at all.

**Patient 35 is the one in the Phase 5 cohort, and it is in BOTH arms.** It is
kept in the table with `has_followup = False` and `in_primary_fit = False`,
recorded in `survival_excluded.tsv` with a reason, and never imputed.
Administrative censoring at the cohort maximum was considered and rejected
(ADR 0025 §2): it would invent an observation AND place it exactly where a
median split is most sensitive, making the fabricated point the most influential
one in the analysis.

Consequence, fixed before any fit: **lung 12 and brain 7 enter a fit**, splits
6/6 and 3/4, and the brain arm lands EXACTLY on `survival.model.min_arm_size`.

WHAT IT DELIBERATELY DOES NOT DO
--------------------------------
  * **Fit, split or test anything.** P5-T3 owns the median split and the
    log-rank. A split computed here would be a threshold applied before the
    table it rests on had been checked.
  * **Filter on `assessable`.** Three of the six signatures failed their Phase 2
    coverage floor (`exhaustion` 1/6 and `tls` 1/5 in brain, `myeloid_m1` 2/10
    at BOTH sites). They are carried with the flag attached, because which
    signatures are measurable is a Phase 2 RESULT and dropping them here would
    hide it. Wherever the flag is False the finding is "not assessable", NEVER a
    prognostic null (ADR 0008) — a score built from genes at background measures
    background, whatever it is regressed against.
  * **Re-derive detection.** `signature_coverage.tsv` already has it, computed
    under ADR 0007's one rule; this script reads that file.
  * **Re-derive the cohort.** P5-T1's `clinical_cohort.tsv` owns `in_phase5`,
    and every patient here is asserted against it.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

FMT = "%.17g"

# aoi_code -> the arm whose endpoint column applies, per ADR 0024 §3. Keyed on
# aoi_code and never on `compartment`: `compartment` is degenerate across sites
# (L, LB and mLN are all `tumour`), so keying on it would pool the two arms.
ARM_BY_AOI_CODE = {"TIME-L": "lung", "TIME-B": "brain"}

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    survival = snakemake.params.survival
    aggregation = survival["aggregation"]
    method_name = aggregation["method"]
    primary_method = aggregation["scoring_method"]
    aoi_codes = survival["aoi_codes"]
    min_arm_size = survival["model"]["min_arm_size"]

    emit("P5-T2 — patient-level signature scores (ADR 0024, ADR 0025)")
    emit("=" * 70)
    emit()
    emit(f"  aggregation.method          {method_name}")
    emit(f"  primary scoring method      {primary_method}")
    emit(f"  sensitivity                 every other method, reported beside it")
    emit(f"  aoi_codes                   {aoi_codes}")
    emit()
    emit("  No Phase 5 result may be a headline claim. Brain is TIME-B n = 8")
    emit("  (hard constraint 8). A null at this n is UNINFORMATIVE, NOT")
    emit("  NEGATIVE (ADR 0024 §7).")
    emit()

    if method_name != "mean":
        raise RuntimeError(
            f"survival.aggregation.method is '{method_name}'. Only 'mean' is "
            "implemented, and it is what ADR 0024 pre-registered (matching "
            "ADR 0021 §1's duplicate_rule). Changing the aggregation is a "
            "stop-and-ask, not a config flip."
        )

    # --- inputs -------------------------------------------------------------
    scores = pd.read_csv(
        snakemake.input.scores, sep="\t", float_precision="round_trip"
    )
    coverage = pd.read_csv(
        snakemake.input.coverage, sep="\t", float_precision="round_trip"
    )
    clinical = pd.read_csv(
        snakemake.input.clinical, sep="\t", float_precision="round_trip"
    )
    emit(f"scores   : {len(scores)} rows  ({snakemake.input.scores})")
    emit(f"coverage : {len(coverage)} rows  ({snakemake.input.coverage})")
    emit(f"clinical : {len(clinical)} rows  ({snakemake.input.clinical})")

    scores = scores.loc[scores["aoi_code"].isin(aoi_codes)].copy()
    if scores.empty:
        raise RuntimeError(
            f"no signature scores carry an aoi_code in {aoi_codes}. Phase 2 "
            "scored TIME-L and TIME-B only; survival.aoi_codes must name those."
        )

    signatures = sorted(scores["signature"].unique())
    methods = sorted(scores["method"].unique())
    if primary_method not in methods:
        raise RuntimeError(
            f"survival.aggregation.scoring_method '{primary_method}' is not "
            f"among the scored methods {methods}."
        )
    emit(f"           {len(signatures)} signatures, methods {methods}")
    emit()

    # --- the collapse -------------------------------------------------------
    # patient_id is "P12"; clinical_cohort.tsv keys on the integer. Derived
    # here rather than parsed twice.
    scores["patient_number"] = (
        scores["patient_id"].astype(str).str.lstrip("P").astype(int)
    )

    emit("COLLAPSE (survival.aggregation, pre-registered)")
    emit("-" * 70)
    grouped = scores.groupby(
        ["patient_number", "site", "aoi_code", "signature", "method"],
        as_index=False,
    )
    patient = grouped.agg(
        score=("score", "mean"),
        n_aoi=("aoi_label", "nunique"),
        aoi_labels=("aoi_label", lambda s: ",".join(sorted(set(s)))),
    )

    combos = patient[["patient_number", "site"]].drop_duplicates()
    emit(f"  {len(scores)} AOI-level rows -> {len(patient)} patient-level rows")
    emit(f"  {len(combos)} patient-site combinations")

    dups = patient.loc[patient["n_aoi"] > 1, ["patient_number", "site", "n_aoi"]]
    dups = dups.drop_duplicates().sort_values("patient_number")
    emit(f"  {len(dups)} patient-sites needed a collapse:")
    for _, row in dups.iterrows():
        emit(f"    P{row['patient_number']:<3d} {row['site']:6s} "
             f"{row['n_aoi']} AOIs")
    emit()

    n_checks = 0

    # ASSERTION — a singleton must pass through BITWISE unchanged. 19 of 21
    # patient-sites are singletons; if the collapse perturbs them, it is wrong
    # everywhere and the two duplicates would hide it.
    # Compared by MERGE rather than by a row-by-row .loc: the source index is
    # non-unique (P12 and P24 each contribute two AOIs), and .loc on a
    # non-unique index returns a Series, which silently turns an equality test
    # into a type error at best and a wrong comparison at worst.
    KEYS = ["patient_number", "site", "signature", "method"]
    singles = patient.loc[patient["n_aoi"] == 1, KEYS + ["score"]]
    # The source is restricted to the same patient-sites before merging, so
    # `validate` can assert one-to-one: it inspects the whole right frame, and
    # the duplicate patient-sites would trip it even though nothing matches them.
    single_sites = set(
        map(tuple, patient.loc[patient["n_aoi"] == 1, ["patient_number", "site"]]
            .drop_duplicates().to_numpy())
    )
    source = scores.loc[
        [
            (pn, st) in single_sites
            for pn, st in zip(scores["patient_number"], scores["site"])
        ],
        KEYS + ["score"],
    ]
    check = singles.merge(
        source, on=KEYS, how="left",
        suffixes=("_agg", "_src"), validate="one_to_one",
    )
    if check["score_src"].isna().any():
        raise RuntimeError(
            "a single-AOI patient-site has no matching source row — the "
            "grouping keys do not round-trip."
        )
    perturbed = check.loc[check["score_agg"] != check["score_src"]]
    if not perturbed.empty:
        raise RuntimeError(
            f"{len(perturbed)} single-AOI patient-sites changed value under "
            f"aggregation:\n{perturbed.head().to_string(index=False)}\n"
            "The mean of one number is that number; a difference here means "
            "the grouping keys are wrong."
        )
    emit(f"  [ok] {len(singles)} single-AOI rows unchanged bitwise")
    n_checks += 1

    # ASSERTION — the duplicates are exactly the ones ADR 0021 §1 names. If
    # this moves, the design changed and the phase stops rather than adapting.
    observed_dups = sorted(set(dups["patient_number"]))
    dup_sites = sorted(set(dups["site"]))
    if observed_dups != [12, 24] or dup_sites != ["lung"]:
        raise RuntimeError(
            f"duplicate patient-sites are {observed_dups} in sites {dup_sites}; "
            "ADR 0021 §1 establishes P12 and P24, lung only, nothing in the "
            "brain set. A change means the design moved."
        )
    emit("  [ok] duplicates are P12 and P24, lung only (ADR 0021 §1)")
    n_checks += 1

    # ASSERTION — a mean must lie between its inputs. Catches a collapse that
    # silently picked one AOI, or summed instead of averaging.
    unbracketed = []
    for _, row in patient.loc[patient["n_aoi"] > 1].iterrows():
        members = scores.loc[
            (scores["patient_number"] == row["patient_number"])
            & (scores["site"] == row["site"])
            & (scores["signature"] == row["signature"])
            & (scores["method"] == row["method"]),
            "score",
        ]
        if not (members.min() <= row["score"] <= members.max()):
            unbracketed.append(
                (row["patient_number"], row["signature"], row["method"])
            )
    if unbracketed:
        raise RuntimeError(
            f"aggregated score outside the range of its inputs for "
            f"{unbracketed}. A mean cannot do that; the collapse is not a mean."
        )
    emit("  [ok] every collapsed score lies between its source AOIs")
    n_checks += 1

    expected_rows = len(combos) * len(signatures) * len(methods)
    if len(patient) != expected_rows:
        raise RuntimeError(
            f"{len(patient)} patient-level rows, expected {expected_rows} "
            f"({len(combos)} patient-sites x {len(signatures)} signatures x "
            f"{len(methods)} methods). The grid is not complete — a missing "
            "cell would silently shrink an arm at fit time."
        )
    if patient.duplicated(
        ["patient_number", "site", "signature", "method"]
    ).any():
        raise RuntimeError("a (patient, site, signature, method) cell appears twice.")
    emit(f"  [ok] complete grid, {expected_rows} rows, no cell twice")
    n_checks += 1

    # --- coverage join ------------------------------------------------------
    emit()
    emit("COVERAGE (Phase 2's floor, carried not applied)")
    emit("-" * 70)
    cov = coverage[["signature", "site", "coverage", "assessable", "genes_absent"]]
    before = len(patient)
    patient = patient.merge(cov, on=["signature", "site"], how="left")
    if len(patient) != before:
        raise RuntimeError(
            "the coverage join changed the row count — signature_coverage.tsv "
            "has more than one row per (signature, site)."
        )
    # A failed join must not read as "assessable unknown". Every cell must have
    # a verdict, because "not assessable" is a reportable FINDING here, not a
    # gap (ADR 0008, Gate 3).
    if patient["assessable"].isna().any():
        missing = patient.loc[
            patient["assessable"].isna(), ["signature", "site"]
        ].drop_duplicates()
        raise RuntimeError(
            f"no coverage verdict for {missing.to_dict('records')}. A null "
            "assessable flag would read as 'unknown' and get reported as a "
            "null rather than as not assessable."
        )
    n_checks += 1

    verdicts = (
        patient[["signature", "site", "assessable"]]
        .drop_duplicates()
        .sort_values(["signature", "site"])
    )
    n_not = int((~verdicts["assessable"]).sum())
    emit(f"  {n_not} of {len(verdicts)} signature x site cells are NOT assessable:")
    for _, row in verdicts.loc[~verdicts["assessable"]].iterrows():
        emit(f"    {row['signature']:22s} {row['site']}")
    emit("  [ok] every cell carries a verdict")
    emit()
    emit("  Carried, not filtered: which signatures are measurable is a Phase 2")
    emit("  RESULT. Wherever assessable is False the finding is 'NOT")
    emit("  ASSESSABLE', never a prognostic null (ADR 0008).")

    # --- outcome join -------------------------------------------------------
    emit()
    emit("OUTCOME (ADR 0024 §3, ADR 0025)")
    emit("-" * 70)
    clin = clinical.set_index("patient_number")

    # Every patient must be one P5-T1 already placed in the Phase 5 cohort.
    absent = sorted(set(patient["patient_number"]) - set(clin.index))
    if absent:
        raise RuntimeError(
            f"patients {absent} have signature scores but no row in "
            "clinical_cohort.tsv. P5-T1's join gate passed, so this means the "
            "two tables were built from different cohorts."
        )
    not_flagged = sorted(
        p for p in set(patient["patient_number"]) if not bool(clin.loc[p, "in_phase5"])
    )
    if not_flagged:
        raise RuntimeError(
            f"patients {not_flagged} carry scores but are not flagged "
            "in_phase5 in clinical_cohort.tsv. P5-T1 derives that flag from "
            "the same aoi_codes; a disagreement means one of them is stale."
        )
    emit(f"  [ok] all {patient['patient_number'].nunique()} patients are "
         "in_phase5 in clinical_cohort.tsv")
    n_checks += 1

    patient["arm"] = patient["aoi_code"].map(ARM_BY_AOI_CODE)
    if patient["arm"].isna().any():
        unmapped = sorted(set(patient.loc[patient["arm"].isna(), "aoi_code"]))
        raise RuntimeError(
            f"aoi_code(s) {unmapped} have no arm in ARM_BY_AOI_CODE. Each "
            "aoi_code must name which endpoint column applies (ADR 0024 §3)."
        )

    patient["months"] = [
        clin.loc[p, f"{arm}_months"]
        for p, arm in zip(patient["patient_number"], patient["arm"])
    ]
    patient["event"] = [
        int(clin.loc[p, f"{arm}_event"])
        for p, arm in zip(patient["patient_number"], patient["arm"])
    ]
    patient["has_followup"] = patient["months"].notna()

    # ASSERTION — the shape of the P5-T1 bug, re-checked AT the join rather
    # than trusted across it. There, `NaN is None` being False counted the
    # censored brain patient as a death with a missing time. A time and an
    # event flag that disagree is the signature of that class of error.
    contradictions = patient.loc[patient["has_followup"] & (patient["event"] == 0)]
    if not contradictions.empty:
        raise RuntimeError(
            f"{len(contradictions)} rows are censored (event = 0) yet carry a "
            "follow-up time. In this cohort censoring is spelled `Alive` or a "
            "blank, and neither carries a duration (ADR 0025) — a time here "
            "means the censoring logic upstream is wrong."
        )
    missing_time_events = patient.loc[
        (~patient["has_followup"]) & (patient["event"] == 1)
    ]
    if not missing_time_events.empty:
        raise RuntimeError(
            f"{len(missing_time_events)} rows record an event (event = 1) with "
            "no time. That is the exact failure ADR 0025 exists to prevent — a "
            "censored patient read as a death."
        )
    emit("  [ok] event flags and follow-up times agree in every row")
    n_checks += 1

    # --- the ADR 0025 exclusion --------------------------------------------
    patient["in_primary_fit"] = patient["has_followup"]

    excluded_rows = []
    excluded = (
        patient.loc[~patient["has_followup"],
                    ["patient_number", "site", "aoi_code", "arm"]]
        .drop_duplicates()
        .sort_values(["patient_number", "site"])
    )
    for _, row in excluded.iterrows():
        excluded_rows.append({
            "patient_id": f"P{row['patient_number']}",
            "patient_number": int(row["patient_number"]),
            "site": row["site"],
            "aoi_code": row["aoi_code"],
            "arm": row["arm"],
            # Phrased PER ARM. The two endpoint columns spell censoring
            # differently (ADR 0024 §3), so one generic sentence would be
            # inaccurate on one of the two rows — and this table is a
            # provenance record, not a note.
            "reason": (
                "censored with no follow-up time — the lung endpoint carries "
                "`Alive`, a status rather than a duration; Supplementary "
                "Data 1 has no last-contact date"
                if row["arm"] == "lung" else
                "censored with no follow-up time — the brain endpoint cell is "
                "blank; Supplementary Data 1 has no last-contact date"
            ),
            "action": "excluded_from_fits_retained_in_table",
            "adr": "ADR 0025",
        })
    excluded_table = pd.DataFrame(
        excluded_rows,
        columns=["patient_id", "patient_number", "site", "aoi_code", "arm",
                 "reason", "action", "adr"],
    )

    emit()
    emit(f"  {len(excluded_table)} patient-arms excluded from fits (ADR 0025):")
    for _, row in excluded_table.iterrows():
        emit(f"    {row['patient_id']:5s} {row['site']:6s} ({row['aoi_code']})")
    emit("    kept in the table with has_followup = False and")
    emit("    in_primary_fit = False — flag-don't-drop, never imputed")

    # --- arm sizes, before and after ---------------------------------------
    emit()
    emit("ARM SIZES (ADR 0025 §3 — stated before any fit)")
    emit("-" * 70)
    arm_stats = {}
    for code in aoi_codes:
        arm = ARM_BY_AOI_CODE[code]
        rows = patient.loc[
            (patient["aoi_code"] == code)
            & (patient["signature"] == signatures[0])
            & (patient["method"] == primary_method)
        ]
        cohort_n = int(len(rows))
        fit_n = int(rows["in_primary_fit"].sum())
        # A median split puts ties on one side; floor/ceil is the balanced case
        # an odd n allows, and it is what ADR 0025 §3 records.
        low, high = fit_n // 2, fit_n - fit_n // 2
        clears = low >= min_arm_size
        arm_stats[arm] = {
            "aoi_code": code,
            "n_cohort": cohort_n,
            "n_in_fit": fit_n,
            "n_excluded": cohort_n - fit_n,
            "split_sizes": [low, high],
            "min_arm_size": min_arm_size,
            "clears_min_arm_size": bool(clears),
        }
        flag = "clears" if clears else "BELOW"
        emit(f"  {arm:6s} ({code:7s})  cohort {cohort_n:2d}  "
             f"in fit {fit_n:2d}  split {low}/{high}  "
             f"{flag} min_arm_size = {min_arm_size}")

    if not arm_stats["brain"]["clears_min_arm_size"]:
        raise RuntimeError(
            f"the brain arm splits {arm_stats['brain']['split_sizes']}, below "
            f"survival.model.min_arm_size = {min_arm_size}. ADR 0025 §3 records "
            "that this arm has no patient to spare; below the floor the brain "
            "analysis is NOT ASSESSABLE rather than null, and P5-T3 must say so "
            "rather than fitting it."
        )
    n_checks += 1
    emit()
    emit("  The brain arm has NO PATIENT TO SPARE: one further exclusion drops")
    emit("  it below the floor, and the brain analysis then becomes not")
    emit("  assessable rather than null. Brain is TIME-B n = 8 in every claim")
    emit("  (hard constraint 8); 7 is what enters a fit. Both belong in any")
    emit("  reported sentence.")

    # --- write --------------------------------------------------------------
    patient["is_primary_method"] = patient["method"] == primary_method
    patient["aggregation"] = method_name
    patient = patient.sort_values(
        ["signature", "site", "method", "patient_number"]
    ).reset_index(drop=True)
    patient = patient[[
        "patient_number", "site", "aoi_code", "arm", "signature", "method",
        "is_primary_method", "score", "aggregation", "n_aoi", "aoi_labels",
        "coverage", "assessable", "genes_absent",
        "months", "event", "has_followup", "in_primary_fit",
    ]]
    patient.to_csv(
        snakemake.output.scores, sep="\t", index=False, float_format=FMT
    )
    excluded_table.to_csv(snakemake.output.excluded, sep="\t", index=False)

    emit()
    emit(f"assertions: {n_checks} of {n_checks} passed")

    summary = {
        "task": "P5-T2",
        "phase": 5,
        "stretch": True,
        "pre_registered": "ADR 0024 (aggregation), ADR 0025 (the exclusion)",
        "aggregation": {
            "method": method_name,
            "primary_scoring_method": primary_method,
            "sensitivity_scoring_methods": [
                m for m in methods if m != primary_method
            ],
            "n_patient_sites_collapsed": int(len(dups)),
            "collapsed_patients": [int(p) for p in observed_dups],
        },
        "n_rows": int(len(patient)),
        "n_patient_sites": int(len(combos)),
        "n_patients": int(patient["patient_number"].nunique()),
        "n_signatures": len(signatures),
        "signatures": signatures,
        "methods": methods,
        "arms": arm_stats,
        "excluded": {
            "n_patient_arms": int(len(excluded_table)),
            "patients": sorted({int(r["patient_number"]) for r in excluded_rows}),
            "adr": "ADR 0025",
            "action": "excluded from fits, retained in the table, never imputed",
        },
        "not_assessable": [
            {"signature": row["signature"], "site": row["site"]}
            for _, row in verdicts.loc[~verdicts["assessable"]].iterrows()
        ],
        "n_assertions": n_checks,
        "reporting_rule": (
            "No Phase 5 result may be a headline claim. Brain is TIME-B n = 8 "
            "inline in every brain claim (hard constraint 8); 7 is what enters "
            "a fit and both belong in the sentence. A null at 6/6 and 3/4 is "
            "UNINFORMATIVE, NOT NEGATIVE (ADR 0024 §7), and a separation at 3 "
            "versus 4 patients describes this cohort rather than supporting an "
            "inference. Every p-value carries n, effect size and a confidence "
            "interval (hard constraint 7). Wherever assessable is False the "
            "finding is 'not assessable', never a prognostic null (ADR 0008)."
        ),
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")

    emit(
        f"wrote {snakemake.output.scores}, {snakemake.output.excluded}, "
        f"{snakemake.output.summary}"
    )
