"""P5-T5 — correct the Phase 5 family ONCE, and account for everything tested.

Owner task: P5-T5. Driven by rule p5t5_multiplicity.

PROJECT_PLAN §6 P5-T5, in full: *"Multiplicity honesty. You tested 5 signatures
x 2 cohorts. Say so, with correction."* This is that task, and it is the only
place in Phase 5 where a q-value is computed.

THE FAMILY IS COUNTED HERE, NEVER QUOTED
-----------------------------------------
ADR 0026 exists because ADR 0024 §8's family size was a number someone wrote
down rather than counted, and it was **wrong**: §8 said "6 signatures x 2
cohorts = 12", which forgot that ADR 0024 §2 itself splits this study into two
arms that are never pooled. An undercounted denominator makes every q too
small -- the direction that manufactures significance.

So this script **derives** the family by counting rows in the two upstream
tables and asserts the total, rather than trusting any recorded number:

    survival_km_models.tsv   primary method, status == modelled   -> 8
    tcga_km_models.tsv       primary method                       -> 6
                                                             family = 14

Both upstream tables assert they carry NO q-column precisely so the correction
cannot happen twice over overlapping families (ADR 0026 §3).

WHAT IS IN THE FAMILY AND WHAT IS NOT
--------------------------------------
IN    the 8 assessable GeoMx signature x arm cells (P5-T3) and the 6 TCGA
      signatures (P5-T4), both on the PRIMARY scoring method. One family, not
      two: "does the same signature stratify an independent cohort" is a single
      question asked twice (ADR 0024 §8, reaffirmed by ADR 0026 §2).

OUT   * `zscore` fits, and TCGA's continuous Cox -- declared SENSITIVITIES
        (ADR 0012). P2-T3 and P3-T3 keep theirs out of their families too;
        running a second estimator does not double the search, it measures
        whether the primary is estimator-dependent. They are reported beside
        the corrected rows and labelled out-of-family.
      * The 4 GeoMx cells below Phase 2's coverage floor. ADR 0008 forbids
        reporting them as prognostic nulls AT ALL, so correcting over them
        would put un-interpretable rows in the denominator and DEFLATE every q.
        They appear in the accounting with no p and no q.

TIES ARE THE FAILURE MODE HERE, AND THIS FAMILY HAS THEM
---------------------------------------------------------
Phase 4 shipped a wrong FDR step-down that was correct-looking and wrong under
ties, and only a monotonicity assertion caught it (NEXT_STEPS lesson 2). This
family contains **two tied pairs**. So the correction is computed twice, by
`statsmodels.multipletests` and by an independent hand-rolled step-up, and the
two must agree exactly -- plus monotonicity is asserted directly.

WHAT IT DELIBERATELY DOES NOT DO
--------------------------------
  * **Re-fit anything.** Every estimate, interval and raw p is lifted from the
    upstream tables unchanged. A number that changed between P5-T3/T4 and here
    would mean one of them is not reproducible.
  * **Choose alpha or the method.** Both are pre-registered in
    `config.yaml -> survival.multiplicity` and asserted against what arrives.
  * **Decide the gate.** Gate 5 is "timebox respected", not "a result was
    found". This script reports; ADR 0028 records the verdict.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

FMT = "%.17g"

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    survival = snakemake.params.survival
    multiplicity = survival["multiplicity"]
    primary_method = survival["aggregation"]["scoring_method"]
    method = multiplicity["method"]
    alpha = multiplicity["alpha"]

    emit("P5-T5 — multiplicity honesty (ADR 0024 §8, ADR 0026, ADR 0027)")
    emit("=" * 72)
    emit()
    emit(f"  correction       {method}")
    emit(f"  alpha            {alpha}")
    emit(f"  primary method   {primary_method}")
    emit()
    emit("  THE FAMILY IS COUNTED HERE, NEVER QUOTED. ADR 0026 exists because")
    emit("  ADR 0024 §8's family size was written down rather than counted, and")
    emit("  was wrong -- an undercounted denominator makes every q too small.")
    emit()

    if method != "fdr_bh":
        raise RuntimeError(
            f"survival.multiplicity.method is '{method}'; only 'fdr_bh' is "
            "implemented and it is what ADR 0024 §8 pre-registered. Changing a "
            "multiplicity correction is a stop-and-ask."
        )

    n_checks = 0

    # --- assemble the family by COUNTING -----------------------------------
    geomx = pd.read_csv(
        snakemake.input.geomx, sep="\t", float_precision="round_trip"
    )
    tcga = pd.read_csv(
        snakemake.input.tcga, sep="\t", float_precision="round_trip"
    )
    emit("FAMILY (counted from the upstream tables)")
    emit("-" * 72)

    # Neither upstream table may already carry a q -- that is what makes "once"
    # enforceable rather than hoped for (ADR 0026 §3).
    for name, frame in (("survival_km_models.tsv", geomx),
                        ("tcga_km_models.tsv", tcga)):
        offending = [c for c in frame.columns if c.startswith("q_")]
        if offending:
            raise RuntimeError(
                f"{name} already carries {offending}. P5-T5 is the only place a "
                "q-value is computed (ADR 0026 §3); a q upstream would mean the "
                "family was corrected twice over overlapping sets."
            )
    emit("  [ok] neither upstream table carries a q-column")
    n_checks += 1

    g_primary = geomx.loc[
        (geomx["method"] == primary_method) & (geomx["status"] == "modelled")
    ].copy()
    g_excluded = geomx.loc[
        (geomx["method"] == primary_method) & (geomx["status"] != "modelled")
    ].copy()
    t_primary = tcga.loc[tcga["method"] == primary_method].copy()

    g_primary["cohort"] = "this study (GeoMx)"
    g_primary["stratum"] = g_primary["arm"]
    t_primary["cohort"] = t_primary["cohort"].fillna("TCGA-LUAD")
    t_primary["stratum"] = "lung (bulk tumour)"

    emit(f"  this study (P5-T3)  {len(g_primary):2d} assessable signature x arm cells")
    emit(f"  TCGA LUAD (P5-T4)   {len(t_primary):2d} signatures")
    family_n = len(g_primary) + len(t_primary)
    emit(f"  FAMILY              {family_n}")

    # Asserted against the pre-registered structure, but DERIVED from the data.
    expected = multiplicity["n_signatures"] * multiplicity["n_cohorts"]
    emit()
    emit(f"  config declares n_signatures x n_cohorts = "
         f"{multiplicity['n_signatures']} x {multiplicity['n_cohorts']} = {expected}")
    emit(f"  the counted family is {family_n}, and the difference is EXPECTED:")
    emit("    +6  this study has TWO ARMS, never pooled (ADR 0024 §2) — the")
    emit("        config keys describe signatures and cohorts, not the family")
    emit("    -4  four GeoMx cells are below Phase 2's coverage floor and are")
    emit("        NOT assessable, so they are not tested (ADR 0008)")
    if family_n != (multiplicity["n_signatures"] * 2 - len(g_excluded)
                    + multiplicity["n_signatures"]):
        raise RuntimeError(
            f"counted family {family_n} does not reconcile: "
            f"{multiplicity['n_signatures']} signatures x 2 arms "
            f"- {len(g_excluded)} not assessable + {multiplicity['n_signatures']} "
            "TCGA. The arithmetic ADR 0026 corrected has drifted again."
        )
    emit("  [ok] the counted family reconciles with the declared structure")
    n_checks += 1

    if len(t_primary) != multiplicity["n_signatures"]:
        raise RuntimeError(
            f"TCGA contributes {len(t_primary)} cells, expected "
            f"{multiplicity['n_signatures']} (ADR 0027 §2)."
        )
    n_checks += 1

    keep = [
        "cohort", "stratum", "signature", "n_patients_col", "n_events_col",
        "hazard_ratio", "ci_low", "ci_high", "p_raw",
    ]
    g_primary["n_patients_col"] = g_primary["n_in_fit"]
    g_primary["n_events_col"] = (
        g_primary["n_events_high"] + g_primary["n_events_low"]
    )
    t_primary["n_patients_col"] = t_primary["n_patients"]
    t_primary["n_events_col"] = t_primary["n_events"]

    family = pd.concat(
        [g_primary[keep], t_primary[keep]], ignore_index=True
    ).rename(columns={"n_patients_col": "n_patients",
                      "n_events_col": "n_events"})
    family["in_family"] = True
    family["comparator"] = np.where(
        (family["cohort"] == "TCGA-LUAD") & (family["signature"] == "myeloid_m1"),
        "TCGA-only — not assessable in GeoMx",
        "paired across cohorts",
    )

    if family["p_raw"].isna().any():
        raise RuntimeError("a family member has no raw p-value.")

    # --- the correction, computed TWICE ------------------------------------
    emit()
    emit(f"CORRECTION — {method} at alpha = {alpha}, over {family_n} tests")
    emit("-" * 72)
    p = family["p_raw"].to_numpy(dtype=float)

    reject, q_sm, _, _ = multipletests(p, alpha=alpha, method=method)

    # INDEPENDENT hand-rolled Benjamini-Hochberg step-up. Phase 4 shipped a
    # step-down that was wrong UNDER TIES and looked right; this family has two
    # tied pairs, so the two implementations must agree exactly.
    order = np.argsort(p, kind="mergesort")
    m = len(p)
    ranked = p[order] * m / np.arange(1, m + 1)
    q_manual = np.minimum.accumulate(ranked[::-1])[::-1]
    q_manual = np.minimum(q_manual, 1.0)
    q_hand = np.empty(m)
    q_hand[order] = q_manual

    delta = float(np.max(np.abs(q_sm - q_hand)))
    if delta > 1e-12:
        raise RuntimeError(
            f"statsmodels and the hand-rolled BH disagree by {delta:.3g}. "
            "This family contains tied p-values, which is exactly where "
            "Phase 4's step-down was wrong and looked right."
        )
    emit(f"  [ok] two implementations agree to {delta:.3g}")
    n_checks += 1

    family["q_bh"] = q_sm
    family["survives_fdr"] = reject

    # Monotonicity: q must be non-decreasing in p. The assertion that caught
    # Phase 4's defect.
    check = family.sort_values("p_raw")["q_bh"].to_numpy()
    if np.any(np.diff(check) < -1e-12):
        raise RuntimeError(
            "q-values are not monotone in p. A BH step-up cannot produce that; "
            "the correction is wrong (NEXT_STEPS lesson 2)."
        )
    emit("  [ok] q is monotone in p")
    n_checks += 1

    n_ties = int(len(p) - len(np.unique(np.round(p, 12))))
    emit(f"  {n_ties} tied p-value(s) in the family — handled identically by both")

    # A q may never be smaller than its raw p under BH.
    if (family["q_bh"] < family["p_raw"] - 1e-12).any():
        raise RuntimeError("a q-value is below its raw p-value.")
    emit("  [ok] no q below its raw p")
    n_checks += 1

    # --- results ------------------------------------------------------------
    family = family.sort_values("p_raw").reset_index(drop=True)
    emit()
    emit("RESULTS — the whole family, corrected once")
    emit("-" * 72)
    emit(f"  {'cohort':20s} {'stratum':18s} {'signature':22s} "
         f"{'n':>4s} {'HR':>7s} {'95% CI':>18s} {'p_raw':>8s} {'q_bh':>8s}")
    for _, row in family.iterrows():
        emit(f"  {row['cohort']:20s} {row['stratum']:18s} "
             f"{row['signature']:22s} {int(row['n_patients']):4d} "
             f"{row['hazard_ratio']:7.3f} "
             f"[{row['ci_low']:6.3f},{row['ci_high']:8.3f}] "
             f"{row['p_raw']:8.4f} {row['q_bh']:8.4f}"
             + ("  SURVIVES" if row["survives_fdr"] else ""))

    n_survive = int(family["survives_fdr"].sum())
    emit()
    emit(f"  {n_survive} of {family_n} tests survive {method} at FDR {alpha}.")
    emit(f"  smallest q = {family['q_bh'].min():.4f} "
         f"(raw p = {family['p_raw'].min():.4f})")

    # --- the honest accounting ---------------------------------------------
    emit()
    emit("ACCOUNTING — everything considered, not only everything tested")
    emit("-" * 72)
    n_cells_geomx = int(len(g_primary) + len(g_excluded))
    emit(f"  this study: {n_cells_geomx} signature x arm cells")
    emit(f"    {len(g_primary):2d} assessable and tested")
    emit(f"    {len(g_excluded):2d} below Phase 2's coverage floor — NOT tested, and")
    emit("       reported as 'not assessable', NEVER as a prognostic null")
    for _, row in g_excluded.iterrows():
        emit(f"         {row['arm']:6s} {row['signature']}")
    emit(f"  TCGA LUAD: {len(t_primary)} signatures, all tested")
    emit(f"  TOTAL CONSIDERED {n_cells_geomx + len(t_primary)}  "
         f"TOTAL TESTED {family_n}")
    emit()
    emit("  OUT OF FAMILY, reported beside it and labelled (ADR 0012):")
    emit(f"    zscore fits in both cohorts, and TCGA's continuous Cox.")
    emit("    Running a second estimator does not double the search; it")
    emit("    measures whether the primary is estimator-dependent.")

    accounting = pd.DataFrame([
        {"cohort": "this study (GeoMx)", "cells_considered": n_cells_geomx,
         "cells_tested": len(g_primary),
         "cells_not_assessable": len(g_excluded),
         "basis": "Phase 2 coverage floor, detected_in_aoi_fraction 0.5 (ADR 0008)"},
        {"cohort": "TCGA-LUAD", "cells_considered": len(t_primary),
         "cells_tested": len(t_primary), "cells_not_assessable": 0,
         "basis": "no detection floor applies — bulk RNA-seq has no negative-probe channel (ADR 0027 §2)"},
    ])

    not_assessable = g_excluded[
        ["arm", "signature", "coverage", "status", "note"]
    ].copy()
    not_assessable["cohort"] = "this study (GeoMx)"
    not_assessable["in_family"] = False

    # --- write ---------------------------------------------------------------
    family.to_csv(snakemake.output.family, sep="\t", index=False, float_format=FMT)
    accounting.to_csv(snakemake.output.accounting, sep="\t", index=False)
    not_assessable.to_csv(
        snakemake.output.not_assessable, sep="\t", index=False, float_format=FMT
    )

    emit()
    emit(f"assertions: {n_checks} of {n_checks} passed")

    summary = {
        "task": "P5-T5",
        "phase": 5,
        "pre_registered": "ADR 0024 §8, corrected by ADR 0026, completed by ADR 0027",
        "correction": {"method": method, "alpha": alpha,
                       "implementations": ["statsmodels.multipletests",
                                           "hand-rolled BH step-up"],
                       "max_disagreement": delta,
                       "n_tied_p_values": n_ties},
        "family": {
            "n": family_n,
            "this_study": int(len(g_primary)),
            "tcga": int(len(t_primary)),
            "counted_not_quoted": True,
            "note": (
                "Derived by counting rows in survival_km_models.tsv and "
                "tcga_km_models.tsv. ADR 0026 exists because ADR 0024 §8's "
                "family size was written down rather than counted, and was "
                "wrong — an undercounted denominator makes every q too small."
            ),
        },
        "out_of_family": {
            "sensitivities": ["zscore (both cohorts)", "TCGA continuous Cox"],
            "not_assessable_cells": int(len(g_excluded)),
            "reason": (
                "Sensitivities do not double the search (ADR 0012). "
                "Not-assessable cells may not be reported as prognostic nulls "
                "at all (ADR 0008), so correcting over them would deflate every q."
            ),
        },
        "result": {
            "n_survive_fdr": n_survive,
            "min_q": float(family["q_bh"].min()),
            "min_p_raw": float(family["p_raw"].min()),
            "smallest_q_test": {
                "cohort": family.iloc[0]["cohort"],
                "stratum": family.iloc[0]["stratum"],
                "signature": family.iloc[0]["signature"],
                "hazard_ratio": float(family.iloc[0]["hazard_ratio"]),
                "ci": [float(family.iloc[0]["ci_low"]),
                       float(family.iloc[0]["ci_high"])],
                "p_raw": float(family.iloc[0]["p_raw"]),
                "q_bh": float(family.iloc[0]["q_bh"]),
            },
        },
        "accounting": {
            "cells_considered": int(n_cells_geomx + len(t_primary)),
            "cells_tested": family_n,
            "cells_not_assessable": int(len(g_excluded)),
        },
        "n_assertions": n_checks,
        "reporting_rule": (
            "This is the ONLY Phase 5 table carrying a q-value, and only q is "
            "reportable — a raw p from P5-T3 or P5-T4 is not. No Phase 5 result "
            "may be a headline claim. Every brain claim states TIME-B n = 8 "
            "inline (hard constraint 8); every p carries n, effect size and a "
            "confidence interval (hard constraint 7). A null remains "
            "UNINFORMATIVE, NOT NEGATIVE in this study (ADR 0024 §7) — but NOT "
            "in TCGA, where n = 502 with 182 events, and that asymmetry is the "
            "deliverable. Cross-cohort comparison is bounded by ADR 0027 §5: "
            "GeoMx measures the PanCK-negative segment and TCGA whole bulk "
            "tumour, and there is no brain comparator at all."
        ),
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")

    emit(f"wrote {snakemake.output.family}, {snakemake.output.accounting}, "
         f"{snakemake.output.not_assessable}, {snakemake.output.summary}")
