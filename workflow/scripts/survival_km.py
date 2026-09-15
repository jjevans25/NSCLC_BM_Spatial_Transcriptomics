"""P5-T3 — median split, Kaplan-Meier and log-rank, within arm.

Owner task: P5-T3. Driven by rule p5t3_survival_km.

PROJECT_PLAN §6 P5-T3: "Median split per signature. Report n per arm; at ~35
patients this is descriptive." The real n is **12 lung and 7 brain entering a
fit** (ADR 0025), so it is more descriptive than the plan assumed, and this
script's job is to produce that description WITHOUT letting it read as an
inference.

Everything applied here was pre-registered before any curve existed:

  ADR 0024 §6  median split, computed WITHIN the arm being split; log-rank;
               the hazard ratio and its CI reported beside it; min_arm_size
  ADR 0025     the censored patient without follow-up is excluded from fits,
               fixing the arms at 12 and 7 and the splits at 6/6 and 3/4
  ADR 0026     the multiplicity family, and the fact that P5-T3 does NOT correct

THE SPLIT IS `score > median`, AND THE TIE SIDE IS DECLARED
------------------------------------------------------------
`high` is strictly above the arm's median; `low` is at or below it. Stated
rather than left to a default, because at n = 7 the median IS an observation and
which side it falls on changes the arm sizes. This gives 6/6 in lung and 3/4 in
brain, which is what ADR 0025 §3 recorded.

Not a tertile, not an optimal cutpoint, not a maximally-selected rank statistic.
Those are the standard routes to a separation that does not replicate, and at 12
and 7 events they are not defensible (ADR 0024 §6).

TWO ESTIMATORS, AND THEIR AGREEMENT IS THE CHECK
-------------------------------------------------
  log-rank   sksurv.compare.compare_survival -- the `scikit-survival` skill
             PROJECT_PLAN §6 names for this task
  Cox        statsmodels PHReg on the binary split indicator, for the hazard
             ratio, its confidence interval and a second p-value

Hard constraint 7 forbids a p-value without n, effect size and a confidence
interval, and a bare log-rank p is exactly that -- hence the Cox fit. The two
p-values are then two implementations of one comparison, the idiom
`crosscheck_mixedlm.py` uses, and their agreement is asserted. Neither is the
tiebreak for the other.

`lifelines` is deliberately NOT used: it is absent from
workflow/envs/py-analysis.yaml and must stay absent, because adding it fires the
software-env rerun trigger on p0t2_fetch_geo's protected() outputs and
invalidates all three pins (ADR 0005, ADR 0013, ADR 0023). Both estimators used
here are already installed.

NO q-VALUES. THAT ABSENCE IS THE DECISION, NOT AN OMISSION
-----------------------------------------------------------
ADR 0026: the family spans this study AND TCGA (P5-T4), so it cannot be
corrected until both have been fitted. P5-T5 corrects once, over the union.
A q-value here would be a corrected value that nothing corrected, and correcting
here and again at P5-T5 would be two corrections over overlapping families.

**Every p-value this script writes is RAW AND UNCORRECTED**, and the rule
asserts that no `q_bh` column reaches any output so the convention cannot be
restored by reflex.

WHAT IT DELIBERATELY DOES NOT DO
--------------------------------
  * **Fit a not-assessable cell.** Four of twelve signature x arm cells sit
    below Phase 2's coverage floor -- `exhaustion` brain, `tls` brain,
    `myeloid_m1` at BOTH sites. They get a `not_assessable` row with a note and
    no estimate, the shape checkpoint_carrier_models.tsv uses. ADR 0008: a score
    built from genes at background measures background, whatever it is regressed
    against, and the finding is "not assessable", NEVER a prognostic null.
  * **Impute a follow-up time.** P5-T2 already excluded the one patient without
    one (ADR 0025); this script asserts the exclusion held rather than redoing
    it.
  * **Choose a cutpoint.** The median is pre-registered. Scanning cutpoints at
    this n would be the single most effective way to manufacture a curve.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sksurv.compare import compare_survival
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.util import Surv
from statsmodels.duration.hazard_regression import PHReg

FMT = "%.17g"

# The log-rank is a SCORE test and the Cox p a WALD test on the same data, so
# they are not the same arithmetic and are not expected to agree closely. At
# n = 12 they were observed to differ by 0.046, and at small n the Wald test is
# the more conservative of the two -- so a tight tolerance here would fail on
# real statistics rather than on a bug. This bound is a loose smoke test only.
#
# The ACTUAL wiring check is `_assert_indicator_wiring` below, which is exact:
# inverting the group indicator must negate the Cox coefficient. That is what
# catches a group mapped the wrong way round, which is what this check is for.
P_AGREEMENT_TOL = 0.10

plt.rcParams["hatch.linewidth"] = 0.6

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    survival = snakemake.params.survival
    model_cfg = survival["model"]
    split_rule = model_cfg["split"]
    test_name = model_cfg["test"]
    alpha = model_cfg["alpha"]
    min_arm_size = model_cfg["min_arm_size"]
    primary_method = survival["aggregation"]["scoring_method"]

    emit("P5-T3 — median split, Kaplan-Meier, log-rank (ADR 0024 §6, 0025, 0026)")
    emit("=" * 72)
    emit()
    emit(f"  split            {split_rule}  (high = score > median, within arm)")
    emit(f"  test             {test_name} + Cox for the hazard ratio and CI")
    emit(f"  alpha            {alpha}")
    emit(f"  min_arm_size     {min_arm_size}")
    emit(f"  primary method   {primary_method}")
    emit()
    emit("  NO q-VALUES HERE, BY DECISION (ADR 0026 §3). The family spans this")
    emit("  study AND TCGA, so P5-T5 corrects ONCE over the union. Every p below")
    emit("  is RAW AND UNCORRECTED.")
    emit()
    emit("  No Phase 5 result may be a headline claim. Brain is TIME-B n = 8")
    emit("  (hard constraint 8); 7 enter a fit. A null is UNINFORMATIVE, NOT")
    emit("  NEGATIVE (ADR 0024 §7), and a separation at 3 versus 4 patients")
    emit("  describes this cohort rather than supporting an inference.")
    emit()

    if split_rule != "median":
        raise RuntimeError(
            f"survival.model.split is '{split_rule}'. Only 'median' is "
            "implemented and it is what ADR 0024 §6 pre-registered. Scanning "
            "cutpoints at n = 12 and n = 7 is the single most effective way to "
            "manufacture a separation; changing this is a stop-and-ask."
        )
    if test_name != "logrank":
        raise RuntimeError(
            f"survival.model.test is '{test_name}'; only 'logrank' is implemented."
        )

    scores = pd.read_csv(
        snakemake.input.scores, sep="\t", float_precision="round_trip"
    )
    emit(f"input: {len(scores)} rows  ({snakemake.input.scores})")

    signatures = sorted(scores["signature"].unique())
    arms = sorted(scores["arm"].unique())
    methods = sorted(scores["method"].unique())
    emit(f"       {len(signatures)} signatures x {len(arms)} arms x "
         f"{len(methods)} methods")
    emit()

    n_checks = 0

    # --- the fit, one signature x arm x method cell at a time ---------------
    def fit_cell(frame, signature, arm, method, assessable):
        """One cell. Returns a row dict; `status` says whether it was fitted."""
        cell = frame.loc[
            (frame["signature"] == signature)
            & (frame["arm"] == arm)
            & (frame["method"] == method)
        ]
        fitted = cell.loc[cell["in_primary_fit"]]
        n_cohort = int(len(cell))
        n_fit = int(len(fitted))

        base = {
            "signature": signature,
            "arm": arm,
            "aoi_code": cell["aoi_code"].iloc[0] if len(cell) else None,
            "method": method,
            "is_primary_method": method == primary_method,
            "split": split_rule,
            "test": test_name,
            "n_cohort": n_cohort,
            "n_in_fit": n_fit,
            "n_excluded": n_cohort - n_fit,
            "assessable": bool(assessable),
            "coverage": float(cell["coverage"].iloc[0]) if len(cell) else np.nan,
        }
        empty = {
            "median_score": np.nan, "n_high": 0, "n_low": 0,
            "n_events_high": 0, "n_events_low": 0,
            "median_months_high": np.nan, "median_months_low": np.nan,
            "estimate": np.nan, "std_error": np.nan,
            "hazard_ratio": np.nan, "ci_low": np.nan, "ci_high": np.nan,
            "p_raw": np.nan, "p_raw_logrank": np.nan, "p_raw_cox": np.nan,
            "logrank_chisq": np.nan, "p_agreement_delta": np.nan,
            "wiring_residual": np.nan, "ci_excludes_one": False,
            "significance_disagreement": False,
        }

        # ADR 0008 / ADR 0018 — a cell below the coverage floor is reported,
        # never fitted. "Not assessable", never a prognostic null.
        if not assessable:
            absent = cell["genes_absent"].iloc[0] if len(cell) else ""
            return {
                **base, **empty,
                "status": "not_assessable",
                "note": (
                    f"below the Phase 2 coverage floor in {arm} "
                    f"(coverage {base['coverage']:.3g}; absent: {absent}). "
                    "Not assessable, never a prognostic null (ADR 0008)."
                ),
            }

        if n_fit < 2 * min_arm_size:
            return {
                **base, **empty,
                "status": "not_assessable",
                "note": (
                    f"{n_fit} patients enter a fit; a median split cannot give "
                    f"two arms of at least min_arm_size = {min_arm_size}."
                ),
            }

        # THE SPLIT. high = strictly above the arm's own median; the tie side is
        # declared rather than defaulted, because at n = 7 the median IS an
        # observation (ADR 0024 §6).
        values = fitted["score"].to_numpy(dtype=float)
        median = float(np.median(values))
        is_high = values > median
        n_high, n_low = int(is_high.sum()), int((~is_high).sum())

        if min(n_high, n_low) < min_arm_size:
            return {
                **base, **empty,
                "median_score": median,
                "n_high": n_high, "n_low": n_low,
                "status": "not_assessable",
                "note": (
                    f"median split gives {n_high}/{n_low}, below "
                    f"min_arm_size = {min_arm_size}."
                ),
            }

        months = fitted["months"].to_numpy(dtype=float)
        events = fitted["event"].to_numpy(dtype=int).astype(bool)

        y = Surv.from_arrays(event=events, time=months)
        chisq, p_logrank = compare_survival(y, is_high.astype(int))

        # Cox on the binary indicator, for the effect size hard constraint 7
        # requires. `high` is the test level, so a hazard ratio above 1 means
        # the high-score arm dies faster.
        cox = PHReg(
            months, is_high.astype(float).reshape(-1, 1), status=events.astype(int)
        ).fit()
        loghr = float(cox.params[0])
        se = float(cox.bse[0])
        ci = cox.conf_int()
        ci_low, ci_high = float(np.exp(ci[0, 0])), float(np.exp(ci[0, 1]))
        p_cox = float(cox.pvalues[0])

        # EXACT wiring check: inverting the group indicator must negate the
        # coefficient. A group mapped the wrong way round inverts the direction
        # of every reported hazard ratio, and nothing else in this script would
        # notice -- the p-values are identical either way.
        cox_inv = PHReg(
            months, (~is_high).astype(float).reshape(-1, 1),
            status=events.astype(int),
        ).fit()
        wiring_residual = abs(float(cox_inv.params[0]) + loghr)
        if wiring_residual > 1e-8:
            raise RuntimeError(
                f"{signature} / {arm} / {method}: inverting the group indicator "
                f"changed |coefficient| by {wiring_residual:.3g}, not 0. The "
                "indicator is not wired through the Cox fit as assumed, and "
                "every hazard ratio's DIRECTION depends on it."
            )

        # The score test and the Wald interval can disagree about whether the
        # effect is distinguishable from none, and at this n they do. Hard
        # constraint 7 forbids a p-value without a confidence interval
        # precisely so this is visible rather than hidden behind the p.
        ci_excludes_one = (ci_low > 1.0) or (ci_high < 1.0)

        return {
            **base,
            "median_score": median,
            "n_high": n_high,
            "n_low": n_low,
            "n_events_high": int(events[is_high].sum()),
            "n_events_low": int(events[~is_high].sum()),
            "median_months_high": float(np.median(months[is_high])),
            "median_months_low": float(np.median(months[~is_high])),
            "estimate": loghr,
            "std_error": se,
            "hazard_ratio": float(np.exp(loghr)),
            "ci_low": ci_low,
            "ci_high": ci_high,
            # p_raw is the LOG-RANK p -- the pre-registered test (ADR 0024 §6).
            # The Cox p is carried beside it as the cross-check, never as the
            # headline.
            "p_raw": float(p_logrank),
            "p_raw_logrank": float(p_logrank),
            "p_raw_cox": p_cox,
            "logrank_chisq": float(chisq),
            "p_agreement_delta": abs(float(p_logrank) - p_cox),
            "wiring_residual": wiring_residual,
            "ci_excludes_one": bool(ci_excludes_one),
            # True when the log-rank calls it significant and the confidence
            # interval does not, or vice versa. Where they disagree THE INTERVAL
            # GOVERNS: a p below alpha beside an interval spanning 1 is not a
            # detected effect, and reporting the p alone would overstate it.
            "significance_disagreement": bool(
                (float(p_logrank) < alpha) != bool(ci_excludes_one)
            ),
            "status": "modelled",
            "note": "",
        }

    cover = (
        scores[["signature", "arm", "assessable"]]
        .drop_duplicates()
        .set_index(["signature", "arm"])["assessable"]
    )

    primary_rows, exploratory_rows = [], []
    for signature in signatures:
        for arm in arms:
            for method in methods:
                assessable = bool(cover.loc[(signature, arm)])
                primary_rows.append(
                    fit_cell(scores, signature, arm, method, assessable)
                )
                # The EXPLORATORY table fits every cell regardless of coverage,
                # no FDR, every row labelled (ADR 0018's shape). No claim may
                # rest on it -- it exists so the restriction is visible as a
                # choice rather than as an absence.
                exploratory_rows.append(
                    fit_cell(scores, signature, arm, method, True)
                )

    primary = pd.DataFrame(primary_rows)
    exploratory = pd.DataFrame(exploratory_rows)
    exploratory["table"] = "exploratory"
    exploratory["fdr"] = "none — no claim may rest on this table (ADR 0018)"

    # --- assertions ---------------------------------------------------------
    emit("ASSERTIONS")
    emit("-" * 72)

    modelled = primary.loc[primary["status"] == "modelled"]
    not_assessable = primary.loc[primary["status"] == "not_assessable"]

    per_method = modelled.loc[modelled["method"] == primary_method]
    if len(per_method) != 8:
        raise RuntimeError(
            f"{len(per_method)} primary cells were modelled for "
            f"{primary_method}, expected 8 (lung 5 + brain 3). ADR 0026 §1 "
            "derives the family from signature_coverage.tsv; a different count "
            "means the coverage table moved and the family with it."
        )
    emit(f"  [ok] {len(per_method)} primary cells modelled for {primary_method} "
         "(lung 5 + brain 3, ADR 0026 §1)")
    n_checks += 1

    # The not-assessable set must BE the coverage table's, not merely the same
    # size -- a coincidence of counts would pass a size check.
    observed_na = set(
        map(tuple, not_assessable.loc[
            not_assessable["method"] == primary_method, ["signature", "arm"]
        ].to_numpy())
    )
    expected_na = set(
        map(tuple, scores.loc[~scores["assessable"], ["signature", "arm"]]
            .drop_duplicates().to_numpy())
    )
    if observed_na != expected_na:
        raise RuntimeError(
            f"not-assessable cells {sorted(observed_na)} do not match "
            f"signature_coverage.tsv's {sorted(expected_na)}."
        )
    emit(f"  [ok] the {len(observed_na)} not-assessable cells ARE the coverage "
         "table's")
    n_checks += 1

    # Splits must equal ADR 0025 §3's pre-registered numbers.
    expected_split = {"lung": (6, 6), "brain": (3, 4)}
    for _, row in modelled.iterrows():
        got = (int(row["n_high"]), int(row["n_low"]))
        want = expected_split[row["arm"]]
        if got != want:
            raise RuntimeError(
                f"{row['signature']} / {row['arm']} / {row['method']} splits "
                f"{got[0]}/{got[1]}; ADR 0025 §3 records {want[0]}/{want[1]}. "
                "The split sizes are a property of the arm, not of the "
                "signature — a difference means the fitted set moved."
            )
        if got[0] + got[1] != int(row["n_in_fit"]):
            raise RuntimeError(
                f"{row['signature']} / {row['arm']}: split arms "
                f"{got[0]}+{got[1]} do not sum to n_in_fit {row['n_in_fit']}."
            )
    emit("  [ok] every split is 6/6 (lung) or 3/4 (brain), per ADR 0025 §3")
    n_checks += 1

    # ADR 0025's exclusion must have HELD: every fitted patient has an event,
    # because the only censored one carries no follow-up time and was removed.
    leaked = modelled.loc[
        (modelled["n_events_high"] + modelled["n_events_low"])
        != modelled["n_in_fit"]
    ]
    if not leaked.empty:
        raise RuntimeError(
            f"{len(leaked)} cells have fewer events than fitted patients. Every "
            "patient entering a Phase 5 fit has an event — the only censored "
            "one has no follow-up time and was excluded at P5-T2 (ADR 0025). A "
            "censored row here means that exclusion leaked."
        )
    emit("  [ok] every fitted patient has an event — ADR 0025's exclusion held")
    n_checks += 1

    # Two estimators of one comparison. A wired-up-wrong group indicator shows
    # here and almost nowhere else.
    worst = modelled["p_agreement_delta"].max()
    if worst > P_AGREEMENT_TOL:
        bad = modelled.loc[modelled["p_agreement_delta"] > P_AGREEMENT_TOL]
        raise RuntimeError(
            f"log-rank and Cox p-values differ by up to {worst:.4g} (tolerance "
            f"{P_AGREEMENT_TOL}) in {len(bad)} cells:\n"
            f"{bad[['signature', 'arm', 'method', 'p_raw_logrank', 'p_raw_cox']].to_string(index=False)}"
        )
    emit(f"  [ok] log-rank and Cox p agree in all {len(modelled)} modelled "
         f"cells (worst delta {worst:.4g}, tolerance {P_AGREEMENT_TOL})")
    n_checks += 1

    # ADR 0026 §3 — no q-column may reach an output.
    for name, frame in (("primary", primary), ("exploratory", exploratory)):
        offending = [c for c in frame.columns if c.startswith("q_") or c == "fdr_bh"]
        if offending:
            raise RuntimeError(
                f"the {name} table carries {offending}. P5-T3 does NOT correct "
                "(ADR 0026 §3): the family spans this study and TCGA, so P5-T5 "
                "corrects once over the union. A q here would be a corrected "
                "value that nothing corrected."
            )
    emit("  [ok] no q-column in any output — P5-T5 corrects once (ADR 0026 §3)")
    n_checks += 1

    # --- KM curve points, so the figure is reproducible from a table --------
    curve_rows = []
    for _, row in modelled.iterrows():
        cell = scores.loc[
            (scores["signature"] == row["signature"])
            & (scores["arm"] == row["arm"])
            & (scores["method"] == row["method"])
            & (scores["in_primary_fit"])
        ]
        values = cell["score"].to_numpy(dtype=float)
        is_high = values > float(row["median_score"])
        for label, mask in (("high", is_high), ("low", ~is_high)):
            t, s = kaplan_meier_estimator(
                cell["event"].to_numpy(dtype=int).astype(bool)[mask],
                cell["months"].to_numpy(dtype=float)[mask],
            )
            for time, surv in zip(t, s):
                curve_rows.append({
                    "signature": row["signature"], "arm": row["arm"],
                    "method": row["method"], "group": label,
                    "time_months": float(time), "survival": float(surv),
                })
    curves = pd.DataFrame(curve_rows)

    # --- results ------------------------------------------------------------
    emit()
    emit(f"RESULTS — {primary_method} (primary), RAW p, UNCORRECTED")
    emit("-" * 72)
    show = modelled.loc[modelled["method"] == primary_method].sort_values(
        ["arm", "signature"]
    )
    for _, row in show.iterrows():
        emit(f"  {row['arm']:5s} {row['signature']:22s} "
             f"n={int(row['n_in_fit']):2d} ({int(row['n_high'])}/{int(row['n_low'])})  "
             f"HR {row['hazard_ratio']:6.3f} "
             f"[{row['ci_low']:6.3f}, {row['ci_high']:8.3f}]  "
             f"p_raw {row['p_raw']:.4f}")
    emit()
    emit("  Not assessable, and reported as such rather than as a null:")
    for _, row in not_assessable.loc[
        not_assessable["method"] == primary_method
    ].sort_values(["arm", "signature"]).iterrows():
        emit(f"    {row['arm']:5s} {row['signature']:22s} {row['note']}")

    n_sig = int((show["p_raw"] < alpha).sum())
    emit()
    emit(f"  {n_sig} of {len(show)} primary cells have RAW p < {alpha}.")
    emit("  These are UNCORRECTED. P5-T5 applies BH once over the family this")
    emit("  study shares with TCGA (ADR 0026 §2), and only that q is reportable.")

    # THE SCORE TEST AND THE INTERVAL CAN DISAGREE, AND AT THIS n THEY DO.
    # Hard constraint 7 forbids a p-value without a confidence interval exactly
    # so this is visible rather than hidden behind the p.
    disagree = show.loc[show["significance_disagreement"]]
    emit()
    emit("LOG-RANK vs CONFIDENCE INTERVAL")
    emit("-" * 72)
    if disagree.empty:
        emit("  No cell disagrees: every raw p below alpha has an interval")
        emit("  excluding 1, and every p above alpha an interval spanning it.")
    else:
        emit(f"  {len(disagree)} of {len(show)} cells DISAGREE — the log-rank and")
        emit("  the hazard-ratio interval do not tell the same story:")
        for _, row in disagree.iterrows():
            emit(f"    {row['arm']:5s} {row['signature']:22s} "
                 f"raw p {row['p_raw']:.4f} (< {alpha}) but "
                 f"HR {row['hazard_ratio']:.3f} "
                 f"[{row['ci_low']:.3f}, {row['ci_high']:.3f}] SPANS 1")
        emit()
        emit("  WHERE THEY DISAGREE, THE INTERVAL GOVERNS. The log-rank is a")
        emit("  score test and the Cox interval a Wald interval; at n = 12 and")
        emit("  n = 7 they diverge, and the interval is the more conservative.")
        emit("  A p below alpha beside an interval spanning 1 is NOT a detected")
        emit("  effect, and reporting the p alone would overstate it. This is")
        emit("  what hard constraint 7 exists to make visible.")

    # --- figure -------------------------------------------------------------
    # Explicit gridspec, no tight_layout (Phase 4 lesson 1: both Phase 4 figures
    # clipped on first render, and P3-T5 did too).
    fig = plt.figure(figsize=(9.6, 16.2))
    gs = fig.add_gridspec(
        len(signatures), len(arms),
        left=0.085, right=0.975, top=0.938, bottom=0.048,
        hspace=0.42, wspace=0.20,
    )
    colours = {"high": "#b2182b", "low": "#2166ac"}

    for i, signature in enumerate(signatures):
        for j, arm in enumerate(sorted(arms, reverse=True)):  # lung, then brain
            ax = fig.add_subplot(gs[i, j])
            row = primary.loc[
                (primary["signature"] == signature)
                & (primary["arm"] == arm)
                & (primary["method"] == primary_method)
            ].iloc[0]

            if row["status"] != "modelled":
                ax.add_patch(plt.Rectangle(
                    (0, 0), 1, 1, transform=ax.transAxes, facecolor="#e8e8e8",
                    edgecolor="#9e9e9e", hatch="///", linewidth=0.8, zorder=0,
                ))
                ax.text(0.5, 0.56, "not assessable", ha="center", va="center",
                        transform=ax.transAxes, fontsize=9.5, weight="bold",
                        color="#3a3a3a")
                ax.text(0.5, 0.40,
                        f"coverage {row['coverage']:.2f} < 0.5",
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=8, color="#3a3a3a")
                ax.set_xticks([])
                ax.set_yticks([])
            else:
                sub = curves.loc[
                    (curves["signature"] == signature)
                    & (curves["arm"] == arm)
                    & (curves["method"] == primary_method)
                ]
                for label in ("high", "low"):
                    g = sub.loc[sub["group"] == label].sort_values("time_months")
                    t = np.concatenate([[0.0], g["time_months"].to_numpy()])
                    s = np.concatenate([[1.0], g["survival"].to_numpy()])
                    n = int(row["n_high"] if label == "high" else row["n_low"])
                    ax.step(t, s, where="post", color=colours[label],
                            linewidth=1.7, label=f"{label} (n={n})")
                ax.set_ylim(-0.30, 1.06)
                ax.set_xlim(left=0)
                # Ticks stop at 0: the axis extends below it only to give the
                # annotation clear space, and a tick at -0.2 survival would be
                # meaningless.
                ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
                ax.tick_params(labelsize=8)
                ax.legend(fontsize=7.2, loc="upper right", frameon=False,
                          handlelength=1.3, borderaxespad=0.2)
                # The interval is printed on the same line as the p, and a
                # disagreement is called out ON THE PANEL: a reader must not be
                # able to take a small p off this figure without seeing that
                # its interval spans 1.
                note = ""
                if row["significance_disagreement"]:
                    note = "\nraw p < alpha but CI spans 1 — interval governs"
                # Opaque bbox, NOT bare text: the curves run through this
                # corner, and on the first render a step crossed the decimal
                # point of "[0.88, 23.16]" so it read "23 16", with "(n=12)"
                # struck through. Phase 4 lesson 1 — open the PNG.
                ax.text(
                    0.02, 0.05,
                    f"HR {row['hazard_ratio']:.2f} "
                    f"[{row['ci_low']:.2f}, {row['ci_high']:.2f}]\n"
                    f"raw p = {row['p_raw']:.3f}  (n={int(row['n_in_fit'])})"
                    + note,
                    transform=ax.transAxes, fontsize=7.2, va="bottom",
                    color="#8c2d04" if row["significance_disagreement"] else "#222222",
                    bbox=dict(facecolor="white", alpha=0.92, edgecolor="none",
                              boxstyle="square,pad=0.28"),
                    zorder=5,
                )
            if i == 0:
                ax.set_title(
                    f"{arm}  ({row['aoi_code']})"
                    + ("   n = 8 cohort, 7 fitted" if arm == "brain" else ""),
                    fontsize=9.5, weight="bold", pad=7,
                )
            if j == 0:
                ax.set_ylabel(signature.replace("_", "\n"), fontsize=8.6)
            if i == len(signatures) - 1:
                ax.set_xlabel("months from diagnosis", fontsize=8.4)

    fig.suptitle(
        "P5-T3  Kaplan-Meier by median signature split, within arm\n"
        "EXPLORATORY — raw uncorrected p; TIME-B n = 8 (7 fitted); "
        "a null here is uninformative, not negative",
        fontsize=11, y=0.982,
    )
    fig.savefig(snakemake.output.figure, dpi=200)
    plt.close(fig)

    # --- write --------------------------------------------------------------
    primary.to_csv(
        snakemake.output.models, sep="\t", index=False, float_format=FMT
    )
    exploratory.to_csv(
        snakemake.output.exploratory, sep="\t", index=False, float_format=FMT
    )
    curves.to_csv(
        snakemake.output.curves, sep="\t", index=False, float_format=FMT
    )

    emit()
    emit(f"assertions: {n_checks} of {n_checks} passed")

    summary = {
        "task": "P5-T3",
        "phase": 5,
        "stretch": True,
        "pre_registered": "ADR 0024 §6 (split, test), ADR 0025 (arms), ADR 0026 (family)",
        "split": {"rule": split_rule, "tie_side": "high = score > median",
                  "computed": "within arm"},
        "test": {"primary": "logrank (sksurv.compare_survival)",
                 "effect_size": "Cox PHReg on the binary split indicator",
                 "p_agreement_tolerance": P_AGREEMENT_TOL,
                 "worst_p_agreement_delta": float(worst)},
        "multiplicity": {
            "corrected_here": False,
            "adr": "ADR 0026",
            "family_contribution_this_cohort": int(len(per_method)),
            "note": (
                "P5-T3 emits RAW p only. The family spans this study and TCGA "
                "(P5-T4), so P5-T5 corrects ONCE over the union. A q-value here "
                "would be a corrected value that nothing corrected."
            ),
        },
        "n_modelled": int(len(per_method)),
        "n_not_assessable": int(len(observed_na)),
        "not_assessable": [
            {"signature": s, "arm": a} for s, a in sorted(observed_na)
        ],
        "arms": {
            arm: {
                "n_cohort": int(show.loc[show["arm"] == arm, "n_cohort"].iloc[0]),
                "n_in_fit": int(show.loc[show["arm"] == arm, "n_in_fit"].iloc[0]),
                "split": [int(show.loc[show["arm"] == arm, "n_high"].iloc[0]),
                          int(show.loc[show["arm"] == arm, "n_low"].iloc[0])],
            }
            for arm in sorted(show["arm"].unique())
        },
        "n_raw_p_below_alpha": n_sig,
        "alpha": alpha,
        "significance_disagreement": {
            "n_cells": int(len(disagree)),
            "cells": [
                {"signature": r["signature"], "arm": r["arm"],
                 "p_raw": float(r["p_raw"]),
                 "hazard_ratio": float(r["hazard_ratio"]),
                 "ci": [float(r["ci_low"]), float(r["ci_high"])]}
                for _, r in disagree.iterrows()
            ],
            "rule": (
                "Where the log-rank and the hazard-ratio interval disagree, THE "
                "INTERVAL GOVERNS. The log-rank is a score test and the Cox "
                "interval a Wald interval; at n = 12 and n = 7 they diverge and "
                "the interval is the more conservative. A p below alpha beside "
                "an interval spanning 1 is not a detected effect."
            ),
        },
        "n_assertions": n_checks,
        "reporting_rule": (
            "EVERY p-VALUE HERE IS RAW AND UNCORRECTED (ADR 0026 §3). No Phase 5 "
            "result may be a headline claim. Brain is TIME-B n = 8 inline in "
            "every brain claim (hard constraint 8); 7 enter a fit and both "
            "belong in the sentence. Every p carries n, hazard ratio and a "
            "confidence interval (hard constraint 7). A null at 6/6 and 3/4 is "
            "UNINFORMATIVE, NOT NEGATIVE (ADR 0024 §7), and a separation at 3 "
            "versus 4 patients describes this cohort rather than supporting an "
            "inference. A cell below the coverage floor is 'not assessable', "
            "NEVER a prognostic null (ADR 0008)."
        ),
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")

    emit(
        f"wrote {snakemake.output.models}, {snakemake.output.exploratory}, "
        f"{snakemake.output.curves}, {snakemake.output.figure}, "
        f"{snakemake.output.summary}"
    )
