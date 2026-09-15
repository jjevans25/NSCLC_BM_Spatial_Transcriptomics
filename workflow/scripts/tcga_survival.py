"""P5-T4 — the TCGA LUAD external check on Phase 5's null.

Owner task: P5-T4. Driven by rule p5t4c_tcga_survival.

P5-T3 found no immune signature stratified survival here: across eight
assessable signature x arm cells every hazard-ratio interval spans 1, at lung
n = 12 and brain n = 7. PROJECT_PLAN §6 P5-T4 exists for that outcome --

    "Does the same signature stratify a properly powered independent cohort?
     A NEGATIVE HERE MOSTLY TELLS YOU ABOUT YOUR n, NOT ABOUT THE BIOLOGY."

-- because at ~500 patients TCGA can separate a power limit from an absent
effect, which is the one thing this project's own n cannot do. **A null in TCGA
means something a null here does not.**

Everything applied was pre-registered in ADR 0027 before any TCGA number existed.

TWO LIMITATIONS THAT ARE PART OF THE DELIVERABLE, NOT FOOTNOTES
---------------------------------------------------------------
1. **THIS IS NOT THE SAME MEASUREMENT.** A GeoMx `TIME-L` score comes from the
   PanCK-negative segment of an ROI sited in a CD45-rich region -- and per Q2's
   live caveat that is not a CD45-sorted population either. A TCGA score comes
   from WHOLE BULK TUMOUR: tumour cells, stroma and immune together. "The same
   signature" is computed on two different things, so a disagreement between the
   cohorts is not necessarily a disagreement about biology, and an agreement is
   not necessarily a replication. Same class of caveat as Phase 2's
   reference-matrix confound.
2. **THERE IS NO BRAIN COMPARATOR, AT ALL.** TCGA LUAD is primary lung. The brain
   arm -- the one carrying TIME-B n = 8, the binding constraint on the project --
   gets no external check from this task. P5-T4 contextualises the LUNG null only.

THE FAMILY IS TCGA'S 6, AND THAT DISCHARGES ADR 0026's DEFERRAL
----------------------------------------------------------------
All six signatures are tested. The Phase 2 detection floor is a property of THIS
assay -- q3 above a multiple of that AOI's NegProbe-WTX (ADR 0007) -- and TCGA
is bulk RNA-seq with no negative-probe channel, so nothing in it corresponds to
that threshold. Restricting to the five assessable in the GeoMx lung arm would
import a limitation the external cohort does not have. `myeloid_m1` is therefore
tested and reportable here, labelled "TCGA-only, no GeoMx comparator".

So the Phase 5 family is 8 (P5-T3) + 6 (here) = 14, and **P5-T5 corrects once
over all of it**. This script computes NO q-value and asserts no `q_` column
reaches any output -- the same guard survival_km.py carries (ADR 0026 §3).

WHAT IT DELIBERATELY DOES NOT DO
--------------------------------
  * **Apply a detection floor.** There is none to apply. It computes a GENE
    PRESENCE table instead -- whether a signature's genes are in the Xena matrix
    at all -- and labels it as presence, never detection. "Not in this matrix"
    and "at background in this assay" are unrelated limitations and conflating
    them would merge two different caveats into one wrong one.
  * **Change the estimator to suit the larger n.** The median split is the
    primary because COMPARABILITY WITH P5-T3 is the entire point: a different
    estimator would confound "different cohort" with "different method". The
    continuous Cox, where TCGA's power actually shows, is a declared sensitivity
    (ADR 0012's shape) and is not the tiebreak.
  * **Impute a follow-up time.** Cases with missing or non-positive OS.time are
    excluded and counted -- ADR 0025's shape, recorded and never invented.
"""

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import gseapy as gp
from sksurv.compare import compare_survival
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.util import Surv
from statsmodels.duration.hazard_regression import PHReg

FMT = "%.17g"

# Lifted from survival_km.py, where the reasoning is recorded: the log-rank is a
# SCORE test and the Cox p a WALD test, so at any n they differ, and a tight
# bound would fail on real statistics. The actual wiring check is exact.
P_AGREEMENT_TOL = 0.10

plt.rcParams["hatch.linewidth"] = 0.6

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    tcga = snakemake.params.tcga
    survival_cfg = snakemake.params.survival
    signatures = snakemake.params.signatures
    scoring = snakemake.params.scoring
    seed = snakemake.params.seed

    cols = tcga["columns"]
    expect = tcga["expect"]
    model_cfg = survival_cfg["model"]
    alpha = model_cfg["alpha"]
    min_arm_size = model_cfg["min_arm_size"]
    agg_method = survival_cfg["aggregation"]["method"]
    primary_method = survival_cfg["aggregation"]["scoring_method"]

    emit("P5-T4 — TCGA LUAD external check (ADR 0027, pre-registered)")
    emit("=" * 72)
    emit()
    emit(f"  cohort           {tcga['cohort']}  via {tcga['provider']}")
    emit(f"  endpoint         {cols['event']} / {cols['time']}  "
         f"({tcga['endpoint_source']})")
    emit(f"  sample type      {expect['sample_type_code']} (primary solid tumour)")
    emit(f"  days per month   {expect['days_per_month']}")
    emit()
    emit("  PROJECT_PLAN §6 P5-T4: 'a negative here mostly tells you about your")
    emit("  n, not about the biology'. At ~500 patients this cohort can separate")
    emit("  a power limit from an absent effect. A NULL IN TCGA MEANS SOMETHING")
    emit("  A NULL IN THIS STUDY DOES NOT.")
    emit()
    emit("  THIS IS NOT THE SAME MEASUREMENT. GeoMx TIME-L is the PanCK-negative")
    emit("  SEGMENT of an ROI; TCGA is WHOLE BULK TUMOUR. A disagreement between")
    emit("  cohorts is not necessarily a disagreement about biology.")
    emit("  THERE IS NO BRAIN COMPARATOR: TCGA LUAD is primary lung, so this")
    emit("  contextualises the LUNG null only (ADR 0027 §5).")
    emit()
    emit("  NO q-VALUES HERE (ADR 0026 §3). The family is 8 (P5-T3) + 6 (here)")
    emit("  = 14 and P5-T5 corrects ONCE over all of it. Every p below is RAW.")
    emit()

    n_checks = 0

    # --- expression ---------------------------------------------------------
    expr_path = snakemake.input.expression
    with gzip.open(expr_path, "rt") as handle:
        expr = pd.read_csv(handle, sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    emit(f"expression: {expr.shape[0]} genes x {expr.shape[1]} samples  ({expr_path})")

    # TCGA barcode positions 14-15 are the sample type. `01` primary solid
    # tumour; `11` solid tissue NORMAL and `02` RECURRENT tumour are both present
    # in this matrix, and scoring either as though it were a primary tumour would
    # be a silent category error.
    def sample_type(barcode):
        parts = str(barcode).split("-")
        return parts[3][:2] if len(parts) > 3 else ""

    types = pd.Series([sample_type(s) for s in expr.columns], index=expr.columns)
    emit(f"            sample-type codes present: "
         f"{types.value_counts().to_dict()}")
    keep = types == expect["sample_type_code"]
    dropped_types = types.loc[~keep].value_counts().to_dict()
    expr = expr.loc[:, keep.to_numpy()]
    emit(f"            kept {expr.shape[1]} type-{expect['sample_type_code']} "
         f"samples; dropped {dropped_types}")

    stray = {t for t in (sample_type(s) for s in expr.columns)
             if t != expect["sample_type_code"]}
    if stray:
        raise RuntimeError(
            f"sample types {sorted(stray)} survived the primary-tumour filter."
        )
    emit("  [ok] every retained sample is a primary solid tumour")
    n_checks += 1

    # --- gene PRESENCE, which is not detection ------------------------------
    emit()
    emit("GENE PRESENCE (not detection — ADR 0027 §2)")
    emit("-" * 72)
    measured = set(expr.index)
    presence_rows = []
    gene_sets = {}
    for name in sorted(signatures):
        genes = list(signatures[name]["genes"])
        present = [g for g in genes if g in measured]
        absent = [g for g in genes if g not in measured]
        fraction = len(present) / len(genes)
        scoreable = fraction >= expect["min_gene_presence"]
        presence_rows.append({
            "signature": name,
            "cohort": tcga["cohort"],
            "n_genes": len(genes),
            "n_genes_present": len(present),
            "presence": fraction,
            "scoreable": scoreable,
            "genes_absent": ",".join(absent),
            "basis": "gene presence in the expression matrix, NOT detection above background",
        })
        if scoreable:
            gene_sets[name] = present
        emit(f"  {name:22s} {len(present):2d}/{len(genes):2d} present "
             f"({fraction:.3f})  {'scoreable' if scoreable else 'NOT SCOREABLE'}"
             + (f"  absent: {','.join(absent)}" if absent else ""))
    presence = pd.DataFrame(presence_rows)
    emit()
    emit("  PRESENCE IS NOT DETECTION. Phase 2's floor is q3 above a multiple of")
    emit("  that AOI's NegProbe-WTX (ADR 0007) — a property of the GeoMx assay.")
    emit("  TCGA has no negative-probe channel, so nothing here corresponds to")
    emit("  it. 'Not in this matrix' and 'at background in this assay' are")
    emit("  unrelated limitations and must never be read as the same thing.")

    if len(gene_sets) != expect["n_signatures"]:
        raise RuntimeError(
            f"{len(gene_sets)} signatures are scoreable, expected "
            f"{expect['n_signatures']} (tcga.expect.n_signatures). That count "
            "is TCGA's contribution to the multiplicity family (ADR 0027 §2); "
            "if it moves, P5-T5's family moves with it."
        )
    n_checks += 1

    # --- scoring, the same two methods Phase 2 used -------------------------
    emit()
    emit(f"SCORING ({len(gene_sets)} sets, {expr.shape[1]} samples)")
    emit("-" * 72)
    # Same call shape as score_signatures.py:160 — min_size explicit because
    # gseapy defaults it to 15 and every set here is 4-13 genes.
    ss = gp.ssgsea(
        data=expr,
        gene_sets={k: list(v) for k, v in gene_sets.items()},
        outdir=None,
        min_size=scoring["min_set_size"],
        max_size=max(len(v) for v in gene_sets.values()),
        permutation_num=0,
        no_plot=True,
        seed=seed,
        threads=snakemake.threads,
    )
    res = ss.res2d.copy()
    if res.empty:
        raise RuntimeError(
            "gseapy.ssgsea returned no rows; the usual cause is min_size "
            f"({scoring['min_set_size']}) exceeding a set's length."
        )
    ssgsea_scores = (
        res.rename(columns={"Name": "sample", "Term": "signature"})
        .assign(score=lambda d: d["NES"].astype(float), method="ssgsea")
        .loc[:, ["sample", "signature", "method", "score"]]
    )
    emit(f"  ssgsea: {len(ssgsea_scores)} rows")

    # z-score mean, standardised across THIS cohort's samples only — the same
    # construction score_signatures.py uses within its 23 AOIs.
    z = (expr.T - expr.T.mean(axis=0)) / expr.T.std(axis=0, ddof=1)
    zrows = []
    for name, genes in gene_sets.items():
        zrows.append(pd.DataFrame({
            "sample": z.index, "signature": name, "method": "zscore",
            "score": z[genes].mean(axis=1).to_numpy(),
        }))
    scores = pd.concat([ssgsea_scores, pd.concat(zrows, ignore_index=True)],
                       ignore_index=True)
    emit(f"  zscore: {len(scores) - len(ssgsea_scores)} rows")

    # --- survival -----------------------------------------------------------
    emit()
    emit("SURVIVAL")
    emit("-" * 72)
    surv = pd.read_csv(snakemake.input.survival, sep="\t", dtype=str)
    emit(f"  {len(surv)} rows  ({snakemake.input.survival})")
    for column in (cols["sample"], cols["patient"], cols["event"],
                   cols["time"], cols["redaction"]):
        if column not in surv.columns:
            raise RuntimeError(
                f"column {column!r} is not in the survival file. Present: "
                f"{list(surv.columns)}."
            )

    n_redacted = int(surv[cols["redaction"]].notna().sum())
    surv = surv.loc[surv[cols["redaction"]].isna()]
    emit(f"  redacted upstream and excluded: {n_redacted}")

    surv["time_days"] = pd.to_numeric(surv[cols["time"]], errors="coerce")
    surv["event"] = pd.to_numeric(surv[cols["event"]], errors="coerce")
    n_no_time = int(surv["time_days"].isna().sum())
    n_nonpos = int((surv["time_days"] < expect["min_followup_days"]).sum())
    surv = surv.loc[
        surv["time_days"].notna()
        & (surv["time_days"] >= expect["min_followup_days"])
        & surv["event"].notna()
    ].copy()
    emit(f"  excluded for missing follow-up time: {n_no_time}")
    emit(f"  excluded for time < {expect['min_followup_days']} day: {n_nonpos}")
    emit("    (excluded and counted, never imputed — ADR 0025's shape)")

    # DAYS -> MONTHS. The factor is declared in config/tcga.yaml, not inlined.
    surv["months"] = surv["time_days"] / expect["days_per_month"]
    surv["event"] = surv["event"].astype(int)

    # --- sample -> patient --------------------------------------------------
    scores = scores.merge(
        surv[[cols["sample"], cols["patient"], "months", "event"]],
        left_on="sample", right_on=cols["sample"], how="inner",
    )
    scores = scores.rename(columns={cols["patient"]: "patient"})

    per_patient = scores.groupby(
        ["patient", "signature", "method"], as_index=False
    ).agg(
        score=("score", agg_method),
        n_samples=("sample", "nunique"),
        months=("months", "first"),
        event=("event", "first"),
    )
    n_patients = per_patient["patient"].nunique()
    n_collapsed = int((per_patient["n_samples"] > 1).sum())
    emit()
    emit(f"  {n_patients} patients after the join")
    emit(f"  {n_collapsed} patient-signature-method cells needed a collapse "
         f"(aggregation.method = {agg_method}, the P5-T2 rule)")

    if n_patients < expect["min_patients"]:
        raise RuntimeError(
            f"{n_patients} patients joined, below tcga.expect.min_patients = "
            f"{expect['min_patients']}. P5-T4's whole purpose is a PROPERLY "
            "POWERED comparator; an underpowered one would defeat it, so the "
            "phase stops rather than adapting."
        )
    emit(f"  [ok] cohort clears min_patients = {expect['min_patients']}")
    n_checks += 1

    # One row per patient per signature per method, or an arm silently shrinks.
    if per_patient.duplicated(["patient", "signature", "method"]).any():
        raise RuntimeError("a (patient, signature, method) cell appears twice.")
    expected_rows = n_patients * len(gene_sets) * 2
    if len(per_patient) != expected_rows:
        raise RuntimeError(
            f"{len(per_patient)} patient-level rows, expected {expected_rows} "
            f"({n_patients} patients x {len(gene_sets)} signatures x 2 methods)."
        )
    emit(f"  [ok] complete grid, {expected_rows} rows")
    n_checks += 1

    # A patient's outcome must not depend on which signature row it came from.
    spread = per_patient.groupby("patient")[["months", "event"]].nunique()
    if (spread > 1).any().any():
        raise RuntimeError(
            "a patient carries more than one follow-up time or event flag. The "
            "survival join duplicated rows."
        )
    emit("  [ok] each patient has exactly one follow-up time and event flag")
    n_checks += 1

    n_events = int(per_patient.groupby("patient")["event"].first().sum())
    emit(f"  events: {n_events} of {n_patients} "
         f"({n_events / n_patients:.1%}); censored {n_patients - n_events}")

    # --- fits ---------------------------------------------------------------
    # Machinery lifted from survival_km.py, including its EXACT wiring check.
    emit()
    emit("FITS — median split (primary) + continuous Cox (sensitivity)")
    emit("-" * 72)

    rows, curve_rows = [], []
    for signature in sorted(gene_sets):
        for method in ("ssgsea", "zscore"):
            cell = per_patient.loc[
                (per_patient["signature"] == signature)
                & (per_patient["method"] == method)
            ]
            values = cell["score"].to_numpy(dtype=float)
            months = cell["months"].to_numpy(dtype=float)
            events = cell["event"].to_numpy(dtype=int).astype(bool)

            # Same split as P5-T3: high = strictly above the median, tie side
            # declared. Comparability is the point (ADR 0027 §4).
            median = float(np.median(values))
            is_high = values > median
            n_high, n_low = int(is_high.sum()), int((~is_high).sum())
            if min(n_high, n_low) < min_arm_size:
                raise RuntimeError(
                    f"{signature}/{method}: split {n_high}/{n_low} below "
                    f"min_arm_size {min_arm_size} at n = {len(values)}."
                )

            y = Surv.from_arrays(event=events, time=months)
            chisq, p_logrank = compare_survival(y, is_high.astype(int))

            cox = PHReg(months, is_high.astype(float).reshape(-1, 1),
                        status=events.astype(int)).fit()
            loghr, se = float(cox.params[0]), float(cox.bse[0])
            ci = cox.conf_int()
            ci_low, ci_high = float(np.exp(ci[0, 0])), float(np.exp(ci[0, 1]))
            p_cox = float(cox.pvalues[0])

            # EXACT wiring check (survival_km.py): inverting the indicator must
            # negate the coefficient, or every reported DIRECTION is suspect.
            cox_inv = PHReg(months, (~is_high).astype(float).reshape(-1, 1),
                            status=events.astype(int)).fit()
            wiring_residual = abs(float(cox_inv.params[0]) + loghr)
            if wiring_residual > 1e-8:
                raise RuntimeError(
                    f"{signature}/{method}: inverting the group indicator "
                    f"changed |coefficient| by {wiring_residual:.3g}, not 0."
                )

            # SENSITIVITY: continuous Cox on the score. This is where ~500
            # patients actually buy power; a primary null beside a significant
            # continuous fit is a statement about the SPLIT, not the biology.
            cont = PHReg(months, values.reshape(-1, 1),
                         status=events.astype(int)).fit()
            cont_ci = cont.conf_int()
            ci_excludes_one = (ci_low > 1.0) or (ci_high < 1.0)

            rows.append({
                "cohort": tcga["cohort"], "signature": signature,
                "method": method, "is_primary_method": method == primary_method,
                "split": model_cfg["split"], "test": model_cfg["test"],
                "n_patients": int(len(values)),
                "n_events": int(events.sum()),
                "median_score": median, "n_high": n_high, "n_low": n_low,
                "n_events_high": int(events[is_high].sum()),
                "n_events_low": int(events[~is_high].sum()),
                "estimate": loghr, "std_error": se,
                "hazard_ratio": float(np.exp(loghr)),
                "ci_low": ci_low, "ci_high": ci_high,
                "p_raw": float(p_logrank),
                "p_raw_logrank": float(p_logrank), "p_raw_cox": p_cox,
                "logrank_chisq": float(chisq),
                "p_agreement_delta": abs(float(p_logrank) - p_cox),
                "wiring_residual": wiring_residual,
                "ci_excludes_one": bool(ci_excludes_one),
                "significance_disagreement": bool(
                    (float(p_logrank) < alpha) != bool(ci_excludes_one)
                ),
                "continuous_hr_per_unit": float(np.exp(cont.params[0])),
                "continuous_ci_low": float(np.exp(cont_ci[0, 0])),
                "continuous_ci_high": float(np.exp(cont_ci[0, 1])),
                "continuous_p_raw": float(cont.pvalues[0]),
                "continuous_ci_excludes_one": bool(
                    float(np.exp(cont_ci[0, 0])) > 1.0
                    or float(np.exp(cont_ci[0, 1])) < 1.0
                ),
                # THE POINT OF CARRYING THE SENSITIVITY (ADR 0027 §4). Where the
                # median split and the continuous fit disagree, the disagreement
                # is a statement about THE SPLIT, not about the biology: a
                # median split discards information, and at this n the
                # continuous fit is the better-powered of the two. Neither is
                # the tiebreak, and this flag exists so the divergence is
                # machine-readable rather than only visible in the log.
                "split_vs_continuous_disagreement": bool(
                    bool(ci_excludes_one)
                    != bool(float(np.exp(cont_ci[0, 0])) > 1.0
                            or float(np.exp(cont_ci[0, 1])) < 1.0)
                ),
                "status": "modelled",
                "geomx_comparator": (
                    "none — not assessable in the GeoMx lung arm"
                    if signature == "myeloid_m1" else "GeoMx lung arm (P5-T3)"
                ),
                "note": "raw uncorrected p; P5-T5 corrects once over 14 (ADR 0026)",
            })

            for label, mask in (("high", is_high), ("low", ~is_high)):
                t_km, s_km = kaplan_meier_estimator(events[mask], months[mask])
                for time, surv_p in zip(t_km, s_km):
                    curve_rows.append({
                        "cohort": tcga["cohort"], "signature": signature,
                        "method": method, "group": label,
                        "time_months": float(time), "survival": float(surv_p),
                    })

    models = pd.DataFrame(rows)
    curves = pd.DataFrame(curve_rows)

    worst = models["p_agreement_delta"].max()
    if worst > P_AGREEMENT_TOL:
        raise RuntimeError(
            f"log-rank and Cox p differ by up to {worst:.4g} "
            f"(tolerance {P_AGREEMENT_TOL})."
        )
    emit(f"  [ok] log-rank and Cox p agree in all {len(models)} cells "
         f"(worst {worst:.4g})")
    n_checks += 1

    for name, frame in (("models", models), ("presence", presence),
                        ("curves", curves), ("patient scores", per_patient)):
        offending = [c for c in frame.columns if c.startswith("q_") or c == "fdr_bh"]
        if offending:
            raise RuntimeError(
                f"the {name} table carries {offending}. P5-T4 does NOT correct "
                "(ADR 0026 §3): the family is 8 + 6 = 14 and P5-T5 corrects "
                "once over all of it."
            )
    emit("  [ok] no q-column in any output — P5-T5 corrects once over 14")
    n_checks += 1

    # --- report -------------------------------------------------------------
    show = models.loc[models["method"] == primary_method].sort_values("signature")
    emit()
    emit(f"RESULTS — {primary_method} (primary), median split, RAW p")
    emit("-" * 72)
    for _, row in show.iterrows():
        flag = " *" if row["ci_excludes_one"] else "  "
        emit(f"  {row['signature']:22s} n={int(row['n_patients'])} "
             f"({int(row['n_high'])}/{int(row['n_low'])}) "
             f"ev={int(row['n_events'])}  "
             f"HR {row['hazard_ratio']:5.3f} "
             f"[{row['ci_low']:5.3f}, {row['ci_high']:5.3f}]  "
             f"p_raw {row['p_raw']:.4f}{flag}")
    emit("  (* = the 95% interval excludes 1)")

    emit()
    emit("  SENSITIVITY — continuous Cox per unit score (where the n buys power):")
    for _, row in show.iterrows():
        emit(f"    {row['signature']:22s} HR/unit "
             f"{row['continuous_hr_per_unit']:.4g} "
             f"[{row['continuous_ci_low']:.4g}, {row['continuous_ci_high']:.4g}]"
             f"  p_raw {row['continuous_p_raw']:.4f}")

    disagree = show.loc[show["significance_disagreement"]]
    if not disagree.empty:
        emit()
        emit(f"  {len(disagree)} cell(s) where the log-rank and the interval")
        emit("  DISAGREE. WHERE THEY DISAGREE THE INTERVAL GOVERNS:")
        for _, row in disagree.iterrows():
            emit(f"    {row['signature']:22s} p {row['p_raw']:.4f} vs "
                 f"CI [{row['ci_low']:.3f}, {row['ci_high']:.3f}]")

    split_cont = show.loc[show["split_vs_continuous_disagreement"]]
    if not split_cont.empty:
        emit()
        emit(f"  {len(split_cont)} of {len(show)} cells: THE SPLIT AND THE")
        emit("  CONTINUOUS FIT DISAGREE about whether the interval excludes 1.")
        for _, row in split_cont.iterrows():
            emit(f"    {row['signature']:22s} "
                 f"split HR {row['hazard_ratio']:.3f} "
                 f"[{row['ci_low']:.3f}, {row['ci_high']:.3f}]  vs  "
                 f"continuous HR/unit {row['continuous_hr_per_unit']:.4g} "
                 f"[{row['continuous_ci_low']:.4g}, "
                 f"{row['continuous_ci_high']:.4g}]")
        emit()
        emit("  A DISAGREEMENT HERE IS A STATEMENT ABOUT THE SPLIT, NOT ABOUT")
        emit("  THE BIOLOGY (ADR 0027 §4). A median split discards information;")
        emit("  at this n the continuous fit is the better powered of the two.")
        emit("  Neither is the tiebreak, and the split stays primary because")
        emit("  COMPARABILITY WITH P5-T3 is what this task exists to provide.")

    n_excl = int(show["ci_excludes_one"].sum())
    emit()
    emit(f"  {n_excl} of {len(show)} primary cells have an interval excluding 1;")
    emit(f"  {int((show['p_raw'] < alpha).sum())} have RAW p < {alpha}.")
    emit("  ALL RAW. P5-T5 applies BH once over the family of 14 (ADR 0026 §2),")
    emit("  and only that q is reportable.")
    emit()
    emit("  INTERPRETATION IS BOUNDED (ADR 0027 §5): TCGA is WHOLE BULK TUMOUR")
    emit("  and GeoMx TIME-L is the PanCK-negative SEGMENT, so this is not the")
    emit("  same measurement — a disagreement is not necessarily about biology.")
    emit("  And there is NO BRAIN COMPARATOR: this bounds the LUNG null only.")

    # --- figure -------------------------------------------------------------
    order = sorted(gene_sets)
    fig = plt.figure(figsize=(11.2, 13.4))
    gs = fig.add_gridspec(
        3, 2, left=0.075, right=0.975, top=0.888, bottom=0.055,
        hspace=0.36, wspace=0.18,
    )
    colours = {"high": "#b2182b", "low": "#2166ac"}
    for idx, signature in enumerate(order):
        ax = fig.add_subplot(gs[idx // 2, idx % 2])
        row = show.loc[show["signature"] == signature].iloc[0]
        sub = curves.loc[
            (curves["signature"] == signature)
            & (curves["method"] == primary_method)
        ]
        for label in ("high", "low"):
            g = sub.loc[sub["group"] == label].sort_values("time_months")
            t_km = np.concatenate([[0.0], g["time_months"].to_numpy()])
            s_km = np.concatenate([[1.0], g["survival"].to_numpy()])
            n = int(row["n_high"] if label == "high" else row["n_low"])
            ax.step(t_km, s_km, where="post", color=colours[label],
                    linewidth=1.7, label=f"{label} (n={n})")
        ax.set_ylim(-0.30, 1.06)
        ax.set_xlim(left=0)
        ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=7.4, loc="upper right", frameon=False,
                  handlelength=1.3, borderaxespad=0.2)
        title = signature.replace("_", " ")
        if signature == "myeloid_m1":
            title += "  (TCGA-only — not assessable in GeoMx)"
        ax.set_title(title, fontsize=9.6, weight="bold", pad=6)
        # Opaque bbox: on P5-T3 a curve crossed the decimal point of an interval
        # and it read as a different number. Phase 4 lesson 1.
        note = ""
        if row["significance_disagreement"]:
            note = "\nraw p < alpha but CI spans 1 — interval governs"
        ax.text(
            0.02, 0.05,
            f"HR {row['hazard_ratio']:.2f} "
            f"[{row['ci_low']:.2f}, {row['ci_high']:.2f}]\n"
            f"raw p = {row['p_raw']:.4f}  (n={int(row['n_patients'])}, "
            f"{int(row['n_events'])} events)" + note,
            transform=ax.transAxes, fontsize=7.4, va="bottom",
            color="#8c2d04" if row["significance_disagreement"] else "#222222",
            bbox=dict(facecolor="white", alpha=0.92, edgecolor="none",
                      boxstyle="square,pad=0.28"),
            zorder=5,
        )
        if idx % 2 == 0:
            ax.set_ylabel("overall survival", fontsize=8.6)
        if idx // 2 == 2:
            ax.set_xlabel("months from diagnosis", fontsize=8.4)

    # THREE SHORT LINES, not two long ones: on the first render the second line
    # was ~150 characters and ran off BOTH edges of the canvas, losing "media"
    # at the left and "exists" at the right. A suptitle is not wrapped or
    # clipped by matplotlib — it is simply drawn past the figure. Phase 4
    # lesson 1, and `top` below is set to leave room for all three.
    fig.suptitle(
        f"P5-T4  TCGA LUAD external check — {tcga['cohort']}, "
        f"n = {n_patients} patients, {n_events} events\n"
        "median split, RAW uncorrected p — P5-T5 corrects once over 14\n"
        "WHOLE BULK TUMOUR, not the PanCK-negative segment GeoMx measures"
        "  ·  no brain comparator exists",
        fontsize=9.4, y=0.985, va="top",
    )
    fig.savefig(snakemake.output.figure, dpi=200)
    plt.close(fig)

    # --- write --------------------------------------------------------------
    per_patient.to_csv(snakemake.output.scores, sep="\t", index=False,
                       float_format=FMT)
    presence.to_csv(snakemake.output.presence, sep="\t", index=False,
                    float_format=FMT)
    models.to_csv(snakemake.output.models, sep="\t", index=False,
                  float_format=FMT)
    curves.to_csv(snakemake.output.curves, sep="\t", index=False,
                  float_format=FMT)

    emit()
    emit(f"assertions: {n_checks} of {n_checks} passed")

    summary = {
        "task": "P5-T4",
        "phase": 5,
        "pre_registered": "ADR 0027",
        "cohort": tcga["cohort"],
        "provider": tcga["provider"],
        "source_note": (
            "UCSC Xena, not cBioPortal: the datahub bulk files return 403/404 "
            "and only the live REST API answers, which cannot be digest-pinned "
            "under ADR 0005. ADR 0024 §1 pre-authorised this fallback."
        ),
        "endpoint": {"event": cols["event"], "time": cols["time"],
                     "units_in": "days", "days_per_month": expect["days_per_month"],
                     "source": tcga["endpoint_source"]},
        "cohort_construction": {
            "sample_type_kept": expect["sample_type_code"],
            "sample_types_dropped": dropped_types,
            "n_redacted_excluded": n_redacted,
            "n_missing_time_excluded": n_no_time,
            "n_nonpositive_time_excluded": n_nonpos,
            "n_patients": int(n_patients),
            "n_events": int(n_events),
            "aggregation": agg_method,
            "n_cells_collapsed": n_collapsed,
        },
        "gene_presence": {
            "basis": "presence in the expression matrix, NOT detection above background",
            "n_scoreable": len(gene_sets),
            "all_genes_present": bool((presence["presence"] == 1.0).all()),
        },
        "multiplicity": {
            "corrected_here": False,
            "tcga_contribution": len(gene_sets),
            "family_total": 8 + len(gene_sets),
            "adr": "ADR 0026 §2, ADR 0027 §2",
            "note": (
                "TCGA contributes 6 — all signatures, because the Phase 2 "
                "detection floor is a property of the GeoMx assay and does not "
                "transfer. Family = 8 (P5-T3) + 6 = 14, corrected once at P5-T5."
            ),
        },
        "n_ci_excluding_one": n_excl,
        "split_vs_continuous": {
            "n_disagreeing": int(len(split_cont)),
            "cells": [
                {"signature": r["signature"],
                 "split_hr": float(r["hazard_ratio"]),
                 "split_ci": [float(r["ci_low"]), float(r["ci_high"])],
                 "continuous_hr_per_unit": float(r["continuous_hr_per_unit"]),
                 "continuous_ci": [float(r["continuous_ci_low"]),
                                   float(r["continuous_ci_high"])]}
                for _, r in split_cont.iterrows()
            ],
            "rule": (
                "A disagreement between the median split and the continuous fit "
                "is a statement about THE SPLIT, not about the biology "
                "(ADR 0027 §4). A median split discards information and at this "
                "n the continuous fit is the better powered; neither is the "
                "tiebreak, and the split stays primary because comparability "
                "with P5-T3 is what this task exists to provide."
            ),
        },
        "n_raw_p_below_alpha": int((show["p_raw"] < alpha).sum()),
        "alpha": alpha,
        "n_assertions": n_checks,
        "limitations": {
            "not_the_same_measurement": (
                "GeoMx TIME-L is the PanCK-negative SEGMENT of an ROI sited in "
                "a CD45-rich region (and not a CD45-sorted population, per Q2); "
                "TCGA is WHOLE BULK TUMOUR. A disagreement between cohorts is "
                "not necessarily a disagreement about biology, and an agreement "
                "is not necessarily a replication."
            ),
            "no_brain_comparator": (
                "TCGA LUAD is primary lung. The brain arm — TIME-B n = 8, the "
                "binding constraint on the project — gets NO external check "
                "from this task. P5-T4 contextualises the lung null only."
            ),
        },
        "reporting_rule": (
            "EVERY p HERE IS RAW AND UNCORRECTED (ADR 0026 §3). Every p carries "
            "n, hazard ratio and a confidence interval (hard constraint 7); "
            "where the log-rank and the interval disagree, THE INTERVAL GOVERNS. "
            "A null in TCGA at this n means something a null in this study does "
            "not — that is the point of the task — but it is bounded by the two "
            "limitations above. myeloid_m1 is TCGA-only and has no GeoMx "
            "comparator."
        ),
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")

    emit(f"wrote {snakemake.output.scores}, {snakemake.output.presence}, "
         f"{snakemake.output.models}, {snakemake.output.curves}, "
         f"{snakemake.output.figure}, {snakemake.output.summary}")
