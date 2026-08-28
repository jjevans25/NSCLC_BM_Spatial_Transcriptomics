"""P2-T3b — refit the Gate 2 signature in statsmodels and compare to lme4.

Owner task: P2-T3. Driven by rule p2t3_model_crosscheck.

PROJECT_PLAN §6 asks for one signature cross-checked in `statsmodels.MixedLM`.
This is an IMPLEMENTATION check, not a statistical one: it asks whether two
independent implementations of the same REML fit land on the same number, so
that a Gate 2 verdict does not rest on one library's handling of a 23-row
design with 18 patients.

`antigen_presentation` is the set checked because it is the set Gate 2 turns
on.

What is comparable and what is not:

  * **Point estimate and standard error are comparable.** Both fit the same
    REML mixed model with a patient random intercept, so they should agree to
    numerical tolerance. A disagreement is a bug in one of them.

  * **Degrees of freedom and p-values are NOT comparable, and are not
    compared.** lmerTest reports Satterthwaite-approximated denominator df;
    statsmodels reports a Wald z with no finite-sample df correction at all. At
    n = 23 the two give visibly different tail probabilities *by construction*.
    Printing them side by side as though a small difference were reassuring —
    or a large one alarming — would be reading noise. The p-value of record is
    lmerTest's.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    target = snakemake.params.crosscheck_signature
    tolerance = snakemake.params.tolerance

    scores = pd.read_csv(snakemake.input.scores, sep="\t")
    lme4 = pd.read_csv(snakemake.input.models, sep="\t")

    emit(f"cross-checking signature: {target}")
    emit(f"tolerance on |estimate| difference: {tolerance}")
    emit()

    comparisons = []
    for method in sorted(scores["method"].unique()):
        d = scores.loc[
            (scores["signature"] == target) & (scores["method"] == method)
        ].copy()

        # lung first, so the coefficient is brain minus lung — the same
        # orientation the R script asserts.
        d["site"] = pd.Categorical(d["site"], categories=["lung", "brain"])

        model = smf.mixedlm("score ~ site", data=d, groups=d["patient_id"])
        fit = model.fit(reml=True, method="lbfgs")

        term = [t for t in fit.params.index if t.startswith("site")]
        if len(term) != 1:
            raise RuntimeError(f"expected one site term, found {term}")
        term = term[0]

        sm_est = float(fit.params[term])
        sm_se = float(fit.bse[term])

        row = lme4.loc[
            (lme4["signature"] == target)
            & (lme4["method"] == method)
            & (lme4["model"] == "primary")
        ]
        if len(row) != 1:
            raise RuntimeError(
                f"expected one lme4 primary row for {target}/{method}, got {len(row)}"
            )
        r_est = float(row["estimate"].iloc[0])
        r_se = float(row["std_error"].iloc[0])

        d_est = abs(sm_est - r_est)
        d_se = abs(sm_se - r_se)
        agree = bool(d_est <= tolerance)

        emit(f"  method = {method}   term = {term}")
        emit(f"    lme4/lmerTest   estimate {r_est:+.6f}   SE {r_se:.6f}")
        emit(f"    statsmodels     estimate {sm_est:+.6f}   SE {sm_se:.6f}")
        emit(f"    |d estimate| = {d_est:.3e}   |d SE| = {d_se:.3e}   agree={agree}")
        emit(
            "    df/p not compared: lmerTest uses Satterthwaite, statsmodels a "
            "Wald z with no df correction. The p-value of record is lmerTest's."
        )
        emit()

        comparisons.append(
            {
                "method": method,
                "lme4_estimate": r_est,
                "lme4_std_error": r_se,
                "statsmodels_estimate": sm_est,
                "statsmodels_std_error": sm_se,
                "abs_diff_estimate": d_est,
                "abs_diff_std_error": d_se,
                "agrees_within_tolerance": agree,
                "converged": bool(fit.converged),
            }
        )

    all_agree = all(c["agrees_within_tolerance"] for c in comparisons)
    if not all_agree:
        emit("CROSS-CHECK FAILED — the two implementations disagree.")
        raise RuntimeError(
            f"lme4 and statsmodels disagree on {target} beyond {tolerance}. "
            "This is an implementation bug, not a finding; do not proceed to "
            "Gate 2 on it."
        )

    emit(f"cross-check PASSED for {target} across {len(comparisons)} methods")

    with open(snakemake.output.crosscheck, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "signature": target,
                "tolerance": tolerance,
                "contrast": "brain_minus_lung",
                "all_agree": all_agree,
                "note": (
                    "Estimates and SEs are compared; df and p-values are not, "
                    "because lmerTest uses Satterthwaite and statsmodels a Wald "
                    "z with no finite-sample correction."
                ),
                "comparisons": comparisons,
            },
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")
