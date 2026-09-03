"""P3-T3b — refit one checkpoint gene in statsmodels and compare to lme4.

Owner task: P3-T3. Driven by rule p3t3_checkpoint_crosscheck.

A NEW FILE rather than an edit to crosscheck_mixedlm.py, deliberately: editing
that file would fire the code rerun-trigger on p2t3_model_crosscheck and re-run
Phase 2 for no numerical gain. Same reason P3-T4's paired check gets its own
file rather than editing paired_check.py.

This is an IMPLEMENTATION check, not a statistical one. Two independent REML
implementations of the same model should land on the same number; a
disagreement is a bug in one of them, and Phase 3's estimates should not rest
on one library's handling of a 45-row design.

It also closes the R/Python boundary END-TO-END. fit_checkpoint_models.R
asserts only STRUCTURE on read, because R_strtod is not correctly rounded and
byte-exact agreement is unattainable Python -> R (ADR 0013 §3). The numeric
guarantee comes from here: this script reads the same TSV with Python's
correctly-rounded parser, refits, and requires agreement with lme4. That check
is not circular, unlike any comparison R can make against its own parse.

`float_precision="round_trip"` is MANDATORY on that read. pandas' default C
parser is not correctly rounded either (ADR 0013 postscript) — reading P3-T2's
own %.17g output back reported 32 of 720 values unequal when the file was
exact.

CD276 is the gene checked, pre-registered in config.yaml as
`checkpoints.model.crosscheck_gene`, because it is the best-detected gene on the
panel (105 of 120 AOIs) and therefore the best-conditioned fit. A cross-check
should test the implementation, not the optimiser's behaviour at the detection
floor.

WHAT IS COMPARED: the point estimate, and only the point estimate. This matches
p2t3_model_crosscheck, which gates on `d_est <= tolerance` and reports the rest.

WHAT IS NOT, AND WHY — all three are reported, none is gated:

  * **Standard error.** The two libraries compute it differently, and this was
    measured here rather than assumed. For CD276 the estimates agree to 7e-9
    (brain) and 1e-6 (lung) while the SEs differ by 4.7e-3 and 2.6e-2 — the lung
    gap is 20% of the SE. statsmodels reaches the same variance components under
    lbfgs, powell and bfgs (group_var 0.1485 lung, 0.1643 brain, all converged),
    so this is not an optimiser artefact: it is two different estimators of the
    same quantity, agreeing only when the variance components are sharply
    determined. They are not sharply determined here — 17 of 30 lung patients
    and 19 of 27 brain patients contribute a single AOI.

    **P3-T4 confirmed the mechanism.** Run on the site contrast, the SE gap is
    17.1% in the immune compartment (23 AOIs, 5 paired patients) and **0.4% in
    the tumour compartment** (57 AOIs, 23 paired). The estimates agree to 4.7e-6
    and 2.4e-6 in both. Where the random intercept is well identified the two
    libraries converge; where it is not, they do not. That is the prediction,
    tested on a contrast it was not derived from.

    **The consequence matters and is not cosmetic: the CI of record is
    lmerTest's**, and a reader who recomputes it in statsmodels will get a wider
    one. Reported in the JSON so that is visible rather than surprising.

  * **Degrees of freedom and p-values.** lmerTest reports Satterthwaite df;
    statsmodels a Wald z with no finite-sample correction. At these n the tail
    probabilities differ by construction. The p-value of record is lmerTest's.

Gating on the SE would fail this check on a property neither library claims to
share, which is a check that tests the wrong thing — and then gets relaxed,
which is worse than never having asserted it.
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

    target = snakemake.params.crosscheck_gene
    tolerance = snakemake.params.tolerance
    contrast_var = snakemake.params.contrast_var
    reference = snakemake.params.reference_level
    test_level = snakemake.params.test_level
    strata_var = snakemake.params.strata_var
    task = snakemake.params.task

    expression = pd.read_csv(
        snakemake.input.expression, sep="\t", float_precision="round_trip"
    )
    lme4 = pd.read_csv(
        snakemake.input.models, sep="\t", float_precision="round_trip"
    )

    emit(f"{task} — cross-checking gene: {target}")
    emit(f"contrast: {contrast_var} '{test_level}' minus '{reference}'")
    emit(f"tolerance on |estimate| difference: {tolerance}")
    emit("std_error, df and p-values are REPORTED but NOT gated — the two")
    emit("libraries use different SE estimators; see the module docstring.")
    emit()

    comparisons = []
    for stratum in sorted(expression[strata_var].unique()):
        d = expression.loc[
            (expression["gene"] == target) & (expression[strata_var] == stratum)
        ].copy()

        # Reference level first, so the coefficient is test minus reference —
        # the same orientation fit_checkpoint_models.R asserts. Getting this
        # backwards would compare a number against its own negation.
        d[contrast_var] = pd.Categorical(
            d[contrast_var], categories=[reference, test_level], ordered=True
        )
        if d[contrast_var].isna().any():
            raise RuntimeError(
                f"{stratum}: {contrast_var} holds levels outside "
                f"[{reference}, {test_level}]"
            )

        model = smf.mixedlm(
            f"expression ~ C({contrast_var})", d, groups=d["patient_id"]
        )
        fit = model.fit(reml=True, method="lbfgs")

        term = [t for t in fit.params.index if t.startswith(f"C({contrast_var})")]
        if len(term) != 1:
            raise RuntimeError(f"{stratum}: expected one contrast term, got {term}")
        term = term[0]
        sm_estimate = float(fit.params[term])
        sm_se = float(fit.bse[term])

        row = lme4.loc[
            (lme4["gene"] == target)
            & (lme4["stratum"] == stratum)
            & (lme4["model"] == "primary")
        ]
        if len(row) != 1:
            # A gene the primary declined to model has no estimate to compare.
            # That is not a failure of the cross-check — but the pre-registered
            # crosscheck_gene must be modelled in every stratum, or the check
            # silently tests nothing.
            raise RuntimeError(
                f"{stratum}: expected exactly one primary row for {target}, "
                f"got {len(row)}. If {target} is not assessable everywhere it "
                "is the wrong choice for checkpoints.model.crosscheck_gene — "
                "that is a stop-and-ask, not something to work around here."
            )
        r_estimate = float(row["estimate"].iloc[0])
        r_se = float(row["std_error"].iloc[0])

        d_est = abs(sm_estimate - r_estimate)
        d_se = abs(sm_se - r_se)
        # Estimate only, as p2t3_model_crosscheck does. See the docstring: the
        # SE is a different estimator in each library, not the same number
        # computed twice.
        agrees = bool(d_est <= tolerance)

        emit(f"  {stratum}:")
        emit(f"    lme4        estimate {r_estimate:+.9f}  se {r_se:.9f}")
        emit(f"    statsmodels estimate {sm_estimate:+.9f}  se {sm_se:.9f}")
        emit(f"    |diff|      estimate {d_est:.3e}  {'ok' if agrees else 'MISMATCH'}"
             f"   (se {d_se:.3e}, {100 * d_se / r_se:.1f}% — reported, not gated)")
        emit(f"    n = {len(d)} AOIs, {d['patient_id'].nunique()} patients")

        comparisons.append(
            {
                "stratum": stratum,
                "gene": target,
                "lme4_estimate": r_estimate,
                "lme4_std_error": r_se,
                "statsmodels_estimate": sm_estimate,
                "statsmodels_std_error": sm_se,
                "abs_diff_estimate": d_est,
                "abs_diff_std_error": d_se,
                "rel_diff_std_error": d_se / r_se,
                "n_obs": int(len(d)),
                "n_patient": int(d["patient_id"].nunique()),
                "agrees": agrees,
            }
        )

    failures = [c for c in comparisons if not c["agrees"]]
    emit()
    if failures:
        raise RuntimeError(
            f"lme4 and statsmodels disagree beyond {tolerance} on "
            f"{[c['stratum'] for c in failures]}. Two independent REML "
            "implementations of the same model do not agree, so one of them is "
            "wrong. Resolve before any Phase 3 number is reported."
        )
    worst = max(c["abs_diff_estimate"] for c in comparisons)
    emit(f"all {len(comparisons)} strata agree; worst |estimate| diff {worst:.3e}")
    emit()
    emit("This also closes the R/Python boundary end-to-end: R asserted only")
    emit("structure on read (R_strtod is not correctly rounded, ADR 0013 §3),")
    emit("and this refit from the same TSV through a correctly-rounded parser")
    emit("is the numeric guarantee.")

    with open(snakemake.output.crosscheck, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "task": task,
                "gene": target,
                "contrast_var": contrast_var,
                "contrast": f"{test_level}_minus_{reference}",
                "tolerance": tolerance,
                "compared": ["estimate"],
                "not_compared": ["std_error", "df", "p_value"],
                "not_compared_reason": (
                    "Different estimators, not the same number computed twice. "
                    "statsmodels' MixedLM fixed-effect SEs and lmerTest's agree "
                    "only when the variance components are sharply determined, "
                    "which they are not here (17/30 lung and 19/27 brain "
                    "patients contribute one AOI); statsmodels reaches the same "
                    "variance components under lbfgs, powell and bfgs, so this "
                    "is not an optimiser artefact. lmerTest reports "
                    "Satterthwaite df, statsmodels a Wald z with no "
                    "finite-sample correction. The SE, CI and p-value of record "
                    "are lmerTest's; recomputing the CI in statsmodels gives a "
                    "wider one."
                ),
                "worst_rel_diff_std_error": max(
                    c["rel_diff_std_error"] for c in comparisons
                ),
                "worst_abs_diff_estimate": worst,
                "all_agree": True,
                "comparisons": comparisons,
            },
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")
