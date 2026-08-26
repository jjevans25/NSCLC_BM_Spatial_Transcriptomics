"""P1-T3 — quantify the variance split. This is Gate 1.

Owner task: P1-T3. Driven by rule p1t3_variance_partition.

Aim A2 asks what drives variance here: patient, site, or compartment. P1-T1/T2
made the answer visible; this puts a number on it, and Gate 1 turns on that
number. The deliverable is one sentence — "X% of variance is
<component>-attributable" — with enough around it that the sentence is
defensible.

Choices worth stating, because each could reasonably go the other way:

  * **Variance components, not eta-squared.** On PC1 the crude numbers are
    compartment 0.754, site 0.163, patient 0.498 — but `patient_id` has 42
    levels across 120 AOIs, whose eta-squared null expectation is already
    (42-1)/(120-1) = 0.345. A naive read of "patient explains half the
    variance" is close to an artefact of the level count. A crossed random-
    effects fit shrinks a factor toward zero in proportion to how little each
    of its levels contributes, so it does not reward level count the way a
    fixed-effect sum of squares does. Both the analytic eta-squared null and an
    empirical permutation null are reported beside every share, so the reader
    can see the floor rather than take the penalisation on trust.

  * **`dsp_run` is a term.** The plan's four-way split is patient / site /
    compartment / residual. P1-T2 found PC1 shifting positive in run B inside
    every compartment that spans both runs (L +19.8, LB +9.7, TIME-L +4.7,
    TIME-B +7.6), so omitting batch pushes that variance into `compartment` —
    the very number the gate turns on. It is only *estimable* in those four
    compartments; `mLN`, `TBME` and `BC` sit wholly inside one run (Q3), so the
    partition is conditional on that. It is also a **two-level** effect with
    only 1 of 42 patients spanning both runs — near-nested in patient, and
    poorly estimated for it. Sensitivity S2 refits without it and reports how
    much actually moves.

  * **All 18,694 genes, HVGs flagged.** Restricting to the 2000 HVGs would let
    feature selection choose the answer; the highest-variance genes here are
    immunoglobulins, which carry a specific and unrepresentative structure. The
    HVG flag is carried as a column so the HVG-only medians come out of the
    same table.

  * **`BC` controls are in the primary fit.** All 120 AOIs, matching the design
    table and the T1/T2 PCA. `BC` is the only `normal_control` compartment and
    is brain-only and run-B-only, so sensitivity S3 refits without it.

  * **`VCSpec` is built explicitly rather than via `vc_formula=`.** The formula
    route returns `vcomp` in **alphabetical** name order regardless of the dict
    order passed in, which silently mislabels every component. Building the
    spec directly keeps the reported order the configured order.

Non-converged genes stay in the table with `converged = False` and are excluded
from every median. All stochastic steps take config["seed"] explicitly.
"""

import json
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from statsmodels.regression.mixed_linear_model import MixedLM, VCSpec

COMPONENT_COLOUR = {
    "patient": "#7B4B94",
    "site": "#C1666B",
    "compartment": "#B0763F",
    "dsp_run": "#4F6D7A",
    "aoi_code": "#B0763F",
    "residual": "#8C8C8C",
}

# Factor column -> the short name used in every output. Keeps `patient_id` from
# becoming a column called `var_patient_id`.
SHORT_NAME = {"patient_id": "patient", "aoi_code": "aoi_code"}


def short(factor):
    return SHORT_NAME.get(factor, factor)


def build_vcspec(obs, factors):
    """A crossed-random-effects VCSpec over `factors`, in the order given.

    statsmodels expresses crossed variance components as a single group
    containing every observation, with one dummy-coded design matrix per
    component. Level order inside a component is irrelevant to the variance;
    component order is not, because `vcomp` is returned positionally.
    """
    names, colnames, mats = [], [], []
    for factor in factors:
        dummies = pd.get_dummies(obs[factor].astype(str), drop_first=False).astype(float)
        names.append(short(factor))
        colnames.append(list(dummies.columns))
        mats.append([dummies.to_numpy()])
    return VCSpec(names, colnames, mats)


def fit_shares(y, vcs, n_obs):
    """Fractional variance per component plus residual, for one response.

    Returns (shares, converged, method). `shares` has one entry per component
    followed by the residual, and sums to 1 by construction.

    statsmodels' default `lbfgs` fails to converge on roughly a third of genes
    here — the likelihood is ridged, not the optimiser broken: `powell`
    converges everywhere and, on genes where both converge, agrees to ~3e-4 at
    the median with a median log-likelihood gap of ~1e-7. So `lbfgs` stays the
    primary (it is the tested default and uses gradients), `powell` is the
    fallback, and the method used is recorded per gene rather than hidden.
    """
    n_comp = len(vcs.names)
    failed = (np.full(n_comp + 1, np.nan), False, "none")
    if not np.isfinite(y).all() or np.ptp(y) == 0:
        return failed
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit, method = None, "none"
        for candidate, kwargs in (("lbfgs", {}), ("powell", {"maxiter": 2000})):
            try:
                attempt = MixedLM(
                    y, np.ones((n_obs, 1)), groups=np.ones(n_obs), exog_vc=vcs
                ).fit(reml=True, method=candidate, **kwargs)
            except (np.linalg.LinAlgError, ValueError):
                continue
            fit, method = attempt, candidate
            if attempt.converged:
                break
    if fit is None:
        return failed
    parts = np.append(np.asarray(fit.vcomp, dtype=float), float(fit.scale))
    total = parts.sum()
    if not np.isfinite(total) or total <= 0:
        return failed
    return parts / total, bool(fit.converged), method


def partition(matrix, obs, factors, gene_names):
    """Per-gene variance partition. `matrix` is AOIs x genes."""
    vcs = build_vcspec(obs, factors)
    n_obs = matrix.shape[0]
    n_genes = matrix.shape[1]
    shares = np.empty((n_genes, len(factors) + 1), dtype=float)
    converged = np.empty(n_genes, dtype=bool)
    methods = np.empty(n_genes, dtype=object)
    for j in range(n_genes):
        shares[j], converged[j], methods[j] = fit_shares(matrix[:, j], vcs, n_obs)
    frame = pd.DataFrame(
        shares,
        index=pd.Index(gene_names, name="gene"),
        columns=[f"var_{short(f)}" for f in factors] + ["var_residual"],
    )
    frame["converged"] = converged
    frame["method"] = methods
    return frame


def medians(frame):
    """Median share per component over converged genes only."""
    ok = frame.loc[frame["converged"]]
    cols = [c for c in frame.columns if c.startswith("var_")]
    return {c.removeprefix("var_"): round(float(ok[c].median()), 6) for c in cols}


def naive_eta2(matrix, obs, factor):
    """Fixed-effect eta-squared per gene, vectorised — the number *not* to use.

    Reported only so the penalisation the mixed model applies is visible rather
    than asserted.
    """
    codes = pd.Categorical(obs[factor].astype(str))
    centred = matrix - matrix.mean(axis=0, keepdims=True)
    ss_total = (centred ** 2).sum(axis=0)
    ss_between = np.zeros(matrix.shape[1])
    for level in range(len(codes.categories)):
        mask = codes.codes == level
        if not mask.any():
            continue
        ss_between += mask.sum() * (centred[mask].mean(axis=0) ** 2)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(ss_total > 0, ss_between / ss_total, np.nan)


log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    seed = snakemake.params.seed
    factors = list(snakemake.params.factors)
    pvca_min_variance = snakemake.params.pvca_min_variance
    n_permutations = snakemake.params.n_permutations
    n_permutation_genes = snakemake.params.n_permutation_genes
    n_top = snakemake.params.n_top_genes
    flavour = snakemake.params.hvg_flavour
    sc.settings.verbosity = 0
    emit(f"seed: {seed} (config['seed'], passed explicitly)")
    emit(f"factors: {factors}")

    adata = sc.read_h5ad(snakemake.input.h5ad)
    obs = adata.obs.copy()
    X = np.asarray(adata.X, dtype=float)
    emit(f"loaded {adata.n_obs} AOIs x {adata.n_vars} genes")

    for factor in factors:
        levels = obs[factor].astype(str).nunique()
        emit(f"  {factor}: {levels} levels across {adata.n_obs} AOIs "
             f"(eta-squared null {(levels - 1) / (adata.n_obs - 1):.3f})")

    # `dsp_run` is nearly nested in patient, which is why S2 exists.
    spanning = int((obs.groupby("patient_id", observed=True)["dsp_run"].nunique() > 1).sum())
    emit(f"  patients spanning both DSP runs: {spanning} of "
         f"{obs['patient_id'].nunique()} — dsp_run is near-nested in patient")

    # HVG flag, recomputed exactly as P1-T1 does it so the two tasks agree on
    # which genes are "the HVGs".
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor=flavour)
    hvg = adata.var["highly_variable"].to_numpy()
    emit(f"HVGs: {int(hvg.sum())} of {adata.n_vars} (flavour={flavour}, target {n_top})")

    # ------------------------------------------------------------ primary fit
    emit("\n--- primary fit: all genes, all 120 AOIs ---")
    result = partition(X, obs, factors, adata.var_names)
    result.insert(len(factors) + 1, "highly_variable", hvg)
    result["mean_log2"] = adata.var["mean_log2"].to_numpy()
    result["detected_in_n_aoi"] = adata.var["detected_in_n_aoi"].to_numpy()

    n_converged = int(result["converged"].sum())
    method_counts = {k: int(v) for k, v in result["method"].value_counts().items()}
    emit(f"converged: {n_converged} of {len(result)} genes "
         f"({100 * n_converged / len(result):.2f}%)")
    emit(f"optimiser used: {method_counts} "
         "(powell is the fallback where lbfgs did not converge)")

    share_cols = [c for c in result.columns if c.startswith("var_")]
    row_sums = result.loc[result["converged"], share_cols].sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-9), "variance shares must sum to 1"
    emit(f"shares sum to 1 within 1e-9 on all {len(row_sums)} converged genes")

    all_genes = medians(result)
    hvg_only = medians(result.loc[result["highly_variable"]])
    emit("\nmedian share, all genes:")
    for name, value in all_genes.items():
        emit(f"  {name:<12} {value:.4f}")
    emit("median share, HVGs only:")
    for name, value in hvg_only.items():
        emit(f"  {name:<12} {value:.4f}")

    result.to_csv(snakemake.output.table, sep="\t", float_format="%.8g")
    emit(f"\nwrote {snakemake.output.table}")

    # ------------------------------------------------------------------ PVCA
    # Bushel's PVCA: eigendecompose the sample-sample correlation matrix, fit
    # the same variance components on each retained eigenvector, and average
    # the shares weighted by eigenvalue fraction. This is the number the
    # headline sentence quotes, because it is one number for the dataset rather
    # than a summary over a gene-wise distribution.
    emit("\n--- PVCA ---")
    corr = np.corrcoef(X - X.mean(axis=0, keepdims=True))
    eigenvalues, eigenvectors = np.linalg.eigh(corr)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
    fraction = eigenvalues / eigenvalues.sum()
    n_keep = int(np.searchsorted(np.cumsum(fraction), pvca_min_variance) + 1)
    emit(f"retained {n_keep} PCs for cumulative variance >= {pvca_min_variance} "
         f"({100 * fraction[:n_keep].sum():.1f}%)")

    pvca_vcs = build_vcspec(obs, factors)
    pvca_rows, pvca_weights = [], []
    for k in range(n_keep):
        shares, converged, _ = fit_shares(eigenvectors[:, k], pvca_vcs, adata.n_obs)
        if converged:
            pvca_rows.append(shares)
            pvca_weights.append(fraction[k])
        else:
            emit(f"  PC{k + 1} did not converge — excluded from the PVCA average")
    weights = np.asarray(pvca_weights) / np.sum(pvca_weights)
    pvca_shares = np.average(np.vstack(pvca_rows), axis=0, weights=weights)
    component_names = [short(f) for f in factors] + ["residual"]
    pvca = {n: round(float(v), 6) for n, v in zip(component_names, pvca_shares)}
    for name, value in pvca.items():
        emit(f"  {name:<12} {value:.4f}")

    # ------------------------------------------------------------ null floors
    emit("\n--- null baselines ---")
    analytic_null = {
        short(f): round((obs[f].astype(str).nunique() - 1) / (adata.n_obs - 1), 6)
        for f in factors
    }
    emit(f"analytic eta-squared null: {analytic_null}")

    naive = {
        short(f): round(float(np.nanmedian(naive_eta2(X, obs, f))), 6) for f in factors
    }
    emit(f"naive eta-squared, median over genes (do not report as a result): {naive}")

    permutation_null, permutation_null_range = None, None
    if n_permutations > 0:
        emit(f"permutation null: {n_permutations} shuffles x {n_permutation_genes} genes")
        per_perm = []
        for p in range(n_permutations):
            rng = np.random.default_rng(seed + p)
            shuffled = obs.copy()
            for factor in factors:
                shuffled[factor] = obs[factor].to_numpy()[rng.permutation(adata.n_obs)]
            genes = rng.choice(adata.n_vars, size=n_permutation_genes, replace=False)
            frame = partition(X[:, genes], shuffled, factors, adata.var_names[genes])
            per_perm.append(medians(frame))
        permutation_null = {
            name: round(float(np.median([m[name] for m in per_perm])), 6)
            for name in component_names
        }
        permutation_null_range = {
            name: [
                round(float(np.min([m[name] for m in per_perm])), 6),
                round(float(np.max([m[name] for m in per_perm])), 6),
            ]
            for name in component_names
        }
        emit(f"permutation null (median of per-shuffle medians): {permutation_null}")
        emit(f"permutation null (min-max across shuffles):       {permutation_null_range}")

    # ----------------------------------------------------------- sensitivity
    # All three run on the HVG subset so they are mutually comparable and the
    # rule stays inside a few minutes. Stated as such in the summary.
    emit(f"\n--- sensitivity fits (HVG subset, n={int(hvg.sum())} genes) ---")
    hvg_X, hvg_names = X[:, hvg], adata.var_names[hvg]

    s1_factors = ["patient_id", "aoi_code", "dsp_run"]
    s1 = medians(partition(hvg_X, obs, s1_factors, hvg_names))
    emit(f"S1 aoi_code (site x compartment as one 7-level cell): {s1}")

    s2_factors = ["patient_id", "site", "compartment"]
    s2 = medians(partition(hvg_X, obs, s2_factors, hvg_names))
    emit(f"S2 no dsp_run: {s2}")

    keep = (obs["aoi_code"].astype(str) != "BC").to_numpy()
    s3 = medians(partition(hvg_X[keep], obs.loc[keep], factors, hvg_names))
    emit(f"S3 no BC (n={int(keep.sum())} AOIs): {s3}")

    # The comparison S2 exists for: how much of `compartment` is really batch?
    baseline = hvg_only
    dsp_absorption = {
        name: round(s2[name] - baseline[name], 6) for name in s2 if name in baseline
    }
    emit(f"S2 - primary (HVG), i.e. what dropping dsp_run absorbs: {dsp_absorption}")

    # -------------------------------------------------------------- headline
    dominant = max(
        (n for n in component_names if n != "residual"), key=lambda n: pvca[n]
    )
    headline = (
        f"{100 * pvca[dominant]:.0f}% of variance is {dominant}-attributable "
        f"(PVCA over {n_keep} PCs spanning {100 * fraction[:n_keep].sum():.0f}% of "
        f"variance; per-gene median {100 * all_genes[dominant]:.0f}% over "
        f"{n_converged} genes"
        + (
            f"; permutation null {100 * permutation_null[dominant]:.1f}%"
            if permutation_null
            else ""
        )
        + "). TIME-B n = 8."
    )
    emit(f"\nGATE 1: {headline}")

    # ---------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.0))

    ax = axes[0]
    ok = result.loc[result["converged"]]
    data = [ok[f"var_{n}"].to_numpy() for n in component_names]
    box = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.62,
                     medianprops={"color": "black", "lw": 1.4})
    for patch, name in zip(box["boxes"], component_names):
        patch.set_facecolor(COMPONENT_COLOUR[name])
        patch.set_alpha(0.85)
    positions = np.arange(1, len(component_names) + 1)
    if permutation_null:
        ax.scatter(positions, [permutation_null[n] for n in component_names],
                   marker="x", s=70, color="crimson", zorder=6, lw=1.8,
                   label=f"permutation null ({n_permutations} shuffles)")
    ax.scatter(
        positions[: len(factors)],
        [analytic_null[short(f)] for f in factors],
        marker="_", s=260, color="black", zorder=6, lw=1.6,
        label=r"analytic $\eta^2$ null",
    )
    ax.set_xticks(positions, component_names, fontsize=9)
    ax.set_ylabel("fraction of variance", fontsize=9)
    ax.set_title(f"Per-gene partition ({n_converged} genes)", fontsize=11, fontweight="bold")
    ax.grid(axis="y", alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    ax = axes[1]
    ax.bar(positions, [pvca[n] for n in component_names],
           color=[COMPONENT_COLOUR[n] for n in component_names], alpha=0.88,
           edgecolor="black", lw=0.5)
    for x, name in zip(positions, component_names):
        ax.text(x, pvca[name] + 0.012, f"{100 * pvca[name]:.1f}%",
                ha="center", fontsize=8.5)
    ax.set_xticks(positions, component_names, fontsize=9)
    ax.set_ylabel("weighted fraction of variance", fontsize=9)
    ax.set_ylim(0, max(pvca.values()) * 1.22)
    ax.set_title(f"PVCA — {n_keep} PCs, "
                 f"{100 * fraction[:n_keep].sum():.0f}% of variance",
                 fontsize=11, fontweight="bold")
    ax.grid(axis="y", alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)

    ax = axes[2]
    compare = [n for n in ["patient", "site", "compartment", "residual"] if n in s2]
    idx = np.arange(len(compare))
    ax.bar(idx - 0.19, [baseline[n] for n in compare], width=0.38,
           color="#4F6D7A", alpha=0.9, edgecolor="black", lw=0.5, label="with dsp_run")
    ax.bar(idx + 0.19, [s2[n] for n in compare], width=0.38,
           color="#C1666B", alpha=0.9, edgecolor="black", lw=0.5, label="without dsp_run")
    for x, name in zip(idx, compare):
        delta = dsp_absorption[name]
        ax.text(x, max(baseline[name], s2[name]) + 0.008,
                f"{delta:+.3f}", ha="center", fontsize=8.5,
                color="crimson" if abs(delta) > 0.01 else "black")
    ax.set_xticks(idx, compare, fontsize=9)
    ax.set_ylabel("median fraction of variance", fontsize=9)
    ax.set_title(f"What dropping dsp_run absorbs (HVG, n={int(hvg.sum())})",
                 fontsize=11, fontweight="bold")
    ax.grid(axis="y", alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc="best")

    fig.suptitle(
        f"P1-T3 — variance partition (A2, Gate 1): {headline}",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(snakemake.output.figure, dpi=200, bbox_inches="tight")
    plt.close(fig)
    emit(f"wrote {snakemake.output.figure}")

    # --------------------------------------------------------------- summary
    Path(snakemake.output.summary).write_text(
        json.dumps(
            {
                "headline": headline,
                "dominant_component": dominant,
                "seed": seed,
                "n_aoi": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "n_genes_converged": n_converged,
                "n_hvg": int(hvg.sum()),
                "factors": factors,
                "model": (
                    "MixedLM, REML, crossed random intercepts via VCSpec "
                    "(single constant group); shares = component / (sum components "
                    "+ residual)"
                ),
                "pvca": pvca,
                "pvca_n_pcs": n_keep,
                "pvca_min_variance": pvca_min_variance,
                "pvca_variance_covered": round(float(fraction[:n_keep].sum()), 6),
                "median_share_all_genes": all_genes,
                "median_share_hvg_only": hvg_only,
                "null_analytic_eta2": analytic_null,
                "null_permutation": permutation_null,
                "null_permutation_range_across_shuffles": permutation_null_range,
                "optimiser_counts": method_counts,
                "n_permutations": n_permutations,
                "n_permutation_genes": n_permutation_genes,
                "naive_eta2_median_do_not_report": naive,
                "sensitivity_scope": (
                    f"S1-S3 are fitted on the {int(hvg.sum())} HVGs only, so they "
                    "are comparable to each other and to median_share_hvg_only, "
                    "not to median_share_all_genes."
                ),
                "sensitivity": {
                    "S1_aoi_code": {"factors": s1_factors, "median_share": s1},
                    "S2_no_dsp_run": {"factors": s2_factors, "median_share": s2},
                    "S3_no_bc": {
                        "factors": factors,
                        "n_aoi": int(keep.sum()),
                        "median_share": s3,
                    },
                },
                "dsp_run_absorption_hvg": dsp_absorption,
                "patients_spanning_both_runs": spanning,
                "caveats": [
                    "dsp_run is estimable only in L, LB, TIME-L and TIME-B; mLN, "
                    "TBME and BC sit wholly within one DSP run (Q3), so the "
                    "partition is conditional on that.",
                    "dsp_run has 2 levels and is near-nested in patient "
                    f"({spanning} of {obs['patient_id'].nunique()} patients span "
                    "both runs), so its component is poorly estimated.",
                    "site and compartment are structurally confounded: lymph_node "
                    "is tumour-only, glial_stroma and normal_control are "
                    "brain-only. S1 reports the alternative in which they are one "
                    "7-level factor.",
                    "TIME-B n = 8.",
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    emit(f"wrote {snakemake.output.summary}")
