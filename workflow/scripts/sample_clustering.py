"""P1-T4 — hierarchical clustering of AOIs, annotated by all three factors.

Owner task: P1-T4. Driven by rule p1t4_sample_clustering.

PROJECT_PLAN §6 asks for a sample-sample correlation heatmap annotated with
compartment / site / patient, and to "note whether patients cluster together
across compartments (*they usually do*)". **In this dataset they do not**, and
that inversion is the finding rather than a footnote — P1-T2's patient panel
showed a patient's AOIs criss-crossing the whole ordination, and P1-T3 put the
number on it: compartment 41% by PVCA against patient 22%.

So this task is written to *test* that expectation rather than illustrate it.
Eyeballing a heatmap is not evidence, so every claim the figure makes is also
computed:

  * **Adjusted Rand index** between the dendrogram cut and each factor, swept
    over k, with a permutation null. ARI is implemented here in numpy rather
    than imported from scikit-learn, which is only *transitively* present in
    this env (via scanpy) and not declared in py-analysis.yaml. The env file
    already records what depending on an undeclared transitive import cost
    once; adding the dependency would also rebuild the env and re-run the
    whole DAG. Twenty lines of contingency-table arithmetic is the cheaper
    correct answer.
  * **Same-cluster rate**: of the AOI pairs sharing a factor level, what
    fraction land in the same cluster — against the baseline rate over all
    pairs. This is the direct form of the plan's question.
  * **Within- vs between-level mean correlation** per factor.

Choices worth stating, because each could reasonably go the other way:

  * **Spearman on the 2000 HVGs.** Over all 18,694 genes every AOI correlates
    with every other at ~0.99 — the shared transcriptome dominates and the
    structure washes out. Restricting to the same HVGs P1-T1/T2 used keeps this
    figure comparable to the ordination. Spearman rather than Pearson because
    the highest-variance genes here are immunoglobulins, whose dynamic range
    would otherwise drive the correlation.
  * **Average linkage on 1 - rho.** UPGMA is the linkage that pairs naturally
    with a correlation distance; Ward's method assumes Euclidean geometry the
    correlation distance does not have.
  * **`BC` controls are kept**, matching the P1-T3 primary fit and the design
    table. They are 7 of 120 AOIs and are labelled in the annotation strip, so
    a reader can see them cluster rather than have them silently removed.

All stochastic steps take config["seed"] explicitly.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata

COMPARTMENT_COLOUR = {
    "tumour": "#B0763F", "immune": "#3E7CB1",
    "glial_stroma": "#5B8C5A", "normal_control": "#8C8C8C",
}
SITE_COLOUR = {"lung": "#C1666B", "brain": "#4F6D7A", "lymph_node": "#D4B483"}
FACTOR_COLOUR = {
    "compartment": "#B0763F", "site": "#C1666B", "patient": "#7B4B94",
    "patient\n(cross-compartment)": "#9B72B0",
}


def adjusted_rand_index(a, b):
    """ARI between two labellings, from the contingency table.

    Implemented rather than imported: scikit-learn is not a declared dependency
    of py-analysis (see the module docstring).
    """
    table = pd.crosstab(pd.Series(a), pd.Series(b)).to_numpy(dtype=float)
    n = table.sum()
    if n < 2:
        return np.nan

    def comb2(x):
        return (x * (x - 1.0) / 2.0).sum()

    sum_ij = comb2(table)
    sum_i = comb2(table.sum(axis=1))
    sum_j = comb2(table.sum(axis=0))
    total = n * (n - 1.0) / 2.0
    expected = sum_i * sum_j / total
    maximum = 0.5 * (sum_i + sum_j)
    if maximum == expected:
        return 0.0
    return float((sum_ij - expected) / (maximum - expected))


def same_cluster_rate(labels, clusters):
    """Fraction of same-level AOI pairs that land in the same cluster."""
    labels = np.asarray(labels)
    clusters = np.asarray(clusters)
    same_level = labels[:, None] == labels[None, :]
    same_cluster = clusters[:, None] == clusters[None, :]
    upper = np.triu(np.ones_like(same_level, dtype=bool), k=1)
    pairs = same_level & upper
    if not pairs.any():
        return np.nan, 0
    return float(same_cluster[pairs].mean()), int(pairs.sum())


def within_between(corr, labels):
    """Mean correlation within a factor's levels vs between them."""
    labels = np.asarray(labels)
    same = labels[:, None] == labels[None, :]
    upper = np.triu(np.ones_like(same, dtype=bool), k=1)
    within = corr[same & upper]
    between = corr[~same & upper]
    return (
        float(within.mean()) if within.size else np.nan,
        float(between.mean()) if between.size else np.nan,
        int(within.size),
    )


def cross_compartment_patient(corr, patient, compartment):
    """Patient similarity with compartment held off the table.

    PROJECT_PLAN §6 asks whether patients cluster together *across
    compartments*. A raw within-patient mean cannot answer that: three patients
    contribute two AOIs of the SAME compartment (P12 and P24 two `TIME-L`, P15
    two `TBME`), so part of any within-patient signal is compartment agreement
    wearing a patient label. Restricting to pairs whose compartments differ
    removes that route entirely.
    """
    patient = np.asarray(patient)
    compartment = np.asarray(compartment)
    upper = np.triu(np.ones((len(patient), len(patient)), dtype=bool), k=1)
    differs = (compartment[:, None] != compartment[None, :]) & upper
    same_patient = (patient[:, None] == patient[None, :]) & differs
    other_patient = (patient[:, None] != patient[None, :]) & differs
    return (
        float(corr[same_patient].mean()) if same_patient.any() else np.nan,
        float(corr[other_patient].mean()) if other_patient.any() else np.nan,
        int(same_patient.sum()),
    )


log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    seed = snakemake.params.seed
    metric = snakemake.params.metric
    linkage_method = snakemake.params.linkage
    k_min, k_max = snakemake.params.k_range
    n_permutations = snakemake.params.n_permutations
    n_top = snakemake.params.n_top_genes
    flavour = snakemake.params.hvg_flavour
    sc.settings.verbosity = 0
    emit(f"seed: {seed} (config['seed'], passed explicitly)")
    emit(f"metric: {metric}, linkage: {linkage_method}, k: {k_min}-{k_max}")

    adata = sc.read_h5ad(snakemake.input.h5ad)
    obs = adata.obs.copy()
    emit(f"loaded {adata.n_obs} AOIs x {adata.n_vars} genes")

    # Same HVG set as P1-T1/T2/T3, so this figure is comparable to the ordination.
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor=flavour)
    hvg = adata.var["highly_variable"].to_numpy()
    X = np.asarray(adata.X, dtype=float)[:, hvg]
    emit(f"HVGs: {int(hvg.sum())} of {adata.n_vars} (flavour={flavour}, target {n_top})")

    # ------------------------------------------------------------- correlation
    if metric == "spearman":
        ranked = np.apply_along_axis(rankdata, 1, X)
    elif metric == "pearson":
        ranked = X
    else:
        raise ValueError(f"unsupported metric: {metric}")
    corr = np.corrcoef(ranked)
    np.fill_diagonal(corr, 1.0)
    off = corr[np.triu(np.ones_like(corr, dtype=bool), k=1)]
    emit(f"AOI-AOI {metric} correlation: median {np.median(off):.3f}, "
         f"range {off.min():.3f}-{off.max():.3f}")

    # -------------------------------------------------------------- clustering
    distance = 1.0 - corr
    np.fill_diagonal(distance, 0.0)
    distance = np.clip((distance + distance.T) / 2.0, 0.0, None)
    tree = linkage(squareform(distance, checks=False), method=linkage_method)
    order = dendrogram(tree, no_plot=True)["leaves"]

    factors = {
        "compartment": obs["compartment"].astype(str).to_numpy(),
        "site": obs["site"].astype(str).to_numpy(),
        "patient": obs["patient_id"].astype(str).to_numpy(),
    }

    # ----------------------------------------------------- agreement across k
    rng = np.random.default_rng(seed)
    rows = []
    for k in range(k_min, k_max + 1):
        clusters = fcluster(tree, t=k, criterion="maxclust")
        baseline_rate, _ = same_cluster_rate(np.zeros(adata.n_obs), clusters)
        for name, labels in factors.items():
            ari = adjusted_rand_index(labels, clusters)
            null = [
                adjusted_rand_index(rng.permutation(labels), clusters)
                for _ in range(n_permutations)
            ]
            rate, n_pairs = same_cluster_rate(labels, clusters)
            rows.append(
                {
                    "k": k,
                    "factor": name,
                    "ari": ari,
                    "ari_null_mean": float(np.mean(null)),
                    "ari_null_p95": float(np.percentile(null, 95)),
                    "same_cluster_rate": rate,
                    "same_cluster_rate_baseline": baseline_rate,
                    "n_same_level_pairs": n_pairs,
                }
            )
    agreement = pd.DataFrame(rows)
    agreement.to_csv(snakemake.output.table, sep="\t", index=False, float_format="%.8g")
    emit(f"\nwrote {snakemake.output.table}")

    best = {
        name: agreement.loc[agreement.query("factor == @name")["ari"].idxmax()]
        for name in factors
    }
    emit("\npeak ARI across k:")
    for name, row in best.items():
        emit(f"  {name:<12} ARI {row['ari']:.3f} at k={int(row['k'])} "
             f"(null mean {row['ari_null_mean']:.3f}, p95 {row['ari_null_p95']:.3f}); "
             f"same-cluster rate {row['same_cluster_rate']:.3f} vs baseline "
             f"{row['same_cluster_rate_baseline']:.3f} over "
             f"{int(row['n_same_level_pairs'])} same-level pairs")

    # ---------------------------------------------- within vs between levels
    emit("\nmean correlation within vs between levels:")
    contrast = {}
    for name, labels in factors.items():
        within, between, n_pairs = within_between(corr, labels)
        contrast[name] = {
            "within": round(within, 6),
            "between": round(between, 6),
            "difference": round(within - between, 6),
            "n_within_pairs": n_pairs,
        }
        emit(f"  {name:<12} within {within:.3f}  between {between:.3f}  "
             f"difference {within - between:+.3f}  ({n_pairs} within-level pairs)")

    within_cc, between_cc, n_cc = cross_compartment_patient(
        corr, factors["patient"], factors["compartment"]
    )
    contrast["patient\n(cross-compartment)"] = {
        "within": round(within_cc, 6),
        "between": round(between_cc, 6),
        "difference": round(within_cc - between_cc, 6),
        "n_within_pairs": n_cc,
    }
    emit(f"  {'patient (cross-compartment)':<12} within {within_cc:.3f}  "
         f"between {between_cc:.3f}  difference {within_cc - between_cc:+.3f}  "
         f"({n_cc} same-patient pairs spanning different compartments)")

    # The plan's question, in its most direct form: at the k where compartment
    # agrees best, do a patient's AOIs land together?
    k_star = int(best["compartment"]["k"])
    clusters_star = fcluster(tree, t=k_star, criterion="maxclust")
    verdict_rate, verdict_pairs = same_cluster_rate(factors["patient"], clusters_star)
    baseline_star, _ = same_cluster_rate(np.zeros(adata.n_obs), clusters_star)
    multi = obs["patient_id"].value_counts()
    n_multi = int((multi > 1).sum())
    emit(f"\nAt k={k_star} (compartment's best): {verdict_pairs} AOI pairs share a "
         f"patient; {100 * verdict_rate:.0f}% of them fall in the same cluster, "
         f"against a {100 * baseline_star:.0f}% baseline over all pairs. "
         f"{n_multi} patients contribute more than one AOI.")

    # ------------------------------------------------------------------ figure
    fig = plt.figure(figsize=(17.5, 9.6))
    outer = fig.add_gridspec(1, 2, width_ratios=[3.3, 1.25], wspace=0.22)
    left = outer[0].subgridspec(4, 2, height_ratios=[0.9, 0.42, 4.4, 0.0],
                                width_ratios=[24, 1], hspace=0.06, wspace=0.03)

    ax_tree = fig.add_subplot(left[0, 0])
    dendrogram(tree, ax=ax_tree, color_threshold=0, above_threshold_color="#444444",
               no_labels=True)
    ax_tree.set_xticks([])
    ax_tree.set_yticks([])
    for spine in ax_tree.spines.values():
        spine.set_visible(False)
    ax_tree.set_title(
        f"P1-T4 — AOI-AOI {metric} correlation over {int(hvg.sum())} HVGs, "
        f"{linkage_method} linkage on 1-rho ({adata.n_obs} AOIs, TIME-B n=8)",
        fontsize=12, fontweight="bold", pad=10,
    )

    # Annotation strips, in the dendrogram's leaf order.
    ax_ann = fig.add_subplot(left[1, 0])
    patient_levels = sorted(set(factors["patient"]))
    patient_cmap = plt.get_cmap("hsv")
    patient_colour = {
        p: patient_cmap(i / max(len(patient_levels) - 1, 1))
        for i, p in enumerate(patient_levels)
    }
    strips = [
        ("compartment", [COMPARTMENT_COLOUR[v] for v in factors["compartment"][order]]),
        ("site", [SITE_COLOUR[v] for v in factors["site"][order]]),
        ("patient", [patient_colour[v] for v in factors["patient"][order]]),
    ]
    for row, (name, colours) in enumerate(strips):
        for col, colour in enumerate(colours):
            ax_ann.add_patch(plt.Rectangle((col, len(strips) - row - 1), 1, 1,
                                           facecolor=colour, edgecolor="none"))
    ax_ann.set_xlim(0, adata.n_obs)
    ax_ann.set_ylim(0, len(strips))
    ax_ann.set_xticks([])
    ax_ann.set_yticks([len(strips) - i - 0.5 for i in range(len(strips))])
    ax_ann.set_yticklabels([name for name, _ in strips], fontsize=8.5)
    for spine in ax_ann.spines.values():
        spine.set_visible(False)

    ax_heat = fig.add_subplot(left[2, 0])
    ordered = corr[np.ix_(order, order)]
    image = ax_heat.imshow(ordered, cmap="RdYlBu_r", aspect="auto",
                           vmin=np.percentile(off, 1), vmax=1.0)
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])
    ax_heat.set_xlabel(f"AOIs, ordered by the dendrogram above", fontsize=9)

    ax_cbar = fig.add_subplot(left[2, 1])
    fig.colorbar(image, cax=ax_cbar).set_label(f"{metric} rho", fontsize=8.5)
    ax_cbar.tick_params(labelsize=7.5)

    legend_handles = (
        [Patch(facecolor=c, label=v) for v, c in COMPARTMENT_COLOUR.items()]
        + [Patch(facecolor=c, label=v) for v, c in SITE_COLOUR.items()]
        + [Patch(facecolor="white", edgecolor="#999999",
                 label=f"patient: {len(patient_levels)} levels, no legend")]
    )
    ax_heat.legend(handles=legend_handles, frameon=False, fontsize=7.5,
                   loc="upper left", bbox_to_anchor=(0.0, -0.03), ncol=4)

    right = outer[1].subgridspec(2, 1, height_ratios=[1, 1], hspace=0.34)

    ax_ari = fig.add_subplot(right[0])
    for name in factors:
        sub = agreement.query("factor == @name")
        ax_ari.plot(sub["k"], sub["ari"], marker="o", ms=4, lw=1.6,
                    color=FACTOR_COLOUR[name], label=name)
    null_band = agreement.groupby("k")["ari_null_p95"].max()
    ax_ari.plot(null_band.index, null_band.to_numpy(), ls="--", lw=1.1,
                color="crimson", label=f"null p95 ({n_permutations} shuffles)")
    ax_ari.axvline(k_star, color="#999999", lw=0.9, ls=":")
    ax_ari.set_xlabel("clusters (k)", fontsize=9)
    ax_ari.set_ylabel("adjusted Rand index", fontsize=9)
    ax_ari.set_title("Cluster agreement with each factor", fontsize=10.5,
                     fontweight="bold")
    ax_ari.grid(alpha=0.25, lw=0.6)
    ax_ari.set_axisbelow(True)
    ax_ari.legend(frameon=False, fontsize=7.5, loc="best")

    ax_wb = fig.add_subplot(right[1])
    names = list(contrast)
    idx = np.arange(len(names))
    ax_wb.bar(idx - 0.19, [contrast[n]["within"] for n in names], width=0.38,
              color=[FACTOR_COLOUR[n] for n in names], alpha=0.95,
              edgecolor="black", lw=0.5, label="within level")
    ax_wb.bar(idx + 0.19, [contrast[n]["between"] for n in names], width=0.38,
              color=[FACTOR_COLOUR[n] for n in names], alpha=0.42,
              edgecolor="black", lw=0.5, label="between levels")
    for x, name in zip(idx, names):
        ax_wb.text(x, max(contrast[name]["within"], contrast[name]["between"]) + 0.004,
                   f"{contrast[name]['difference']:+.3f}", ha="center", fontsize=8)
    ax_wb.set_xticks(idx, names, fontsize=7.6)
    ax_wb.set_ylabel(f"mean {metric} rho", fontsize=9)
    ax_wb.set_title("Correlation within vs between levels", fontsize=10.5,
                    fontweight="bold")
    ax_wb.grid(axis="y", alpha=0.25, lw=0.6)
    ax_wb.set_axisbelow(True)
    ax_wb.legend(frameon=False, fontsize=7.5, loc="best")

    fig.savefig(snakemake.output.figure, dpi=200, bbox_inches="tight")
    plt.close(fig)
    emit(f"\nwrote {snakemake.output.figure}")

    # ----------------------------------------------------------------- summary
    # Both halves of this are true and they are not in tension: compartment owns
    # the top-level partition, while patient is a real second-order effect that
    # never becomes cluster structure. Stating only the first half would
    # understate the non-independence that P1-T6 has to argue from.
    verdict = (
        f"Compartment structures the dendrogram and patients do not cluster "
        f"together: at k={k_star}, {100 * verdict_rate:.0f}% of the "
        f"{verdict_pairs} AOI pairs sharing a patient fall in the same cluster, "
        f"against a {100 * baseline_star:.0f}% baseline over all pairs, and peak "
        f"patient ARI is {best['patient']['ari']:.3f} against compartment's "
        f"{best['compartment']['ari']:.3f}. PROJECT_PLAN §6 expected the "
        f"opposite. Patient similarity is nonetheless real and is not compartment "
        f"in disguise: across {n_cc} same-patient AOI pairs whose compartments "
        f"differ, mean rho is {within_cc:.3f} against {between_cc:.3f} for "
        f"different-patient pairs ({within_cc - between_cc:+.3f}). That is the "
        f"non-independence (1|patient) rests on. TIME-B n = 8."
    )
    emit(f"\nVERDICT: {verdict}")

    Path(snakemake.output.summary).write_text(
        json.dumps(
            {
                "verdict": verdict,
                "seed": seed,
                "n_aoi": int(adata.n_obs),
                "n_hvg": int(hvg.sum()),
                "metric": metric,
                "linkage": linkage_method,
                "distance": "1 - rho",
                "k_range": [k_min, k_max],
                "n_permutations": n_permutations,
                "correlation_median": round(float(np.median(off)), 6),
                "correlation_min": round(float(off.min()), 6),
                "correlation_max": round(float(off.max()), 6),
                "peak_ari": {
                    name: {
                        "k": int(row["k"]),
                        "ari": round(float(row["ari"]), 6),
                        "ari_null_mean": round(float(row["ari_null_mean"]), 6),
                        "ari_null_p95": round(float(row["ari_null_p95"]), 6),
                        "same_cluster_rate": round(float(row["same_cluster_rate"]), 6),
                        "same_cluster_rate_baseline": round(
                            float(row["same_cluster_rate_baseline"]), 6
                        ),
                        "n_same_level_pairs": int(row["n_same_level_pairs"]),
                    }
                    for name, row in best.items()
                },
                "within_vs_between": contrast,
                "cross_compartment_patient": {
                    "within_patient": round(within_cc, 6),
                    "between_patient": round(between_cc, 6),
                    "difference": round(within_cc - between_cc, 6),
                    "n_pairs": n_cc,
                    "note": (
                        "Restricted to AOI pairs whose compartments differ, so "
                        "no part of this contrast can be compartment agreement "
                        "wearing a patient label."
                    ),
                },
                "k_star": k_star,
                "patients_with_multiple_aoi": n_multi,
                "caveats": [
                    "Correlation is computed on the 2000 HVGs, not all 18,694 "
                    "genes: over the full transcriptome every AOI pair "
                    "correlates at ~0.99 and the structure is not visible.",
                    "site and compartment are structurally confounded "
                    "(lymph_node is tumour-only, glial_stroma and "
                    "normal_control are brain-only), so their ARI curves are "
                    "not independent of each other.",
                    "TIME-B n = 8.",
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    emit(f"wrote {snakemake.output.summary}")
