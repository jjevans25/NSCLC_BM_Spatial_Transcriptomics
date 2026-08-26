"""P1-T1/T2 — feature selection, PCA, and the tri-coloured ordination.

Owner tasks: P1-T1, P1-T2. Driven by rule p1t1_pca_landscape.

Aim A2: what actually drives variance here — patient, site, or compartment?
This produces the ordination that makes the answer visible; P1-T3 puts a number
on it. The two tasks share one PCA, so they share one rule rather than
recomputing it.

Choices worth stating, because each could reasonably go the other way:

  * **HVG selection on log2(Q3+1)**, scanpy's `seurat` flavour, which expects
    log input. X is exactly that (P0-T7).
  * **PCA is zero-centred but not unit-scaled.** Scaling every HVG to unit
    variance is the single-cell convention; it would give a gene detected in
    three AOIs the same leverage as one detected in all 120. Given P0-T5's
    finding that detection varies systematically by compartment, that would
    convert a detection artefact into a principal component. Bulk-style
    zero-centred PCA is the more conservative reading.
  * **QC-flagged AOIs are retained and marked**, never dropped (`qc.flag_only`).
    The ordination marks them so a reviewer can see whether they sit apart —
    which is the question flagging exists to raise.

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

CODE_ORDER = ["L", "LB", "mLN", "TIME-L", "TIME-B", "TBME", "BC"]
COMPARTMENT_COLOUR = {
    "tumour": "#B0763F", "immune": "#3E7CB1",
    "glial_stroma": "#5B8C5A", "normal_control": "#8C8C8C",
}
SITE_COLOUR = {"lung": "#C1666B", "brain": "#4F6D7A", "lymph_node": "#D4B483"}

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    seed = snakemake.params.seed
    n_top = snakemake.params.n_top_genes
    flavour = snakemake.params.hvg_flavour
    n_pcs = snakemake.params.n_pcs
    sc.settings.verbosity = 0
    np.random.seed(seed)
    emit(f"seed: {seed} (config['seed'], passed explicitly)")

    adata = sc.read_h5ad(snakemake.input.h5ad)
    emit(f"loaded {adata.n_obs} AOIs x {adata.n_vars} genes")

    # ------------------------------------------------------- feature selection
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor=flavour)
    hvg = adata.var["highly_variable"]
    emit(f"HVGs: {int(hvg.sum())} of {adata.n_vars} (flavour={flavour}, target {n_top})")

    sub = adata[:, hvg].copy()

    # -------------------------------------------------------------------- PCA
    sc.tl.pca(sub, n_comps=n_pcs, svd_solver="arpack", random_state=seed, zero_center=True)
    variance_ratio = sub.uns["pca"]["variance_ratio"]
    emit(f"PCA: {n_pcs} components, PC1-PC4 explain "
         f"{100 * variance_ratio[:4].sum():.1f}% of variance")
    for i in range(4):
        emit(f"  PC{i + 1}: {100 * variance_ratio[i]:5.2f}%")

    # ------------------------------------------------------------ PC coordinates
    coords = pd.DataFrame(
        sub.obsm["X_pca"][:, :n_pcs],
        index=sub.obs_names,
        columns=[f"PC{i + 1}" for i in range(n_pcs)],
    )
    meta_columns = ["aoi_code", "patient_id", "site", "compartment", "dsp_run", "qc_flag"]
    out = sub.obs[meta_columns].join(coords)
    out.index.name = "aoi_label"
    out.to_csv(snakemake.output.coords, sep="\t", float_format="%.8g")
    emit(f"\nwrote {snakemake.output.coords}")

    # ------------------------------------------------------- PC1-PC4 loadings
    loadings = pd.DataFrame(
        sub.varm["PCs"][:, :4],
        index=sub.var_names,
        columns=[f"PC{i + 1}" for i in range(4)],
    )
    loadings.index.name = "gene"
    loadings.to_csv(snakemake.output.loadings, sep="\t", float_format="%.8g")
    emit(f"wrote {snakemake.output.loadings}")
    for pc in ["PC1", "PC2", "PC3", "PC4"]:
        top = loadings[pc].abs().sort_values(ascending=False).head(8).index.tolist()
        emit(f"  {pc} top |loading|: {', '.join(top)}")

    # ------------------------------------------------------------ scree figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    idx = np.arange(1, n_pcs + 1)
    ax1.bar(idx, 100 * variance_ratio, color="#3E7CB1", alpha=0.85)
    ax1.set_xlabel("principal component", fontsize=9)
    ax1.set_ylabel("% variance explained", fontsize=9)
    ax1.set_title("Scree", fontsize=11, fontweight="bold")
    ax1.grid(axis="y", alpha=0.25, lw=0.6)
    ax1.set_axisbelow(True)

    ax2.plot(idx, 100 * np.cumsum(variance_ratio), marker="o", ms=3.5, lw=1.5, color="#B0763F")
    ax2.axhline(80, ls="--", lw=1.0, color="crimson")
    n80 = int(np.searchsorted(np.cumsum(variance_ratio), 0.80) + 1)
    ax2.text(n_pcs * 0.55, 82, f"80% at PC{n80}", color="crimson", fontsize=9)
    ax2.set_xlabel("principal component", fontsize=9)
    ax2.set_ylabel("cumulative % variance", fontsize=9)
    ax2.set_title("Cumulative", fontsize=11, fontweight="bold")
    ax2.grid(alpha=0.25, lw=0.6)
    ax2.set_axisbelow(True)

    fig.suptitle(
        f"P1-T1 — PCA on {int(hvg.sum())} highly variable genes, {adata.n_obs} AOIs "
        f"(TIME-B n=8)", fontsize=12, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(snakemake.output.scree, dpi=200, bbox_inches="tight")
    plt.close(fig)
    emit(f"wrote {snakemake.output.scree}")

    # --------------------------------------------------------------- UMAP
    sc.pp.neighbors(sub, n_neighbors=10, n_pcs=min(n_pcs, 15), random_state=seed)
    sc.tl.umap(sub, random_state=seed)
    umap = pd.DataFrame(sub.obsm["X_umap"], index=sub.obs_names, columns=["UMAP1", "UMAP2"])
    umap.index.name = "aoi_label"
    umap.to_csv(snakemake.output.umap, sep="\t", float_format="%.8g")
    emit(f"wrote {snakemake.output.umap}")

    # ----------------------------------------------- tri-coloured ordination
    pc1, pc2 = coords["PC1"], coords["PC2"]
    v1, v2 = 100 * variance_ratio[0], 100 * variance_ratio[1]
    flagged = sub.obs["qc_flag"].astype(str).isin(["True", "true"]).to_numpy()

    u1, u2 = umap["UMAP1"], umap["UMAP2"]
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 10.4))

    def scatter(ax, x, y, colour_map, key, title, xlab, ylab, legend):
        for value, colour in colour_map.items():
            mask = (sub.obs[key] == value).to_numpy()
            if not mask.any():
                continue
            ax.scatter(x[mask], y[mask], s=46, c=colour, edgecolor="black",
                       linewidth=0.5, alpha=0.88, label=value)
        # Mark, never drop: flagged AOIs get a ring so a reviewer can see whether
        # they sit apart. That is the question flagging exists to raise.
        ax.scatter(x[flagged], y[flagged], s=170, facecolors="none",
                   edgecolors="crimson", linewidth=1.4, zorder=5,
                   label=f"QC-flagged (n={int(flagged.sum())})")
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel(xlab, fontsize=9)
        ax.set_ylabel(ylab, fontsize=9)
        ax.grid(alpha=0.22, lw=0.6)
        ax.set_axisbelow(True)
        if legend:
            ax.legend(frameon=False, fontsize=8, loc="best")

    def by_patient(ax, x, y, title, xlab, ylab):
        patients = sorted(sub.obs["patient_id"].unique())
        cmap = plt.get_cmap("hsv")
        for i, patient in enumerate(patients):
            mask = (sub.obs["patient_id"] == patient).to_numpy()
            colour = cmap(i / max(len(patients) - 1, 1))
            if mask.sum() > 1:
                order = np.argsort(x[mask].to_numpy())
                ax.plot(x[mask].to_numpy()[order], y[mask].to_numpy()[order],
                        color=colour, lw=0.8, alpha=0.55, zorder=1)
            ax.scatter(x[mask], y[mask], s=46, color=colour, edgecolor="black",
                       linewidth=0.5, alpha=0.9, zorder=2)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel(xlab, fontsize=9)
        ax.set_ylabel(ylab, fontsize=9)
        ax.grid(alpha=0.22, lw=0.6)
        ax.set_axisbelow(True)

    px, py = f"PC1 ({v1:.1f}%)", f"PC2 ({v2:.1f}%)"
    scatter(axes[0, 0], pc1, pc2, COMPARTMENT_COLOUR, "compartment",
            "PCA — by compartment", px, py, True)
    scatter(axes[0, 1], pc1, pc2, SITE_COLOUR, "site", "PCA — by site", px, py, True)
    scatter(axes[1, 0], u1, u2, COMPARTMENT_COLOUR, "compartment",
            "UMAP — by compartment", "UMAP1", "UMAP2", False)
    scatter(axes[1, 1], u1, u2, SITE_COLOUR, "site", "UMAP — by site",
            "UMAP1", "UMAP2", False)
    by_patient(axes[1, 2], u1, u2,
               f"UMAP — by patient ({sub.obs['patient_id'].nunique()} subjects)",
               "UMAP1", "UMAP2")

    # Patient: 42 subjects, so a categorical legend is useless. Colour by hue and
    # join each patient's AOIs with a line — the visual question is whether a
    # patient's AOIs stay together, not which patient is which.
    by_patient(axes[0, 2], pc1, pc2,
               f"PCA — by patient ({sub.obs['patient_id'].nunique()} subjects, AOIs joined)",
               px, py)

    fig.suptitle(
        "P1-T2 — PCA (top) and UMAP (bottom), each coloured three ways: "
        f"what drives the variance? (A2; {adata.n_obs} AOIs, TIME-B n=8)",
        fontsize=12.5, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    fig.savefig(snakemake.output.ordination, dpi=200, bbox_inches="tight")
    plt.close(fig)
    emit(f"wrote {snakemake.output.ordination}")

    Path(snakemake.output.summary).write_text(
        json.dumps(
            {
                "seed": seed,
                "n_aoi": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "n_hvg": int(hvg.sum()),
                "hvg_flavour": flavour,
                "n_pcs": n_pcs,
                "zero_centered_not_scaled": True,
                # Full vector, not just the first four: P1-T5's explorer lets a
                # reader put any PC on either axis and needs the % for the axis
                # label. An app-tier notebook may only read declared rule
                # outputs, so the number has to be here rather than recomputed.
                "variance_ratio": [round(float(v), 8) for v in variance_ratio],
                "variance_ratio_pc1_4": [round(float(v), 6) for v in variance_ratio[:4]],
                "cumulative_pc1_4": round(float(variance_ratio[:4].sum()), 6),
                "n_pcs_for_80pct": n80,
                "n_qc_flagged_retained": int(flagged.sum()),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    emit(f"wrote {snakemake.output.summary}")
