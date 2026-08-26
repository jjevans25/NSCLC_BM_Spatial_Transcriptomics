# ---------------------------------------------------------------------------
# Phase 1 — variance landscape: PCA, neighbours, clustering, UMAP over AOIs.
#
# Owner task: Phase 1 (see Markdowns/PROJECT_PLAN.md §6).
# STUB — no rules yet. TARGETS_LANDSCAPE stays empty so `snakemake -n` reports an
# empty DAG rather than a missing-input error. Populate both the rules and
# this list in the same commit.
# ---------------------------------------------------------------------------

rule p1t1_pca_landscape:
    """Feature selection, PCA, scree, and the tri-coloured ordination (A2)."""
    input:
        h5ad=f"{PATHS['interim']}/aoi_normalised.h5ad",
    output:
        scree=report(
            f"{PATHS['figures']}/pca_scree.png",
            caption="../report/pca_scree.rst",
            category="Phase 1 — variance landscape",
            labels={"task": "P1-T1", "figure": "scree"},
        ),
        ordination=report(
            f"{PATHS['figures']}/pca_ordination.png",
            caption="../report/pca_ordination.rst",
            category="Phase 1 — variance landscape",
            labels={"task": "P1-T2", "figure": "tri-coloured ordination"},
        ),
        coords=f"{PATHS['tables']}/pca_coords.tsv",
        loadings=f"{PATHS['tables']}/pca_loadings.tsv",
        umap=f"{PATHS['tables']}/umap_coords.tsv",
        summary=f"{PATHS['tables']}/pca_summary.json",
    params:
        seed=SEED,
        n_top_genes=config["landscape"]["n_top_genes"],
        hvg_flavour=config["landscape"]["hvg_flavour"],
        n_pcs=config["landscape"]["n_pcs"],
    log:
        f"{PATHS['logs']}/p1t1_pca_landscape.log",
    benchmark:
        f"{PATHS['benchmarks']}/p1t1_pca_landscape.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 2
    script:
        "../scripts/pca_landscape.py"


rule p1t3_variance_partition:
    """Variance components per gene + PVCA: patient vs site vs compartment (A2).

    This is Gate 1. `dsp_run` is a term because P1-T2 showed a consistent run-B
    shift on PC1, and leaving it out would inflate `compartment` — the number
    the gate turns on.
    """
    input:
        h5ad=f"{PATHS['interim']}/aoi_normalised.h5ad",
    output:
        table=f"{PATHS['tables']}/variance_partition.tsv",
        summary=f"{PATHS['tables']}/variance_partition_summary.json",
        figure=report(
            f"{PATHS['figures']}/variance_partition.png",
            caption="../report/variance_partition.rst",
            category="Phase 1 — variance landscape",
            labels={"task": "P1-T3", "figure": "variance partition"},
        ),
    params:
        seed=SEED,
        factors=config["landscape"]["variance"]["factors"],
        pvca_min_variance=config["landscape"]["variance"]["pvca_min_variance"],
        n_permutations=config["landscape"]["variance"]["n_permutations"],
        n_permutation_genes=config["landscape"]["variance"]["n_permutation_genes"],
        n_top_genes=config["landscape"]["n_top_genes"],
        hvg_flavour=config["landscape"]["hvg_flavour"],
    log:
        f"{PATHS['logs']}/p1t3_variance_partition.log",
    benchmark:
        f"{PATHS['benchmarks']}/p1t3_variance_partition.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 2
    script:
        "../scripts/variance_partition.py"


rule p1t4_sample_clustering:
    """AOI-AOI correlation heatmap + cluster agreement with each factor (A2).

    PROJECT_PLAN §6 expects patients to cluster together across compartments.
    P1-T2 and P1-T3 both say they do not, so this rule computes the agreement
    (ARI, same-cluster rate, within-vs-between correlation) rather than leaving
    the claim to the eye.
    """
    input:
        h5ad=f"{PATHS['interim']}/aoi_normalised.h5ad",
    output:
        table=f"{PATHS['tables']}/cluster_agreement.tsv",
        summary=f"{PATHS['tables']}/sample_clustering_summary.json",
        figure=report(
            f"{PATHS['figures']}/sample_correlation_heatmap.png",
            caption="../report/sample_clustering.rst",
            category="Phase 1 — variance landscape",
            labels={"task": "P1-T4", "figure": "correlation heatmap"},
        ),
    params:
        seed=SEED,
        metric=config["landscape"]["clustering"]["metric"],
        linkage=config["landscape"]["clustering"]["linkage"],
        k_range=config["landscape"]["clustering"]["k_range"],
        n_permutations=config["landscape"]["clustering"]["n_permutations"],
        n_top_genes=config["landscape"]["n_top_genes"],
        hvg_flavour=config["landscape"]["hvg_flavour"],
    log:
        f"{PATHS['logs']}/p1t4_sample_clustering.log",
    benchmark:
        f"{PATHS['benchmarks']}/p1t4_sample_clustering.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 2
    script:
        "../scripts/sample_clustering.py"


rule p1t5_landscape_explorer:
    """Export the app-tier landscape explorer to self-contained WASM HTML.

    Every path the notebook reads is declared here. An undeclared read is a
    silent provenance hole (PROJECT_PLAN §A.5).

    The notebook is staged into a build directory beside a `public/` folder
    holding the three pipeline outputs, because `marimo export html-wasm`
    copies `public/` into the export and `mo.notebook_location()` resolves to
    the served URL in the browser. That is what lets the exported app fetch its
    data over HTTP, which is what keeps `mo.ui.table` interactive: bundling the
    frames through `cache_cells` instead would cache each table cell's output
    as an `UnhashableStub` and render stub text in place of the tables.

    Staged under results/ rather than next to the notebook so the source tree
    stays free of generated files.
    """
    input:
        notebook="notebooks/apps/landscape_explorer.py",
        coords=f"{PATHS['tables']}/pca_coords.tsv",
        umap=f"{PATHS['tables']}/umap_coords.tsv",
        summary=f"{PATHS['tables']}/pca_summary.json",
    output:
        app=report(
            directory(f"{PATHS['reports']}/landscape_explorer"),
            htmlindex="index.html",
            caption="../report/landscape_explorer.rst",
            category="Interactive",
            labels={"task": "P1-T5", "app": "landscape explorer"},
        ),
    params:
        build=f"{PATHS['interim']}/explorer_build",
    log:
        f"{PATHS['logs']}/p1t5_landscape_explorer.log",
    benchmark:
        f"{PATHS['benchmarks']}/p1t5_landscape_explorer.tsv"
    conda:
        "../envs/py-analysis.yaml"
    threads: 1
    shell:
        "rm -rf {params.build} && mkdir -p {params.build}/public && "
        "cp {input.notebook} {params.build}/landscape_explorer.py && "
        "cp {input.coords} {input.umap} {input.summary} {params.build}/public/ && "
        "marimo export html-wasm --execute --mode run -f "
        "{params.build}/landscape_explorer.py -o {output.app} > {log} 2>&1"


TARGETS_LANDSCAPE = [
    f"{PATHS['figures']}/pca_scree.png",
    f"{PATHS['figures']}/pca_ordination.png",
    f"{PATHS['tables']}/pca_coords.tsv",
    f"{PATHS['tables']}/pca_loadings.tsv",
    f"{PATHS['tables']}/umap_coords.tsv",
    f"{PATHS['tables']}/pca_summary.json",
    f"{PATHS['tables']}/variance_partition.tsv",
    f"{PATHS['tables']}/variance_partition_summary.json",
    f"{PATHS['figures']}/variance_partition.png",
    f"{PATHS['tables']}/cluster_agreement.tsv",
    f"{PATHS['tables']}/sample_clustering_summary.json",
    f"{PATHS['figures']}/sample_correlation_heatmap.png",
    f"{PATHS['reports']}/landscape_explorer",
]
