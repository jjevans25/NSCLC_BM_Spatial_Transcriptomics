"""P3-T1 — resolve the pre-registered checkpoint panel against the transcriptome.

Owner task: P3-T1. Driven by rule p3t1_resolve_checkpoints.

The Phase 3 analogue of resolve_signatures.py, and deliberately the same shape:
`config/checkpoints.yaml` is hand-written, and hard constraint 5 says no output
is a result until a rule produces it. This is that rule. It turns the panel YAML
into a table stating, per gene, whether the symbol exists in the 18,694 genes
the assay measured, so the panel can never quietly shrink between the file and
the detection audit.

**A symbol that does not resolve stops the workflow.** Recovery is a
stop-and-ask (adding or removing a checkpoint gene is a scientific decision,
CLAUDE.md), not an edit to this script.

What it deliberately does NOT do:

  * **Filter by detection.** Whether a gene clears background is a per-AOI,
    per-compartment property and it is P3-T2's job. Conflating "absent from the
    panel" with "below background in brain" is the exact error ADR 0008 exists
    to prevent — brain background is HIGHER, so near-background genes would
    read as missing ones. Which genes turn out to be measurable is the Phase 3
    RESULT, not a Phase 3 input.

  * **Rank or select.** The panel is pre-registered (ADR 0015). Every gene in
    the file reaches every downstream table, including the ones that will be
    reported as not assessable everywhere. A panel that shed its undetected
    members would be a panel chosen after seeing the answer.

Overlap with the Phase 2 signatures is reported rather than rejected. Five of
the nine panel genes (PDCD1, CTLA4, LAG3, HAVCR2, TIGIT) are members of the
`exhaustion` set, which Phase 2 measured at 1 of 6 genes above background in
brain. A reader comparing the two phases' tables needs to see that they are not
independent draws — and that overlap is also the best available prior on what
P3-T2 will find.
"""

import json
from pathlib import Path

import anndata as ad
import pandas as pd

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    checkpoints = snakemake.params.checkpoints
    overlap = snakemake.params.signature_overlap
    emit(f"checkpoint panel: {len(checkpoints)} genes (PRE-REGISTERED, ADR 0015)")
    emit()

    adata = ad.read_h5ad(snakemake.input.h5ad)
    measured = set(adata.var_names)
    emit(f"genes measured: {len(measured)} (from {snakemake.input.h5ad})")
    emit()

    rows = []
    unresolved = []
    for gene in sorted(checkpoints):
        spec = checkpoints[gene]
        in_matrix = gene in measured
        if not in_matrix:
            unresolved.append(gene)
        rows.append(
            {
                "gene": gene,
                "alias": spec["alias"],
                "trial_stage": spec["trial_stage"],
                "in_matrix": in_matrix,
                "source": spec["source"],
                "rationale": " ".join(spec["rationale"].split()),
            }
        )
        mark = "" if in_matrix else "   <-- NOT IN MATRIX"
        emit(
            f"  {gene:8s} {spec['alias']:8s} {spec['trial_stage']:16s} "
            f"{spec['source']}{mark}"
        )

    if unresolved:
        emit()
        emit(f"FAILED — unresolved symbols: {unresolved}")
        raise RuntimeError(
            f"Checkpoint symbols absent from the expression matrix: {unresolved}. "
            "A panel that silently shrinks is a different panel, so this stops "
            "rather than drops. Changing panel membership is a stop-and-ask "
            "decision (CLAUDE.md); if a symbol was renamed, resolve it against "
            "HGNC and record the change in an ADR before editing "
            "config/checkpoints.yaml. Note that a gene being BELOW BACKGROUND "
            "is not this error — that is P3-T2's finding and is reported, never "
            "dropped."
        )

    membership = pd.DataFrame(rows)
    membership.to_csv(snakemake.output.membership, sep="\t", index=False)

    # --- overlap with the Phase 2 signatures, recorded ---------------------
    emit()
    if overlap:
        for set_name, shared in sorted(overlap.items()):
            emit(f"overlap with signature `{set_name}`: {shared}")
        emit(
            "  Recorded, not rejected. Phase 2 measured `exhaustion` at 1 of 6 "
            "genes above background in TIME-B; the overlapping panel genes "
            "should be expected to behave the same way (docs/limitations.md §10)."
        )
    else:
        emit("no panel gene appears in any Phase 2 signature")

    summary = {
        "n_genes": len(checkpoints),
        "all_resolved": True,
        "n_genes_measured": len(measured),
        "genes": sorted(checkpoints),
        "aliases": {g: checkpoints[g]["alias"] for g in sorted(checkpoints)},
        "trial_stages": {g: checkpoints[g]["trial_stage"] for g in sorted(checkpoints)},
        "sources": {g: checkpoints[g]["source"] for g in sorted(checkpoints)},
        "signature_overlap": {k: v for k, v in sorted(overlap.items())},
        "pre_registered": "ADR 0015",
        "h5ad": str(snakemake.input.h5ad),
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    emit()
    emit(f"resolved {len(checkpoints)}/{len(checkpoints)} panel genes")
