"""P2-T1 — resolve the versioned gene sets against the measured transcriptome.

Owner task: P2-T1. Driven by rule p2t1_resolve_signatures.

`config/signatures.yaml` is a hand-written file, and hard constraint 5 says no
output is a result until a rule produces it. This is that rule. It turns the
YAML into a table that states, per gene, whether the symbol actually exists in
the 18,694 genes the assay measured — so a signature can never quietly shrink
between the file and the score.

The failure mode this exists to prevent: a symbol is renamed or mistyped, the
scoring step drops it without comment, and a four-gene signature is silently
scored on three. That is not a smaller signature, it is a different one.
**A symbol that does not resolve stops the workflow.** Recovery is a
stop-and-ask (adding or removing a gene from a signature is a scientific
decision, CLAUDE.md), not an edit to this script.

What it deliberately does NOT do: filter by detection. Whether a gene clears
background is a per-AOI, per-site property and it is P2-T2's job; conflating
"absent from the panel" with "below background in brain" is exactly the error
ADR 0008 warns about, since brain background is higher and would make
near-background genes look like missing ones.

Set-to-set overlap is reported rather than rejected. TIGIT and CTLA4 sit in
both `exhaustion` and the Phase 3 checkpoint panel, so the signature scores are
not independent draws and a reader should be able to see that from the table.
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

    signatures = snakemake.params.signatures
    emit(f"signature sets: {len(signatures)}")

    adata = ad.read_h5ad(snakemake.input.h5ad)
    measured = set(adata.var_names)
    emit(f"genes measured: {len(measured)} (from {snakemake.input.h5ad})")
    emit()

    rows = []
    unresolved = {}
    for name in sorted(signatures):
        spec = signatures[name]
        genes = list(spec["genes"])
        missing = [g for g in genes if g not in measured]
        if missing:
            unresolved[name] = missing
        for gene in genes:
            rows.append(
                {
                    "signature": name,
                    "gene": gene,
                    "in_matrix": gene in measured,
                    "n_genes_in_set": len(genes),
                    "source": spec["source"],
                }
            )
        emit(f"{name:22s} {len(genes) - len(missing)}/{len(genes)} resolved | {spec['source']}")

    if unresolved:
        detail = "; ".join(f"{k}: {v}" for k, v in unresolved.items())
        emit()
        emit(f"FAILED — unresolved symbols: {detail}")
        raise RuntimeError(
            "Signature symbols absent from the expression matrix: "
            f"{detail}. A signature that silently shrinks is a different "
            "signature, so this stops rather than drops. Changing set "
            "membership is a stop-and-ask decision (CLAUDE.md); if a symbol "
            "was renamed, resolve it against HGNC and record the change in "
            "an ADR before editing config/signatures.yaml."
        )

    membership = pd.DataFrame(rows)

    # --- overlap between sets, recorded ------------------------------------
    emit()
    names = sorted(signatures)
    overlaps = []
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            shared = sorted(set(signatures[a]["genes"]) & set(signatures[b]["genes"]))
            if shared:
                overlaps.append({"set_a": a, "set_b": b, "shared_genes": ",".join(shared)})
                emit(f"overlap {a} <-> {b}: {shared}")
    if not overlaps:
        emit("no gene is shared between two signature sets")

    membership.to_csv(snakemake.output.membership, sep="\t", index=False)

    summary = {
        "n_sets": len(signatures),
        "n_gene_slots": int(len(membership)),
        "n_unique_genes": int(membership["gene"].nunique()),
        "all_resolved": True,
        "n_genes_measured": len(measured),
        "set_sizes": {k: len(v["genes"]) for k, v in sorted(signatures.items())},
        "sources": {k: v["source"] for k, v in sorted(signatures.items())},
        "overlaps": overlaps,
        "h5ad": str(snakemake.input.h5ad),
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    emit()
    emit(
        f"resolved {summary['n_gene_slots']} gene slots "
        f"({summary['n_unique_genes']} unique) across {summary['n_sets']} sets"
    )
