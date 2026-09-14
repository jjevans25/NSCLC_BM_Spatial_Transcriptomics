"""P4-T7 — the inferred-crosstalk network figure.

Owner task: P4-T7. Driven by rule p4t7_network_figure.

Nodes are genes, coloured by the compartment they are measured in; edges are
inferred ligand–receptor relationships, coloured by rho. One panel per
adjacency, laid out bipartite — tumour-side partners on the left, immune-side on
the right — because that is what the design actually is, and a spring layout
would invent a topology the data does not have.

**Never "colocalisation" and never "spatially adjacent."** There are no
coordinates. An edge here means *the two compartments' expression correlated
across patients*, which is inference from a database plus a correlation, not a
measurement of proximity. `p4t7_language_audit` enforces the wording over the
whole repository rather than trusting anyone to remember it.

**Edge colour is `RdBu_r`, and that is a deliberate reuse rather than a
borrow.** The project's colour language is already fixed: `RdBu_r` means a
signed effect (P2-T7) and `PRGn` centred on the detection multiple means a
background ratio (P3-T5, ADR 0016 §5). A Spearman rho *is* a signed effect on
[-1, 1] — the same quantity class `RdBu_r` was established for — so speaking
that language is consistent, and inventing a third one here would imply rho is
a third kind of quantity. The caption states the choice.

**Edge STYLE, not edge presence, carries the FDR.** Solid means the edge clears
the pre-registered empirical FDR; dashed means it does not. Drawing only the
survivors would make an empty network for an outcome Gate 4 explicitly licenses
— "no LR pair exceeded chance expectation at n=13" — and a blank panel looks
like a broken rule rather than a finding. The top edges by |rho| are always
drawn; whether any of them is real is what the style says.

**Complex interactions are not on this figure.** They live in a declared
exploratory table with no FDR of any kind (ADR 0021 §4), so they have no style
to draw.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

# Reused verbatim from pca_landscape.py / checkpoint_dotplot.py so a compartment
# is the same colour everywhere in the project.
SITE_COLOUR = {"lung": "#C1666B", "brain": "#4F6D7A", "lymph_node": "#D4B483"}
TUMOUR_NODE = "#8C6A5D"
IMMUNE_NODE = "#4F6D7A"
CMAP = "RdBu_r"

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    alpha = snakemake.params.null_calibration["fdr_alpha"]
    top_n = snakemake.params.top_n_edges
    adjacencies = snakemake.params.adjacencies

    emit("P4-T7 — inferred crosstalk network")
    emit("=" * 70)
    emit()
    emit("  NEVER 'colocalisation' and never 'spatially adjacent'. There are no")
    emit("  coordinates. An edge is a CROSS-PATIENT CORRELATION between two")
    emit("  compartments, inferred from a database — not a measured proximity.")
    emit()
    emit(f"  edge colour   RdBu_r on rho — a signed effect, the quantity class")
    emit(f"                P2-T7 established that language for. Consistent")
    emit(f"                reuse, not a borrow; PRGn stays reserved for a")
    emit(f"                background ratio (ADR 0016 §5).")
    emit(f"  edge style    solid = clears empirical FDR {alpha}; dashed = does not")
    emit(f"  edges drawn   top {top_n} by |rho| per site, survivors or not")
    emit()

    fdr = pd.read_csv(snakemake.input.fdr, sep="\t", float_precision="round_trip")
    emit(f"primary direction-rows: {len(fdr)}")

    # Explicit gridspec with a dedicated colourbar column, NOT
    # fig.colorbar(ax=axes) + tight_layout. Attaching a colourbar to the axes
    # steals width from them, tight_layout then warns that it cannot handle the
    # result, and the right-hand node labels get clipped — CEACAM1 and CEACAM5
    # both rendered as "CEACAM" on the first pass. This is the `main | key |
    # colourbar` arrangement checkpoint_dotplot.py already uses, for the same
    # reason.
    n_panels = len(adjacencies)
    fig = plt.figure(figsize=(7.8 * n_panels, 8.4))
    gs = fig.add_gridspec(
        1, n_panels + 1,
        width_ratios=[1] * n_panels + [0.035],
        # right leaves room for the colourbar LABEL, not just the bar —
        # at 0.97 the label rendered off-canvas. top leaves the panel
        # titles clear of the four-line suptitle.
        wspace=0.30, left=0.02, right=0.93, top=0.80, bottom=0.03,
    )
    axes = np.array([fig.add_subplot(gs[0, i]) for i in range(n_panels)])
    cax = fig.add_subplot(gs[0, n_panels])
    norm = Normalize(vmin=-1.0, vmax=1.0)
    cmap = plt.get_cmap(CMAP)

    stats = {}
    for ax, site in zip(axes, sorted(adjacencies)):
        sel = (
            fdr[fdr["site"] == site]
            .sort_values("abs_rho", ascending=False)
            .head(top_n)
            .copy()
        )
        n_pat = int(sel["n"].iloc[0]) if len(sel) else 0
        n_clear = int(sel["clears_empirical_fdr"].sum()) if len(sel) else 0
        emit(f"{site}: {len(sel)} edges drawn, {n_clear} clear the FDR, "
             f"n = {n_pat} patients")

        G = nx.DiGraph()
        for r in sel.itertuples(index=False):
            # Node identity is (gene, side): the same gene measured in the two
            # compartments is two different measurements, and merging them
            # would draw an edge that crosses a compartment it never touched.
            t_side = (
                r.ligand if r.ligand_aoi_code == adjacencies[site]["tumour"]
                else r.receptor
            )
            i_side = (
                r.receptor if r.ligand_aoi_code == adjacencies[site]["tumour"]
                else r.ligand
            )
            u, v = f"{t_side}\n({adjacencies[site]['tumour']})", \
                   f"{i_side}\n({adjacencies[site]['immune']})"
            G.add_node(u, side="tumour")
            G.add_node(v, side="immune")
            G.add_edge(u, v, rho=float(r.rho),
                       clears=bool(r.clears_empirical_fdr))

        left = sorted(n for n, d in G.nodes(data=True) if d["side"] == "tumour")
        right = sorted(n for n, d in G.nodes(data=True) if d["side"] == "immune")
        pos = {}
        for i, n in enumerate(left):
            pos[n] = (0.0, -i * (max(len(right), 1) / max(len(left), 1)))
        for i, n in enumerate(right):
            pos[n] = (1.0, -i * (max(len(left), 1) / max(len(right), 1)))

        for u, v, d in G.edges(data=True):
            ax.annotate(
                "",
                xy=pos[v], xytext=pos[u],
                arrowprops=dict(
                    arrowstyle="-|>", color=cmap(norm(d["rho"])),
                    lw=1.0 + 2.0 * abs(d["rho"]),
                    linestyle="-" if d["clears"] else (0, (3, 3)),
                    alpha=0.95 if d["clears"] else 0.42,
                    shrinkA=16, shrinkB=16,
                ),
            )
        for nodes, colour in ((left, TUMOUR_NODE), (right, IMMUNE_NODE)):
            for n in nodes:
                x, y = pos[n]
                ax.scatter([x], [y], s=210, color=colour, zorder=3,
                           edgecolors="white", linewidths=1.1)
                ax.text(x + (-0.045 if colour == TUMOUR_NODE else 0.045), y,
                        n.split("\n")[0], fontsize=7.4, zorder=4,
                        ha="right" if colour == TUMOUR_NODE else "left",
                        va="center")

        ax.set_title(
            f"{site}   {adjacencies[site]['tumour']} <-> "
            f"{adjacencies[site]['immune']}\n"
            f"{adjacencies[site]['immune']} n = {n_pat} patients\n"
            f"{len(sel)} strongest |rho|; {n_clear} clear the empirical FDR",
            fontsize=9.5, color=SITE_COLOUR.get(site, "black"),
        )
        # Widen for the LABELS, not for the nodes. The nodes sit at x = 0 and
        # x = 1; the text runs outward from them, and the longest symbol decides
        # how far. P3-T5's title clipped at both canvas edges on first render
        # despite the plan warning about it, and this figure clipped CEACAM1 and
        # CEACAM5 to "CEACAM" for the same reason: a margin guessed rather than
        # measured. Scale it by the longest symbol actually drawn.
        longest = max(
            [len(n.split("\n")[0]) for n in list(left) + list(right)] or [8]
        )
        margin = 0.06 + 0.055 * longest
        ax.set_xlim(-margin, 1.0 + margin)
        ax.axis("off")
        stats[site] = {
            "n_edges_drawn": int(len(sel)),
            "n_clears_fdr": n_clear,
            "n_patients": n_pat,
            "n_nodes": int(G.number_of_nodes()),
        }

    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("Spearman rho (signed effect — RdBu_r, P2-T7)", fontsize=8.5)
    cbar.ax.tick_params(labelsize=7.5)

    # KEEP EVERY TITLE LINE SHORT. P3-T5's clipped at both canvas edges on first
    # render despite the plan warning about it; verify by opening the PNG.
    fig.suptitle(
        "P4-T7  inferred crosstalk between adjacent compartments\n"
        "left = tumour compartment, right = immune; arrow = ligand to receptor\n"
        "solid clears the empirical FDR, dashed does not — exploratory, ADR 0014\n"
        "NOT colocalisation: this assay carries no coordinates",
        fontsize=10, y=0.985, va="top",
    )
    # No tight_layout: the gridspec above already fixes the geometry, and
    # tight_layout would undo the margins the labels need.
    fig.savefig(snakemake.output.figure, dpi=200)
    plt.close(fig)
    emit()
    emit(f"wrote {snakemake.output.figure}")

    summary = {
        "task": "P4-T7",
        "aim": "A5",
        "exploratory": True,
        "exploratory_adr": "ADR 0014",
        "fdr_alpha": alpha,
        "top_n_edges_per_site": top_n,
        "by_site": stats,
        "edge_colour": (
            "RdBu_r on rho. A deliberate REUSE, not a borrow: rho is a signed "
            "effect on [-1, 1], the quantity class P2-T7 established RdBu_r "
            "for. PRGn centred on the detection multiple stays reserved for a "
            "background ratio (ADR 0016 §5)."
        ),
        "edge_style": (
            "Solid clears the empirical FDR, dashed does not. Style rather than "
            "presence, because drawing only survivors would make a blank panel "
            "for an outcome Gate 4 explicitly licenses, and a blank panel reads "
            "as a broken rule rather than a finding."
        ),
        "node_identity": (
            "(gene, side). The same gene measured in the two compartments is "
            "two measurements; merging them would draw an edge across a "
            "compartment it never touched."
        ),
        "excluded": (
            "Complex interactions are not on this figure — they carry no FDR of "
            "any kind (ADR 0021 §4), so they have no style to draw."
        ),
        "wording": (
            "inferred crosstalk between adjacent compartments. NEVER "
            "colocalisation, never spatially adjacent — this assay carries no "
            "coordinates (hard constraint 6)."
        ),
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    emit(f"wrote {snakemake.output.summary}")
