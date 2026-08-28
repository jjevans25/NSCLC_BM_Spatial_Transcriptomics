"""P2-T5 figure — stacked cell composition, one panel per reference arm.

Owner task: P2-T5. Driven by rule p2t5_decon_figure.

Three panels, and the split is the argument:

  A  safeTME across BOTH sites. One reference, so this is the only panel whose
     brain-vs-lung comparison means anything.
  B  Lung_HCA on TIME-L      } PROJECT_PLAN §6's two-matrix approach. B and C
  C  Brain_Darmanis on TIME-B} ARE NOT COMPARABLE WITH EACH OTHER.

Drawing B and C as one panel with a shared legend would invite exactly the
comparison the assay cannot support, so they are separated and the figure says
why on its face: Lung_HCA resolves 9 lymphoid types, Brain_Darmanis resolves 0,
so a lymphoid difference between B and C would be the reference rather than the
biology. Those two counts are read from decon_summary.json rather than
hardcoded, so the caption cannot drift from what was actually fitted.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    aoi_codes = list(snakemake.params.aoi_codes)
    comp = pd.read_csv(snakemake.input.table, sep="\t")
    with open(snakemake.input.summary, encoding="utf-8") as handle:
        summary = json.load(handle)

    refs = summary["references"]
    lym = {k: v["n_lymphoid"] for k, v in refs.items()}
    emit(f"lymphoid type counts per reference: {lym}")

    panels = [("safetme", None, "A")]
    panels += [("tissue", c, l) for c, l in zip(aoi_codes, "BC")]

    fig, axes = plt.subplots(
        1, 3, figsize=(17.5, 6.4),
        gridspec_kw={"width_ratios": [23, 15, 8], "wspace": 0.02},
    )

    for ax, (arm, code, letter) in zip(axes, panels):
        d = comp.loc[comp["arm"] == arm]
        if code is not None:
            d = d.loc[d["aoi_code"] == code]

        order = (
            d[["aoi_label", "aoi_code", "patient_id"]]
            .drop_duplicates()
            .assign(_k=lambda x: x["aoi_code"].map({c: i for i, c in enumerate(aoi_codes)}))
            .sort_values(["_k", "patient_id", "aoi_label"])
        )
        cols = order["aoi_label"].tolist()

        mat = d.pivot_table(index="cell_type", columns="aoi_label", values="proportion")
        # Largest mean share at the bottom of the stack, so the eye compares the
        # dominant populations against a common baseline rather than a ragged one.
        mat = mat.loc[mat.mean(axis=1).sort_values(ascending=False).index, cols]

        cmap = plt.get_cmap("tab20")
        colours = [cmap(i % 20) for i in range(mat.shape[0])]

        bottom = np.zeros(len(cols))
        for (name, row), colour in zip(mat.iterrows(), colours):
            ax.bar(range(len(cols)), row.to_numpy(), bottom=bottom, width=0.9,
                   label=name, color=colour, linewidth=0)
            bottom += row.to_numpy()

        ax.set_xlim(-0.5, len(cols) - 0.5)
        ax.set_ylim(0, 1)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(order["patient_id"].tolist(), rotation=90, fontsize=7)
        ax.set_xlabel("patient", fontsize=9)

        if arm == "safetme":
            ref = refs["safetme"]["reference"]
            n = order["aoi_code"].value_counts().to_dict()
            title = (
                f"{letter}. {ref} — both sites, ONE reference\n"
                f"TIME-L n = {n[aoi_codes[0]]}, TIME-B n = {n[aoi_codes[1]]}"
                "  ·  brain vs lung IS comparable here"
            )
            b = np.cumsum([n[c] for c in aoi_codes])[:-1]
            for x in b:
                ax.axvline(x - 0.5, color="black", lw=2.0)
            ax.set_ylabel("proportion of deconvolved abundance", fontsize=10)
        else:
            meta = refs[f"tissue_{code}"]
            # Wrapped onto three lines rather than two: the single-line form
            # ran past the right edge of panel C and clipped the very words
            # that stop the panels being compared.
            title = (
                f"{letter}. {meta['reference']} — {code} only (n = {len(cols)})\n"
                f"{meta['n_types']} types, {meta['n_lymphoid']} lymphoid\n"
                "NOT comparable with the other tissue panel"
            )
            ax.set_yticklabels([])

        ax.set_title(title, fontsize=9.5, loc="left")
        ax.legend(fontsize=6.2, loc="upper left", bbox_to_anchor=(0, -0.16),
                  ncol=2, frameon=False, handlelength=1.0, columnspacing=0.8)

    lung_lym = refs[f"tissue_{aoi_codes[0]}"]["n_lymphoid"]
    brain_lym = refs[f"tissue_{aoi_codes[1]}"]["n_lymphoid"]
    fig.suptitle(
        "P2-T5 cell composition of the TIME AOIs, renormalised within each AOI"
        "   ·   TIME-B n = 8\n"
        f"Panels B and C use different references ({lung_lym} vs {brain_lym} "
        "lymphoid cell types) and must not be compared with each other — "
        "a lymphoid gap between them would be the reference, not the biology",
        fontsize=10.5, y=0.995,
    )
    fig.subplots_adjust(top=0.80, bottom=0.42, left=0.055, right=0.99)
    fig.savefig(snakemake.output.figure, dpi=200)
    plt.close(fig)
    emit(f"wrote {snakemake.output.figure}")
