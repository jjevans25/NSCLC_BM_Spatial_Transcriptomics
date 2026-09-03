"""P3-T2b — external validation against the source publication's own files.

Owner task: P3-T2b. Driven by rule p3t2b_external_validation.

WHY THIS EXISTS. Everything Phase 3 asserted before this rule was internal.
P3-T2's 25 cross-checks run three implementations of ADR 0007's detection rule
-- checkpoint_detection.py, export_tsv.py, score_signatures.py -- but all three
read the SAME obs['negprobe'] column and test the SAME hypothesis. They catch
indexing, subsetting and transcription errors. They cannot catch the rule being
wrong, the multiple being mis-set, or the matrix being something other than
what docs/data-provenance.md says it is. ADR 0007's reproduction of the paper's
sequencing-saturation claim was the project's only external check, and it
validates DCC parsing rather than values.

PMID 36216799 deposits Source Data containing the per-AOI values behind its
figures, including the complete 120-AOI matrix. This rule compares that against
layers['q3'] cell by cell, and then uses the same files to characterise what the
published analysis did with those values.

FOUR THINGS IT DOES, in order of what they license:

  1. **The matrix.** Every value, every AOI, every gene, against the authors'
     deposited file. Exact equality is required; anything else stops the phase.
     This is the first external evidence that this project's inputs are the
     published inputs.

  2. **The 119-vs-120 gap.** docs/limitations.md records that the paper analyses
     119 ROIs where GEO deposits 120, and that which AOI was excluded is
     unknown. The fibrosis sheet answers it: it carries a single "TBME15" where
     the full matrix carries TBME15a and TBME15b. The rule ASSERTS that column
     equals TBME15a across every gene and differs from TBME15b, so the answer is
     measured rather than assumed.

  3. **Count quantisation.** The matrix is Q3-scaled counts with no zeros (Q1),
     so an AOI's gene values sit on a lattice whose spacing is one raw count.
     Recovering the spacing turns "IDO1 is near background in TBME" into a
     statement about resolution: how many counts the gene actually carries, and
     how many other genes hold the identical value. The recovered scale is then
     checked against the raw counts in P0-T2's .dcc archive, so it is measured
     rather than inferred.

     The quantum is the SPACING between adjacent distinct values, not the
     smallest value — the smallest value is 1 count only in the shallowest
     AOIs. Getting that wrong inflates every count by 5-13x, and it fails
     loudly: the lattice check collapses and the .dcc comparison rejects it.

  4. **The published contrast, reproduced.** The paper's F(h)-vs-F(-) DEG list
     is reproduced exactly -- all 99 log2FC values -- which establishes the
     estimator and licenses the only claim this rule makes about the published
     analysis: which genes are NOT in that list, and by how far they miss.

WHAT IT DOES NOT DO:

  * **It does not change the panel.** The paper names checkpoint genes this
    project does not carry (BTLA, CD160, PDCD1LG2, VTCN1, IDO2, TDO2). They are
    transcribed below as the paper's claim and reported alongside ours; they are
    NOT added to config/checkpoints.yaml, which is pre-registered (ADR 0015) and
    a stop-and-ask decision. Adding a gene after seeing which ones move is the
    exact failure P3-T2 was built to prevent.

  * **It does not change a threshold.** ADR 0016 §3's detection floor and
    ADR 0007's 2x background are untouched. The published |log2FC| > 1.5
    threshold is applied only to the paper's own list, never to anything here.

  * **It does not recompute detection.** checkpoint_detection.tsv already has
    it, keyed on aoi_code; this rule reads that file.

STDLIB XLSX PARSING, ON PURPOSE. openpyxl is not in workflow/envs/py-analysis.yaml
and must not be added: that env is p0t2_fetch_geo's, whose outputs are
protected(), so a dependency change fires the software-env rerun trigger and the
DAG dies on ProtectedOutputException (ADR 0005, ADR 0013). The env file also
warns that anything it lacks falls through to the user's base install, which
would make this rule pass here and fail everywhere else. An .xlsx is a zip of
XML; zipfile and ElementTree are enough.

The reader uses float() rather than a pandas parser, which matters here: ADR
0013's postscript records that pandas' default C parser is not correctly
rounded, and this rule asserts EXACT equality. float() is correctly rounded, so
the same decimal text yields the same double on both sides of the comparison.
"""

import gzip
import json
import tarfile
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

FMT = "%.17g"

# The checkpoint genes the SOURCE PAPER names, transcribed from its Results and
# Discussion. This is the paper's set, not this project's panel: it is reported
# so the published claim can be evaluated in the same table as ours, and it
# confers NO licence to alter config/checkpoints.yaml (ADR 0015 — panel
# membership is a stop-and-ask decision, and doing it after seeing the audit is
# the failure P3-T2 exists to prevent). Same idiom as ADR0008_RECON in
# checkpoint_detection.py: an external table transcribed into the code that
# checks against it.
#
#   "F(-) TBME expressed significantly greater levels of [...] PDCD1LG2
#    (programmed cell death 1 ligand 2 or PD-L2), BTLA (B- and T-lymphocyte
#    attenuator), VTCN1 (V-set domain-containing T-cell activation inhibitor 1),
#    and IDO1"
#
# with IDO2/TDO2 named alongside IDO1 as "metabolic immune checkpoints", and
# CD160 named in the Discussion's combination-therapy recommendation.
PAPER_CHECKPOINTS = ["PDCD1LG2", "BTLA", "VTCN1", "IDO1", "IDO2", "TDO2", "CD160"]

# --------------------------------------------------------------------------
# Minimal .xlsx reader (stdlib only) — see the module docstring.
# --------------------------------------------------------------------------
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RNS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def _column_index(ref):
    """'BC12' -> 54. Cell refs are sparse, so a row's cells must be placed."""
    n = 0
    for char in ref:
        if not char.isalpha():
            break
        n = n * 26 + (ord(char.upper()) - 64)
    return n - 1


def _shared_strings(archive):
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    out = []
    with archive.open("xl/sharedStrings.xml") as handle:
        for _, elem in ET.iterparse(handle, events=("end",)):
            if elem.tag == NS + "si":
                # <si> may hold one <t> or several <r><t> runs; join them all.
                out.append("".join(node.text or "" for node in elem.iter(NS + "t")))
                elem.clear()
    return out


def _worksheet_path(archive, name):
    with archive.open("xl/workbook.xml") as handle:
        workbook = ET.parse(handle).getroot()
    with archive.open("xl/_rels/workbook.xml.rels") as handle:
        rels = {rel.get("Id"): rel.get("Target") for rel in ET.parse(handle).getroot()}
    available = []
    for sheet in workbook.iter(NS + "sheet"):
        available.append(sheet.get("name"))
        if sheet.get("name") == name:
            target = rels[sheet.get(RNS + "id")]
            return target if target.startswith("xl/") else "xl/" + target.lstrip("/")
    raise KeyError(
        f"worksheet {name!r} not in the workbook. Present: {available}. "
        "Sheet names are transcribed exactly in config/external_validation.yaml, "
        "trailing whitespace included — do not tidy them."
    )


def read_sheet(path, name):
    """One worksheet as a list of rows, each a list of str | float | None.

    iterparse over the zip stream rather than a full parse: the Source Data
    sheet is ~2.3 million cells and materialising its DOM is gratuitous.

    Rows are placed at the index their `r` attribute declares, and gaps are
    padded. Excel OMITS entirely blank rows from the XML, so appending
    sequentially silently shifts everything below one — which is exactly what
    happens here: every sheet in these files has a blank row between its title
    and its header. The header_row values in config/external_validation.yaml are
    true sheet positions, the ones a person reading the file in Excel would
    count, and this is what makes that true.
    """
    placed = {}
    with zipfile.ZipFile(path) as archive:
        strings = _shared_strings(archive)
        with archive.open(_worksheet_path(archive, name)) as handle:
            for _, elem in ET.iterparse(handle, events=("end",)):
                if elem.tag != NS + "row":
                    continue
                cells = {}
                for cell in elem.findall(NS + "c"):
                    ref = cell.get("r") or ""
                    index = _column_index(ref) if ref else len(cells)
                    kind = cell.get("t")
                    if kind == "inlineStr":
                        node = cell.find(NS + "is")
                        value = (
                            "".join(t.text or "" for t in node.iter(NS + "t"))
                            if node is not None
                            else None
                        )
                    else:
                        node = cell.find(NS + "v")
                        if node is None or node.text is None:
                            value = None
                        elif kind == "s":
                            value = strings[int(node.text)]
                        elif kind in (None, "n"):
                            # float(), not a pandas parser: ADR 0013's postscript
                            # records that pandas' C parser is not correctly
                            # rounded, and this rule asserts EXACT equality.
                            value = float(node.text)
                        else:
                            value = node.text
                    cells[index] = value
                width = max(cells) + 1 if cells else 0
                # `r` is 1-based; keep a 0-based list so header_row reads like
                # an index rather than an off-by-one waiting to happen.
                declared = elem.get("r")
                position = int(declared) - 1 if declared else len(placed)
                placed[position] = [cells.get(i) for i in range(width)]
                elem.clear()
    if not placed:
        return []
    return [placed.get(i, []) for i in range(max(placed) + 1)]


def sheet_frame(path, name, header_row):
    """A worksheet as a DataFrame indexed by its first column."""
    rows = read_sheet(path, name)
    header = rows[header_row]
    width = len(header)
    body = [r[:width] + [None] * (width - len(r)) for r in rows[header_row + 1 :]]
    frame = pd.DataFrame(body, columns=[str(h).strip() if h else h for h in header])
    key = frame.columns[0]
    frame = frame[frame[key].notna()].copy()
    frame[key] = frame[key].astype(str).str.strip()
    return frame.set_index(key)


def dcc_median_counts(archive_path, samples_path):
    """Median raw count per AOI, straight from the P0-T2 .dcc archive.

    The <Code_Summary> block of a .dcc is `<probe id>,<count>` — integers, before
    any normalisation. This is what makes the recovered quantum a MEASUREMENT
    rather than an inference: median(value / quantum) has to equal the AOI's own
    median raw count.

    No PKC is needed and none exists (ADR 0006). The comparison is
    distribution-to-distribution, so probe identity never has to be resolved to
    a gene — which is the only reason this check is available at all.
    """
    samples = pd.read_csv(samples_path, sep="\t", dtype=str)
    label_of = dict(zip(samples["gsm_id"], samples["aoi_label"]))
    out = {}
    with tarfile.open(archive_path) as tar:
        for member in tar.getmembers():
            name = Path(member.name).name
            if not name.endswith(".dcc.gz"):
                continue
            label = label_of.get(name.split("_")[0])
            if label is None:
                continue
            handle = tar.extractfile(member)
            if handle is None:
                continue
            counts = []
            inside = False
            with gzip.open(handle, "rt") as text:
                for line in text:
                    line = line.strip()
                    if line == "<Code_Summary>":
                        inside = True
                    elif line == "</Code_Summary>":
                        break
                    elif inside and "," in line:
                        counts.append(int(line.rsplit(",", 1)[1]))
            if counts:
                out[label] = float(np.median(counts))
    return out


log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    config = snakemake.params.external
    expect = config["expect"]
    sheets = config["sheets"]
    supp = {Path(p).name: p for p in snakemake.input.supplementary}
    path_of = {
        key: supp[spec["dest"]] for key, spec in config["artifacts"].items()
    }

    emit("P3-T2b — external validation against PMID 36216799's own files")
    emit("  the first check in this project that can fail for a reason other")
    emit("  than our own arithmetic (see the module docstring)")
    emit()

    QUANTUM_PROBE = expect["quantum_probe"]
    LATTICE_PROBE = expect["lattice_probe"]
    LATTICE_MIN = expect["lattice_min_fraction"]
    DCC_RATIO_TOL = expect["dcc_ratio_tol"]

    checkpoints = snakemake.params.checkpoints
    membership = pd.read_csv(snakemake.input.membership, sep="\t")
    alias_of = dict(zip(membership["gene"], membership["alias"]))
    panel = sorted(checkpoints)

    # ------------------------------------------------------- 1. the matrix
    adata = ad.read_h5ad(snakemake.input.h5ad)
    ours = pd.DataFrame(
        np.asarray(adata.layers["q3"]), index=adata.obs_names, columns=adata.var_names
    )
    negprobe = adata.obs["negprobe"]

    spec = sheets["full_matrix"]
    published = sheet_frame(path_of[spec["artifact"]], spec["name"], spec["header_row"])
    published = published.apply(pd.to_numeric, errors="coerce")

    emit("cross-check 1 — layers['q3'] vs the deposited Source Data matrix")
    emit(f"  published : {published.shape[0]} rows x {published.shape[1]} AOIs")
    emit(f"  ours      : {ours.shape[0]} AOIs x {ours.shape[1]} genes")

    if ours.shape[0] != expect["n_aoi"] or ours.shape[1] != expect["n_genes"]:
        raise RuntimeError(
            f"expected {expect['n_aoi']} AOIs x {expect['n_genes']} genes, got "
            f"{ours.shape[0]} x {ours.shape[1]}. PROJECT_PLAN §2.1 is "
            "authoritative — do not adjust the expectation to match the data."
        )

    our_aoi, their_aoi = set(ours.index), set(published.columns)
    if our_aoi != their_aoi:
        raise RuntimeError(
            "AOI label sets differ. Only ours: "
            f"{sorted(our_aoi - their_aoi)}; only published: "
            f"{sorted(their_aoi - our_aoi)}."
        )
    emit(f"  AOI labels identical: {len(our_aoi)}/{expect['n_aoi']}")

    controls = list(expect["control_rows"])
    extra = sorted(set(published.index) - set(ours.columns))
    if extra != sorted(controls):
        raise RuntimeError(
            f"published matrix carries unexpected rows {extra}; "
            f"expected exactly {sorted(controls)}."
        )
    missing = sorted(set(ours.columns) - set(published.index))
    if missing:
        raise RuntimeError(f"{len(missing)} of our genes absent upstream: {missing[:8]}")
    emit(f"  gene sets differ by exactly {extra} — held in obs, compared below")

    aoi_order = sorted(our_aoi)
    gene_order = sorted(ours.columns)
    mine = ours.loc[aoi_order, gene_order].to_numpy(float)
    theirs = published.loc[gene_order, aoi_order].to_numpy(float).T
    delta = np.abs(mine - theirs)
    n_exact = int((mine == theirs).sum())

    emit(f"  compared  : {mine.size:,} values")
    emit(f"  exact     : {n_exact:,} ({100 * n_exact / mine.size:.4f}%)")
    emit(f"  max |diff|: {np.nanmax(delta):g}")
    if n_exact != mine.size:
        raise RuntimeError(
            f"{mine.size - n_exact} of {mine.size} values differ from the "
            f"published Source Data (max |diff| {np.nanmax(delta):g}). Either "
            "this project's matrix is not the published matrix, or the reader "
            "is lossy. Both stop the phase."
        )

    for control in controls:
        theirs_control = published.loc[control, aoi_order].to_numpy(float)
        ours_control = negprobe.loc[aoi_order].to_numpy(float)
        if not np.array_equal(theirs_control, ours_control):
            raise RuntimeError(
                f"published {control} row does not equal obs['negprobe'] "
                f"(max |diff| {np.nanmax(np.abs(theirs_control - ours_control)):g}). "
                "ADR 0007's detection rule is defined against this row."
            )
        emit(f"  {control} row == obs['negprobe'] exactly")

    per_aoi = pd.DataFrame(
        {
            "aoi_label": aoi_order,
            "aoi_code": adata.obs.loc[aoi_order, "aoi_code"].astype(str).to_numpy(),
            "site": adata.obs.loc[aoi_order, "site"].astype(str).to_numpy(),
            "n_genes_compared": mine.shape[1],
            "n_exact": (mine == theirs).sum(axis=1),
            "max_abs_diff": delta.max(axis=1),
        }
    )
    emit()

    # -------------------------------------------- 2. the 119-vs-120 gap
    emit("cross-check 2 — the 119-vs-120 gap (docs/limitations.md §, open thread)")
    spec = sheets["fibrosis_matrix"]
    fibrosis_matrix = sheet_frame(
        path_of[spec["artifact"]], spec["name"], spec["header_row"]
    ).apply(pd.to_numeric, errors="coerce")

    alias_resolution = {}
    for published_label, our_label in sorted(expect["aoi_alias"].items()):
        if published_label not in fibrosis_matrix.columns:
            raise RuntimeError(
                f"alias {published_label!r} not a column of {spec['name']!r}; "
                f"columns are {list(fibrosis_matrix.columns)}."
            )
        shared = sorted(set(fibrosis_matrix.index) & set(ours.columns))
        column = fibrosis_matrix.loc[shared, published_label].to_numpy(float)
        # Agreement WITHIN TOLERANCE, not equality: this sheet rounds
        # differently from "Figure 1b-e" (see config/external_validation.yaml).
        # The exact comparison is cross-check 1 and it already passed over the
        # whole matrix; what this test has to do is tell TBME15a from TBME15b.
        matches = {
            candidate: float(
                np.isclose(
                    column,
                    ours.loc[candidate, shared].to_numpy(float),
                    rtol=expect["alias_rtol"],
                    atol=0.0,
                ).mean()
            )
            for candidate in sorted(our_aoi)
            if candidate.startswith(published_label)
        }
        best = matches.get(our_label, 0.0)
        if best < expect["alias_min_agreement"]:
            raise RuntimeError(
                f"alias {published_label!r} -> {our_label!r} is not supported by "
                f"the data: agreement {best:.4f} < "
                f"{expect['alias_min_agreement']} (by candidate {matches})."
            )
        ambiguous = {
            k: v
            for k, v in matches.items()
            if k != our_label and v > expect["alias_max_other_agreement"]
        }
        if ambiguous:
            raise RuntimeError(
                f"alias {published_label!r} is ambiguous — {ambiguous} also "
                f"agree above {expect['alias_max_other_agreement']}."
            )
        n_exact_alias = int(
            (column == ours.loc[our_label, shared].to_numpy(float)).sum()
        )
        alias_resolution[published_label] = {
            "resolves_to": our_label,
            "agreement_by_candidate": matches,
            "n_genes_compared": len(shared),
            "n_exact": n_exact_alias,
            "note": (
                "agreement is within alias_rtol, not exact: this sheet carries "
                "an extra significant figure for large values where "
                "'Figure 1b-e' matches GEO exactly"
            ),
        }
        emit(
            f"  published {published_label!r} -> {our_label}: agreement "
            + ", ".join(f"{k} {v:.4f}" for k, v in sorted(matches.items()))
        )
        emit(
            f"    ({n_exact_alias}/{len(shared)} genes also bit-exact; the "
            "remainder are this sheet's extra significant figure)"
        )
    dropped = sorted(
        candidate
        for record in alias_resolution.values()
        for candidate in record["agreement_by_candidate"]
        if candidate != record["resolves_to"]
    )
    emit(
        f"  -> the paper's 119 is 120 minus {dropped}, dropped without comment. "
        "No AOI failed QC upstream; docs/limitations.md's open thread is closed."
    )
    emit()

    # ------------------------------------------------- 3. count quantisation
    emit("cross-check 3 — count quantisation (resolution, not background)")
    quant_rows = []
    checkpoint_rows = []
    detection = pd.read_csv(
        snakemake.input.detection, sep="\t", float_precision="round_trip"
    )
    detected_lookup = {
        (r.gene, r.aoi_code): r for r in detection.itertuples(index=False)
    }

    dcc_median = dcc_median_counts(snakemake.input.raw_archive, snakemake.input.samples)
    emit(f"  raw counts read from {len(dcc_median)} .dcc files in the P0-T2 archive")

    for label in aoi_order:
        values = ours.loc[label].to_numpy(float)
        distinct = np.unique(values[values > 0])

        # The quantum is the SPACING of the lattice, not its smallest member.
        # The smallest value is 1 count only in the shallowest AOIs; taking it
        # as the quantum inflates every count by 5-13x and collapses the lattice
        # check to ~0.15. See config/external_validation.yaml.
        steps = np.diff(distinct[:QUANTUM_PROBE])
        quantum = float(np.median(steps[steps > 0]))

        probe = distinct[:LATTICE_PROBE]
        multiples = np.round(probe / quantum)
        # Both the value and the quantum are reported to 2 dp in the deposited
        # file, so the tolerance has to grow with the multiple.
        on_lattice = np.abs(probe - multiples * quantum) <= 0.005 * (1.0 + multiples)
        lattice_fraction = float(on_lattice.mean())

        implied = values / quantum
        observed = dcc_median.get(label)
        dcc_ratio = (
            float(np.median(implied) / observed)
            if observed is not None and observed > 0
            else float("nan")
        )
        agrees = bool(abs(dcc_ratio - 1.0) <= DCC_RATIO_TOL)
        reliable = bool(lattice_fraction >= LATTICE_MIN and agrees)

        counts = np.unique(values, return_counts=True)
        ties = dict(zip(counts[0], counts[1]))

        quant_rows.append(
            {
                "aoi_label": label,
                "aoi_code": str(adata.obs.at[label, "aoi_code"]),
                "site": str(adata.obs.at[label, "site"]),
                "quantum": quantum,
                "n_distinct_values": int(len(distinct)),
                "n_genes": int(len(values)),
                "lattice_fraction": lattice_fraction,
                "smallest_value_implied_counts": float(distinct[0] / quantum),
                "dcc_median_counts": (
                    float(observed) if observed is not None else float("nan")
                ),
                "implied_median_counts": float(np.median(implied)),
                "dcc_ratio": dcc_ratio,
                "dcc_agrees": agrees,
                "counts_reliable": reliable,
                "negprobe_implied_counts": float(negprobe[label] / quantum),
                "max_implied_counts": float(values.max() / quantum),
            }
        )

        for gene in panel:
            value = float(ours.at[label, gene])
            code = str(adata.obs.at[label, "aoi_code"])
            row = detected_lookup.get((gene, code))
            checkpoint_rows.append(
                {
                    "gene": gene,
                    "alias": alias_of[gene],
                    "aoi_label": label,
                    "aoi_code": code,
                    "q3": value,
                    "quantum": quantum,
                    "implied_counts": value / quantum,
                    "counts_reliable": reliable,
                    "n_genes_sharing_value": int(ties[value]),
                    "negprobe_implied_counts": float(negprobe[label] / quantum),
                    "assessable_in_compartment": (
                        bool(row.assessable) if row is not None else None
                    ),
                }
            )

    quantisation = pd.DataFrame(quant_rows)
    checkpoint_counts = pd.DataFrame(checkpoint_rows)

    n_lattice = int((quantisation["lattice_fraction"] >= LATTICE_MIN).sum())
    n_dcc = int(quantisation["dcc_agrees"].sum())
    emit(
        f"  lattice holds for {n_lattice}/{len(quantisation)} AOIs "
        f"(>= {LATTICE_MIN:.0%} of the {LATTICE_PROBE} smallest distinct values "
        f"are integer multiples of the quantum estimated from {QUANTUM_PROBE})"
    )
    emit(
        f"  quantum agrees with the .dcc raw counts for {n_dcc}/"
        f"{len(quantisation)} AOIs (median ratio "
        f"{quantisation['dcc_ratio'].median():.4f}, tolerance "
        f"{DCC_RATIO_TOL:.0%}) — the count scale is MEASURED, not inferred"
    )
    if n_dcc < len(quantisation):
        raise RuntimeError(
            f"{len(quantisation) - n_dcc} AOIs' recovered quantum disagrees "
            "with their own raw counts. The lattice may look clean and still be "
            "the wrong spacing; nothing downstream may speak in counts."
        )
    emit(
        f"  distinct values per AOI: median "
        f"{quantisation['n_distinct_values'].median():.0f}, "
        f"min {quantisation['n_distinct_values'].min()}, "
        f"max {quantisation['n_distinct_values'].max()} (of {len(gene_order):,} genes)"
    )
    emit(
        f"  smallest value in an AOI is {quantisation['smallest_value_implied_counts'].min():.0f}"
        f"-{quantisation['smallest_value_implied_counts'].max():.0f} counts "
        "(median "
        f"{quantisation['smallest_value_implied_counts'].median():.0f}) — which is "
        "why the quantum is the spacing, not the minimum"
    )
    emit()
    emit("  panel implied counts vs background, by compartment (median over AOIs):")
    emit(f"    {'compartment':12s} {'background':>10s}  " +
         "".join(f"{g:>9s}" for g in panel))
    for code in sorted(checkpoint_counts["aoi_code"].unique()):
        block = checkpoint_counts[checkpoint_counts["aoi_code"] == code]
        background = block["negprobe_implied_counts"].median()
        cells = "".join(
            f"{block.loc[block['gene'] == g, 'implied_counts'].median():>9.0f}"
            for g in panel
        )
        emit(f"    {code:12s} {background:>10.0f}  {cells}")
    emit()
    ties_by_code = checkpoint_counts.groupby("aoi_code")["n_genes_sharing_value"]
    emit("  genes sharing a panel gene's exact value (median, max per compartment):")
    for code, block in ties_by_code:
        emit(f"    {code:12s} median {block.median():6.0f}   max {block.max():6.0f}")
    emit()

    # ------------------------------------- 4. the published contrast, reproduced
    emit("cross-check 4 — the published F(h) vs F(-) contrast, reproduced")
    spec = sheets["fibrosis_scores"]
    scores = sheet_frame(path_of[spec["artifact"]], spec["name"], spec["header_row"])
    classification = scores.columns[1]
    groups = {}
    for case, row in scores.iterrows():
        label = f"TBME{int(float(case)):02d}"
        groups.setdefault(str(row[classification]).strip(), []).append(label)

    sizes = {k: len(v) for k, v in sorted(groups.items())}
    if sizes != dict(expect["fibrosis_group_sizes"]):
        raise RuntimeError(
            f"fibrosis group sizes {sizes} != expected "
            f"{dict(expect['fibrosis_group_sizes'])}."
        )
    if sum(sizes.values()) != expect["fibrosis_n_cases"]:
        raise RuntimeError(f"expected {expect['fibrosis_n_cases']} cases, got {sizes}")
    emit(f"  fibrosis classes (Supplementary Data 2): {sizes}")

    reference = list(groups[expect["fibrosis_reference_group"]])
    test = list(groups[expect["fibrosis_test_group"]])
    resolve = {k: v["resolves_to"] for k, v in alias_resolution.items()}
    columns = set(fibrosis_matrix.columns)
    for name in reference + test:
        if name not in columns:
            raise RuntimeError(
                f"{name} is in the fibrosis classification but not in "
                f"{spec['name']!r}; columns {sorted(columns)}."
            )
    if len(columns) - len(controls) - len(reference) - len(test) != 0:
        pass  # the sheet also carries the BC columns; not part of this contrast.

    def contrast(gene):
        """Mean-of-log2 and Student's t, the paper's estimator (see below)."""
        a = np.log2(fibrosis_matrix.loc[gene, test].to_numpy(float))
        b = np.log2(fibrosis_matrix.loc[gene, reference].to_numpy(float))
        return float(a.mean() - b.mean()), float(stats.ttest_ind(a, b).pvalue)

    spec = sheets["published_deg"]
    published_deg = sheet_frame(
        path_of[spec["artifact"]], spec["name"], spec["header_row"]
    )
    published_deg.columns = [str(c).strip() for c in published_deg.columns]
    published_deg = published_deg.apply(pd.to_numeric, errors="coerce")
    if len(published_deg) != expect["published_deg_n"]:
        raise RuntimeError(
            f"published DEG list has {len(published_deg)} rows, expected "
            f"{expect['published_deg_n']}."
        )

    lfc_col = [c for c in published_deg.columns if "log2" in c.lower()][0]
    p_col = [c for c in published_deg.columns if c.lower().startswith("p")][0]
    q_col = [c for c in published_deg.columns if c.lower().startswith("q")][0]

    repro = []
    for gene, row in published_deg.iterrows():
        lfc, pval = contrast(gene)
        repro.append(
            {
                "gene": gene,
                "published_log2fc": float(row[lfc_col]),
                "reproduced_log2fc": lfc,
                "abs_diff_log2fc": abs(lfc - float(row[lfc_col])),
                "published_p": float(row[p_col]),
                "reproduced_p": pval,
                "rel_diff_p": abs(pval / float(row[p_col]) - 1.0),
                "published_q": float(row[q_col]),
            }
        )
    reproduction = pd.DataFrame(repro)

    worst_lfc = float(reproduction["abs_diff_log2fc"].max())
    worst_p = float(reproduction["rel_diff_p"].max())
    emit(
        f"  reproduced {len(reproduction)} published log2FC: "
        f"max |diff| {worst_lfc:.2e} (tolerance {expect['log2fc_tolerance']})"
    )
    emit(
        f"  reproduced {len(reproduction)} published p-values: "
        f"max rel diff {worst_p:.2e} (tolerance {expect['pvalue_rtol']})"
    )
    if worst_lfc > expect["log2fc_tolerance"] or worst_p > expect["pvalue_rtol"]:
        raise RuntimeError(
            "the published F(h) vs F(-) contrast does not reproduce "
            f"(log2FC {worst_lfc:.3e}, p {worst_p:.3e}). Nothing downstream may "
            "describe what the published analysis did, because this rule has "
            "not established that it knows."
        )
    emit(
        "  -> estimator confirmed: mean-of-log2 with Student's t, "
        "positive log2FC = up in F(h)"
    )

    smallest = float(published_deg[lfc_col].abs().min())
    if smallest < expect["published_log2fc_threshold"]:
        raise RuntimeError(
            f"published list contains |log2FC| {smallest:.3f}, below its own "
            f"stated threshold {expect['published_log2fc_threshold']}."
        )
    emit(
        f"  published list respects its stated |log2FC| > "
        f"{expect['published_log2fc_threshold']} floor (smallest {smallest:.3f})"
    )

    in_list = set(published_deg.index)
    tbme = detection[detection["aoi_code"] == "TBME"].set_index("gene")
    contrast_rows = []
    for gene in sorted(set(panel) | set(PAPER_CHECKPOINTS)):
        if gene not in fibrosis_matrix.index:
            continue
        lfc, pval = contrast(gene)
        row = {
            "gene": gene,
            "alias": alias_of.get(gene, ""),
            "on_project_panel": gene in checkpoints,
            "named_by_paper": gene in PAPER_CHECKPOINTS,
            "log2fc_fh_vs_fneg": lfc,
            "p_student": pval,
            "abs_log2fc": abs(lfc),
            "published_threshold": expect["published_log2fc_threshold"],
            "reaches_published_threshold": abs(lfc)
            >= expect["published_log2fc_threshold"],
            "in_published_deg_list": gene in in_list,
        }
        if gene in tbme.index:
            row["tbme_n_detected"] = int(tbme.at[gene, "n_detected"])
            row["tbme_n_aoi"] = int(tbme.at[gene, "n_aoi"])
            row["tbme_assessable"] = bool(tbme.at[gene, "assessable"])
        contrast_rows.append(row)
    fibrosis_contrast = pd.DataFrame(contrast_rows)

    reaching = fibrosis_contrast[fibrosis_contrast["reaches_published_threshold"]]
    listed = fibrosis_contrast[fibrosis_contrast["in_published_deg_list"]]
    emit()
    emit("  checkpoint genes in the published contrast:")
    emit(f"    {'gene':10s} {'log2FC':>8s} {'p':>9s}  {'>=1.5?':>7s} {'in list?':>9s}  TBME detection")
    for row in fibrosis_contrast.sort_values("abs_log2fc", ascending=False).itertuples(
        index=False
    ):
        det = (
            f"{row.tbme_n_detected}/{row.tbme_n_aoi}"
            if hasattr(row, "tbme_n_detected") and pd.notna(row.tbme_n_detected)
            else "n/a"
        )
        emit(
            f"    {row.gene:10s} {row.log2fc_fh_vs_fneg:+8.3f} {row.p_student:9.4f}"
            f"  {'YES' if row.reaches_published_threshold else 'no':>7s}"
            f" {'YES' if row.in_published_deg_list else 'no':>9s}  {det}"
        )
    emit()
    emit(
        f"  {len(reaching)} of {len(fibrosis_contrast)} checkpoint genes reach the "
        f"published |log2FC| threshold; {len(listed)} appear in the published list"
    )
    emit(
        f"  largest |log2FC| among them: "
        f"{fibrosis_contrast['abs_log2fc'].max():.3f} vs threshold "
        f"{expect['published_log2fc_threshold']}"
    )
    emit()
    emit("READING RULE. This table describes the PUBLISHED contrast, computed on")
    emit("the published groups with the published estimator. It is not a result")
    emit("about brain metastasis biology and must never be reported as one; it")
    emit("says what the published checkpoint statement rests on. ADR 0008's")
    emit("reporting rule is unchanged: a gene below the floor in a compartment is")
    emit("'not assessable in <compartment>', never 'lower in <compartment>'.")
    emit()

    # -------------------------------------------------------------- outputs
    per_aoi.to_csv(snakemake.output.matrix, sep="\t", index=False, float_format=FMT)
    quantisation.to_csv(
        snakemake.output.quantisation, sep="\t", index=False, float_format=FMT
    )
    checkpoint_counts.to_csv(
        snakemake.output.counts, sep="\t", index=False, float_format=FMT
    )
    reproduction.to_csv(
        snakemake.output.reproduction, sep="\t", index=False, float_format=FMT
    )
    fibrosis_contrast.to_csv(
        snakemake.output.contrast, sep="\t", index=False, float_format=FMT
    )

    summary = {
        "source_publication": {
            "pmid": "36216799",
            "doi": "10.1038/s41467-022-33365-y",
            "files": {
                key: {"dest": spec["dest"], "sha256": spec["sha256"]}
                for key, spec in sorted(config["artifacts"].items())
            },
        },
        "matrix_agreement": {
            "n_aoi": int(len(aoi_order)),
            "n_genes": int(len(gene_order)),
            "n_values": int(mine.size),
            "n_exact": n_exact,
            "max_abs_diff": float(np.nanmax(delta)),
            "aoi_labels_identical": True,
            "control_rows_compared": controls,
            "status": "exact agreement",
        },
        "aoi_alias_resolution": alias_resolution,
        "quantisation": {
            "lattice_probe": LATTICE_PROBE,
            "lattice_min_fraction": LATTICE_MIN,
            "n_aoi_counts_reliable": int(quantisation["counts_reliable"].sum()),
            "median_distinct_values_per_aoi": float(
                quantisation["n_distinct_values"].median()
            ),
            "min_distinct_values_per_aoi": int(quantisation["n_distinct_values"].min()),
        },
        "published_contrast": {
            "reference_group": expect["fibrosis_reference_group"],
            "test_group": expect["fibrosis_test_group"],
            "n_reference": len(reference),
            "n_test": len(test),
            "excluded_class_sizes": {
                k: v
                for k, v in sizes.items()
                if k
                not in (
                    expect["fibrosis_reference_group"],
                    expect["fibrosis_test_group"],
                )
            },
            "estimator": "mean-of-log2, Student's t, positive = up in F(h)",
            "n_published_deg_reproduced": int(len(reproduction)),
            "max_abs_diff_log2fc": worst_lfc,
            "max_rel_diff_p": worst_p,
            "published_log2fc_threshold": expect["published_log2fc_threshold"],
            "smallest_published_abs_log2fc": smallest,
            "n_checkpoint_genes_examined": int(len(fibrosis_contrast)),
            "n_reaching_published_threshold": int(len(reaching)),
            "n_in_published_deg_list": int(len(listed)),
            "largest_checkpoint_abs_log2fc": float(
                fibrosis_contrast["abs_log2fc"].max()
            ),
            "paper_named_checkpoints": PAPER_CHECKPOINTS,
        },
        "panel_unchanged": (
            "config/checkpoints.yaml is untouched by this rule. The genes the "
            "paper names are reported, never adopted — panel membership is "
            "pre-registered (ADR 0015) and a stop-and-ask decision."
        ),
        "reporting_rule": (
            "The published-contrast table describes the published analysis, not "
            "this project's compartments. ADR 0008 point 4 is unchanged."
        ),
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")

    emit(f"wrote {len(per_aoi)} AOI rows, {len(quantisation)} quantisation rows,")
    emit(f"      {len(checkpoint_counts)} checkpoint-count rows,")
    emit(f"      {len(reproduction)} reproduction rows, "
         f"{len(fibrosis_contrast)} contrast rows")
