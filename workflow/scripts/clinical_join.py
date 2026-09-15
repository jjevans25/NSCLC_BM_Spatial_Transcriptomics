"""P5-T1 — parse Supplementary Data 1 and join it to the GEO design, or stop the phase.

Owner task: P5-T1. Driven by rule p5t1_clinical_join.

PROJECT_PLAN §6 makes this a HARD GATE: "if it doesn't join cleanly on patient
ID, stop the phase. Do not spend three days on fuzzy ID matching for a stretch
goal." Its acceptance criterion is "a joined table **or** a documented
abandonment", and both are legitimate outcomes. So this script has **no fuzzy
matching, no fallback path and no second threshold** — the same shape
`lr_nominations.py` has, and for the same reason: an honest failure must be
REACHABLE rather than escapable, or the gate is decorative.

Everything it applies was pre-registered in ADR 0024 and committed before any
Phase 5 number existed. The parameters live in `config/clinical.yaml`
(what the FILE contains) and `config/config.yaml → survival` (what the ANALYSIS
uses), which are deliberately redundant: common.smk asserts they agree at parse
time, because an endpoint assigned to the wrong arm inverts the phase's result
and nothing downstream would notice.

WHAT THE DEPOSITED FILE DOES THAT A NAIVE READ GETS WRONG
---------------------------------------------------------
Four things, all measured from the file rather than assumed, all recorded in
ADR 0024 §3-§4 and `docs/data-provenance.md` §Q5:

  1. **Header cells carry embedded NEWLINES.** `Age at \nNSCLC diagnosis`,
     `Primary lung cancer \ndiagnosis to death (Months)`; and
     ` Location of BrM ` has a leading AND trailing space.
     `docs/data-provenance.md` §Q5 transcribed the column list
     already-normalised, so it must not be matched verbatim. Normalisation —
     collapse every whitespace run to one space, then strip — is applied to
     headers and string values alike.

  2. **The sheet ends in a legend, not in data.** 53 rows: title, blank,
     header, 44 patients, three blanks, then three treatment-abbreviation lines
     (`SR: Surgical resection`, …). A read to end-of-sheet ingests those as
     three patients with a null ID. The manifest bounds the data rows.

  3. **The two endpoint columns spell censoring DIFFERENTLY.** The lung column
     carries the literal `Alive`; the brain column carries an EMPTY CELL for the
     same patients. The lung `Alive` set and the brain blank set must be
     IDENTICAL, and that identity is what licenses reading a blank as CENSORED
     rather than MISSING. It is asserted here, not assumed — reading those
     blanks as missing would silently drop every censored observation in the
     study, and there are only three.

  4. **Missingness has five spellings** (`N/A`, `NA`, `n/a`, `Unspecified`, the
     empty cell) and `Gender` carries a lowercase `m`. A token outside the
     declared set is a HARD FAILURE, because a sixth spelling that reads as data
     is the failure the declaration exists to prevent.

STDLIB XLSX PARSING, LIFTED NOT IMPORTED
----------------------------------------
The reader below is lifted from `workflow/scripts/external_validation.py`, the
way `checkpoint_detection.py` lifts its detection rule from
`score_signatures.py`. This project has no cross-script imports and Snakemake
script-mode sibling imports are not reliable across versions; editing
`external_validation.py` to extract a shared module would also fire its code
rerun-trigger and re-parse a 45 MB workbook.

`openpyxl` is deliberately absent from `workflow/envs/py-analysis.yaml` and must
stay absent: that env is `p0t2_fetch_geo`'s, whose outputs are `protected()`, so
a dependency change fires the software-env rerun trigger and the DAG dies on
ProtectedOutputException (ADR 0005, ADR 0013, ADR 0023). An .xlsx is a zip of
XML; `zipfile` and `ElementTree` are enough.

WHAT IT DELIBERATELY DOES NOT DO
--------------------------------
  * **Fit anything.** This task joins and asserts. P5-T2 aggregates, P5-T3
    models. A survival estimate produced here would be a number computed before
    the table it rests on had been checked.
  * **Drop a patient.** Flag-don't-drop (CLAUDE.md): the nine supplementary rows
    with no expression data are RECORDED as `in_geo = False`, not removed, and
    the 35/44 shortfall is expected (Q5 constraint 3) rather than a failure.
  * **Discover the cohort.** Every count in `clinical.expect` is asserted
    against the file. If one moves, the design changed and the phase stops
    rather than adapting to it.
"""

import json
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

FMT = "%.17g"

# --------------------------------------------------------------------------
# Minimal .xlsx reader (stdlib only) — lifted from external_validation.py.
# See the module docstring for why it is lifted rather than imported.
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
        "Sheet names are transcribed exactly in config/clinical.yaml, trailing "
        "whitespace included — do not tidy them."
    )


def read_sheet(path, name):
    """One worksheet as a list of rows, each a list of str | float | None.

    Rows are placed at the index their `r` attribute declares, and gaps are
    padded. Excel OMITS entirely blank rows from the XML, so appending
    sequentially silently shifts everything below one — which is exactly what
    happens here: this sheet has a blank row between its title and its header.
    The `header_row` value in config/clinical.yaml is a true sheet position, the
    one a person reading the file in Excel would count, and this is what makes
    that true.

    float() rather than a pandas parser: ADR 0013's postscript records that
    pandas' default C parser is not correctly rounded, and float() is.
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
                        elif kind == "str":
                            value = node.text
                        else:
                            try:
                                value = float(node.text)
                            except ValueError:
                                value = node.text
                    cells[index] = value
                row_index = int(elem.get("r")) - 1
                width = (max(cells) + 1) if cells else 0
                placed[row_index] = [cells.get(i) for i in range(width)]
                elem.clear()
    if not placed:
        return []
    return [placed.get(i, []) for i in range(max(placed) + 1)]


log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg=""):
        print(msg, file=log, flush=True)

    clinical = snakemake.params.clinical
    survival = snakemake.params.survival
    sheet_spec = clinical["sheets"]["cohort"]
    norm_spec = clinical["normalisation"]
    columns = clinical["columns"]
    expect = clinical["expect"]

    emit("P5-T1 — clinical join (ADR 0024, pre-registered)")
    emit("=" * 70)
    emit()
    emit("  HARD GATE. PROJECT_PLAN §6: if it does not join cleanly on patient")
    emit("  ID, STOP THE PHASE. No fuzzy matching, no fallback path. P5-T1's")
    emit("  acceptance criterion is a joined table OR a documented abandonment,")
    emit("  and both are legitimate outcomes.")
    emit()
    emit("  Phase 5 is a STRETCH goal and Gate 5 is 'timebox respected'. No")
    emit("  Phase 5 result may be a headline claim.")
    emit()

    # --- the normalisation rule, applied everywhere -------------------------
    def normalise(value):
        """The one rule, per config/clinical.yaml → normalisation.

        Applied to headers and to string values alike. Non-strings pass
        through untouched so a float endpoint stays a float.
        """
        if not isinstance(value, str):
            return value
        out = value
        if norm_spec["collapse_whitespace"]:
            out = re.sub(r"\s+", " ", out)
        if norm_spec["strip"]:
            out = out.strip()
        return out

    missing_tokens = set(norm_spec["missing_tokens"])

    def is_missing(value):
        if value is None:
            return True
        if isinstance(value, str):
            return normalise(value) in missing_tokens
        return False

    # --- read the sheet -----------------------------------------------------
    # input.clinical is a LIST (the manifest may hold more than one artifact),
    # and this sheet lives in exactly one of them. Resolved by the `artifact`
    # key the sheet spec names rather than by position, so adding a second
    # clinical artifact later cannot silently repoint this read.
    wanted = clinical["artifacts"][sheet_spec["artifact"]]["dest"]
    candidates = [
        path for path in snakemake.input.clinical if Path(path).name == wanted
    ]
    if len(candidates) != 1:
        raise RuntimeError(
            f"sheet 'cohort' names artifact "
            f"'{sheet_spec['artifact']}' (dest {wanted!r}), which matched "
            f"{len(candidates)} of the rule's inputs {list(snakemake.input.clinical)}."
        )
    workbook = candidates[0]

    rows = read_sheet(workbook, sheet_spec["name"])
    emit(f"workbook : {workbook}")
    emit(f"sheet    : {sheet_spec['name']!r}  ({len(rows)} rows)")

    header_raw = rows[sheet_spec["header_row"]]
    header = [normalise(h) for h in header_raw]
    emit(f"header   : row {sheet_spec['header_row']} (0-based), "
         f"{len(header)} columns")

    # The embedded newlines are the reason docs/data-provenance.md §Q5's column
    # list is normalised. Logged verbatim so the quirk stays visible rather than
    # being silently absorbed by the rule that handles it.
    quirky = [
        (j, raw) for j, raw in enumerate(header_raw)
        if isinstance(raw, str) and normalise(raw) != raw
    ]
    emit(f"           {len(quirky)} header cells needed normalising:")
    for j, raw in quirky:
        emit(f"             [{j:2d}] {raw!r} -> {normalise(raw)!r}")
    emit()

    if len(header) != expect["n_columns"]:
        raise RuntimeError(
            f"Supplementary Data 1 has {len(header)} columns, expected "
            f"{expect['n_columns']} (clinical.expect.n_columns). The deposited "
            "file changed shape; stop rather than adapt to it."
        )

    index_of = {}
    for name in [columns["patient_id"], columns["age_band"],
                 columns["endpoints"]["lung"], columns["endpoints"]["brain"]]:
        if name not in header:
            raise RuntimeError(
                f"required column {name!r} is not in the normalised header.\n"
                f"  header: {header}\n"
                "Names in config/clinical.yaml → columns are NORMALISED "
                "(whitespace collapsed, stripped). The deposited header cells "
                "carry embedded newlines, so a verbatim name will never match."
            )
        index_of[name] = header.index(name)

    # --- bounded data rows --------------------------------------------------
    # NOT read to end-of-sheet: rows 47-49 are blank and 50-52 are the
    # treatment-abbreviation legend, which an unbounded read ingests as three
    # patients with a null ID.
    first, last = sheet_spec["first_data_row"], sheet_spec["last_data_row"]
    data_rows = rows[first:last + 1]
    emit(f"data rows: {first}..{last} inclusive ({len(data_rows)} rows)")

    trailing = rows[last + 1:]
    non_empty_trailing = [
        [v for v in r if v is not None] for r in trailing
        if any(v is not None for v in r)
    ]
    emit(f"           {len(trailing)} rows below the bound, "
         f"{len(non_empty_trailing)} of them non-empty and EXCLUDED:")
    for r in non_empty_trailing:
        emit(f"             {r}")
    emit()

    # --- build the frame ----------------------------------------------------
    records = []
    for offset, raw_row in enumerate(data_rows):
        padded = list(raw_row) + [None] * (len(header) - len(raw_row))
        record = {}
        for j, name in enumerate(header):
            value = padded[j]
            record[name] = None if is_missing(value) else normalise(value)
        record["_sheet_row"] = first + offset
        records.append(record)

    frame = pd.DataFrame.from_records(records)

    # Patient ID arrives as a float (Excel stores it numerically). Coerce to a
    # nullable integer so the join key is an integer on both sides — a float
    # key would never match samples.tsv's parsed "P12" -> 12.
    pid_col = columns["patient_id"]
    if frame[pid_col].isna().any():
        bad = frame.loc[frame[pid_col].isna(), "_sheet_row"].tolist()
        raise RuntimeError(
            f"sheet rows {bad} have no Patient ID. The data-row bounds in "
            "config/clinical.yaml are wrong, or the sheet changed — a null ID "
            "here is normally the treatment-abbreviation legend being read as "
            "data."
        )
    frame["patient_number"] = frame[pid_col].astype(float).astype(int)
    frame = frame.drop(columns=[pid_col])

    # --- assertions on the cohort ------------------------------------------
    n_checks = 0
    emit("COHORT ASSERTIONS (clinical.expect, pre-registered)")
    emit("-" * 70)

    if len(frame) != expect["n_cohort"]:
        raise RuntimeError(
            f"read {len(frame)} patients, expected {expect['n_cohort']} "
            "(clinical.expect.n_cohort)."
        )
    emit(f"  [ok] cohort size                 {len(frame)}")
    n_checks += 1

    numbers = sorted(frame["patient_number"])
    span = list(range(expect["patient_id_min"], expect["patient_id_max"] + 1))
    if numbers != span:
        raise RuntimeError(
            f"patient IDs are not contiguous {expect['patient_id_min']}.."
            f"{expect['patient_id_max']}. Missing: "
            f"{sorted(set(span) - set(numbers))}; unexpected: "
            f"{sorted(set(numbers) - set(span))}."
        )
    emit(f"  [ok] IDs contiguous              "
         f"{expect['patient_id_min']}..{expect['patient_id_max']}")
    n_checks += 1

    age_col = columns["age_band"]
    observed_bands = sorted({v for v in frame[age_col] if v is not None})
    if observed_bands != sorted(expect["age_bands"]):
        raise RuntimeError(
            f"age bands {observed_bands} do not match the declared "
            f"{sorted(expect['age_bands'])}. Declared rather than inferred so a "
            "band that appears or vanishes is caught (ADR 0024 §4)."
        )
    # Ordered factor, NEVER continuous (Q5 constraint 1). Recorded as an ordered
    # categorical here so the type carries the rule downstream rather than a
    # comment asking P5-T3 to remember it.
    frame[age_col] = pd.Categorical(
        frame[age_col], categories=expect["age_bands"], ordered=True
    )
    emit(f"  [ok] age bands (ordered factor)  {expect['age_bands']}")
    n_checks += 1

    # --- censoring: the two columns spell it differently --------------------
    emit()
    emit("CENSORING (ADR 0024 §3 — measured, not assumed)")
    emit("-" * 70)
    lung_col = columns["endpoints"]["lung"]
    brain_col = columns["endpoints"]["brain"]
    token = expect["censoring_token_lung"]

    lung_alive = sorted(
        frame.loc[frame[lung_col] == token, "patient_number"]
    )
    brain_blank = sorted(
        frame.loc[frame[brain_col].isna(), "patient_number"]
    )
    emit(f"  lung column  {lung_col!r}")
    emit(f"    == {token!r} for patients {lung_alive}")
    emit(f"  brain column {brain_col!r}")
    emit(f"    EMPTY        for patients {brain_blank}")

    if lung_alive != brain_blank:
        raise RuntimeError(
            f"the lung '{token}' set {lung_alive} and the brain blank set "
            f"{brain_blank} are NOT identical.\n"
            "That identity is the whole licence for reading a blank brain "
            "value as CENSORED rather than MISSING (ADR 0024 §3). Without it, "
            "one of the two columns means something this project has not "
            "established, and the phase stops rather than guessing which."
        )
    if lung_alive != sorted(expect["censored_patients"]):
        raise RuntimeError(
            f"censored patients {lung_alive} do not match the pre-registered "
            f"{sorted(expect['censored_patients'])} (clinical.expect)."
        )
    emit(f"  [ok] the two sets are IDENTICAL  {lung_alive}")
    emit("       -> a blank brain value is CENSORED, not missing")
    n_checks += 1

    if not expect["censoring_is_blank_brain"]:
        raise RuntimeError(
            "clinical.expect.censoring_is_blank_brain is false, but this "
            "script's censoring logic depends on it being true."
        )

    # Build the (time, event) pair per arm. event = 1 death, 0 censored.
    #
    # `pd.isna`, NOT `value is None`. The brain column is numeric, so
    # DataFrame.from_records turns its empty cells into NaN — and `NaN is None`
    # is False, which silently counts the censored patient as an EVENT with a
    # missing time. The lung column escapes that only because its censoring is
    # a string. The two columns spelling censoring differently is exactly what
    # makes this asymmetric, and the per-arm assertion below is what catches it.
    censored_set = set(lung_alive)
    for arm, col in (("lung", lung_col), ("brain", brain_col)):
        times, events = [], []
        for value in frame[col]:
            censored = pd.isna(value) or (value == token)
            if censored:
                times.append(pd.NA)
                events.append(0)
            else:
                times.append(float(value))
                events.append(1)
        frame[f"{arm}_months"] = times
        frame[f"{arm}_event"] = events

        # The event flag must agree with the censoring set derived above, per
        # patient, in BOTH arms. Without this the NaN-vs-None asymmetry above
        # produces a plausible table with one extra death in it.
        flagged = set(frame.loc[frame[f"{arm}_event"] == 0, "patient_number"])
        if flagged != censored_set:
            raise RuntimeError(
                f"{arm} arm: patients flagged censored {sorted(flagged)} do not "
                f"match the censoring set {sorted(censored_set)}. The two "
                "endpoint columns spell censoring differently (ADR 0024 §3), so "
                "a censoring test that works on one column can silently fail on "
                "the other."
            )
        n_checks += 1
    emit()
    emit(f"  [ok] event flags agree with the censoring set in BOTH arms")

    # A censored patient has no time in EITHER column — the lung cell says
    # `Alive` and the brain cell is blank — so neither arm can supply a
    # follow-up duration. Recorded rather than imputed: an invented censoring
    # time is a fabricated observation, and P5-T3 must decide what to do with
    # three patients it cannot place on a time axis.
    emit()
    emit("  NOTE: a censored patient has NO follow-up time in either column —")
    emit("        `Alive` carries no duration and the brain cell is blank. The")
    emit("        time is left missing rather than imputed; P5-T3 decides how")
    emit("        to handle three patients with an event flag and no time.")

    # --- the join to the GEO design ----------------------------------------
    emit()
    emit("JOIN ASSERTIONS (docs/data-provenance.md §Q5 — P5-T1's HARD GATE)")
    emit("-" * 70)

    samples = pd.read_csv(snakemake.input.samples, sep="\t", dtype=str)
    geo = samples.loc[samples["is_control"] == "False"].copy()
    emit(f"  samples.tsv: {len(samples)} AOIs, {len(geo)} non-control")

    pid_pattern = re.compile(r"P(\d+)$")
    unparseable = sorted(
        {p for p in geo["patient_id"] if not pid_pattern.fullmatch(str(p))}
    )
    if unparseable:
        raise RuntimeError(
            f"patient_id values {unparseable} do not match 'P<number>'. "
            "P0-T3 writes this column; the join key cannot be derived without "
            "it."
        )
    geo["patient_number"] = [
        int(pid_pattern.fullmatch(str(p)).group(1)) for p in geo["patient_id"]
    ]

    geo_numbers = sorted(set(geo["patient_number"]))
    emit(f"  distinct GEO patients: {len(geo_numbers)}")

    # ASSERTION 1 — every GEO patient appears in Supplementary Data 1.
    if expect["assert_every_geo_patient_present"]:
        absent = sorted(set(geo_numbers) - set(frame["patient_number"]))
        if absent:
            raise RuntimeError(
                f"GEO patients {absent} are ABSENT from Supplementary Data 1.\n"
                "This is P5-T1's hard gate and its failure action is to STOP "
                "THE PHASE (PROJECT_PLAN §6, ADR 0024 §5). Do NOT add fuzzy "
                "matching: P5-T1's acceptance criterion is a joined table OR a "
                "documented abandonment, and this is the abandonment."
            )
        emit(f"  [ok] ASSERTION 1  all {len(geo_numbers)} GEO patients present "
             "in Supplementary Data 1")
        n_checks += 1

    # ASSERTION 2 — the AOI-code numeric suffix equals the patient number.
    if expect["assert_aoi_suffix_equals_patient"]:
        aoi_pattern = re.compile(r"(?P<code>[A-Za-z-]+)(?P<num>\d+)(?P<rep>[a-z]?)$")
        mismatches = []
        for _, row in geo.iterrows():
            match = aoi_pattern.fullmatch(str(row["aoi_label"]))
            if match is None:
                mismatches.append((row["aoi_label"], row["patient_id"], "unparseable"))
            elif int(match.group("num")) != row["patient_number"]:
                mismatches.append((row["aoi_label"], row["patient_id"], "mismatch"))
        if mismatches:
            raise RuntimeError(
                f"{len(mismatches)} AOI labels disagree with their patient_id: "
                f"{mismatches[:10]}.\n"
                "This is P5-T1's hard gate (docs/data-provenance.md §Q5, "
                "ADR 0024 §5). STOP THE PHASE rather than matching loosely."
            )
        emit(f"  [ok] ASSERTION 2  all {len(geo)} AOI-code suffixes equal their "
             "patient number")
        n_checks += 1

    # Coverage is 35 of 44 and that is EXPECTED, not a failure (Q5 constraint 3).
    frame["in_geo"] = frame["patient_number"].isin(geo_numbers)
    n_in_geo = int(frame["in_geo"].sum())
    if n_in_geo != expect["n_geo_patients"]:
        raise RuntimeError(
            f"{n_in_geo} supplementary rows have expression data, expected "
            f"{expect['n_geo_patients']} (clinical.expect.n_geo_patients)."
        )
    emit(f"  [ok] coverage                    {n_in_geo} of {len(frame)} "
         "supplementary rows have expression data")
    emit("       the other "
         f"{len(frame) - n_in_geo} are RECORDED as in_geo = False, not dropped "
         "(flag-don't-drop)")
    n_checks += 1

    # --- the TIME cohort Phase 5 actually models ---------------------------
    emit()
    emit("THE PHASE 5 COHORT (ADR 0024 §2)")
    emit("-" * 70)
    aoi_codes = survival["aoi_codes"]
    time_aois = geo.loc[geo["aoi_code"].isin(aoi_codes)]
    by_code = {
        code: sorted(set(time_aois.loc[time_aois["aoi_code"] == code,
                                       "patient_number"]))
        for code in aoi_codes
    }
    for code in aoi_codes:
        emit(f"  {code:8s} {len(by_code[code]):2d} patients  {by_code[code]}")
    union = sorted(set().union(*by_code.values()))
    shared = sorted(set.intersection(*(set(v) for v in by_code.values())))
    emit(f"  shared   {len(shared):2d} patients  {shared}")
    emit(f"  DISTINCT {len(union):2d} patients — this is the n, not 35")
    emit()
    emit("  PROJECT_PLAN §6 P5-T3's '~35 patients' counts patients with a")
    emit("  TUMOUR AOI, none of which carries a signature score. Phase 5 does")
    emit("  not chase that number: reaching it means scoring L/LB, which is a")
    emit("  stop-and-ask scope change and re-opens the detection floor where")
    emit("  the panel clears far less (ADR 0024 §2).")

    frame["in_phase5"] = frame["patient_number"].isin(union)

    # Event counts per arm, over the patients each arm actually models. These
    # are what make a Phase 5 null uninformative rather than negative, so they
    # are computed rather than quoted.
    emit()
    emit("  events per arm (over the patients that arm models):")
    arm_stats = {}
    code_to_arm = {"TIME-L": "lung", "TIME-B": "brain"}
    for code in aoi_codes:
        arm = code_to_arm[code]
        subset = frame.loc[frame["patient_number"].isin(by_code[code])]
        n_events = int(subset[f"{arm}_event"].sum())
        n_cens = int(len(subset) - n_events)
        arm_stats[arm] = {
            "aoi_code": code,
            "n_patients": int(len(subset)),
            "n_events": n_events,
            "n_censored": n_cens,
            "patients": by_code[code],
        }
        emit(f"    {arm:5s} ({code})  n = {len(subset):2d}  "
             f"events = {n_events:2d}  censored = {n_cens}")

    emit()
    emit("  With these event counts the analysis detects only very large hazard")
    emit("  ratios. A Phase 5 null is UNINFORMATIVE, NOT NEGATIVE (ADR 0024 §7),")
    emit("  and a separation at this n describes this cohort rather than")
    emit("  supporting an inferential claim. Brain is TIME-B n = 8 (hard")
    emit("  constraint 8).")

    # --- write ---------------------------------------------------------------
    frame = frame.sort_values("patient_number").reset_index(drop=True)
    ordered = (
        ["patient_number", "in_geo", "in_phase5"]
        + [c for c in frame.columns
           if c not in {"patient_number", "in_geo", "in_phase5", "_sheet_row"}]
        + ["_sheet_row"]
    )
    frame = frame[ordered]
    frame.to_csv(snakemake.output.clinical, sep="\t", index=False, float_format=FMT)

    emit()
    emit(f"assertions: {n_checks} of {n_checks} passed")

    summary = {
        "task": "P5-T1",
        "phase": 5,
        "stretch": True,
        "gate": "Gate 5 is 'timebox respected', not 'a result was found'",
        "pre_registered": "ADR 0024",
        "source": {
            "artifact": str(workbook),
            "sheet": sheet_spec["name"],
            "header_row": sheet_spec["header_row"],
            "data_rows": [first, last],
            "n_header_cells_normalised": len(quirky),
            "n_trailing_rows_excluded": len(non_empty_trailing),
        },
        "n_cohort": int(len(frame)),
        "n_in_geo": n_in_geo,
        "n_phase5_patients": len(union),
        "phase5_patients": union,
        "phase5_shared_patients": shared,
        "patients_by_aoi_code": by_code,
        "censoring": {
            "lung_token": token,
            "brain_is_blank": True,
            "censored_patients": lung_alive,
            "sets_identical": True,
            "note": (
                "The two endpoint columns spell censoring differently and the "
                "sets are asserted identical. A censored patient carries no "
                "follow-up time in either column; the time is left missing "
                "rather than imputed."
            ),
        },
        "arms": arm_stats,
        "join": {
            "every_geo_patient_present": True,
            "aoi_suffix_equals_patient": True,
            "n_geo_aoi_checked": int(len(geo)),
            "fuzzy_matching": "none — forbidden by PROJECT_PLAN §6 and ADR 0024 §5",
        },
        "n_assertions": n_checks,
        "reporting_rule": (
            "No Phase 5 result may be a headline claim. Brain is TIME-B n = 8 "
            "(hard constraint 8). With 12 events in lung and 7 in brain a null "
            "is an assay- and cohort-size limit — UNINFORMATIVE, NOT NEGATIVE "
            "(ADR 0024 §7) — and a separation at 6 vs 7 patients describes this "
            "cohort rather than supporting an inference. Every p-value carries "
            "n, effect size and a confidence interval (hard constraint 7). "
            "Three of the six signatures failed their Phase 2 detection floor "
            "and are 'not assessable' wherever they did, never a prognostic "
            "null (ADR 0008)."
        ),
    }
    with open(snakemake.output.summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")

    emit(f"wrote {snakemake.output.clinical}, {snakemake.output.summary}")
