# P4-T2 — export the pinned CellChatDB .rda across the R/Python boundary.
#
# Owner task: P4-T2. Driven by rule p4t2c_export_lr_database.
#
# The narrowest possible R rule, and it exists for one reason: CellChatDB is
# distributed as a base-R .rda and nothing in the Python env can read one.
# `load()` is a base function, so this needs no CellChat installation and no
# package beyond what r-geomx already has — the same env that reads .RData
# deconvolution references in deconvolve.R. No conda env YAML is edited.
#
# It does NO science. It resolves no symbol, applies no detection floor, drops
# no interaction except the pre-registered annotation classes, and ranks
# nothing. Everything downstream of the parse is Python's (p4t2d), because the
# alternative is a second place for the detection rule to live.
#
# The handoff is plain TSV, never zellkonverter (ADR 0001), and — as everywhere
# else on this boundary — the reader asserts STRUCTURE, not byte-exact floats
# (ADR 0013 §3). Nothing here is a float: the database is symbols and labels,
# which is the one case where a TSV round-trip is exact by construction.

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

emit <- function(...) cat(..., "\n", sep = "")

emit("P4-T2 — export the pinned ligand-receptor database")
emit(strrep("=", 70))
emit("")

expect <- snakemake@params[["expect"]]
exclude <- unlist(snakemake@params[["exclude_annotations"]])
database <- snakemake@params[["database"]]
version <- snakemake@params[["version"]]

emit("database   : ", database, " ", version)
emit("pinned at  : ", snakemake@input[["rda"]])
emit("excluding  : ", paste(exclude, collapse = ", "),
     "   (pre-registered, ADR 0021 section 2)")
emit("")

# --- load ------------------------------------------------------------------
# Into a private environment rather than the global one: an .rda names its own
# objects, and load() into globalenv() would let a renamed object silently
# shadow something. Asserting the name is how the pin is checked to be the file
# it claims to be.
env <- new.env(parent = emptyenv())
loaded <- load(snakemake@input[["rda"]], envir = env)
emit("objects in the .rda: ", paste(loaded, collapse = ", "))

if (length(loaded) != 1L) {
  stop("expected exactly one object in the .rda, found ", length(loaded),
       " (", paste(loaded, collapse = ", "), "). The pin names a single ",
       "CellChatDB object; anything else means the pinned commit moved.")
}
db <- get(loaded[1], envir = env)

required <- c("interaction", "complex")
missing <- setdiff(required, names(db))
if (length(missing) > 0L) {
  stop("CellChatDB object is missing component(s): ",
       paste(missing, collapse = ", "), ". Present: ",
       paste(names(db), collapse = ", "), ".")
}
emit("components: ", paste(names(db), collapse = ", "))
emit("")

interaction <- db$interaction
complex <- db$complex

# --- assert the manifest's expectations ------------------------------------
# config/ligand_receptor.yaml declares these BEFORE the parse runs. A database
# that parses to different numbers means the pin moved or the parse is wrong,
# and either stops the phase rather than proceeding against a database nobody
# checked.
check <- function(label, got, want) {
  ok <- identical(as.integer(got), as.integer(want))
  emit(sprintf("  %-28s %6d   expected %6d   %s",
               label, as.integer(got), as.integer(want),
               if (ok) "ok" else "MISMATCH"))
  ok
}

emit("assertions against config/ligand_receptor.yaml expect:")
results <- c(
  check("n_interactions", nrow(interaction), expect[["n_interactions"]]),
  check("n_complex", nrow(complex), expect[["n_complex"]])
)

ann_counts <- table(interaction$annotation)
for (cls in sort(names(expect[["annotations"]]))) {
  got <- if (cls %in% names(ann_counts)) ann_counts[[cls]] else 0L
  results <- c(results, check(paste0("annotation: ", cls), got,
                              expect[["annotations"]][[cls]]))
}

excluded_mask <- interaction$annotation %in% exclude
usable <- interaction[!excluded_mask, , drop = FALSE]
results <- c(
  results,
  check("n_excluded_by_annotation", sum(excluded_mask),
        expect[["n_excluded_by_annotation"]]),
  check("n_usable", nrow(usable), expect[["n_usable"]])
)

if (!all(results)) {
  stop("The pinned database does not match config/ligand_receptor.yaml's ",
       "`expect` block. Either the pin moved (base_url must name a COMMIT, ",
       "never a branch) or this parse is wrong. Both are stop-and-ask: the ",
       "interaction panel is pre-registered (ADR 0021 section 2), and editing ",
       "`expect` to match what was read would be recording an answer rather ",
       "than checking one.")
}
emit("")
emit("all ", length(results), " manifest assertions passed")
emit("")

# --- classify 1:1 vs complex, and expand subunits --------------------------
# ADR 0021 section 4: the PRIMARY admits only interactions whose ligand and
# receptor each resolve to exactly one gene symbol. Complexes go to a declared
# EXPLORATORY table, admitted only when EVERY subunit clears the floor, and
# carrying no FDR of any kind. The split is made here, in one place, so no
# downstream table can disagree about which row is which.
#
# A complex is named by a rowname of db$complex; anything else is a bare symbol.
complex_names <- rownames(complex)
is_complex <- function(x) x %in% complex_names

subunits_of <- function(x) {
  if (!is_complex(x)) return(x)
  row <- complex[x, , drop = FALSE]
  vals <- unlist(row, use.names = FALSE)
  vals <- vals[!is.na(vals)]
  vals <- trimws(vals)
  vals[nzchar(vals)]
}

one_to_one <- !is_complex(usable$ligand) & !is_complex(usable$receptor)
emit("interaction split (ADR 0021 section 4):")
emit(sprintf("  one-to-one (primary candidates)   %5d", sum(one_to_one)))
emit(sprintf("  complex    (exploratory only)     %5d", sum(!one_to_one)))
results <- check("n_one_to_one", sum(one_to_one), expect[["n_one_to_one"]])
if (!results) {
  stop("one-to-one count disagrees with the manifest; see above.")
}
emit("")

out <- data.frame(
  interaction_name = as.character(usable$interaction_name),
  pathway_name = as.character(usable$pathway_name),
  ligand = as.character(usable$ligand),
  receptor = as.character(usable$receptor),
  annotation = as.character(usable$annotation),
  interaction_name_2 = as.character(usable$interaction_name_2),
  evidence = as.character(usable$evidence),
  is_one_to_one = one_to_one,
  ligand_is_complex = is_complex(usable$ligand),
  receptor_is_complex = is_complex(usable$receptor),
  ligand_subunits = vapply(usable$ligand,
                           function(x) paste(subunits_of(x), collapse = ";"),
                           character(1), USE.NAMES = FALSE),
  receptor_subunits = vapply(usable$receptor,
                             function(x) paste(subunits_of(x), collapse = ";"),
                             character(1), USE.NAMES = FALSE),
  stringsAsFactors = FALSE
)

# v2 carries per-partner molecular annotation v1 lacks. Exported because it
# lets Phase 4's likely central finding — that membrane and matrix ligands
# survive detection while soluble ones do not — be stated in the database's own
# terms rather than inferred from the three-way class column alone.
#
# DESCRIPTIVE ONLY, NEVER A FILTER. Nothing enters or leaves an analysis table
# on the strength of these columns (ADR 0021 section 2).
optional <- c("ligand.secreted_type", "ligand.transmembrane",
              "receptor.transmembrane", "receptor.surfaceome_main",
              "receptor.adhesome", "is_neurotransmitter")
for (col in optional) {
  target <- gsub(".", "_", col, fixed = TRUE)
  out[[target]] <- if (col %in% colnames(usable)) {
    as.character(usable[[col]])
  } else {
    NA_character_
  }
}
emit("descriptive columns carried (never used as a filter): ",
     paste(gsub(".", "_", optional, fixed = TRUE), collapse = ", "))

# A subunit string must never be empty — that would silently admit a complex
# with no genes to the exploratory table.
bad <- which(!nzchar(out$ligand_subunits) | !nzchar(out$receptor_subunits))
if (length(bad) > 0L) {
  stop("interaction(s) resolved to an empty subunit list: ",
       paste(out$interaction_name[bad], collapse = ", "),
       ". A complex with no subunits cannot be detection-filtered.")
}

cx_out <- data.frame(
  complex_name = complex_names,
  subunits = vapply(complex_names,
                    function(x) paste(subunits_of(x), collapse = ";"),
                    character(1), USE.NAMES = FALSE),
  n_subunits = vapply(complex_names,
                      function(x) length(subunits_of(x)),
                      integer(1), USE.NAMES = FALSE),
  stringsAsFactors = FALSE
)

write.table(out, snakemake@output[["interactions"]], sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")
write.table(cx_out, snakemake@output[["complexes"]], sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")

emit("")
emit("wrote ", nrow(out), " usable interactions -> ",
     snakemake@output[["interactions"]])
emit("wrote ", nrow(cx_out), " complex definitions -> ",
     snakemake@output[["complexes"]])
emit("")
emit("This rule resolves NO symbol and applies NO detection floor. Which")
emit("interactions are measurable is the Phase 4 RESULT (p4t2d), not an input.")

sink(type = "message")
sink(type = "output")
