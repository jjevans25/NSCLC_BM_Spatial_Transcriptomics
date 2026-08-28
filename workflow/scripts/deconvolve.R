# P2-T5 -- cell-composition deconvolution of the TIME AOIs, two arms.
#
# Owner task: P2-T5. Driven by rule p2t5_deconvolution.
#
# Choices worth stating, because each could reasonably go the other way:
#
#   * TWO ARMS, and only one of them supports a cross-site comparison.
#
#       safetme  safeTME (906 x 18, collapsed to 14 via safeTME.matches),
#                applied to BOTH sites. One reference for both, so its
#                composition IS comparable brain vs lung.
#       tissue   PROJECT_PLAN §6's two-matrix approach: Lung_HCA for TIME-L,
#                Brain_Darmanis for TIME-B. Kept because the plan specifies it,
#                and reported with the caveat below.
#
#     The caveat is not generic and it is not a footnote. Lung_HCA resolves
#     CD8+ T, CD4+ T, Treg, NK, B, plasma and pDC; Brain_Darmanis -- a GBM
#     infiltrating-front reference, not normal brain -- carries Neoplastic,
#     Myeloid, Astrocyte, OPC, Oligodendrocyte, Neuron and Vascular types and
#     NO LYMPHOID TYPE AT ALL. The tissue arm therefore cannot structurally
#     show a lymphoid difference between sites: a lymphoid gap there would be
#     the reference, not the biology. The script counts the lymphoid types in
#     each reference and writes the count out, so that sentence rests on a
#     measured number rather than on prose.
#
#   * `raw` IS DELIBERATELY OMITTED from the spatialdecon() call. It would give
#     SpatialDecon its Poisson error model and per-observation weights, but raw
#     gene-level counts do not exist for this dataset: Q1 established the DCCs
#     are probe-level and GEO carries no PKC, so gene counts are not
#     recoverable even in principle. The fit is unweighted and a reader should
#     know that rather than assume the error model was used.
#
#   * Background is each AOI's NegProbe-WTX level broadcast across genes. This
#     is the same background quantity the project's detection rule uses
#     (qc.detection_background_multiple x negprobe), so deconvolution and
#     detection are anchored to one definition rather than two.
#
#   * Composition is renormalised to proportions WITHIN each AOI, and every
#     comparison is made within compartment. That is P2-T5's "renormalise
#     within compartment" requirement, made explicit here rather than taken
#     from SpatialDecon's own prop_of_all so the arithmetic is visible.
#
#   * Input arrives as plain TSV from p0t7_export_tsv -- expr_q3.tsv (linear
#     Q3, which is what spatialdecon's `norm` wants) plus obs.tsv and var.tsv.
#     That is exactly the boundary ADR 0001 designed, and it means no
#     zellkonverter and no second export.
#
# TIME-B n = 8 throughout.

suppressPackageStartupMessages(library(SpatialDecon))

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

emit <- function(...) cat(..., "\n", sep = "")

aoi_codes <- unlist(snakemake@params[["aoi_codes"]])
ref_by_code <- snakemake@params[["ref_by_code"]]
primary_reference <- snakemake@params[["primary_reference"]]

# Anything a lymphoid label could plausibly be called across two independently
# curated references. Deliberately broad: over-matching would understate the
# asymmetry the caveat rests on, so a false positive here is the safe error.
LYMPHOID <- "T\\.cell|T\\.CD|CD4|CD8|B\\.cell|B\\.naive|B\\.memory|^B$|NK|natural.killer|plasma|Treg|regulatory\\.T|lymph"

# ------------------------------------------------------------------- inputs
obs <- read.delim(snakemake@input[["obs"]], check.names = FALSE)
var <- read.delim(snakemake@input[["var"]], check.names = FALSE)
q3 <- as.matrix(read.delim(snakemake@input[["expr_q3"]], header = FALSE))

if (nrow(q3) != nrow(obs) || ncol(q3) != nrow(var)) {
  stop("expr_q3.tsv is ", nrow(q3), " x ", ncol(q3), " but obs/var say ",
       nrow(obs), " x ", nrow(var), " -- the TSV handoff is misaligned.")
}

norm_all <- t(q3)
rownames(norm_all) <- var$gene
colnames(norm_all) <- obs$aoi_label

keep <- obs$aoi_code %in% aoi_codes
o_all <- obs[keep, ]
norm_all <- norm_all[, keep, drop = FALSE]

emit("AOIs: ", ncol(norm_all), " across ", paste(aoi_codes, collapse = ", "))
print(table(o_all$aoi_code))
if (ncol(norm_all) != 23L) {
  stop("expected 23 TIME AOIs (TIME-L 15 + TIME-B 8), got ", ncol(norm_all))
}
emit("")

make_bg <- function(n, negprobe) {
  matrix(rep(negprobe, each = nrow(n)), nrow = nrow(n), dimnames = dimnames(n))
}

load_profile <- function(path) {
  e <- new.env()
  load(path, envir = e)
  # CellProfileLibrary stores profile_matrix as a data.frame; spatialdecon
  # requires a matrix and errors with "X should be a matrix" otherwise.
  as.matrix(e$profile_matrix)
}

to_long <- function(beta, o, arm, reference) {
  prop <- sweep(beta, 2, colSums(beta), "/")
  stopifnot(all(abs(colSums(prop) - 1) < 1e-8))
  idx <- match(colnames(beta), o$aoi_label)
  data.frame(
    arm = arm,
    reference = reference,
    aoi_label = rep(colnames(beta), each = nrow(beta)),
    aoi_code = rep(o$aoi_code[idx], each = nrow(beta)),
    site = rep(o$site[idx], each = nrow(beta)),
    patient_id = rep(o$patient_id[idx], each = nrow(beta)),
    cell_type = rep(rownames(beta), times = ncol(beta)),
    abundance = as.vector(beta),
    proportion = as.vector(prop),
    stringsAsFactors = FALSE
  )
}

lymphoid_types <- function(x) grep(LYMPHOID, rownames(x), ignore.case = TRUE, value = TRUE)

parts <- list()
ref_meta <- list()

# ------------------------------------------------- arm 1: safeTME, both sites
data(safeTME, envir = environment())
data(safeTME.matches, envir = environment())

emit("=== arm 'safetme' -- ", primary_reference, ", applied to BOTH sites ===")
bg <- make_bg(norm_all, o_all$negprobe)
fit <- spatialdecon(norm = norm_all, bg = bg, X = safeTME,
                    cellmerges = safeTME.matches)
beta <- fit$beta
emit("  ", nrow(beta), " cell types x ", ncol(beta), " AOIs")
emit("  types: ", paste(rownames(beta), collapse = ", "))
lym <- lymphoid_types(beta)
emit("  lymphoid types: ", length(lym), " (", paste(lym, collapse = ", "), ")")
emit("  -> ONE reference for both sites, so this arm's composition IS ",
     "comparable brain vs lung.")
emit("")
parts[[length(parts) + 1L]] <- to_long(beta, o_all, "safetme", primary_reference)
ref_meta[["safetme"]] <- list(reference = primary_reference,
                              n_types = nrow(beta), n_lymphoid = length(lym))

# ------------------------------------ arm 2: tissue references, one per site
emit("=== arm 'tissue' -- PROJECT_PLAN §6's two-matrix approach ===")
for (code in aoi_codes) {
  path <- ref_by_code[[code]]
  ref_name <- sub("\\.RData$", "", basename(path))
  sel <- o_all$aoi_code == code
  n <- norm_all[, sel, drop = FALSE]
  o <- o_all[sel, ]
  X <- load_profile(path)
  fit <- spatialdecon(norm = n, bg = make_bg(n, o$negprobe), X = X)
  b <- fit$beta
  lym <- lymphoid_types(b)
  emit("  ", code, " vs ", ref_name, ": ", nrow(b), " types x ", ncol(b), " AOIs")
  emit("    types: ", paste(rownames(b), collapse = ", "))
  emit("    lymphoid types: ", length(lym),
       if (length(lym)) paste0(" (", paste(lym, collapse = ", "), ")") else " -- NONE")
  parts[[length(parts) + 1L]] <- to_long(b, o, "tissue", ref_name)
  ref_meta[[paste0("tissue_", code)]] <- list(reference = ref_name,
                                              n_types = nrow(b),
                                              n_lymphoid = length(lym))
}

tissue_lym <- vapply(ref_meta[grep("^tissue_", names(ref_meta))],
                     function(x) x$n_lymphoid, integer(1))
emit("")
emit("CAVEAT, measured rather than asserted: the two tissue references resolve ",
     paste(tissue_lym, collapse = " and "), " lymphoid types respectively.")
if (min(tissue_lym) == 0L) {
  emit("One of them resolves NONE. The tissue arm therefore CANNOT show a ",
       "lymphoid difference between sites -- a gap there would be the ")
  emit("reference, not the biology. Cross-site comparison of this arm is ",
       "reference-confounded and must not be made.")
}
emit("")

composition <- do.call(rbind, parts)
write.table(composition, snakemake@output[["table"]], sep = "\t",
            row.names = FALSE, quote = FALSE)
emit("wrote ", nrow(composition), " composition rows")

# ------------------------------------------------- mean proportion per arm/site
emit("")
emit("mean proportion by arm x site (top types):")
for (a in unique(composition$arm)) {
  for (s in unique(composition$aoi_code)) {
    d <- composition[composition$arm == a & composition$aoi_code == s, ]
    if (!nrow(d)) next
    m <- tapply(d$proportion, d$cell_type, mean)
    m <- sort(m, decreasing = TRUE)[1:min(5, length(m))]
    emit("  ", a, " / ", s, ": ",
         paste(sprintf("%s %.3f", names(m), m), collapse = "  "))
  }
}

# --------------------------------------------------------------- summary json
esc <- function(x) gsub('"', '\\\\"', x)
meta_json <- paste(vapply(names(ref_meta), function(k) {
  v <- ref_meta[[k]]
  sprintf('"%s": {"reference": "%s", "n_types": %d, "n_lymphoid": %d}',
          k, esc(v$reference), v$n_types, v$n_lymphoid)
}, character(1)), collapse = ",\n    ")

writeLines(sprintf(paste0(
  '{\n  "n_aoi": %d,\n  "n_time_b": 8,\n',
  '  "renormalised_within": "aoi_then_compared_within_compartment",\n',
  '  "error_model_weights_used": false,\n',
  '  "error_model_note": "raw gene-level counts do not exist (Q1: probe-level DCCs, no PKC), so spatialdecon was called without `raw` and the fit is unweighted",\n',
  '  "cross_site_comparable": {"safetme": true, "tissue": false},\n',
  '  "cross_site_note": "the tissue arm uses a different reference per site; Brain_Darmanis resolves no lymphoid type, so a lymphoid difference there would be the reference, not the biology",\n',
  '  "references": {\n    %s\n  }\n}'),
  ncol(norm_all), meta_json), snakemake@output[["summary"]])

sink(type = "message"); sink(type = "output")
