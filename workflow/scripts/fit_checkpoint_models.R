# P3-T3 / P3-T4 -- the checkpoint contrasts, primary + declared exploratory.
#
# Owner task: P3-T3 (carrier: which compartment carries each checkpoint, within
# site). Written parameterised so P3-T4 (shift: lung -> brain, within
# compartment) drives the same script with a different grouping variable -- the
# arithmetic is identical and a second copy would be a second place for the
# reference level to drift.
#
# A4 IS EXPLORATORY (ADR 0008). Everything here is computed; the demotion is
# about what may be CLAIMED. Three rules bind every row:
#
#   1. Every gene reported states its per-group detection.
#   2. A gene below the pre-registered floor is "not assessable in <group>",
#      NEVER "lower in <group>". Background is HIGHER in the immune compartments
#      (median negprobe_log2: L 4.96 -> TIME-L 5.36), so a near-background gene
#      reads as depleted there artefactually.
#   3. Every brain claim states TIME-B n = 8 inline (hard constraint 8).
#
# TWO TABLES, AND THE DIFFERENCE BETWEEN THEM IS THE POINT (ADR 0018):
#
#   PRIMARY      genes clearing the floor in BOTH groups. ADR 0016 s3 exactly as
#                pre-registered, before p3t2_checkpoint_detection existed. The
#                ONLY table in Phase 3 carrying a q-value. Genes that fail get a
#                not_assessable row with their detection counts and NO estimate
#                -- present in the table, never silently absent.
#
#   EXPLORATORY  genes clearing the floor in AT LEAST ONE group, identical
#                formulas, NO FDR OF ANY KIND, every row labelled exploratory
#                and carrying the background-gradient direction. This is a
#                second step down from an already-exploratory aim: no P3-T7
#                sentence may rest on it.
#
# The exploratory table also refits the genes the primary accepted, marked
# in_primary. Those estimates MUST equal the primary's and are asserted to 1e-10
# -- it re-derives the primary through the secondary's code path, so the two
# tables cannot silently disagree, and it is what makes them comparable rather
# than rival.
#
# WHY NOT JUST RELAX THE PRIMARY. Because the restriction is symmetric and the
# artefact is directional: where the background gradient is large it runs
# AGAINST the observed direction (background is higher in immune, yet these
# genes are enriched there -- VSIR 8/30 -> 14/15, HAVCR2 2/30 -> 10/15), so the
# artefact can only shrink those effects. Admitting them to the primary on that
# reasoning would still be a threshold changed after seeing the audit, and a
# reader cannot tell that from relaxing it until something appeared. ADR 0018
# records the decision and the rejected alternatives.
#
# gradient_direction is computed per row rather than argued in prose:
#   protective  gradient runs against the estimate -- artefact can only shrink it
#   permissive  gradient runs with the estimate -- artefact could contribute
#   negligible  |background_delta| < GRADIENT_NEGLIGIBLE
#
# Choices inherited from fit_signature_models.R, each for its stated reason:
#
#   * (1|patient_id) is MANDATORY and justified on NON-INDEPENDENCE, not variance
#     share (ADR 0009 s1). Patient is 12% of per-gene variance and its clustering
#     ARI is 0.003; the load-bearing evidence is the +0.107 rho contrast across
#     the 102 same-patient AOI pairs whose compartments differ. Singular fits are
#     NOT a reason to drop it -- half of Phase 2's fits were singular and it
#     stayed. Changing it is a stop-and-ask.
#   * dsp_run and negprobe_log2 are SENSITIVITIES, never the primary (ADR 0012,
#     ADR 0016 s3). Both are reported with their shift from the primary, which
#     quantifies the confound instead of asserting it away.
#   * Warnings and messages are trapped SEPARATELY. lme4 signals a singular fit
#     with message(), not warning(); a warning-only handler leaves the column
#     empty on every row while "boundary (singular) fit" scrolls past.
#   * Satterthwaite df via lmerTest.
#
# No p-value without n, effect size and a confidence interval (hard constraint 7).

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

emit <- function(...) cat(..., "\n", sep = "")

alpha <- snakemake@params[["fdr_alpha"]]
df_method <- snakemake@params[["df_method"]]
contrast_var <- snakemake@params[["contrast_var"]]   # "compartment" | "site"
reference <- snakemake@params[["reference_level"]]   # "tumour" | "lung"
test_level <- snakemake@params[["test_level"]]       # "immune" | "brain"
strata_var <- snakemake@params[["strata_var"]]       # "site" | "compartment"
gene_floor <- snakemake@params[["detected_in_aoi_fraction"]]
task <- snakemake@params[["task"]]

# |background_delta| below this is reported as negligible rather than as a
# direction. 0.1 log2 is ~7% on the background level; the measured compartment
# gradients are 0.31-0.40 and the site gradients 0.03-0.06, so this cleanly
# separates the two cases it exists to separate.
GRADIENT_NEGLIGIBLE <- 0.1

if (!identical(df_method, "satterthwaite")) {
  stop("checkpoints.model.df_method is '", df_method, "', but this script reads ",
       "df from summary(m)$coefficients, which is Satterthwaite. Implement the ",
       "requested method before changing the config -- do not let the field ",
       "drift out of agreement with the code.")
}

# ---------------------------------------------------------------- read + verify
# ADR 0001's boundary. What is asserted is structural, not byte-exact: R_strtod
# is not correctly rounded, so exactness is unattainable Python -> R (ADR 0013
# s3). The numeric check that matters is end-to-end, in
# p3t3_checkpoint_crosscheck.
raw <- read.delim(snakemake@input[["expression"]], colClasses = "character",
                  check.names = FALSE)
d <- raw
for (col in c("expression", "negprobe_log2")) {
  d[[col]] <- suppressWarnings(as.numeric(raw[[col]]))
}
d$detected <- tolower(raw$detected) %in% c("true", "1")

required <- c("aoi_label", "patient_id", "aoi_code", "compartment", "site",
              "dsp_run", "gene", "expression", "detected", "negprobe_log2")
missing_cols <- setdiff(required, names(d))
if (length(missing_cols)) {
  stop("checkpoint_expression.tsv is missing columns: ",
       paste(missing_cols, collapse = ", "))
}

n_bad <- sum(!is.finite(d$expression) | !is.finite(d$negprobe_log2))
if (n_bad > 0L) {
  stop(n_bad, " of ", nrow(raw), " values did not parse to a finite double. ",
       "This is a structural failure at the R/Python boundary, not parser ",
       "rounding (ADR 0001).")
}
emit("boundary: ", nrow(d), " rows parsed, all finite (ADR 0001; numeric ",
     "agreement is checked end-to-end by p3t3_checkpoint_crosscheck)")

d$patient_id <- factor(d$patient_id)
d$dsp_run <- factor(d$dsp_run)

# The contrast factor is releveled and ASSERTED. fit_signature_models.R asserts
# the same thing for `site`; getting it backwards would invert the sign of the
# phase's central number and every downstream sentence with it.
d[[contrast_var]] <- factor(d[[contrast_var]], levels = c(reference, test_level))
stopifnot(levels(d[[contrast_var]])[1] == reference)
term <- paste0(contrast_var, test_level)
emit(task, ": ", contrast_var, " reference level is '", reference,
     "' -> coefficient '", term, "' is ",
     toupper(test_level), " MINUS ", toupper(reference))

# `compartment` is DEGENERATE across sites -- L, LB and mLN are all `tumour`.
# Subsetting by strata_var first is what makes it a clean two-level factor
# within each stratum. Never fit across strata.
strata <- sort(unique(d[[strata_var]]))
emit("strata (", strata_var, "): ", paste(strata, collapse = ", "),
     " -- fitted separately, NEVER pooled")

genes <- sort(unique(d$gene))
emit("panel: ", length(genes), " genes (PRE-REGISTERED, ADR 0015)")
emit("detection floor: gene present in a group when detected in >= ",
     gene_floor, " of its AOIs (PRE-REGISTERED, ADR 0016)")
emit("")

# ------------------------------------------------------------- assessability
# Recomputed from the `detected` column that p3t2_checkpoint_detection wrote,
# NOT from a second detection rule. checkpoint_detection.tsv already holds the
# same counts; this derives them per (stratum, gene, group) because that is the
# grouping the models need and the audit is keyed on aoi_code.
assess <- list()
for (s in strata) {
  for (g in genes) {
    for (lvl in c(reference, test_level)) {
      i <- d[[strata_var]] == s & d$gene == g & d[[contrast_var]] == lvl
      assess[[paste(s, g, lvl)]] <- list(
        n = sum(i), n_detected = sum(d$detected[i]),
        rate = mean(d$detected[i]),
        code = unique(d$aoi_code[i]),
        background = median(d$negprobe_log2[i])
      )
    }
  }
}
present <- function(s, g, lvl) assess[[paste(s, g, lvl)]]$rate >= gene_floor

# ------------------------------------------------------------------- fit
fit_one <- function(dd, formula_text, label) {
  warns <- character(0)
  msgs <- character(0)
  m <- withCallingHandlers(
    lmerTest::lmer(as.formula(formula_text), data = dd, REML = TRUE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(mm) {
      msgs <<- c(msgs, sub("\n$", "", conditionMessage(mm)))
      invokeRestart("muffleMessage")
    }
  )
  co <- summary(m)$coefficients
  if (!(term %in% rownames(co))) {
    stop("term ", term, " absent from ", label, " -- factor level order changed?")
  }
  est <- co[term, "Estimate"]
  se <- co[term, "Std. Error"]
  df <- co[term, "df"]
  tcrit <- qt(0.975, df)
  list(estimate = est, std_error = se, df = df,
       t_value = co[term, "t value"], p_raw = co[term, "Pr(>|t|)"],
       ci_low = est - tcrit * se, ci_high = est + tcrit * se,
       is_singular = isSingular(m), n_obs = nrow(dd),
       n_patient = nlevels(droplevels(dd$patient_id)),
       warnings = if (length(warns)) paste(unique(warns), collapse = " | ") else "",
       messages = if (length(msgs)) paste(unique(msgs), collapse = " | ") else "")
}

models <- c(primary = snakemake@params[["primary"]],
            batch_sensitivity = snakemake@params[["batch_sensitivity"]],
            background_sensitivity = snakemake@params[["background_sensitivity"]])

blank_row <- function(s, g, status, note) {
  a_ref <- assess[[paste(s, g, reference)]]
  a_test <- assess[[paste(s, g, test_level)]]
  data.frame(
    task = task, stratum = s, gene = g, model = NA_character_,
    formula = NA_character_, term = NA_character_,
    estimate = NA_real_, std_error = NA_real_, df = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_, t_value = NA_real_,
    p_raw = NA_real_, is_singular = NA, n_obs = NA_integer_,
    n_patient = NA_integer_,
    group_reference = reference, group_test = test_level,
    aoi_code_reference = a_ref$code, aoi_code_test = a_test$code,
    n_reference = a_ref$n, n_test = a_test$n,
    n_detected_reference = a_ref$n_detected, n_detected_test = a_test$n_detected,
    detection_rate_reference = a_ref$rate, detection_rate_test = a_test$rate,
    assessable_reference = a_ref$rate >= gene_floor,
    assessable_test = a_test$rate >= gene_floor,
    background_reference = a_ref$background, background_test = a_test$background,
    background_delta = a_test$background - a_ref$background,
    gradient_direction = NA_character_,
    status = status, note = note,
    df_method = df_method, warnings = "", messages = "",
    stringsAsFactors = FALSE
  )
}

fit_rows <- function(s, g, status) {
  a_ref <- assess[[paste(s, g, reference)]]
  a_test <- assess[[paste(s, g, test_level)]]
  bg_delta <- a_test$background - a_ref$background
  dd <- d[d[[strata_var]] == s & d$gene == g, ]
  out <- list()
  for (model in names(models)) {
    r <- fit_one(dd, models[[model]], paste(task, s, g, model))
    # Direction is relative to THIS row's estimate: the artefact pushes the
    # estimate in the direction of the background shift, so a gradient with the
    # same sign as the estimate could have contributed to it (permissive), and
    # one with the opposite sign can only have shrunk it (protective).
    direction <- if (abs(bg_delta) < GRADIENT_NEGLIGIBLE) {
      "negligible"
    } else if (sign(bg_delta) == sign(r$estimate)) {
      "permissive"
    } else {
      "protective"
    }
    row <- blank_row(s, g, status, "")
    row$model <- model
    row$formula <- models[[model]]
    row$term <- paste0(test_level, "_vs_", reference)
    row$estimate <- r$estimate; row$std_error <- r$std_error; row$df <- r$df
    row$ci_low <- r$ci_low; row$ci_high <- r$ci_high
    row$t_value <- r$t_value; row$p_raw <- r$p_raw
    row$is_singular <- r$is_singular
    row$n_obs <- r$n_obs; row$n_patient <- r$n_patient
    row$gradient_direction <- direction
    row$warnings <- r$warnings; row$messages <- r$messages
    out[[length(out) + 1L]] <- row
  }
  do.call(rbind, out)
}

# --------------------------------------------------------------- PRIMARY
# ADR 0016 s3, unchanged: BOTH groups must clear the floor.
emit("PRIMARY -- genes clearing the floor in BOTH groups (ADR 0016 s3,")
emit("pre-registered before the audit existed). The only table with a q-value.")
prim_rows <- list()
for (s in strata) {
  fitted <- character(0)
  for (g in genes) {
    if (present(s, g, reference) && present(s, g, test_level)) {
      prim_rows[[length(prim_rows) + 1L]] <- fit_rows(s, g, "modelled")
      fitted <- c(fitted, g)
    } else {
      a_ref <- assess[[paste(s, g, reference)]]
      a_test <- assess[[paste(s, g, test_level)]]
      where <- if (!present(s, g, reference) && !present(s, g, test_level)) {
        paste0("not assessable in either group (", a_ref$code, " ",
               a_ref$n_detected, "/", a_ref$n, ", ", a_test$code, " ",
               a_test$n_detected, "/", a_test$n, ")")
      } else if (!present(s, g, reference)) {
        paste0("not assessable in ", a_ref$code, " (", a_ref$n_detected, "/",
               a_ref$n, ")")
      } else {
        paste0("not assessable in ", a_test$code, " (", a_test$n_detected, "/",
               a_test$n, ")")
      }
      prim_rows[[length(prim_rows) + 1L]] <- blank_row(s, g, "not_assessable", where)
    }
  }
  emit("  ", s, ": ", length(fitted), " of ", length(genes), " genes modelled",
       if (length(fitted)) paste0(" (", paste(fitted, collapse = ", "), ")") else "")
}
primary <- do.call(rbind, prim_rows)

# BH across GENES, within stratum, within model. T3 and T4 are different
# families (config comment); the sensitivity refits are the same hypothesis
# under a different adjustment, not extra hypotheses, so they are corrected
# within model rather than across models.
primary$q_bh <- NA_real_
for (s in strata) {
  for (model in names(models)) {
    i <- primary$stratum == s & primary$model == model & !is.na(primary$p_raw)
    if (any(i)) primary$q_bh[i] <- p.adjust(primary$p_raw[i], method = "BH")
  }
}

# ----------------------------------------------------------- EXPLORATORY
# ADR 0018. At least one group, identical formulas, NO FDR.
emit("")
emit("EXPLORATORY -- genes clearing the floor in AT LEAST ONE group (ADR 0018).")
emit("NO q-value, NO FDR. Exploratory-within-exploratory: no P3-T7 sentence")
emit("may rest on this table.")
expl_rows <- list()
for (s in strata) {
  fitted <- character(0)
  for (g in genes) {
    if (present(s, g, reference) || present(s, g, test_level)) {
      rows <- fit_rows(s, g, "exploratory")
      rows$in_primary <- present(s, g, reference) && present(s, g, test_level)
      expl_rows[[length(expl_rows) + 1L]] <- rows
      fitted <- c(fitted, g)
    }
  }
  emit("  ", s, ": ", length(fitted), " of ", length(genes), " genes fitted",
       if (length(fitted)) paste0(" (", paste(fitted, collapse = ", "), ")") else "")
}
exploratory <- do.call(rbind, expl_rows)
exploratory$q_bh <- NULL

# The overlap MUST agree with the primary. This re-derives the primary through
# the secondary's code path; a disagreement means the two tables describe
# different fits and a reader comparing them would be misled.
overlap <- exploratory[exploratory$in_primary, ]
key_e <- paste(overlap$stratum, overlap$gene, overlap$model)
key_p <- paste(primary$stratum, primary$gene, primary$model)
mm <- match(key_e, key_p)
if (any(is.na(mm))) {
  stop("exploratory rows marked in_primary have no primary counterpart: ",
       paste(key_e[is.na(mm)], collapse = "; "))
}
worst <- max(abs(overlap$estimate - primary$estimate[mm]))
emit("")
emit("overlap check: ", nrow(overlap), " exploratory rows are also primary; ",
     "max |difference| in estimate = ", sprintf("%.3e", worst))
if (worst > 1e-10) {
  stop("exploratory and primary disagree by ", worst, " on genes present in ",
       "both. The two tables must describe the same fit for the same gene.")
}

# ---------------------------------------------------------------- sensitivity shift
add_delta <- function(df) {
  df$delta_vs_primary <- NA_real_
  k <- paste(df$stratum, df$gene)
  p_est <- df$estimate[df$model == "primary"]
  p_key <- k[df$model == "primary"]
  for (model in c("batch_sensitivity", "background_sensitivity")) {
    i <- which(df$model == model)
    j <- match(k[i], p_key)
    df$delta_vs_primary[i] <- df$estimate[i] - p_est[j]
  }
  df
}
primary <- add_delta(primary)
exploratory <- add_delta(exploratory)

primary <- primary[order(primary$stratum, primary$model, primary$gene), ]
exploratory <- exploratory[order(exploratory$stratum, exploratory$model,
                                 exploratory$gene), ]
write.table(primary, snakemake@output[["primary"]], sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")
write.table(exploratory, snakemake@output[["exploratory"]], sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")

# ------------------------------------------------------------------ report
report <- function(df, title, show_q) {
  emit("")
  emit(title)
  for (s in strata) {
    for (model in names(models)) {
      sub <- df[df$stratum == s & df$model == model & !is.na(df$estimate), ]
      if (!nrow(sub)) next
      emit("")
      emit("  ", s, " / ", model, ":  ", models[[model]])
      for (i in seq_len(nrow(sub))) {
        r <- sub[i, ]
        emit(sprintf(
          "    %-8s %+.3f [%+.3f, %+.3f]  df=%5.1f  p=%.4f%s  %-11s %2d/%-2d vs %2d/%-2d%s",
          r$gene, r$estimate, r$ci_low, r$ci_high, r$df, r$p_raw,
          if (show_q) sprintf("  q=%.4f", r$q_bh) else "",
          r$gradient_direction,
          r$n_detected_reference, r$n_reference,
          r$n_detected_test, r$n_test,
          if (isTRUE(r$is_singular)) "  SINGULAR" else ""))
      }
    }
  }
}

report(primary, "PRIMARY (q-values are BH across genes, within stratum, within model):", TRUE)
report(exploratory, "EXPLORATORY (NO FDR -- ADR 0018; p_raw is uncorrected):", FALSE)

na_rows <- primary[primary$status == "not_assessable", ]
emit("")
emit("not assessable in the primary: ", nrow(na_rows), " gene x stratum cells")
for (i in seq_len(nrow(na_rows))) {
  emit("  ", na_rows$stratum[i], " ", na_rows$gene[i], ": ", na_rows$note[i])
}

n_sing <- sum(primary$is_singular, na.rm = TRUE) +
  sum(exploratory$is_singular, na.rm = TRUE)
emit("")
emit("singular fits: ", n_sing, " (surfaced, not suppressed -- ADR 0012;")
emit("  NOT a reason to drop (1|patient_id): half of Phase 2's fits were")
emit("  singular and it stayed. ADR 0009 s1.)")

grad <- table(exploratory$gradient_direction[exploratory$model == "primary"])
emit("")
emit("background gradient on the exploratory primary fits:")
for (nm in names(grad)) emit("  ", nm, ": ", grad[[nm]])
emit("")
emit("  DIRECTION, STATED CAREFULLY -- detection and expression move OPPOSITE")
emit("  ways under the same gradient, and conflating them inverts the argument:")
emit("    detection  is q3 > 2x negprobe, so higher background raises the bar")
emit("               and the gene is detected LESS -> looks DEPLETED.")
emit("    expression is log2(q3+1), where background adds to the signal, so")
emit("               higher background INFLATES it -> looks ENRICHED.")
emit("  These models are on expression. Background is higher in the immune")
emit("  compartments, every estimate here is positive, so the gradient is")
emit("  PERMISSIVE: it could have contributed to these effects. ADR 0018.")

# The background sensitivity measures the gradient's contribution directly,
# which is worth more than the categorical label: it is the difference the
# adjustment actually makes, per gene.
bg <- exploratory[exploratory$model == "background_sensitivity", ]
bg <- bg[order(bg$delta_vs_primary), ]
emit("")
emit("adjusting for negprobe_log2 -- shift from the primary estimate")
emit("(negative = background was inflating the compartment effect):")
for (i in seq_len(nrow(bg))) {
  r <- bg[i, ]
  emit(sprintf("  %-6s %-8s %+.3f  (%+.3f -> %+.3f)  %s",
               r$stratum, r$gene, r$delta_vs_primary,
               r$estimate - r$delta_vs_primary, r$estimate,
               if (isTRUE(r$in_primary)) "primary" else "exploratory only"))
}
in_p <- bg$in_primary
emit(sprintf("  mean shift: %+.3f exploratory-only, %+.3f primary",
             mean(bg$delta_vs_primary[!in_p]), mean(bg$delta_vs_primary[in_p])))
emit("  The genes the primary EXCLUDES are the ones the adjustment moves. That")
emit("  is the restriction working, not an argument for relaxing it.")
emit("")
emit("REPORTING RULE (ADR 0008 point 4). A gene below the floor in a group is")
emit("'not assessable in <group>', NEVER 'lower in <group>'. The exploratory")
emit("table reports what a fit WOULD say; it does not promote a gene to")
emit("assessable. Every brain claim states TIME-B n = 8 inline, and the power")
emit("floor is ~1.1-1.3 SD at 80% -- a null here is uninformative, not negative.")

# --------------------------------------------------------------- summary json
esc <- function(x) gsub('"', '\\\\"', as.character(x))
jnum <- function(x) ifelse(is.na(x), "null", sprintf("%.17g", as.numeric(x)))
jrow <- function(r) sprintf(
  paste0('{"stratum":"%s","gene":"%s","model":%s,"estimate":%s,"std_error":%s,',
         '"df":%s,"ci_low":%s,"ci_high":%s,"p_raw":%s,"q_bh":%s,',
         '"n_detected_reference":%s,"n_reference":%s,"n_detected_test":%s,',
         '"n_test":%s,"gradient_direction":%s,"status":"%s"}'),
  esc(r[["stratum"]]), esc(r[["gene"]]),
  if (is.na(r[["model"]])) "null" else sprintf('"%s"', esc(r[["model"]])),
  jnum(r[["estimate"]]), jnum(r[["std_error"]]), jnum(r[["df"]]),
  jnum(r[["ci_low"]]), jnum(r[["ci_high"]]), jnum(r[["p_raw"]]),
  if ("q_bh" %in% names(r)) jnum(r[["q_bh"]]) else "null",
  r[["n_detected_reference"]], r[["n_reference"]],
  r[["n_detected_test"]], r[["n_test"]],
  if (is.na(r[["gradient_direction"]])) "null"
  else sprintf('"%s"', esc(r[["gradient_direction"]])),
  esc(r[["status"]]))

prim_json <- vapply(seq_len(nrow(primary)),
                    function(i) jrow(primary[i, ]), character(1))
expl_json <- vapply(seq_len(nrow(exploratory)),
                    function(i) jrow(exploratory[i, ]), character(1))

writeLines(sprintf(paste0(
  '{\n  "task": "%s",\n  "contrast_var": "%s",\n',
  '  "reference_level": "%s",\n  "test_level": "%s",\n',
  '  "contrast": "%s_minus_%s",\n  "strata_var": "%s",\n',
  '  "n_time_b": 8,\n  "power_floor_sd": "1.1-1.3",\n',
  '  "detected_in_aoi_fraction": %s,\n  "fdr_alpha": %s,\n',
  '  "fdr_method": "BH within stratum within model",\n  "df_method": "%s",\n',
  '  "primary_rule": "gene clears the floor in BOTH groups (ADR 0016 s3, pre-registered)",\n',
  '  "exploratory_rule": "gene clears the floor in AT LEAST ONE group (ADR 0018); NO FDR",\n',
  '  "overlap_max_abs_diff": %s,\n  "n_singular": %s,\n',
  '  "reporting_rule": "A gene below the floor in a group is not assessable there, never lower there (ADR 0008 point 4).",\n',
  '  "primary": [\n    %s\n  ],\n  "exploratory": [\n    %s\n  ]\n}'),
  task, contrast_var, reference, test_level, test_level, reference, strata_var,
  gene_floor, alpha, df_method, sprintf("%.17g", worst), n_sing,
  paste(prim_json, collapse = ",\n    "),
  paste(expl_json, collapse = ",\n    ")),
  snakemake@output[["summary"]])

emit("")
emit("wrote ", nrow(primary), " primary rows and ", nrow(exploratory),
     " exploratory rows")
sink(type = "message"); sink(type = "output")
