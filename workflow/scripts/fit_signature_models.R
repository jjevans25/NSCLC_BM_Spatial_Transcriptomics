# P2-T3 -- the lung-vs-brain signature contrast, within the TIME compartment.
#
# Owner task: P2-T3. Driven by rule p2t3_signature_models.
#
# Fits, per signature x scoring method:
#
#   primary      score ~ site + (1|patient_id)
#   sensitivity  score ~ site + dsp_run + (1|patient_id)
#
# Choices worth stating, because each could reasonably go the other way:
#
#   * (1|patient_id) is MANDATORY and is justified on NON-INDEPENDENCE, not on
#     variance share (ADR 0009 s1). Patient is only 12% of per-gene variance and
#     its clustering ARI is 0.003; the load-bearing evidence is the +0.107 rho
#     contrast across the 102 same-patient AOI pairs whose compartments differ.
#     A reader who sees only the ARI will conclude this term can be dropped. It
#     cannot. Changing it is a stop-and-ask.
#
#   * dsp_run is the SENSITIVITY, not the primary. It carries ~10% of variance
#     (P1-T3), which argues for adjusting -- but TIME-B splits 6/2 across the two
#     runs, so the adjusted fit is barely identified. Reporting both quantifies
#     the absorption instead of asserting it, which is the shape P1-T3's S2 used
#     to put +0.014 on batch rather than a claim. ADR 0012.
#
#   * site is releveled to lung, so the reported coefficient is BRAIN MINUS LUNG.
#     Gate 2 expects a NEGATIVE antigen-presentation estimate (reduced antigen
#     presentation in brain). Getting this reference level backwards would
#     invert the gate's verdict, so it is set explicitly and asserted.
#
#   * Singular fits are REPORTED, NEVER SUPPRESSED. With TIME-B run-B n = 2 the
#     sensitivity fit is expected to go singular. That is the evidence for why
#     it is not the primary, not a failure to hide.
#
#   * Satterthwaite df via lmerTest. At n = 23 with 5 repeated patients the
#     denominator df are nowhere near n-2, and reporting a naive t would
#     overstate precision.
#
# No p-value is reported without n, effect size and a confidence interval
# (CLAUDE.md hard constraint 7), so every row here carries all four.

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

emit <- function(...) cat(..., "\n", sep = "")

alpha <- snakemake@params[["fdr_alpha"]]
crosscheck_set <- snakemake@params[["crosscheck_signature"]]
df_method <- snakemake@params[["df_method"]]

# The df below are read straight off summary(m)$coefficients, which lmerTest
# computes by Satterthwaite unless asked otherwise. A config that said
# kenward-roger would therefore be silently ignored -- the field would look
# load-bearing and not be. Fail instead, so the config means what it says.
if (!identical(df_method, "satterthwaite")) {
  stop("contexture.model.df_method is '", df_method, "', but this script reads ",
       "df from summary(m)$coefficients, which is Satterthwaite. Implement the ",
       "requested method (summary(m, ddf = \"Kenward-Roger\")) before changing ",
       "the config -- do not let the field drift out of agreement with the code.")
}

# ---------------------------------------------------------------- read + verify
# ADR 0001 puts the R/Python boundary at plain files with digits = 17 and an
# asserted round-trip. Python wrote these with %.17g. What is asserted on this
# side is NOT byte-exact agreement, and the reason is measured rather than
# assumed:
#
#   R's number parser is not correctly rounded. `as.numeric`, `scan` and
#   `read.table` all go through R_strtod, and all three read Python's
#   "0.44525532065683371" as 0x1.c7f102c25d47bp-2 where the correctly-rounded
#   value -- what Python wrote, and what Python's own parser recovers -- is
#   0x1.c7f102c25d47cp-2. One ULP low, on 34 of 276 values.
#
# So exactness holds R -> Python (Python's parser IS correctly rounded, which is
# why assemble_h5ad.py asserts it and passes) but is unattainable Python -> R.
# Nor is it rescuable by comparing the two 17-digit renderings as strings: %g
# strips trailing zeros, so a 1-ULP difference can change the string LENGTH
# ("0.09206579933276339" vs "0.092065799332763404") and any prefix rule then
# fails on formatting rather than on precision.
#
# Two things guard the boundary instead, and together they are stronger than the
# string comparison was:
#
#   1. Structural checks here -- every value parses, is finite, and the frame has
#      the expected shape. This is what actually catches the realistic failures:
#      a truncated file, a shifted column, a locale decimal comma, text where a
#      double belongs.
#   2. rule p2t3_model_crosscheck, which refits the Gate 2 signature in
#      statsmodels from THIS SAME TSV and requires the estimate to agree with
#      lme4's to 1e-4. That is an end-to-end check across the boundary and it is
#      not circular, unlike any comparison R can make against its own parse.
#
# A 1-ULP difference is ~1e-16 relative, against an lmer convergence tolerance
# of ~1e-8. ADR 0013 records this.
raw <- read.delim(snakemake@input[["scores"]], colClasses = "character",
                  check.names = FALSE)
scores <- raw
scores$score <- suppressWarnings(as.numeric(raw$score))

required <- c("aoi_label", "signature", "method", "score", "site",
              "patient_id", "dsp_run")
missing_cols <- setdiff(required, names(scores))
if (length(missing_cols)) {
  stop("signature_scores.tsv is missing columns: ",
       paste(missing_cols, collapse = ", "))
}

n_bad <- sum(is.na(scores$score) | !is.finite(scores$score))
if (n_bad > 0L) {
  ex <- raw$score[which(is.na(scores$score) | !is.finite(scores$score))[1]]
  stop(n_bad, " of ", nrow(raw), " score values did not parse to a finite ",
       "double (first offending token: '", ex, "'). This is a structural ",
       "failure at the R/Python boundary, not parser rounding (ADR 0001).")
}

n_exact <- sum(raw$score == sprintf("%.17g", scores$score))
emit("boundary: ", nrow(scores), " values parsed, all finite; ", n_exact,
     " of them reprint byte-identically at %.17g.")
emit("The remainder differ by one ULP because R_strtod is not correctly ",
     "rounded -- see the header and ADR 0013. Numeric agreement across the ")
emit("boundary is checked end-to-end by rule p2t3_model_crosscheck.")

scores$patient_id <- factor(scores$patient_id)
scores$dsp_run <- factor(scores$dsp_run)

# Brain minus lung. Asserted rather than assumed -- see the header.
scores$site <- factor(scores$site, levels = c("lung", "brain"))
stopifnot(levels(scores$site)[1] == "lung")
emit("site reference level: ", levels(scores$site)[1],
     " -> coefficient 'sitebrain' is BRAIN MINUS LUNG")

n_by_site <- table(scores$site[!duplicated(scores$aoi_label)])
emit("AOIs: ", paste(names(n_by_site), n_by_site, sep = " = ", collapse = ", "))
emit("TIME-B n = 8 is the binding constraint; power floor ~1.1-1.3 SD at 80%.")
emit("")

signatures <- sort(unique(scores$signature))
methods <- sort(unique(scores$method))

# ------------------------------------------------------------------- fit
fit_one <- function(d, formula_text, label) {
  # Warnings and messages are trapped SEPARATELY and reported in separate
  # columns. lme4 signals a singular fit with message(), not warning() -- so a
  # warning-only handler leaves the `warnings` column empty on every row while
  # "boundary (singular) fit" scrolls past on the console, which reads as "no
  # fit had anything to say". Conflating the two condition classes into one
  # column would hide which mechanism actually fired.
  warns <- character(0)
  msgs <- character(0)
  m <- withCallingHandlers(
    lmerTest::lmer(as.formula(formula_text), data = d, REML = TRUE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      msgs <<- c(msgs, sub("\n$", "", conditionMessage(m)))
      invokeRestart("muffleMessage")
    }
  )
  co <- summary(m)$coefficients
  term <- "sitebrain"
  if (!(term %in% rownames(co))) {
    stop("term ", term, " absent from ", label, " -- site level order changed?")
  }
  est <- co[term, "Estimate"]
  se <- co[term, "Std. Error"]
  df <- co[term, "df"]
  tcrit <- qt(1 - (1 - 0.95) / 2, df)
  list(
    estimate = est,
    std_error = se,
    df = df,
    t_value = co[term, "t value"],
    p_raw = co[term, "Pr(>|t|)"],
    ci_low = est - tcrit * se,
    ci_high = est + tcrit * se,
    is_singular = isSingular(m),
    n_obs = nrow(d),
    n_patient = nlevels(droplevels(d$patient_id)),
    warnings = if (length(warns)) paste(unique(warns), collapse = " | ") else "",
    messages = if (length(msgs)) paste(unique(msgs), collapse = " | ") else ""
  )
}

rows <- list()
for (meth in methods) {
  for (sig in signatures) {
    d <- scores[scores$method == meth & scores$signature == sig, ]
    for (model in c("primary", "sensitivity")) {
      f <- snakemake@params[[model]]
      r <- fit_one(d, f, paste(sig, meth, model))
      rows[[length(rows) + 1L]] <- data.frame(
        signature = sig, method = meth, model = model, formula = f,
        term = "brain_vs_lung",
        estimate = r$estimate, std_error = r$std_error, df = r$df,
        ci_low = r$ci_low, ci_high = r$ci_high,
        t_value = r$t_value, p_raw = r$p_raw,
        is_singular = r$is_singular,
        n_obs = r$n_obs, n_patient = r$n_patient,
        n_lung = sum(d$site == "lung"), n_brain = sum(d$site == "brain"),
        df_method = df_method,
        warnings = r$warnings, messages = r$messages,
        stringsAsFactors = FALSE
      )
    }
  }
}
res <- do.call(rbind, rows)

# BH across signatures, within each method x model family. Correcting across
# models too would treat the sensitivity refit as an extra hypothesis, which it
# is not -- it is the same hypothesis under a different adjustment.
res$q_bh <- NA_real_
for (meth in methods) {
  for (model in c("primary", "sensitivity")) {
    i <- res$method == meth & res$model == model
    res$q_bh[i] <- p.adjust(res$p_raw[i], method = "BH")
  }
}

# How much does the batch adjustment actually move the estimate?
key <- paste(res$signature, res$method)
prim <- res[res$model == "primary", ]
sens <- res[res$model == "sensitivity", ]
mm <- match(paste(sens$signature, sens$method), paste(prim$signature, prim$method))
res$delta_vs_primary <- NA_real_
res$delta_vs_primary[res$model == "sensitivity"] <- sens$estimate - prim$estimate[mm]

res <- res[order(res$method, res$model, res$signature), ]
write.table(res, snakemake@output[["table"]], sep = "\t",
            row.names = FALSE, quote = FALSE)

# ------------------------------------------------------------------ report
emit("brain vs lung, per signature (estimate [95% CI], Satterthwaite df):")
for (meth in methods) {
  emit("")
  emit("  method = ", meth)
  for (model in c("primary", "sensitivity")) {
    emit("    ", model, ":  ", snakemake@params[[model]])
    sub <- res[res$method == meth & res$model == model, ]
    for (i in seq_len(nrow(sub))) {
      r <- sub[i, ]
      emit(sprintf(
        "      %-22s %+.3f [%+.3f, %+.3f]  df=%5.1f  p=%.4f  q=%.4f%s",
        r$signature, r$estimate, r$ci_low, r$ci_high, r$df, r$p_raw, r$q_bh,
        if (r$is_singular) "  SINGULAR" else ""))
      if (nzchar(r$warnings)) emit("        warning: ", r$warnings)
    }
  }
}

n_sing <- sum(res$is_singular)
emit("")
emit("singular fits: ", n_sing, " of ", nrow(res),
     " (surfaced, not suppressed -- ADR 0012)")

sens_rows <- res[res$model == "sensitivity", ]
emit("batch absorption |sensitivity - primary|: median ",
     sprintf("%.4f", median(abs(sens_rows$delta_vs_primary))),
     ", max ", sprintf("%.4f", max(abs(sens_rows$delta_vs_primary))))

emit("")
emit("Gate 2 turns on ", crosscheck_set,
     ": the published direction is REDUCED antigen presentation in brain,")
emit("i.e. a NEGATIVE estimate here. Observed (primary):")
g <- res[res$signature == crosscheck_set & res$model == "primary", ]
for (i in seq_len(nrow(g))) {
  emit(sprintf("  %-8s %+.3f [%+.3f, %+.3f]  p=%.4f  q=%.4f",
               g$method[i], g$estimate[i], g$ci_low[i], g$ci_high[i],
               g$p_raw[i], g$q_bh[i]))
}

# --------------------------------------------------------------- summary json
esc <- function(x) gsub('"', '\\\\"', x)
json_rows <- apply(res, 1, function(r) {
  sprintf(paste0('{"signature":"%s","method":"%s","model":"%s","estimate":%s,',
                 '"std_error":%s,"df":%s,"ci_low":%s,"ci_high":%s,"p_raw":%s,',
                 '"q_bh":%s,"is_singular":%s,"n_lung":%s,"n_brain":%s}'),
          r[["signature"]], r[["method"]], r[["model"]],
          r[["estimate"]], r[["std_error"]], r[["df"]], r[["ci_low"]],
          r[["ci_high"]], r[["p_raw"]], r[["q_bh"]],
          tolower(r[["is_singular"]]), r[["n_lung"]], r[["n_brain"]])
})
writeLines(sprintf(paste0('{\n  "reference_level": "lung",\n',
  '  "contrast": "brain_minus_lung",\n',
  '  "n_time_b": 8,\n  "fdr_alpha": %s,\n  "fdr_method": "BH",\n',
  '  "df_method": "%s",\n  "n_singular": %s,\n',
  '  "gate2_signature": "%s",\n  "results": [\n    %s\n  ]\n}'),
  alpha, df_method, n_sing, crosscheck_set, paste(json_rows, collapse = ",\n    ")),
  snakemake@output[["summary"]])

emit("")
emit("wrote ", nrow(res), " model rows")
sink(type = "message"); sink(type = "output")
