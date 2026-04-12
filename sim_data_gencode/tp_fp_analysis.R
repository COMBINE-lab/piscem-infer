#!/usr/bin/env Rscript
# Compute TP/FP/FN for transcript detection across methods.
# "Expressed" = estimated TPM > 0; "truly expressed" = ground truth TPM > 0
# in that sample's condition.

library(data.table)

simdir <- "sim_data_gencode"

gt <- fread(file.path(simdir, "ground_truth.csv"))
gt[, short_id := sub(" .*", "", transcript_id)]

sample_info <- fread(file.path(simdir, "sample_info.csv"))

read_quant <- function(path) {
  if (!file.exists(path)) { warning("Not found: ", path); return(NULL) }
  df <- fread(path, sep = "\t")
  setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
  df
}

apply_filter_and_renorm <- function(qdf, keep_set = NULL) {
  est_tpm <- qdf$tpm
  if (!is.null(keep_set)) {
    est_tpm <- est_tpm * (qdf$target_name %in% keep_set)
    total <- sum(est_tpm)
    if (total > 0) est_tpm <- est_tpm * 1e6 / total
  }
  est_tpm
}

# NoCond+AV expressed set for EM+Selection
av_expressed <- character(0)
for (i in 1:nrow(sample_info)) {
  sn <- sample_info$sample_name[i]
  qpath <- file.path(simdir, "quant_hier_nocond_av", sn, paste0(sn, ".quant"))
  qdf <- read_quant(qpath)
  if (!is.null(qdf)) av_expressed <- union(av_expressed, qdf$target_name[qdf$tpm > 0])
}

methods <- list(
  list(name = "Single EM",     dir = "quant_em",                filter_set = NULL),
  list(name = "EM+Selection",  dir = "quant_em",                filter_set = av_expressed),
  list(name = "Oracle Sel",    dir = "quant_oracle_em",         filter_set = NULL),
  list(name = "Cons TPM",      dir = "quant_consensus",         filter_set = NULL),
  list(name = "Cons UES",      dir = "quant_consensus_ues",     filter_set = NULL),
  list(name = "Cons Support",  dir = "quant_consensus_support", filter_set = NULL),
  list(name = "Cons Supp CA",  dir = "quant_consensus_support_condaware", filter_set = NULL),
  list(name = "NullSink .25",  dir = "quant_consensus_nullsink_s025", filter_set = NULL),
  list(name = "NullSink .50",  dir = "quant_consensus_nullsink_s050", filter_set = NULL),
  list(name = "Cons Score",    dir = "quant_consensus_score",   filter_set = NULL),
  list(name = "PenEM Support", dir = "quant_consensus_penem",   filter_set = NULL),
  list(name = "Hier Cond",     dir = "quant_multi",             filter_set = NULL),
  list(name = "Hier NoCond",   dir = "quant_hier_nocond",       filter_set = NULL),
  list(name = "NoCond+AV",     dir = "quant_hier_nocond_av",    filter_set = NULL)
)

cat("Total transcripts in reference: ", nrow(gt), "\n")
cat("Truly expressed (either condition): ", sum(gt$expected_tpm_control > 0 | gt$expected_tpm_treatment > 0), "\n")
cat("Truly expressed (control only):     ", sum(gt$expected_tpm_control > 0), "\n")
cat("Truly expressed (treatment only):   ", sum(gt$expected_tpm_treatment > 0), "\n\n")

# ---- Per-sample TP/FP/FN ----
cat("=== Per-sample detection (averaged across 6 samples) ===\n\n")
cat(sprintf("%-14s | %6s %6s %6s %6s | %8s %8s %8s\n",
            "Method", "TP", "FP", "FN", "TN", "Prec", "Recall", "F1"))
cat(paste0(rep("-", 80), collapse = ""), "\n")

for (meth in methods) {
  all_tp <- c(); all_fp <- c(); all_fn <- c(); all_tn <- c()

  for (i in 1:nrow(sample_info)) {
    sn <- sample_info$sample_name[i]
    cond <- sample_info$condition[i]
    tpm_col <- paste0("expected_tpm_", cond)

    qpath <- file.path(simdir, meth$dir, sn, paste0(sn, ".quant"))
    qdf <- read_quant(qpath)
    if (is.null(qdf)) next

    est_tpm <- apply_filter_and_renorm(qdf, meth$filter_set)

    m <- match(qdf$target_name, gt$short_id)
    valid <- !is.na(m)
    true_tpm <- gt[[tpm_col]][m[valid]]
    est <- est_tpm[valid]

    tp <- sum(est > 0 & true_tpm > 0)
    fp <- sum(est > 0 & true_tpm == 0)
    fn <- sum(est == 0 & true_tpm > 0)
    tn <- sum(est == 0 & true_tpm == 0)

    all_tp <- c(all_tp, tp)
    all_fp <- c(all_fp, fp)
    all_fn <- c(all_fn, fn)
    all_tn <- c(all_tn, tn)
  }

  tp <- round(mean(all_tp))
  fp <- round(mean(all_fp))
  fn <- round(mean(all_fn))
  tn <- round(mean(all_tn))
  prec <- mean(all_tp) / (mean(all_tp) + mean(all_fp))
  rec  <- mean(all_tp) / (mean(all_tp) + mean(all_fn))
  f1   <- 2 * prec * rec / (prec + rec)

  cat(sprintf("%-14s | %6d %6d %6d %6d | %8.4f %8.4f %8.4f\n",
              meth$name, tp, fp, fn, tn, prec, rec, f1))
}

# ---- Also break down by DE status ----
cat("\n=== Detection of DE transcripts (averaged across 6 samples) ===\n\n")
cat(sprintf("%-14s | %6s %6s %6s | %8s %8s\n",
            "Method", "TP_DE", "FN_DE", "N_DE", "Recall", "  (non-DE FP)"))
cat(paste0(rep("-", 65), collapse = ""), "\n")

for (meth in methods) {
  all_tp_de <- c(); all_fn_de <- c(); all_n_de <- c(); all_fp_nonde <- c()

  for (i in 1:nrow(sample_info)) {
    sn <- sample_info$sample_name[i]
    cond <- sample_info$condition[i]
    tpm_col <- paste0("expected_tpm_", cond)

    qpath <- file.path(simdir, meth$dir, sn, paste0(sn, ".quant"))
    qdf <- read_quant(qpath)
    if (is.null(qdf)) next

    est_tpm <- apply_filter_and_renorm(qdf, meth$filter_set)

    m <- match(qdf$target_name, gt$short_id)
    valid <- !is.na(m)
    true_tpm <- gt[[tpm_col]][m[valid]]
    est <- est_tpm[valid]
    is_de <- gt$is_de[m[valid]]

    # DE transcripts that are truly expressed
    de_expr <- is_de & true_tpm > 0
    tp_de <- sum(est[de_expr] > 0)
    fn_de <- sum(est[de_expr] == 0)

    # Non-DE transcripts not truly expressed that get called
    nonde_noexpr <- !is_de & true_tpm == 0
    fp_nonde <- sum(est[nonde_noexpr] > 0)

    all_tp_de <- c(all_tp_de, tp_de)
    all_fn_de <- c(all_fn_de, fn_de)
    all_n_de <- c(all_n_de, sum(de_expr))
    all_fp_nonde <- c(all_fp_nonde, fp_nonde)
  }

  cat(sprintf("%-14s | %6d %6d %6d | %8.4f %8d\n",
              meth$name,
              round(mean(all_tp_de)), round(mean(all_fn_de)), round(mean(all_n_de)),
              mean(all_tp_de) / mean(all_n_de),
              round(mean(all_fp_nonde))))
}

cat("\nDone.\n")
