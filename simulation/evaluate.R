#!/usr/bin/env Rscript
#
# Evaluate piscem-infer quantification accuracy against ground truth.
#
# Compares single-sample EM vs multi-sample hierarchical estimates.
#
# Usage:
#   Rscript simulation/evaluate.R [sim_data_dir]

args <- commandArgs(trailingOnly = TRUE)
simdir <- if (length(args) >= 1) args[1] else "sim_data"

cat("=== Evaluation ===\n")
cat("Data directory:", simdir, "\n\n")

# ---- Load ground truth ----

gt_file <- file.path(simdir, "ground_truth.csv")
if (!file.exists(gt_file)) {
  stop("Ground truth file not found: ", gt_file)
}
gt <- read.csv(gt_file, stringsAsFactors = FALSE)
cat("Loaded ground truth for", nrow(gt), "transcripts\n")

# ---- Load sample info ----

sample_info <- read.csv(file.path(simdir, "sample_info.csv"), stringsAsFactors = FALSE)
cat("Found", nrow(sample_info), "samples\n\n")

# ---- Helper: read .quant file ----

read_quant <- function(path) {
  if (!file.exists(path)) {
    warning("Quant file not found: ", path)
    return(NULL)
  }
  df <- read.delim(path, stringsAsFactors = FALSE)
  colnames(df) <- c("target_name", "len", "eelen", "tpm", "ecount")
  df
}

# ---- Compare function ----

compare_to_gt <- function(quant_df, gt, condition) {
  # Get expected TPM for this condition
  tpm_col <- paste0("expected_tpm_", condition)
  if (!tpm_col %in% colnames(gt)) {
    warning("Column not found: ", tpm_col)
    return(NULL)
  }

  # Match by transcript name
  m <- match(quant_df$target_name, gt$transcript_id)
  if (any(is.na(m))) {
    warning(sum(is.na(m)), " transcripts not found in ground truth")
  }

  valid <- !is.na(m)
  est_tpm <- quant_df$tpm[valid]
  true_tpm <- gt[[tpm_col]][m[valid]]

  # Add pseudocount for log correlation
  pseudo <- 0.01
  log_est <- log2(est_tpm + pseudo)
  log_true <- log2(true_tpm + pseudo)

  list(
    pearson = cor(log_est, log_true, method = "pearson"),
    spearman = cor(est_tpm, true_tpm, method = "spearman"),
    mard = median(abs(est_tpm - true_tpm) / (true_tpm + pseudo)),
    n = sum(valid)
  )
}

# ---- Evaluate single-sample ----

cat("=== Single-sample EM results ===\n")
single_dir <- file.path(simdir, "quant_single")
single_results <- list()

for (i in 1:nrow(sample_info)) {
  sname <- sample_info$sample_name[i]
  condition <- sample_info$condition[i]
  quant_path <- file.path(single_dir, sname, paste0(sname, ".quant"))

  quant_df <- read_quant(quant_path)
  if (is.null(quant_df)) next

  metrics <- compare_to_gt(quant_df, gt, condition)
  single_results[[sname]] <- c(condition = condition, metrics)

  cat(sprintf("  %s (%s): Pearson(log2)=%.4f, Spearman=%.4f, MARD=%.4f\n",
              sname, condition, metrics$pearson, metrics$spearman, metrics$mard))
}

# ---- Evaluate multi-sample ----

cat("\n=== Multi-sample hierarchical results ===\n")
multi_dir <- file.path(simdir, "quant_multi")
multi_results <- list()

for (i in 1:nrow(sample_info)) {
  sname <- sample_info$sample_name[i]
  condition <- sample_info$condition[i]
  quant_path <- file.path(multi_dir, sname, paste0(sname, ".quant"))

  quant_df <- read_quant(quant_path)
  if (is.null(quant_df)) next

  metrics <- compare_to_gt(quant_df, gt, condition)
  multi_results[[sname]] <- c(condition = condition, metrics)

  cat(sprintf("  %s (%s): Pearson(log2)=%.4f, Spearman=%.4f, MARD=%.4f\n",
              sname, condition, metrics$pearson, metrics$spearman, metrics$mard))
}

# ---- Summary comparison ----

if (length(single_results) > 0 && length(multi_results) > 0) {
  cat("\n=== Summary ===\n")

  single_pearson <- mean(sapply(single_results, function(x) x$pearson))
  single_spearman <- mean(sapply(single_results, function(x) x$spearman))
  single_mard <- mean(sapply(single_results, function(x) x$mard))

  multi_pearson <- mean(sapply(multi_results, function(x) x$pearson))
  multi_spearman <- mean(sapply(multi_results, function(x) x$spearman))
  multi_mard <- mean(sapply(multi_results, function(x) x$mard))

  cat(sprintf("  Single-sample:  Pearson=%.4f  Spearman=%.4f  MARD=%.4f\n",
              single_pearson, single_spearman, single_mard))
  cat(sprintf("  Multi-sample:   Pearson=%.4f  Spearman=%.4f  MARD=%.4f\n",
              multi_pearson, multi_spearman, multi_mard))

  cat(sprintf("\n  Delta Pearson:  %+.4f (%s)\n",
              multi_pearson - single_pearson,
              ifelse(multi_pearson > single_pearson, "improved", "worse")))
  cat(sprintf("  Delta Spearman: %+.4f (%s)\n",
              multi_spearman - single_spearman,
              ifelse(multi_spearman > single_spearman, "improved", "worse")))
  cat(sprintf("  Delta MARD:     %+.4f (%s)\n",
              multi_mard - single_mard,
              ifelse(multi_mard < single_mard, "improved", "worse")))

  # ---- DE-specific analysis ----
  cat("\n=== DE transcript analysis ===\n")

  de_tx <- gt$transcript_id[gt$is_de]
  non_de_tx <- gt$transcript_id[!gt$is_de]

  # Compare MARD for DE vs non-DE transcripts
  for (method_name in c("single", "multi")) {
    results <- if (method_name == "single") single_results else multi_results
    de_mards <- c()
    non_de_mards <- c()

    for (sname in names(results)) {
      condition <- results[[sname]]$condition
      tpm_col <- paste0("expected_tpm_", condition)
      quant_dir <- if (method_name == "single") single_dir else multi_dir
      quant_path <- file.path(quant_dir, sname, paste0(sname, ".quant"))
      quant_df <- read_quant(quant_path)
      if (is.null(quant_df)) next

      m <- match(quant_df$target_name, gt$transcript_id)
      valid <- !is.na(m)
      est_tpm <- quant_df$tpm[valid]
      true_tpm <- gt[[tpm_col]][m[valid]]
      tx_names <- gt$transcript_id[m[valid]]

      pseudo <- 0.01
      ard <- abs(est_tpm - true_tpm) / (true_tpm + pseudo)

      de_mask <- tx_names %in% de_tx
      de_mards <- c(de_mards, median(ard[de_mask]))
      non_de_mards <- c(non_de_mards, median(ard[!de_mask]))
    }

    cat(sprintf("  %s-sample:  DE MARD=%.4f,  non-DE MARD=%.4f\n",
                method_name, mean(de_mards), mean(non_de_mards)))
  }
}

cat("\n=== Done ===\n")
