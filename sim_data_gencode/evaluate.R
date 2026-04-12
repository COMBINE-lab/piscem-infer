#!/usr/bin/env Rscript
#
# Evaluate piscem-infer quantification accuracy against ground truth.
#
# Compares single-sample EM vs multi-sample hierarchical estimates.
# Handles full transcriptomes where many transcripts have 0 true abundance.
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
# Extract short ID (first field before space) to match piscem's truncated names
gt$short_id <- sub(" .*", "", gt$transcript_id)

# Derive is_expressed from base_reads (may not exist as a column)
if (!"is_expressed" %in% colnames(gt)) {
  gt$is_expressed <- gt$base_reads > 0
}

num_expressed <- sum(gt$is_expressed)
num_total <- nrow(gt)
cat("Loaded ground truth for", num_total, "transcripts")
cat(" (", num_expressed, "expressed,", num_total - num_expressed, "unexpressed )\n")

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

# Compare estimated TPM against ground truth.
# Returns metrics for: all transcripts, expressed only, unexpressed only.
compare_to_gt <- function(quant_df, gt, condition) {
  tpm_col <- paste0("expected_tpm_", condition)
  if (!tpm_col %in% colnames(gt)) {
    warning("Column not found: ", tpm_col)
    return(NULL)
  }

  # Match by short ID (piscem truncates names at first space)
  m <- match(quant_df$target_name, gt$short_id)
  if (any(is.na(m))) {
    warning(sum(is.na(m)), " transcripts not found in ground truth")
  }

  valid <- !is.na(m)
  est_tpm <- quant_df$tpm[valid]
  true_tpm <- gt[[tpm_col]][m[valid]]
  is_expr <- gt$is_expressed[m[valid]]

  pseudo <- 0.01

  # --- All transcripts ---
  log_est <- log2(est_tpm + pseudo)
  log_true <- log2(true_tpm + pseudo)

  all_metrics <- list(
    pearson = cor(log_est, log_true, method = "pearson"),
    spearman = cor(est_tpm, true_tpm, method = "spearman"),
    mard = median(abs(est_tpm - true_tpm) / (true_tpm + pseudo)),
    n = sum(valid)
  )

  # --- Expressed only ---
  expr_metrics <- NULL
  if (sum(is_expr) > 0) {
    e_est <- est_tpm[is_expr]
    e_true <- true_tpm[is_expr]
    expr_metrics <- list(
      pearson = cor(log2(e_est + pseudo), log2(e_true + pseudo), method = "pearson"),
      spearman = cor(e_est, e_true, method = "spearman"),
      mard = median(abs(e_est - e_true) / (e_true + pseudo)),
      n = sum(is_expr)
    )
  }

  # --- Unexpressed: false positive rate ---
  unexpr_metrics <- NULL
  n_unexpr <- sum(!is_expr)
  if (n_unexpr > 0) {
    u_est <- est_tpm[!is_expr]
    # For unexpressed transcripts, any estimated TPM > threshold is a false positive
    fp_thresh <- 1.0  # TPM threshold for "detected"
    n_fp <- sum(u_est > fp_thresh)
    fp_rate <- n_fp / n_unexpr
    unexpr_metrics <- list(
      n = n_unexpr,
      n_fp = n_fp,
      fp_rate = fp_rate,
      median_tpm = median(u_est),
      max_tpm = max(u_est)
    )
  }

  list(all = all_metrics, expressed = expr_metrics, unexpressed = unexpr_metrics)
}

# ---- Evaluate both methods ----

for (method_name in c("single", "multi")) {
  cat(sprintf("=== %s-sample results ===\n", method_name))
  quant_dir <- file.path(simdir, paste0("quant_", method_name))

  results <- list()

  for (i in 1:nrow(sample_info)) {
    sname <- sample_info$sample_name[i]
    condition <- sample_info$condition[i]
    quant_path <- file.path(quant_dir, sname, paste0(sname, ".quant"))

    quant_df <- read_quant(quant_path)
    if (is.null(quant_df)) next

    metrics <- compare_to_gt(quant_df, gt, condition)
    results[[sname]] <- c(condition = condition, metrics)

    a <- metrics$all
    cat(sprintf("  %s (%s): Pearson(log2)=%.4f, Spearman=%.4f, MARD=%.4f",
                sname, condition, a$pearson, a$spearman, a$mard))
    if (!is.null(metrics$expressed)) {
      cat(sprintf(" | Expr: P=%.4f, S=%.4f, M=%.4f",
                  metrics$expressed$pearson, metrics$expressed$spearman, metrics$expressed$mard))
    }
    if (!is.null(metrics$unexpressed)) {
      cat(sprintf(" | FP: %d/%d (%.1f%%)",
                  metrics$unexpressed$n_fp, metrics$unexpressed$n, metrics$unexpressed$fp_rate * 100))
    }
    cat("\n")
  }

  cat("\n")

  # Store for summary
  assign(paste0(method_name, "_results"), results)
}

# ---- Summary comparison ----

if (length(single_results) > 0 && length(multi_results) > 0) {
  cat("=== Summary ===\n")

  for (scope in c("all", "expressed")) {
    scope_label <- if (scope == "all") "All transcripts" else "Expressed only"

    single_pearson <- mean(sapply(single_results, function(x) x[[scope]]$pearson), na.rm = TRUE)
    single_spearman <- mean(sapply(single_results, function(x) x[[scope]]$spearman), na.rm = TRUE)
    single_mard <- mean(sapply(single_results, function(x) x[[scope]]$mard), na.rm = TRUE)

    multi_pearson <- mean(sapply(multi_results, function(x) x[[scope]]$pearson), na.rm = TRUE)
    multi_spearman <- mean(sapply(multi_results, function(x) x[[scope]]$spearman), na.rm = TRUE)
    multi_mard <- mean(sapply(multi_results, function(x) x[[scope]]$mard), na.rm = TRUE)

    cat(sprintf("\n  --- %s ---\n", scope_label))
    cat(sprintf("  Single-sample:  Pearson=%.4f  Spearman=%.4f  MARD=%.4f\n",
                single_pearson, single_spearman, single_mard))
    cat(sprintf("  Multi-sample:   Pearson=%.4f  Spearman=%.4f  MARD=%.4f\n",
                multi_pearson, multi_spearman, multi_mard))

    cat(sprintf("  Delta Pearson:  %+.4f (%s)\n",
                multi_pearson - single_pearson,
                ifelse(multi_pearson > single_pearson, "improved", "worse")))
    cat(sprintf("  Delta Spearman: %+.4f (%s)\n",
                multi_spearman - single_spearman,
                ifelse(multi_spearman > single_spearman, "improved", "worse")))
    cat(sprintf("  Delta MARD:     %+.4f (%s)\n",
                multi_mard - single_mard,
                ifelse(multi_mard < single_mard, "improved", "worse")))
  }

  # Unexpressed false positives (only if there are unexpressed transcripts)
  has_unexpr <- any(sapply(single_results, function(x) !is.null(x$unexpressed)))
  if (has_unexpr) {
    cat("\n  --- Unexpressed (false positives, TPM > 1) ---\n")
    single_fp <- mean(sapply(single_results, function(x) if (!is.null(x$unexpressed)) x$unexpressed$fp_rate else NA), na.rm = TRUE)
    multi_fp <- mean(sapply(multi_results, function(x) if (!is.null(x$unexpressed)) x$unexpressed$fp_rate else NA), na.rm = TRUE)
    single_max <- max(sapply(single_results, function(x) if (!is.null(x$unexpressed)) x$unexpressed$max_tpm else NA), na.rm = TRUE)
    multi_max <- max(sapply(multi_results, function(x) if (!is.null(x$unexpressed)) x$unexpressed$max_tpm else NA), na.rm = TRUE)

    cat(sprintf("  Single-sample:  FP rate=%.2f%%  max spurious TPM=%.1f\n",
                single_fp * 100, single_max))
    cat(sprintf("  Multi-sample:   FP rate=%.2f%%  max spurious TPM=%.1f\n",
                multi_fp * 100, multi_max))
  } else {
    cat("\n  (No unexpressed transcripts in ground truth — FP analysis skipped)\n")
  }

  # ---- DE-specific analysis ----
  cat("\n  --- DE transcript analysis ---\n")

  de_tx <- gt$short_id[gt$is_de]
  non_de_tx <- gt$short_id[!gt$is_de & gt$is_expressed]

  for (method_name in c("single", "multi")) {
    results <- get(paste0(method_name, "_results"))
    de_mards <- c()
    non_de_mards <- c()

    for (sname in names(results)) {
      condition <- results[[sname]]$condition
      tpm_col <- paste0("expected_tpm_", condition)
      quant_dir <- if (method_name == "single") {
        file.path(simdir, "quant_single")
      } else {
        file.path(simdir, "quant_multi")
      }
      quant_path <- file.path(quant_dir, sname, paste0(sname, ".quant"))
      quant_df <- read_quant(quant_path)
      if (is.null(quant_df)) next

      m <- match(quant_df$target_name, gt$short_id)
      valid <- !is.na(m)
      est_tpm <- quant_df$tpm[valid]
      true_tpm <- gt[[tpm_col]][m[valid]]
      tx_ids <- gt$short_id[m[valid]]

      pseudo <- 0.01
      ard <- abs(est_tpm - true_tpm) / (true_tpm + pseudo)

      de_mask <- tx_ids %in% de_tx
      non_de_mask <- tx_ids %in% non_de_tx

      if (sum(de_mask) > 0) de_mards <- c(de_mards, median(ard[de_mask]))
      if (sum(non_de_mask) > 0) non_de_mards <- c(non_de_mards, median(ard[non_de_mask]))
    }

    cat(sprintf("  %s-sample:  DE MARD=%.4f,  non-DE expressed MARD=%.4f\n",
                method_name, mean(de_mards), mean(non_de_mards)))
  }
}

cat("\n=== Done ===\n")
