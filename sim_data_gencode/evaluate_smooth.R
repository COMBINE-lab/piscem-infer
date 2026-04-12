#!/usr/bin/env Rscript
#
# Evaluate coverage smoothness EM vs plain EM vs selection variants.
# Single-sample evaluation on sample_01.

simdir <- "sim_data_gencode"

cat("=== Coverage Smoothness Evaluation ===\n\n")

# ---- Load ground truth ----
gt <- read.csv(file.path(simdir, "ground_truth.csv"), stringsAsFactors = FALSE)
gt$short_id <- sub(" .*", "", gt$transcript_id)
if (!"is_expressed" %in% colnames(gt)) {
  gt$is_expressed <- gt$base_reads > 0
}
num_expressed <- sum(gt$is_expressed)
num_total <- nrow(gt)
cat("Ground truth:", num_total, "transcripts (", num_expressed, "expressed,",
    num_total - num_expressed, "unexpressed)\n\n")

# ---- Helper ----
read_quant <- function(path) {
  if (!file.exists(path)) { warning("Not found: ", path); return(NULL) }
  df <- read.delim(path, stringsAsFactors = FALSE)
  colnames(df) <- c("target_name", "len", "eelen", "tpm", "ecount")
  df
}

compare_to_gt <- function(quant_df, gt, condition = "control") {
  tpm_col <- paste0("expected_tpm_", condition)
  m <- match(quant_df$target_name, gt$short_id)
  valid <- !is.na(m)
  est_tpm <- quant_df$tpm[valid]
  true_tpm <- gt[[tpm_col]][m[valid]]
  is_expr <- gt$is_expressed[m[valid]]

  pseudo <- 0.01

  # All transcripts
  log_est <- log2(est_tpm + pseudo)
  log_true <- log2(true_tpm + pseudo)
  all_m <- list(
    pearson = cor(log_est, log_true, method = "pearson"),
    spearman = cor(est_tpm, true_tpm, method = "spearman"),
    mard = median(abs(est_tpm - true_tpm) / (true_tpm + pseudo))
  )

  # Expressed only
  e_est <- est_tpm[is_expr]
  e_true <- true_tpm[is_expr]
  expr_m <- list(
    pearson = cor(log2(e_est + pseudo), log2(e_true + pseudo), method = "pearson"),
    spearman = cor(e_est, e_true, method = "spearman"),
    mard = median(abs(e_est - e_true) / (e_true + pseudo))
  )

  # False positives (unexpressed with TPM > 1)
  u_est <- est_tpm[!is_expr]
  fp_thresh <- 1.0
  n_fp <- sum(u_est > fp_thresh)
  fp_rate <- n_fp / sum(!is_expr)

  # False negatives (expressed with TPM < 0.01)
  fn_thresh <- 0.01
  n_fn <- sum(e_est < fn_thresh)
  fn_rate <- n_fn / sum(is_expr)

  list(all = all_m, expressed = expr_m,
       fp = n_fp, fp_rate = fp_rate,
       fn = n_fn, fn_rate = fn_rate)
}

# ---- Evaluate all methods ----
methods <- list(
  "Plain EM"         = file.path(simdir, "quant_em",         "sample_01", "sample_01.quant"),
  "EM + Selection"   = file.path(simdir, "quant_sel",        "sample_01", "sample_01.quant"),
  "EM + Smooth"      = file.path(simdir, "quant_smooth",     "sample_01", "sample_01.quant"),
  "EM + Smooth + Sel"= file.path(simdir, "quant_smooth_sel", "sample_01", "sample_01.quant")
)

cat(sprintf("%-20s | %8s %8s %8s | %8s %8s %8s | %6s %6s | %6s %6s\n",
            "Method", "P(all)", "S(all)", "M(all)",
            "P(expr)", "S(expr)", "M(expr)",
            "FP", "FP%", "FN", "FN%"))
cat(paste0(rep("-", 110), collapse = ""), "\n")

for (name in names(methods)) {
  qf <- read_quant(methods[[name]])
  if (is.null(qf)) { cat(name, ": MISSING\n"); next }

  r <- compare_to_gt(qf, gt, "control")

  cat(sprintf("%-20s | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f | %6d %5.1f%% | %6d %5.1f%%\n",
              name,
              r$all$pearson, r$all$spearman, r$all$mard,
              r$expressed$pearson, r$expressed$spearman, r$expressed$mard,
              r$fp, r$fp_rate * 100,
              r$fn, r$fn_rate * 100))
}

cat("\n=== Done ===\n")
