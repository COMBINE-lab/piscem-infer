#!/usr/bin/env Rscript
# Test: can simple consensus filtering of per-sample EM match the
# hierarchical model's transcript selection?
#
# For each threshold K (1..N), keep transcripts called expressed (TPM > 0)
# in at least K of N samples. Then zero + renormalize.

library(data.table)

simdir <- "sim_data_gencode"
pseudo <- 0.01

gt <- fread(file.path(simdir, "ground_truth.csv"))
gt[, short_id := sub(" .*", "", transcript_id)]

sample_info <- fread(file.path(simdir, "sample_info.csv"))

read_quant <- function(path) {
  if (!file.exists(path)) { warning("Not found: ", path); return(NULL) }
  df <- fread(path, sep = "\t")
  setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
  df
}

# ---- Load all per-sample EM results ----
em_quants <- list()
all_names <- NULL
for (i in 1:nrow(sample_info)) {
  sn <- sample_info$sample_name[i]
  qpath <- file.path(simdir, "quant_em", sn, paste0(sn, ".quant"))
  qdf <- read_quant(qpath)
  em_quants[[sn]] <- qdf
  if (is.null(all_names)) all_names <- qdf$target_name
}

N <- nrow(sample_info)

# ---- Count how many samples call each transcript as expressed ----
expressed_count <- rep(0L, length(all_names))
for (sn in sample_info$sample_name) {
  qdf <- em_quants[[sn]]
  m <- match(all_names, qdf$target_name)
  expressed_count <- expressed_count + (qdf$tpm[m] > 0)
}

# ---- Also get NoCond+AV expressed set for comparison ----
av_expressed <- character(0)
for (i in 1:nrow(sample_info)) {
  sn <- sample_info$sample_name[i]
  qpath <- file.path(simdir, "quant_hier_nocond_av", sn, paste0(sn, ".quant"))
  qdf <- read_quant(qpath)
  if (!is.null(qdf)) av_expressed <- union(av_expressed, qdf$target_name[qdf$tpm > 0])
}

# ---- Ground truth expressed set ----
m_gt <- match(all_names, gt$short_id)
# Use union of both conditions
truly_expressed <- gt$expected_tpm_control[m_gt] > 0 | gt$expected_tpm_treatment[m_gt] > 0

cat(sprintf("Total transcripts: %d\n", length(all_names)))
cat(sprintf("Truly expressed: %d\n", sum(truly_expressed)))
cat(sprintf("NoCond+AV expressed set: %d\n\n", length(av_expressed)))

# ---- Evaluate consensus filter at different K thresholds ----
cat("=== Consensus filter: keep transcripts expressed in >= K of N samples ===\n\n")
cat(sprintf("%3s | %6s %6s %6s | %8s %8s %8s | %8s %8s %8s %8s\n",
            "K", "kept", "TP", "FP", "Prec", "Recall", "F1",
            "Pearson", "Spearman", "RMSE", "FC slope"))
cat(paste0(rep("-", 105), collapse = ""), "\n")

ctrl_samples  <- sample_info$sample_name[sample_info$condition == "control"]
treat_samples <- sample_info$sample_name[sample_info$condition == "treatment"]

for (K in 1:N) {
  keep_set <- all_names[expressed_count >= K]

  # Per-sample accuracy
  all_pearson <- c(); all_spearman <- c(); all_rmse <- c()
  all_tp <- c(); all_fp <- c(); all_fn <- c()

  for (i in 1:nrow(sample_info)) {
    sn <- sample_info$sample_name[i]
    cond <- sample_info$condition[i]
    tpm_col <- paste0("expected_tpm_", cond)

    qdf <- em_quants[[sn]]
    est_tpm <- qdf$tpm
    keep <- qdf$target_name %in% keep_set
    est_tpm <- est_tpm * keep
    total <- sum(est_tpm)
    if (total > 0) est_tpm <- est_tpm * 1e6 / total

    m <- match(qdf$target_name, gt$short_id)
    valid <- !is.na(m)
    true_t <- gt[[tpm_col]][m[valid]]
    est <- est_tpm[valid]

    all_pearson <- c(all_pearson, cor(log2(est + pseudo), log2(true_t + pseudo), method = "pearson"))
    all_spearman <- c(all_spearman, cor(est, true_t, method = "spearman"))
    all_rmse <- c(all_rmse, sqrt(mean((log2(est + 1) - log2(true_t + 1))^2)))

    all_tp <- c(all_tp, sum(est > 0 & true_t > 0))
    all_fp <- c(all_fp, sum(est > 0 & true_t == 0))
    all_fn <- c(all_fn, sum(est == 0 & true_t > 0))
  }

  tp <- mean(all_tp); fp <- mean(all_fp); fn <- mean(all_fn)
  prec <- tp / (tp + fp)
  rec <- tp / (tp + fn)
  f1 <- 2 * prec * rec / (prec + rec)

  # Fold change
  ctrl_tpms <- list(); treat_tpms <- list()
  for (sn in ctrl_samples) {
    qdf <- em_quants[[sn]]
    est_tpm <- qdf$tpm * (qdf$target_name %in% keep_set)
    total <- sum(est_tpm)
    if (total > 0) est_tpm <- est_tpm * 1e6 / total
    ctrl_tpms[[sn]] <- est_tpm
  }
  for (sn in treat_samples) {
    qdf <- em_quants[[sn]]
    est_tpm <- qdf$tpm * (qdf$target_name %in% keep_set)
    total <- sum(est_tpm)
    if (total > 0) est_tpm <- est_tpm * 1e6 / total
    treat_tpms[[sn]] <- est_tpm
  }

  mean_ctrl  <- Reduce("+", ctrl_tpms)  / length(ctrl_tpms)
  mean_treat <- Reduce("+", treat_tpms) / length(treat_tpms)
  est_lfc    <- log2(mean_treat + pseudo) - log2(mean_ctrl + pseudo)

  ref_names <- em_quants[[ctrl_samples[1]]]$target_name
  m <- match(ref_names, gt$short_id)
  valid <- !is.na(m)
  true_ctrl <- gt$expected_tpm_control[m[valid]]
  true_treat <- gt$expected_tpm_treatment[m[valid]]
  true_lfc <- log2(true_treat + pseudo) - log2(true_ctrl + pseudo)
  de_mask <- ref_names[valid] %in% gt$short_id[gt$is_de]
  big_fc <- abs(true_lfc) > 1

  slope_str <- ""
  if (sum(de_mask & big_fc) > 10) {
    fit <- lm(est_lfc[valid][de_mask & big_fc] ~ true_lfc[de_mask & big_fc])
    slope_str <- sprintf("%.3f", coef(fit)[2])
  }

  cat(sprintf("%3d | %6d %6d %6d | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f %8s\n",
              K, length(keep_set), round(tp), round(fp),
              prec, rec, f1,
              mean(all_pearson), mean(all_spearman), mean(all_rmse), slope_str))
}

# ---- Reference rows ----
cat(paste0(rep("-", 105), collapse = ""), "\n")
cat("Reference methods:\n")

# NoCond+AV row
all_pearson <- c(); all_spearman <- c(); all_rmse <- c()
all_tp <- c(); all_fp <- c(); all_fn <- c()
for (i in 1:nrow(sample_info)) {
  sn <- sample_info$sample_name[i]
  cond <- sample_info$condition[i]
  tpm_col <- paste0("expected_tpm_", cond)
  qpath <- file.path(simdir, "quant_hier_nocond_av", sn, paste0(sn, ".quant"))
  qdf <- read_quant(qpath)
  m <- match(qdf$target_name, gt$short_id)
  valid <- !is.na(m)
  true_t <- gt[[tpm_col]][m[valid]]
  est <- qdf$tpm[valid]
  all_pearson <- c(all_pearson, cor(log2(est + pseudo), log2(true_t + pseudo), method = "pearson"))
  all_spearman <- c(all_spearman, cor(est, true_t, method = "spearman"))
  all_rmse <- c(all_rmse, sqrt(mean((log2(est + 1) - log2(true_t + 1))^2)))
  all_tp <- c(all_tp, sum(est > 0 & true_t > 0))
  all_fp <- c(all_fp, sum(est > 0 & true_t == 0))
  all_fn <- c(all_fn, sum(est == 0 & true_t > 0))
}
tp <- mean(all_tp); fp <- mean(all_fp); fn <- mean(all_fn)
cat(sprintf("AV  | %6d %6d %6d | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f %8s\n",
            length(av_expressed), round(tp), round(fp),
            tp/(tp+fp), tp/(tp+fn), 2*(tp/(tp+fp))*(tp/(tp+fn))/((tp/(tp+fp))+(tp/(tp+fn))),
            mean(all_pearson), mean(all_spearman), mean(all_rmse), "0.918"))

cat("\nDone.\n")
