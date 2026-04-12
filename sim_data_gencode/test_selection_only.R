#!/usr/bin/env Rscript
# Test: does the hierarchical model's benefit come from selection or shrinkage?
#
# Approach: take single-sample EM estimates, zero out transcripts NOT in the
# NoCond+AV expressed set, renormalize to TPM. Compare this "EM + selection"
# against plain EM, NoCond, and NoCond+AV.

library(data.table)

simdir <- "sim_data_gencode"
pseudo <- 0.01

# ---- Load ground truth ----
gt <- fread(file.path(simdir, "ground_truth.csv"))
gt[, short_id := sub(" .*", "", transcript_id)]

sample_info <- fread(file.path(simdir, "sample_info.csv"))

read_quant <- function(path) {
  if (!file.exists(path)) { warning("Not found: ", path); return(NULL) }
  df <- fread(path, sep = "\t")
  setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
  df
}

# ---- Get expressed set from NoCond+AV (union across all samples) ----
av_expressed <- character(0)
for (i in 1:nrow(sample_info)) {
  sn <- sample_info$sample_name[i]
  qpath <- file.path(simdir, "quant_hier_nocond_av", sn, paste0(sn, ".quant"))
  qdf <- read_quant(qpath)
  if (!is.null(qdf)) {
    av_expressed <- union(av_expressed, qdf$target_name[qdf$tpm > 0])
  }
}
cat(sprintf("NoCond+AV expressed set: %d transcripts (out of %d)\n",
            length(av_expressed), nrow(gt)))

# ---- Build "EM + selection" estimates ----
# For each sample: take EM TPM, zero out non-expressed, renormalize
em_sel_quants <- list()
for (i in 1:nrow(sample_info)) {
  sn <- sample_info$sample_name[i]
  qpath <- file.path(simdir, "quant_em", sn, paste0(sn, ".quant"))
  qdf <- read_quant(qpath)
  if (is.null(qdf)) next

  # Zero out transcripts not in expressed set
  keep <- qdf$target_name %in% av_expressed
  tpm_filtered <- qdf$tpm * keep

  # Renormalize to sum to 1e6
  total <- sum(tpm_filtered)
  if (total > 0) {
    tpm_renorm <- tpm_filtered * 1e6 / total
  } else {
    tpm_renorm <- tpm_filtered
  }

  em_sel_quants[[sn]] <- data.table(
    target_name = qdf$target_name,
    tpm = tpm_renorm,
    eelen = qdf$eelen
  )
}

# ---- Evaluate all methods ----
methods <- list(
  list(name = "Single EM",     dir = "quant_em",             subdir = TRUE),
  list(name = "EM+Selection",  custom = TRUE),
  list(name = "Hier NoCond",   dir = "quant_hier_nocond",    subdir = TRUE),
  list(name = "NoCond+AV",     dir = "quant_hier_nocond_av", subdir = TRUE)
)

ctrl_samples  <- sample_info$sample_name[sample_info$condition == "control"]
treat_samples <- sample_info$sample_name[sample_info$condition == "treatment"]

cat("\n=== Per-sample accuracy ===\n\n")
cat(sprintf("%-14s | %8s %8s %8s | %8s\n",
            "Method", "Pearson", "Spearman", "RMSE", "N>0"))
cat(paste0(rep("-", 60), collapse = ""), "\n")

for (meth in methods) {
  all_pearson <- c()
  all_spearman <- c()
  all_rmse <- c()
  all_npos <- c()

  for (i in 1:nrow(sample_info)) {
    sn <- sample_info$sample_name[i]
    cond <- sample_info$condition[i]
    tpm_col <- paste0("expected_tpm_", cond)

    if (!is.null(meth$custom) && meth$custom) {
      qdf <- em_sel_quants[[sn]]
      if (is.null(qdf)) next
      est_tpm <- qdf$tpm
      tx_names <- qdf$target_name
    } else {
      qpath <- file.path(simdir, meth$dir, sn, paste0(sn, ".quant"))
      qdf <- read_quant(qpath)
      if (is.null(qdf)) next
      est_tpm <- qdf$tpm
      tx_names <- qdf$target_name
    }

    m <- match(tx_names, gt$short_id)
    valid <- !is.na(m)
    est <- est_tpm[valid]
    true_t <- gt[[tpm_col]][m[valid]]

    p <- cor(log2(est + pseudo), log2(true_t + pseudo), method = "pearson")
    s <- cor(est, true_t, method = "spearman")
    rmse <- sqrt(mean((log2(est + 1) - log2(true_t + 1))^2))

    all_pearson <- c(all_pearson, p)
    all_spearman <- c(all_spearman, s)
    all_rmse <- c(all_rmse, rmse)
    all_npos <- c(all_npos, sum(est > 0))
  }

  cat(sprintf("%-14s | %8.4f %8.4f %8.4f | %8.0f\n",
              meth$name, mean(all_pearson), mean(all_spearman),
              mean(all_rmse), mean(all_npos)))
}

# ---- Fold change preservation ----
cat("\n=== Fold change preservation (DE transcripts, |FC| > 1) ===\n\n")
cat(sprintf("%-14s | %8s %8s %8s\n", "Method", "Slope", "R²", "N"))
cat(paste0(rep("-", 45), collapse = ""), "\n")

for (meth in methods) {
  ctrl_tpms <- list()
  treat_tpms <- list()
  ref_names <- NULL

  for (sn in ctrl_samples) {
    if (!is.null(meth$custom) && meth$custom) {
      qdf <- em_sel_quants[[sn]]
    } else {
      qpath <- file.path(simdir, meth$dir, sn, paste0(sn, ".quant"))
      qdf <- read_quant(qpath)
    }
    if (!is.null(qdf)) {
      ctrl_tpms[[sn]] <- qdf$tpm
      if (is.null(ref_names)) ref_names <- qdf$target_name
    }
  }

  for (sn in treat_samples) {
    if (!is.null(meth$custom) && meth$custom) {
      qdf <- em_sel_quants[[sn]]
    } else {
      qpath <- file.path(simdir, meth$dir, sn, paste0(sn, ".quant"))
      qdf <- read_quant(qpath)
    }
    if (!is.null(qdf)) treat_tpms[[sn]] <- qdf$tpm
  }

  if (length(ctrl_tpms) == 0 || length(treat_tpms) == 0) next

  mean_ctrl  <- Reduce("+", ctrl_tpms)  / length(ctrl_tpms)
  mean_treat <- Reduce("+", treat_tpms) / length(treat_tpms)
  est_lfc    <- log2(mean_treat + pseudo) - log2(mean_ctrl + pseudo)

  m <- match(ref_names, gt$short_id)
  valid <- !is.na(m)
  true_ctrl <- gt$expected_tpm_control[m[valid]]
  true_treat <- gt$expected_tpm_treatment[m[valid]]
  true_lfc <- log2(true_treat + pseudo) - log2(true_ctrl + pseudo)

  est_fc <- est_lfc[valid]
  de_mask <- ref_names[valid] %in% gt$short_id[gt$is_de]
  big_fc <- abs(true_lfc) > 1

  if (sum(de_mask & big_fc) > 10) {
    fit <- lm(est_fc[de_mask & big_fc] ~ true_lfc[de_mask & big_fc])
    slope <- coef(fit)[2]
    r2 <- summary(fit)$r.squared
    cat(sprintf("%-14s | %8.3f %8.3f %8d\n",
                meth$name, slope, r2, sum(de_mask & big_fc)))
  }
}

cat("\nDone.\n")
