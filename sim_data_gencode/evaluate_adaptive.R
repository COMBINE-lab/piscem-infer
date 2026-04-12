#!/usr/bin/env Rscript
#
# Compare: single-sample EM, multi-sample nocond, multi-sample nocond+adaptive-variance
# against ground truth for the gencode simulation data.

simdir <- "sim_data_gencode"

# ---- Load ground truth ----
gt <- read.csv(file.path(simdir, "ground_truth.csv"), stringsAsFactors = FALSE)
gt$short_id <- sub(" .*", "", gt$transcript_id)
if (!"is_expressed" %in% colnames(gt)) gt$is_expressed <- gt$base_reads > 0

sample_info <- read.csv(file.path(simdir, "sample_info.csv"), stringsAsFactors = FALSE)

read_quant <- function(path) {
  if (!file.exists(path)) { warning("Not found: ", path); return(NULL) }
  df <- read.delim(path, stringsAsFactors = FALSE)
  colnames(df) <- c("target_name", "len", "eelen", "tpm", "ecount")
  df
}

pseudo <- 0.01

cat("=== Gencode Simulation: EM vs NoCond vs NoCond+AdaptiveVar ===\n\n")

methods <- list(
  list(name = "Single EM",   dir = "quant_em",             subdir = TRUE),
  list(name = "Hier NoCond", dir = "quant_hier_nocond",    subdir = TRUE),
  list(name = "NoCond+AV",   dir = "quant_hier_nocond_av", subdir = TRUE)
)

# Also check if condition-aware multi exists
if (dir.exists(file.path(simdir, "quant_multi"))) {
  methods <- c(list(list(name = "Hier Cond",  dir = "quant_multi", subdir = TRUE)), methods)
}

for (meth in methods) {
  cat(sprintf("=== %s ===\n", meth$name))

  all_pearson <- c()
  all_spearman <- c()
  all_mard <- c()
  de_mards <- c()
  non_de_mards <- c()

  for (i in 1:nrow(sample_info)) {
    sname <- sample_info$sample_name[i]
    condition <- sample_info$condition[i]
    tpm_col <- paste0("expected_tpm_", condition)

    qpath <- file.path(simdir, meth$dir, sname, paste0(sname, ".quant"))
    qdf <- read_quant(qpath)
    if (is.null(qdf)) next

    m <- match(qdf$target_name, gt$short_id)
    valid <- !is.na(m)
    est <- qdf$tpm[valid]
    true_t <- gt[[tpm_col]][m[valid]]

    log_est <- log2(est + pseudo)
    log_true <- log2(true_t + pseudo)

    p <- cor(log_est, log_true, method = "pearson")
    s <- cor(est, true_t, method = "spearman")
    md <- median(abs(est - true_t) / (true_t + pseudo))

    all_pearson <- c(all_pearson, p)
    all_spearman <- c(all_spearman, s)
    all_mard <- c(all_mard, md)

    # DE-specific
    tx_ids <- gt$short_id[m[valid]]
    ard <- abs(est - true_t) / (true_t + pseudo)

    de_tx <- gt$short_id[gt$is_de]
    non_de_tx <- gt$short_id[!gt$is_de & gt$is_expressed]

    de_mask <- tx_ids %in% de_tx
    non_de_mask <- tx_ids %in% non_de_tx

    if (sum(de_mask) > 0) de_mards <- c(de_mards, median(ard[de_mask]))
    if (sum(non_de_mask) > 0) non_de_mards <- c(non_de_mards, median(ard[non_de_mask]))

    cat(sprintf("  %s (%s): Pearson=%.4f, Spearman=%.4f, MARD=%.4f\n",
                sname, condition, p, s, md))
  }

  cat(sprintf("  MEAN:          Pearson=%.4f, Spearman=%.4f, MARD=%.4f\n",
              mean(all_pearson), mean(all_spearman), mean(all_mard)))
  if (length(de_mards) > 0) {
    cat(sprintf("  DE MARD=%.4f, non-DE MARD=%.4f\n",
                mean(de_mards), mean(non_de_mards)))
  }
  cat("\n")
}

# ---- Fold-change analysis ----
cat("=== Fold Change Preservation ===\n\n")

# Compare control vs treatment fold changes for each method
ctrl_samples <- sample_info$sample_name[sample_info$condition == "control"]
treat_samples <- sample_info$sample_name[sample_info$condition == "treatment"]

for (meth in methods) {
  # Average TPM across control replicates
  ctrl_tpms <- list()
  for (sname in ctrl_samples) {
    qpath <- file.path(simdir, meth$dir, sname, paste0(sname, ".quant"))
    qdf <- read_quant(qpath)
    if (!is.null(qdf)) ctrl_tpms[[sname]] <- qdf$tpm
  }

  treat_tpms <- list()
  for (sname in treat_samples) {
    qpath <- file.path(simdir, meth$dir, sname, paste0(sname, ".quant"))
    qdf <- read_quant(qpath)
    if (!is.null(qdf)) treat_tpms[[sname]] <- qdf$tpm
  }

  if (length(ctrl_tpms) == 0 || length(treat_tpms) == 0) next

  # Use first sample's names as reference
  qpath_ref <- file.path(simdir, meth$dir, ctrl_samples[1], paste0(ctrl_samples[1], ".quant"))
  ref_df <- read_quant(qpath_ref)

  mean_ctrl <- Reduce("+", ctrl_tpms) / length(ctrl_tpms)
  mean_treat <- Reduce("+", treat_tpms) / length(treat_tpms)

  est_log2fc <- log2(mean_treat + pseudo) - log2(mean_ctrl + pseudo)

  # Ground truth fold change
  m <- match(ref_df$target_name, gt$short_id)
  valid <- !is.na(m)
  true_ctrl <- gt$expected_tpm_control[m[valid]]
  true_treat <- gt$expected_tpm_treatment[m[valid]]
  true_log2fc <- log2(true_treat + pseudo) - log2(true_ctrl + pseudo)

  est_fc <- est_log2fc[valid]

  # Only look at DE transcripts with substantial true FC
  de_mask <- ref_df$target_name[valid] %in% gt$short_id[gt$is_de]
  big_fc <- abs(true_log2fc) > 1

  if (sum(de_mask & big_fc) > 10) {
    fit <- lm(est_fc[de_mask & big_fc] ~ true_log2fc[de_mask & big_fc])
    slope <- coef(fit)[2]
    r2 <- summary(fit)$r.squared
    cat(sprintf("  %-15s: FC attenuation slope=%.3f, R²=%.3f (n=%d DE transcripts with |FC|>1)\n",
                meth$name, slope, r2, sum(de_mask & big_fc)))
  }
}

cat("\nDone.\n")
