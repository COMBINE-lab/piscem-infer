#!/usr/bin/env Rscript
# Evaluate structural selection + consensus methods vs baselines.

library(data.table)

simdir <- "sim_data_gencode"
pseudo <- 0.01

gt <- fread(file.path(simdir, "ground_truth.csv"))
gt[, short_id := sub(" .*", "", transcript_id)]
sample_info <- fread(file.path(simdir, "sample_info.csv"))
ctrl_samples  <- sample_info$sample_name[sample_info$condition == "control"]
treat_samples <- sample_info$sample_name[sample_info$condition == "treatment"]

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

methods <- list(
  list(name = "Single EM",       dir = "quant_em",                      filter_set = NULL),
  list(name = "Cons TPM",        dir = "quant_consensus",               filter_set = NULL),
  list(name = "Cons UES",        dir = "quant_consensus_ues",           filter_set = NULL),
  list(name = "Cons Support",    dir = "quant_consensus_support",       filter_set = NULL),
  list(name = "Sel+Cons TPM",    dir = "quant_consensus_sel_tpm",       filter_set = NULL),
  list(name = "Sel+Cons UES",    dir = "quant_consensus_sel_ues",       filter_set = NULL),
  list(name = "Sel+Cons Supp",   dir = "quant_consensus_sel_support",   filter_set = NULL),
  list(name = "NoCond+AV",       dir = "quant_hier_nocond_av",          filter_set = NULL)
)

cat(sprintf("%-16s | %6s %6s %6s | %8s %8s %8s | %8s %8s %8s %8s\n",
            "Method", "TP", "FP", "FN", "Prec", "Recall", "F1",
            "Pearson", "Spearman", "RMSE", "FC slope"))
cat(paste0(rep("-", 120), collapse = ""), "\n")

for (meth in methods) {
  all_pearson <- c(); all_spearman <- c(); all_rmse <- c()
  all_tp <- c(); all_fp <- c(); all_fn <- c()

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
    true_t <- gt[[tpm_col]][m[valid]]
    est <- est_tpm[valid]

    all_pearson <- c(all_pearson, cor(log2(est + pseudo), log2(true_t + pseudo), method = "pearson"))
    all_spearman <- c(all_spearman, cor(est, true_t, method = "spearman"))
    all_rmse <- c(all_rmse, sqrt(mean((log2(est + 1) - log2(true_t + 1))^2)))

    all_tp <- c(all_tp, sum(est > 0 & true_t > 0))
    all_fp <- c(all_fp, sum(est > 0 & true_t == 0))
    all_fn <- c(all_fn, sum(est == 0 & true_t > 0))
  }

  tp <- mean(all_tp); fp <- mean(all_fp); fn_ <- mean(all_fn)
  prec <- tp / (tp + fp)
  rec <- tp / (tp + fn_)
  f1 <- 2 * prec * rec / (prec + rec)

  # Fold change
  ctrl_tpms <- list(); treat_tpms <- list()
  for (sn in ctrl_samples) {
    qdf <- read_quant(file.path(simdir, meth$dir, sn, paste0(sn, ".quant")))
    if (is.null(qdf)) next
    est_tpm <- apply_filter_and_renorm(qdf, meth$filter_set)
    ctrl_tpms[[sn]] <- est_tpm
  }
  for (sn in treat_samples) {
    qdf <- read_quant(file.path(simdir, meth$dir, sn, paste0(sn, ".quant")))
    if (is.null(qdf)) next
    est_tpm <- apply_filter_and_renorm(qdf, meth$filter_set)
    treat_tpms[[sn]] <- est_tpm
  }

  slope_str <- "   N/A"
  if (length(ctrl_tpms) > 0 && length(treat_tpms) > 0) {
    mean_ctrl  <- Reduce("+", ctrl_tpms)  / length(ctrl_tpms)
    mean_treat <- Reduce("+", treat_tpms) / length(treat_tpms)
    est_lfc    <- log2(mean_treat + pseudo) - log2(mean_ctrl + pseudo)

    ref_names <- read_quant(file.path(simdir, meth$dir, ctrl_samples[1],
                                       paste0(ctrl_samples[1], ".quant")))$target_name
    m <- match(ref_names, gt$short_id)
    valid <- !is.na(m)
    true_ctrl <- gt$expected_tpm_control[m[valid]]
    true_treat <- gt$expected_tpm_treatment[m[valid]]
    true_lfc <- log2(true_treat + pseudo) - log2(true_ctrl + pseudo)
    de_mask <- ref_names[valid] %in% gt$short_id[gt$is_de]
    big_fc <- abs(true_lfc) > 1

    if (sum(de_mask & big_fc) > 10) {
      fit <- lm(est_lfc[valid][de_mask & big_fc] ~ true_lfc[de_mask & big_fc])
      slope_str <- sprintf("%6.3f", coef(fit)[2])
    }
  }

  cat(sprintf("%-16s | %6d %6d %6d | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f %8s\n",
              meth$name, round(tp), round(fp), round(fn_),
              prec, rec, f1,
              mean(all_pearson), mean(all_spearman), mean(all_rmse), slope_str))
}
