#!/usr/bin/env Rscript
# FC calibration experiment: use Phase 1 (pre-consensus) vs Phase 2 (post-consensus)
# estimates to calibrate fold change attenuation.
#
# Approach: For the Sel+Cons Support method, Phase 2 attenuates fold changes because
# reads from filtered transcripts redistribute. We estimate the attenuation by
# regressing Phase 2 FC against plain EM FC (which has no masking-induced attenuation)
# for high-confidence consensus transcripts, then apply the inverse correction.

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
  df <- fread(path, sep = "\t"); setnames(df, c("target_name","len","eelen","tpm","ecount")); df
}

# Load both methods: plain EM (no attenuation) and Sel+Cons Support (attenuated)
em_dir <- "quant_em"
cons_dir <- "quant_consensus_sel_support"

# Compute mean TPM per condition for each method
load_condition_means <- function(dir) {
  ctrl_list <- list(); treat_list <- list()
  for (sn in ctrl_samples) {
    qdf <- read_quant(file.path(simdir, dir, sn, paste0(sn, ".quant")))
    if (!is.null(qdf)) ctrl_list[[sn]] <- qdf$tpm
  }
  for (sn in treat_samples) {
    qdf <- read_quant(file.path(simdir, dir, sn, paste0(sn, ".quant")))
    if (!is.null(qdf)) treat_list[[sn]] <- qdf$tpm
  }
  list(
    ctrl = Reduce("+", ctrl_list) / length(ctrl_list),
    treat = Reduce("+", treat_list) / length(treat_list),
    names = read_quant(file.path(simdir, dir, ctrl_samples[1], paste0(ctrl_samples[1], ".quant")))$target_name
  )
}

em_means <- load_condition_means(em_dir)
cons_means <- load_condition_means(cons_dir)

# Compute LFCs
em_lfc <- log2(em_means$treat + pseudo) - log2(em_means$ctrl + pseudo)
cons_lfc <- log2(cons_means$treat + pseudo) - log2(cons_means$ctrl + pseudo)

# Match transcript names
m <- match(cons_means$names, em_means$names)
valid <- !is.na(m)

# Ground truth
m_gt <- match(cons_means$names[valid], gt$short_id)
valid_gt <- !is.na(m_gt)

true_ctrl <- gt$expected_tpm_control[m_gt[valid_gt]]
true_treat <- gt$expected_tpm_treatment[m_gt[valid_gt]]
true_lfc <- log2(true_treat + pseudo) - log2(true_ctrl + pseudo)
is_de <- cons_means$names[valid][valid_gt] %in% gt$short_id[gt$is_de]

# Focus on transcripts expressed in both conditions for both methods
both_expr <- cons_means$ctrl[valid] >= 1 & cons_means$treat[valid] >= 1 &
             em_means$ctrl[m[valid]] >= 1 & em_means$treat[m[valid]] >= 1

cat("=== FC Calibration Experiment ===\n\n")

# Method 1: Direct regression of consensus FC on EM FC
# The idea: EM FC has no masking attenuation, so it's a better reference for FC direction.
# Consensus FC is biased toward zero due to read redistribution.
cons_fc_both <- cons_lfc[valid][both_expr]
em_fc_both <- em_lfc[m[valid]][both_expr]

# Only use high-confidence transcripts for calibration (large |FC|)
high_fc <- abs(em_fc_both) > 0.5
if (sum(high_fc) > 100) {
  fit <- lm(cons_fc_both[high_fc] ~ em_fc_both[high_fc])
  slope <- coef(fit)[2]
  intercept <- coef(fit)[1]
  cat(sprintf("Calibration regression (cons_fc ~ em_fc): slope=%.4f, intercept=%.4f (n=%d)\n",
              slope, intercept, sum(high_fc)))

  # Apply inverse calibration: corrected_fc = (cons_fc - intercept) / slope
  corrected_fc <- (cons_lfc[valid] - intercept) / slope
} else {
  cat("Not enough high-FC transcripts for calibration\n")
  corrected_fc <- cons_lfc[valid]
}

# Evaluate all variants against ground truth
cat("\n=== FC slope (estimated vs true, DE transcripts with |true FC| > 1) ===\n\n")

eval_fc_slope <- function(est_lfc, label) {
  # Match to ground truth
  m2 <- match(cons_means$names[valid], gt$short_id)
  v2 <- !is.na(m2)
  true_c <- gt$expected_tpm_control[m2[v2]]
  true_t <- gt$expected_tpm_treatment[m2[v2]]
  t_lfc <- log2(true_t + pseudo) - log2(true_c + pseudo)
  de <- cons_means$names[valid][v2] %in% gt$short_id[gt$is_de]
  big <- abs(t_lfc) > 1

  if (sum(de & big) > 10) {
    fit <- lm(est_lfc[v2][de & big] ~ t_lfc[de & big])
    cat(sprintf("  %-25s: slope = %.4f\n", label, coef(fit)[2]))
  } else {
    cat(sprintf("  %-25s: insufficient data\n", label))
  }
}

eval_fc_slope(em_lfc[m[valid]], "Plain EM")
eval_fc_slope(cons_lfc[valid], "Sel+Cons Support")
eval_fc_slope(corrected_fc, "Sel+Cons Calibrated")

# Also test: what if we just use Phase 1 EM estimates (pre-masking) for FC?
# Phase 1 runs on the structurally filtered set but without consensus masking.
# This avoids masking-induced attenuation entirely.
cat("\n=== Alternative: Use Phase 1 estimates for FC ===\n")
cat("(Phase 1 = EM on structurally-filtered set, no consensus masking)\n")
cat("This could be implemented by outputting Phase 1 TPMs alongside Phase 2.\n")

# Method 2: Weighted average of EM FC and consensus FC
# Weight by EC support — transcripts with high support use consensus FC,
# transcripts with low support use EM FC (less masking effect on them).
cat("\n=== Method 2: Weighted blend ===\n")
blend_fc <- 0.5 * cons_lfc[valid] + 0.5 * em_lfc[m[valid]]
eval_fc_slope(blend_fc, "50/50 blend (EM + Cons)")

blend_fc2 <- 0.3 * cons_lfc[valid] + 0.7 * em_lfc[m[valid]]
eval_fc_slope(blend_fc2, "30/70 blend (EM-heavy)")

cat("\nDone.\n")
