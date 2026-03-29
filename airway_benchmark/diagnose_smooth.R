#!/usr/bin/env Rscript
#
# Diagnose WHY coverage smoothness is increasing CV.
# Look at the specific changes smoothing makes to isoform abundances.

library(data.table)

benchdir <- "airway_benchmark"
samples <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")

load_quant <- function(method_dir) {
  dfs <- list()
  for (s in samples) {
    path <- file.path(benchdir, method_dir, paste0(s, ".quant"))
    df <- fread(path, sep = "\t")
    setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
    dfs[[s]] <- df[, .(target_name, tpm)]
    setnames(dfs[[s]], "tpm", s)
  }
  Reduce(function(a, b) merge(a, b, by = "target_name", all = TRUE), dfs)
}

parse_gene <- function(target_name) {
  sapply(strsplit(target_name, "\\|"), function(x) x[6])
}

em <- load_quant("quant/em")
sm <- load_quant("quant/smooth_ma")

# Add gene info
em$gene_name <- parse_gene(em$target_name)
sm$gene_name <- parse_gene(sm$target_name)

# ---- 1. Mass redistribution: how does smoothing change total TPM per gene? ----
cat("=== 1. TPM conservation check ===\n")
cat("Total TPM should be ~1M per sample for both methods.\n\n")

for (s in samples) {
  cat(sprintf("  %s: EM total=%.0f, Smooth total=%.0f\n",
              s, sum(em[[s]]), sum(sm[[s]])))
}

# ---- 2. How many isoforms gain/lose TPM? ----
cat("\n=== 2. Isoform TPM changes (averaged across samples) ===\n")

em_mean <- rowMeans(as.matrix(em[, ..samples]))
sm_mean <- rowMeans(as.matrix(sm[, ..samples]))

# Pair them
comp <- data.table(
  target_name = em$target_name,
  gene_name = em$gene_name,
  em_mean = em_mean,
  sm_mean = sm_mean,
  diff = sm_mean - em_mean,
  ratio = sm_mean / (em_mean + 1e-6)
)

# Expressed in EM (TPM >= 1)
expr <- comp[em_mean >= 1]
cat(sprintf("\nAmong %d isoforms expressed in EM (mean TPM >= 1):\n", nrow(expr)))
cat(sprintf("  Gained TPM (smooth > EM):  %d (%.1f%%)\n",
            sum(expr$diff > 0.1), sum(expr$diff > 0.1) / nrow(expr) * 100))
cat(sprintf("  Lost TPM (smooth < EM):    %d (%.1f%%)\n",
            sum(expr$diff < -0.1), sum(expr$diff < -0.1) / nrow(expr) * 100))
cat(sprintf("  Unchanged (|diff| < 0.1):  %d (%.1f%%)\n",
            sum(abs(expr$diff) < 0.1), sum(abs(expr$diff) < 0.1) / nrow(expr) * 100))

# Newly expressed (EM < 0.1, smooth >= 1)
newly_expr <- comp[em_mean < 0.1 & sm_mean >= 1]
cat(sprintf("\nNewly expressed (EM < 0.1, Smooth >= 1 TPM): %d isoforms\n", nrow(newly_expr)))

# Lost expression (EM >= 1, smooth < 0.1)
lost_expr <- comp[em_mean >= 1 & sm_mean < 0.1]
cat(sprintf("Lost expression (EM >= 1, Smooth < 0.1 TPM): %d isoforms\n", nrow(lost_expr)))

# ---- 3. Per-sample variability of the CHANGES ----
cat("\n=== 3. Are the changes consistent across replicates? ===\n")
cat("For each isoform, compute the per-sample diff (smooth - EM).\n")
cat("If changes are consistent, the SD of diffs should be small.\n\n")

em_mat <- as.matrix(em[, ..samples])
sm_mat <- as.matrix(sm[, ..samples])
diff_mat <- sm_mat - em_mat

# For expressed isoforms
expr_mask <- em_mean >= 1
diff_expr <- diff_mat[expr_mask, ]
mean_diff <- rowMeans(diff_expr)
sd_diff <- apply(diff_expr, 1, sd)
cv_diff <- sd_diff / (abs(mean_diff) + 1e-6)

cat(sprintf("Among %d expressed isoforms:\n", sum(expr_mask)))
cat(sprintf("  Median |mean_diff|: %.4f TPM\n", median(abs(mean_diff))))
cat(sprintf("  Median SD of diff:  %.4f TPM\n", median(sd_diff)))
cat(sprintf("  Median CV of diff:  %.4f\n", median(cv_diff[abs(mean_diff) > 0.1])))

# ---- 4. Direction consistency: does smoothing move the same isoform the same way? ----
cat("\n=== 4. Direction consistency ===\n")
cat("For isoforms where |mean_diff| > 1 TPM, how often does the sign\n")
cat("of (smooth - EM) agree across all 4 replicates?\n\n")

big_movers <- which(expr_mask & abs(mean_diff) > 1)
n_all_agree <- 0
n_some_disagree <- 0
for (i in big_movers) {
  diffs <- diff_mat[i, ]
  if (all(diffs > 0) || all(diffs < 0)) {
    n_all_agree <- n_all_agree + 1
  } else {
    n_some_disagree <- n_some_disagree + 1
  }
}
cat(sprintf("  Isoforms with |mean diff| > 1 TPM: %d\n", length(big_movers)))
cat(sprintf("  All 4 replicates agree on direction: %d (%.1f%%)\n",
            n_all_agree, n_all_agree / length(big_movers) * 100))
cat(sprintf("  Disagreement across replicates:      %d (%.1f%%)\n",
            n_some_disagree, n_some_disagree / length(big_movers) * 100))

# ---- 5. Look at specific genes where smoothing hurts most ----
cat("\n=== 5. Case studies: genes where smoothing increases CV most ===\n\n")

em_cv <- apply(em_mat, 1, function(x) sd(x) / (mean(x) + 1e-6))
sm_cv <- apply(sm_mat, 1, function(x) sd(x) / (mean(x) + 1e-6))

comp$cv_em <- em_cv
comp$cv_sm <- sm_cv
comp$cv_diff <- sm_cv - em_cv

# Per-gene aggregation
gene_iso_count <- comp[, .N, by = gene_name]
setnames(gene_iso_count, "N", "n_isoforms")
comp <- merge(comp, gene_iso_count, by = "gene_name")

# Look at complex genes where smoothing hurts
bad_genes <- comp[n_isoforms >= 10 & em_mean >= 1,
                  .(mean_cv_em = mean(cv_em), mean_cv_sm = mean(cv_sm),
                    n_expr = .N, n_isoforms = n_isoforms[1]),
                  by = gene_name]
bad_genes$cv_increase <- bad_genes$mean_cv_sm - bad_genes$mean_cv_em
bad_genes <- bad_genes[order(-cv_increase)]

for (g in head(bad_genes$gene_name, 5)) {
  cat(sprintf("--- %s (%d isoforms) ---\n", g,
              bad_genes[gene_name == g]$n_isoforms))
  gene_data <- comp[gene_name == g & em_mean >= 0.5]
  gene_data <- gene_data[order(-em_mean)]

  # Show per-sample TPM for EM and smooth
  cat(sprintf("  %-20s | %8s %8s %8s %8s | %8s %8s %8s %8s | %6s %6s\n",
              "Isoform", "EM_S1", "EM_S2", "EM_S3", "EM_S4",
              "SM_S1", "SM_S2", "SM_S3", "SM_S4", "CV_EM", "CV_SM"))

  em_sub <- em[em$target_name %in% gene_data$target_name]
  sm_sub <- sm[sm$target_name %in% gene_data$target_name]

  for (i in seq_len(min(8, nrow(gene_data)))) {
    tn <- gene_data$target_name[i]
    short <- sub("\\|.*", "", tn)
    em_row <- as.numeric(em_sub[target_name == tn, ..samples])
    sm_row <- as.numeric(sm_sub[target_name == tn, ..samples])
    cat(sprintf("  %-20s | %8.1f %8.1f %8.1f %8.1f | %8.1f %8.1f %8.1f %8.1f | %6.3f %6.3f\n",
                short,
                em_row[1], em_row[2], em_row[3], em_row[4],
                sm_row[1], sm_row[2], sm_row[3], sm_row[4],
                gene_data$cv_em[i], gene_data$cv_sm[i]))
  }
  cat("\n")
}

# ---- 6. Coverage profile comparison: is the smoothing penalty actually high? ----
cat("=== 6. Magnitude of TPM redistribution within genes ===\n")
cat("For each gene, compute: total absolute TPM change / total gene TPM\n\n")

gene_stats <- comp[, .(
  total_em_tpm = sum(em_mean),
  total_sm_tpm = sum(sm_mean),
  total_abs_diff = sum(abs(diff)),
  n_isoforms = .N
), by = gene_name]

gene_stats$frac_redistributed <- gene_stats$total_abs_diff / (gene_stats$total_em_tpm + 1e-6)
gene_stats <- gene_stats[total_em_tpm >= 10]  # only genes with substantial expression

cat(sprintf("Among %d genes with total TPM >= 10:\n", nrow(gene_stats)))
cat(sprintf("  Median fraction redistributed: %.4f (%.2f%%)\n",
            median(gene_stats$frac_redistributed),
            median(gene_stats$frac_redistributed) * 100))

# By complexity
for (tier in list(c(1,1), c(2,5), c(6,10), c(11,20), c(21,50), c(51,9999))) {
  sub <- gene_stats[n_isoforms >= tier[1] & n_isoforms <= tier[2]]
  if (nrow(sub) > 0) {
    cat(sprintf("  %3d-%3d isoforms (%4d genes): median redistrib = %.4f (%.2f%%)\n",
                tier[1], tier[2], nrow(sub),
                median(sub$frac_redistributed),
                median(sub$frac_redistributed) * 100))
  }
}

cat("\n=== Done ===\n")
