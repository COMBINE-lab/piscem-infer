#!/usr/bin/env Rscript
#
# Evaluate replicate concordance across quantification methods.
#
# For each method, loads TPM estimates from 4 untreated biological replicates
# and computes per-isoform coefficient of variation (CV) across replicates.
# Stratifies by gene complexity (number of isoforms per gene).
#
# A method that better resolves isoform ambiguity should produce lower CV
# across replicates, especially at complex loci.

library(data.table)

benchdir <- "airway_benchmark"
samples <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")

methods <- list(
  "Plain EM"     = "quant/em",
  "NullSink .50e3" = "quant/nullsink_u4",
  "NullSink .50e3 P2Sq" = "quant/nullsink_u4_phase2squarem",
  "NullSink .50e3 Fast" = "quant/nullsink_u4_fast",
  "EM+Sel"       = "quant/sel",
  "Smooth(MA)"   = "quant/smooth_ma",
  "Smooth(SG)"   = "quant/smooth_sg",
  "MA+Sel"       = "quant/smooth_ma_sel",
  "SG+Sel"       = "quant/smooth_sg_sel"
)

# ---- Load all quant files ----
cat("Loading quantification results...\n")

load_method <- function(method_dir) {
  dfs <- list()
  for (s in samples) {
    path <- file.path(benchdir, method_dir, paste0(s, ".quant"))
    if (!file.exists(path)) {
      cat("  MISSING:", path, "\n")
      return(NULL)
    }
    df <- fread(path, sep = "\t")
    setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
    dfs[[s]] <- df[, .(target_name, tpm)]
    setnames(dfs[[s]], "tpm", s)
  }
  # Merge all samples
  merged <- Reduce(function(a, b) merge(a, b, by = "target_name", all = TRUE), dfs)
  merged
}

all_results <- list()
for (name in names(methods)) {
  cat("  Loading", name, "...\n")
  all_results[[name]] <- load_method(methods[[name]])
}

# ---- Compute gene complexity ----
# Parse gene name from GENCODE header: field 6 (pipe-delimited)
parse_gene <- function(target_name) {
  sapply(strsplit(target_name, "\\|"), function(x) x[6])
}

# Use the first method's target list
ref <- all_results[[1]]
ref$gene_name <- parse_gene(ref$target_name)
gene_iso_count <- ref[, .N, by = gene_name]
setnames(gene_iso_count, "N", "n_isoforms")

ref <- merge(ref, gene_iso_count, by = "gene_name")

# ---- Concordance analysis ----
cat("\n=== Replicate Concordance (CV of TPM across 4 replicates) ===\n\n")

compute_cv <- function(merged, min_mean_tpm = 1.0) {
  # Only consider transcripts with mean TPM >= threshold (expressed)
  tpm_mat <- as.matrix(merged[, ..samples])
  row_mean <- rowMeans(tpm_mat)
  row_sd <- apply(tpm_mat, 1, sd)
  cv <- row_sd / (row_mean + 1e-6)  # add small constant to avoid /0

  data.table(
    target_name = merged$target_name,
    mean_tpm = row_mean,
    sd_tpm = row_sd,
    cv = cv
  )
}

# Global metrics
cat(sprintf("%-14s | %8s %8s %8s | %8s %8s | %8s\n",
            "Method", "Med CV", "Mean CV", "CV>1",
            "Med CV>5", "Med CV>50", "N expr"))
cat(paste0(rep("-", 85), collapse = ""), "\n")

for (name in names(all_results)) {
  m <- all_results[[name]]
  if (is.null(m)) { cat(sprintf("%-14s | MISSING\n", name)); next }

  cv_df <- compute_cv(m)

  # Filter to expressed transcripts (mean TPM >= 1)
  expr <- cv_df[mean_tpm >= 1]
  expr5 <- cv_df[mean_tpm >= 5]
  expr50 <- cv_df[mean_tpm >= 50]

  cat(sprintf("%-14s | %8.4f %8.4f %8d | %8.4f %8.4f | %8d\n",
              name,
              median(expr$cv), mean(expr$cv), sum(expr$cv > 1),
              median(expr5$cv), median(expr50$cv),
              nrow(expr)))
}

# ---- Stratified by gene complexity ----
cat("\n=== Stratified by gene complexity (isoforms per gene) ===\n")
cat("(Median CV of expressed isoforms, TPM >= 1)\n\n")

# Complexity tiers
tiers <- list(
  "1 iso"     = c(1, 1),
  "2-5 iso"   = c(2, 5),
  "6-10 iso"  = c(6, 10),
  "11-20 iso" = c(11, 20),
  "21-50 iso" = c(21, 50),
  "51+ iso"   = c(51, 9999)
)

cat(sprintf("%-14s", "Method"))
for (t in names(tiers)) cat(sprintf(" | %10s", t))
cat("\n")
cat(paste0(rep("-", 16 + length(tiers) * 13), collapse = ""), "\n")

for (name in names(all_results)) {
  m <- all_results[[name]]
  if (is.null(m)) next

  cv_df <- compute_cv(m)
  cv_df$gene_name <- parse_gene(cv_df$target_name)
  cv_df <- merge(cv_df, gene_iso_count, by = "gene_name")

  cat(sprintf("%-14s", name))
  for (tier_name in names(tiers)) {
    bounds <- tiers[[tier_name]]
    tier <- cv_df[n_isoforms >= bounds[1] & n_isoforms <= bounds[2] & mean_tpm >= 1]
    if (nrow(tier) > 0) {
      cat(sprintf(" | %10.4f", median(tier$cv)))
    } else {
      cat(sprintf(" | %10s", "n/a"))
    }
  }
  cat("\n")
}

# ---- Per-gene CV comparison ----
cat("\n=== Per-gene mean CV comparison (genes with 10+ isoforms, >=3 expressed) ===\n")
cat("Showing genes where smoothing changes CV most\n\n")

# For the plain EM baseline and best smooth method
em <- all_results[["Plain EM"]]
smooth <- all_results[["Smooth(MA)"]]
if (!is.null(em) && !is.null(smooth)) {
  cv_em <- compute_cv(em)
  cv_sm <- compute_cv(smooth)
  cv_em$gene_name <- parse_gene(cv_em$target_name)
  cv_sm$gene_name <- parse_gene(cv_sm$target_name)

  cv_em <- merge(cv_em, gene_iso_count, by = "gene_name")
  cv_sm <- merge(cv_sm, gene_iso_count, by = "gene_name")

  # Per-gene mean CV for expressed isoforms
  gene_cv_em <- cv_em[mean_tpm >= 1, .(mean_cv_em = mean(cv), n_expr = .N), by = .(gene_name, n_isoforms)]
  gene_cv_sm <- cv_sm[mean_tpm >= 1, .(mean_cv_sm = mean(cv)), by = gene_name]

  gene_comp <- merge(gene_cv_em, gene_cv_sm, by = "gene_name")
  gene_comp$cv_diff <- gene_comp$mean_cv_em - gene_comp$mean_cv_sm
  gene_comp$cv_ratio <- gene_comp$mean_cv_em / (gene_comp$mean_cv_sm + 1e-6)

  # Filter to complex genes
  complex <- gene_comp[n_isoforms >= 10 & n_expr >= 3]
  complex <- complex[order(-abs(cv_diff))]

  cat(sprintf("%-15s %5s %5s | %8s %8s | %8s\n",
              "Gene", "#Iso", "#Expr", "CV(EM)", "CV(Smooth)", "Diff"))
  cat(paste0(rep("-", 65), collapse = ""), "\n")

  for (i in seq_len(min(20, nrow(complex)))) {
    r <- complex[i]
    cat(sprintf("%-15s %5d %5d | %8.4f %8.4f | %+8.4f\n",
                r$gene_name, r$n_isoforms, r$n_expr,
                r$mean_cv_em, r$mean_cv_sm, r$cv_diff))
  }
}

# ---- Summary statistics ----
cat("\n=== Summary: Paired Wilcoxon test (per-isoform CV) ===\n")
cat("Is smooth CV significantly different from plain EM CV?\n\n")

em <- all_results[["Plain EM"]]
for (smooth_name in c("Smooth(MA)", "Smooth(SG)", "MA+Sel", "SG+Sel")) {
  smooth <- all_results[[smooth_name]]
  if (is.null(em) || is.null(smooth)) next

  cv_em <- compute_cv(em)
  cv_sm <- compute_cv(smooth)

  # Match by target name
  both <- merge(cv_em, cv_sm, by = "target_name", suffixes = c(".em", ".sm"))
  both <- both[mean_tpm.em >= 1 & mean_tpm.sm >= 1]

  wt <- wilcox.test(both$cv.em, both$cv.sm, paired = TRUE, alternative = "greater")
  n_better <- sum(both$cv.sm < both$cv.em)
  n_worse <- sum(both$cv.sm > both$cv.em)

  cat(sprintf("  %-14s: p=%.2e, %d improved / %d worsened / %d tied (n=%d)\n",
              smooth_name, wt$p.value, n_better, n_worse,
              nrow(both) - n_better - n_worse, nrow(both)))
}

cat("\n=== Done ===\n")
