#!/usr/bin/env Rscript
#
# Evaluate replicate concordance across correction clamp levels.

library(data.table)

benchdir <- "airway_benchmark"
samples <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")

methods <- list(
  "Plain EM"     = "quant/em",
  "Struct c=0.02"= "quant/smooth_structural",
  "Struct c=0.10"= "quant/structural_c10",
  "Struct LBFGS" = "quant/structural_lbfgs"
)

load_method <- function(method_dir) {
  dfs <- list()
  for (s in samples) {
    path <- file.path(benchdir, method_dir, paste0(s, ".quant"))
    if (!file.exists(path)) return(NULL)
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

cat("Loading...\n")
all_results <- list()
for (name in names(methods)) {
  all_results[[name]] <- load_method(methods[[name]])
}

# Gene complexity
ref <- all_results[[1]]
ref$gene_name <- parse_gene(ref$target_name)
gene_iso_count <- ref[, .N, by = gene_name]
setnames(gene_iso_count, "N", "n_isoforms")

compute_cv <- function(merged) {
  tpm_mat <- as.matrix(merged[, ..samples])
  row_mean <- rowMeans(tpm_mat)
  row_sd <- apply(tpm_mat, 1, sd)
  cv <- row_sd / (row_mean + 1e-6)
  data.table(target_name = merged$target_name, mean_tpm = row_mean, cv = cv)
}

# ---- Global CV ----
cat("\n=== Replicate Concordance (CV of TPM, isoforms with mean TPM >= 1) ===\n\n")
cat(sprintf("%-14s | %8s %8s %8s | %8s\n",
            "Method", "Med CV", "Mean CV", "CV>1", "N expr"))
cat(paste0(rep("-", 60), collapse = ""), "\n")

for (name in names(all_results)) {
  m <- all_results[[name]]
  if (is.null(m)) { cat(sprintf("%-14s | MISSING\n", name)); next }
  cv_df <- compute_cv(m)
  expr <- cv_df[mean_tpm >= 1]
  cat(sprintf("%-14s | %8.4f %8.4f %8d | %8d\n",
              name, median(expr$cv), mean(expr$cv), sum(expr$cv > 1), nrow(expr)))
}

# ---- Stratified ----
cat("\n=== Stratified by gene complexity ===\n\n")
tiers <- list("1"=c(1,1), "2-5"=c(2,5), "6-10"=c(6,10),
              "11-20"=c(11,20), "21-50"=c(21,50), "51+"=c(51,9999))

cat(sprintf("%-14s", "Method"))
for (t in names(tiers)) cat(sprintf(" | %8s", t))
cat("\n")
cat(paste0(rep("-", 16 + length(tiers) * 11), collapse = ""), "\n")

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
    cat(sprintf(" | %8.4f", if (nrow(tier) > 0) median(tier$cv) else NA))
  }
  cat("\n")
}

# ---- Direction consistency ----
cat("\n=== Direction consistency (isoforms with |mean diff| > 1 TPM) ===\n")

em <- all_results[["Plain EM"]]
em_mat <- as.matrix(em[, ..samples])

for (name in names(all_results)) {
  if (name == "Plain EM") next
  m <- all_results[[name]]
  if (is.null(m)) next

  sm_mat <- as.matrix(m[, ..samples])
  diff_mat <- sm_mat - em_mat
  em_mean <- rowMeans(em_mat)
  mean_diff <- rowMeans(diff_mat)

  big <- which(em_mean >= 1 & abs(mean_diff) > 1)
  if (length(big) == 0) {
    cat(sprintf("  %-14s: no isoforms with |mean diff| > 1 TPM\n", name))
    next
  }
  n_agree <- sum(apply(diff_mat[big, , drop=FALSE], 1, function(x) all(x > 0) || all(x < 0)))
  cat(sprintf("  %-14s: %d/%d agree (%.1f%%), median |diff|=%.2f TPM\n",
              name, n_agree, length(big), n_agree/length(big)*100,
              median(abs(mean_diff[big]))))
}

# ---- Redistribution magnitude ----
cat("\n=== TPM redistribution (mean across samples) ===\n")

em_mean <- rowMeans(em_mat)

for (name in names(all_results)) {
  if (name == "Plain EM") next
  m <- all_results[[name]]
  if (is.null(m)) next

  sm_mat <- as.matrix(m[, ..samples])
  sm_mean <- rowMeans(sm_mat)

  abs_diff <- abs(sm_mean - em_mean)
  expr <- em_mean >= 1
  n_newly_expr <- sum(em_mean < 0.1 & sm_mean >= 1)
  cat(sprintf("  %-14s: median |diff|=%.4f, newly expressed=%d, total redistrib=%.0f TPM\n",
              name, median(abs_diff[expr]), n_newly_expr, sum(abs_diff)))
}

# ---- Wilcoxon ----
cat("\n=== Paired Wilcoxon (smooth CV < EM CV?) ===\n")
cv_em <- compute_cv(em)
for (name in names(all_results)) {
  if (name == "Plain EM") next
  m <- all_results[[name]]
  if (is.null(m)) next
  cv_sm <- compute_cv(m)
  both <- merge(cv_em, cv_sm, by = "target_name", suffixes = c(".em", ".sm"))
  both <- both[mean_tpm.em >= 1 & mean_tpm.sm >= 1]
  wt <- wilcox.test(both$cv.sm, both$cv.em, paired = TRUE, alternative = "less")
  n_better <- sum(both$cv.sm < both$cv.em)
  n_worse <- sum(both$cv.sm > both$cv.em)
  cat(sprintf("  %-14s: p=%.2e, %d better / %d worse (n=%d)\n",
              name, wt$p.value, n_better, n_worse, nrow(both)))
}

cat("\n=== Done ===\n")
