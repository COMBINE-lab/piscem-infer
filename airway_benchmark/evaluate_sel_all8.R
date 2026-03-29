#!/usr/bin/env Rscript
# Evaluate structural selection + consensus on all 8 airway samples.
# Within-condition replicate concordance + between-condition fold change.

library(data.table)

benchdir <- "airway_benchmark"

samples_untreated <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")
samples_treated   <- c("SRR1039509", "SRR1039513", "SRR1039517", "SRR1039521")
all_samples <- c(samples_untreated, samples_treated)

methods <- list(
  "Plain EM"         = list(dir = "quant/em",                    subdir = FALSE),
  "Cons Supp(8)"     = list(dir = "quant/consensus_support_all8", subdir = TRUE),
  "Sel+Supp(8)"      = list(dir = "quant/sel_support_all8",      subdir = TRUE),
  "Sel+Supp(8)CA"    = list(dir = "quant/sel_support_all8_ca",   subdir = TRUE)
)

load_method <- function(method_info, sample_list) {
  dfs <- list()
  for (s in sample_list) {
    if (method_info$subdir) {
      path <- file.path(benchdir, method_info$dir, s, paste0(s, ".quant"))
    } else {
      path <- file.path(benchdir, method_info$dir, paste0(s, ".quant"))
    }
    if (!file.exists(path)) {
      cat("  MISSING:", path, "\n")
      return(NULL)
    }
    df <- fread(path, sep = "\t")
    setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
    dfs[[s]] <- df[, .(target_name, tpm)]
    setnames(dfs[[s]], "tpm", s)
  }
  Reduce(function(a, b) merge(a, b, by = "target_name", all = TRUE), dfs)
}

compute_cv <- function(merged, sample_subset) {
  tpm_mat <- as.matrix(merged[, ..sample_subset])
  row_mean <- rowMeans(tpm_mat)
  row_sd <- apply(tpm_mat, 1, sd)
  cv <- row_sd / (row_mean + 1e-6)
  data.table(target_name = merged$target_name, mean_tpm = row_mean, cv = cv)
}

# ---- Load ----
cat("Loading...\n")
all_results <- list()
for (name in names(methods)) {
  cat("  ", name, "...\n")
  all_results[[name]] <- load_method(methods[[name]], all_samples)
}

# Gene complexity for stratification
parse_gene <- function(target_name) {
  sapply(strsplit(target_name, "\\|"), function(x) x[6])
}
ref <- all_results[[1]]
if (is.null(ref)) stop("Plain EM not found")
ref$gene_name <- parse_gene(ref$target_name)
gene_iso_count <- ref[, .N, by = gene_name]
setnames(gene_iso_count, "N", "n_isoforms")

# ---- Within-condition CV ----
for (cond_name in c("Untreated", "Treated")) {
  subs <- if (cond_name == "Untreated") samples_untreated else samples_treated
  cat(sprintf("\n=== %s Replicate Concordance (CV, mean TPM >= 1) ===\n\n", cond_name))
  cat(sprintf("%-18s | %8s %8s %8s | %8s\n",
              "Method", "Med CV", "Mean CV", "CV>1", "N expr"))
  cat(paste0(rep("-", 65), collapse = ""), "\n")

  for (name in names(all_results)) {
    m <- all_results[[name]]
    if (is.null(m)) { cat(sprintf("%-18s | MISSING\n", name)); next }
    cv_df <- compute_cv(m, subs)
    expr <- cv_df[mean_tpm >= 1]
    cat(sprintf("%-18s | %8.4f %8.4f %8d | %8d\n",
                name, median(expr$cv), mean(expr$cv), sum(expr$cv > 1), nrow(expr)))
  }
}

# ---- Stratified by gene complexity (untreated only) ----
cat("\n=== Untreated: Stratified by gene complexity ===\n")
cat("(Median CV, TPM >= 1)\n\n")

tiers <- list(
  "1 iso"     = c(1, 1),
  "2-5 iso"   = c(2, 5),
  "6-10 iso"  = c(6, 10),
  "11-20 iso" = c(11, 20),
  "21-50 iso" = c(21, 50),
  "51+ iso"   = c(51, 9999)
)

cat(sprintf("%-18s", "Method"))
for (t in names(tiers)) cat(sprintf(" | %10s", t))
cat("\n")
cat(paste0(rep("-", 20 + length(tiers) * 13), collapse = ""), "\n")

for (name in names(all_results)) {
  m <- all_results[[name]]
  if (is.null(m)) next

  cv_df <- compute_cv(m, samples_untreated)
  cv_df$gene_name <- parse_gene(cv_df$target_name)
  cv_df <- merge(cv_df, gene_iso_count, by = "gene_name")

  cat(sprintf("%-18s", name))
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

# ---- Between-condition fold change ----
cat("\n=== Between-condition fold change (dex vs untreated) ===\n")
cat("(Isoforms with mean TPM >= 1 in both conditions)\n\n")

pseudo <- 0.01
cat(sprintf("%-18s | %10s %10s %10s | %8s\n",
            "Method", "Mean|LFC|", "Med|LFC|", "SD(LFC)", "N both"))
cat(paste0(rep("-", 75), collapse = ""), "\n")

for (name in names(all_results)) {
  m <- all_results[[name]]
  if (is.null(m)) next
  untrt_mean <- rowMeans(as.matrix(m[, ..samples_untreated]))
  trt_mean   <- rowMeans(as.matrix(m[, ..samples_treated]))
  both_expr <- untrt_mean >= 1 & trt_mean >= 1
  if (sum(both_expr) > 0) {
    lfc <- log2(trt_mean[both_expr] + pseudo) - log2(untrt_mean[both_expr] + pseudo)
    cat(sprintf("%-18s | %10.4f %10.4f %10.4f | %8d\n",
                name, mean(abs(lfc)), median(abs(lfc)), sd(lfc), sum(both_expr)))
  }
}

# ---- Pairwise correlation of fold changes ----
cat("\n=== Fold change correlation (Spearman) between methods ===\n")
cat("(Based on log2FC for isoforms expressed in both conditions across all methods)\n\n")

# Compute LFC for all methods
all_lfc <- list()
common_names <- NULL
for (name in names(all_results)) {
  m <- all_results[[name]]
  if (is.null(m)) next
  untrt_mean <- rowMeans(as.matrix(m[, ..samples_untreated]))
  trt_mean   <- rowMeans(as.matrix(m[, ..samples_treated]))
  both_expr <- untrt_mean >= 1 & trt_mean >= 1
  lfc_dt <- data.table(
    target_name = m$target_name[both_expr],
    lfc = log2(trt_mean[both_expr] + pseudo) - log2(untrt_mean[both_expr] + pseudo)
  )
  setnames(lfc_dt, "lfc", name)
  all_lfc[[name]] <- lfc_dt
  if (is.null(common_names)) {
    common_names <- lfc_dt$target_name
  } else {
    common_names <- intersect(common_names, lfc_dt$target_name)
  }
}

# Print correlation of LFC between Plain EM and each method
cat(sprintf("%-18s | %10s\n", "Method", "Cor w/ EM"))
cat(paste0(rep("-", 35), collapse = ""), "\n")
em_lfc <- all_lfc[["Plain EM"]]
if (!is.null(em_lfc)) {
  for (name in names(all_lfc)) {
    m_lfc <- all_lfc[[name]]
    merged <- merge(em_lfc, m_lfc, by = "target_name")
    r <- cor(merged[["Plain EM"]], merged[[name]], method = "spearman")
    cat(sprintf("%-18s | %10.4f\n", name, r))
  }
}

cat("\nDone.\n")
