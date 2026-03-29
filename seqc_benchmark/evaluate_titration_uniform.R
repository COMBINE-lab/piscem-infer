#!/usr/bin/env Rscript
#
# SEQC Titration evaluation on UNIFORM transcript ground set.
# All methods evaluated on the same set of transcripts — the union of
# transcripts with mean TPM >= 1 in A or B across ANY method.
# Undetected transcripts get TPM=0.

library(data.table)

benchdir <- "seqc_benchmark"
pseudo <- 0.01

read_piscem <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- fread(path, sep = "\t")
  setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
  df[, .(target_name, tpm)]
}

read_salmon <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- fread(path, sep = "\t")
  df[, .(target_name = Name, tpm = TPM)]
}

read_kallisto <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- fread(path, sep = "\t")
  df[, .(target_name = target_id, tpm)]
}

samples <- c("A", "B", "C", "D")
reps <- 1:4

methods <- list(
  "Plain EM" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_em", paste0(s, "_", r, ".quant"))
  ),
  "Sel+Supp" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_support", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Sel+Adapt" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_adapt", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Adpt+Gene10" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_adapt_generescue", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Cond+Gene10" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_adapt_condgene", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Pos5+CG" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_pos5_condgene", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Adpt+Pos5" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_adapt_pos5", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Adpt+GeneRsc" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_adapt_generescue", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Adpt+Rescue" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_adapt_rescue", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Salmon" = list(
    reader = read_salmon,
    path_fn = function(s, r) file.path(benchdir, "quant_salmon", paste0(s, "_", r), "quant.sf")
  ),
  "Salmon gcB" = list(
    reader = read_salmon,
    path_fn = function(s, r) file.path(benchdir, "quant_salmon_gcbias", paste0(s, "_", r), "quant.sf")
  ),
  "Kallisto" = list(
    reader = read_kallisto,
    path_fn = function(s, r) file.path(benchdir, "quant_kallisto", paste0(s, "_", r), "abundance.tsv")
  )
)

# ---- Helper: load sample mean TPMs keyed by target_name ----
load_sample_means <- function(meth, sample_type) {
  tpm_list <- list()
  for (r in reps) {
    df <- meth$reader(meth$path_fn(sample_type, r))
    if (!is.null(df)) {
      setnames(df, "tpm", paste0("tpm_", r))
      tpm_list[[as.character(r)]] <- df
    }
  }
  if (length(tpm_list) == 0) return(NULL)
  merged <- Reduce(function(a, b) merge(a, b, by = "target_name", all = TRUE), tpm_list)
  tpm_cols <- setdiff(names(merged), "target_name")
  merged$mean_tpm <- rowMeans(as.matrix(merged[, ..tpm_cols]), na.rm = TRUE)
  merged[, .(target_name, mean_tpm)]
}

# ============================================================
# Build uniform ground set: union of transcripts with mean TPM >= 1
# in A or B from the most permissive method (plain EM).
# ============================================================

cat("Building uniform transcript ground set...\n")

# Use plain EM's full transcript list (all ~245K transcripts).
# Every method outputs all transcripts (TPM=0 for filtered ones),
# so the full set is the unbiased ground for comparison.
em_A <- load_sample_means(methods[["Plain EM"]], "A")
ground_set <- em_A$target_name
cat(sprintf("Ground set: %d transcripts (full annotation)\n", length(ground_set)))

# Helper: get mean TPMs on the ground set (0 for undetected)
get_on_ground <- function(meth, sample_type) {
  sm <- load_sample_means(meth, sample_type)
  if (is.null(sm)) return(rep(0, length(ground_set)))
  m <- match(ground_set, sm$target_name)
  result <- rep(0.0, length(ground_set))
  found <- !is.na(m)
  result[found] <- sm$mean_tpm[m[found]]
  result
}

# ============================================================
# Titration accuracy on uniform ground set
# ============================================================

cat("\n=== Titration Accuracy (UNIFORM ground set) ===\n")
cat("Expected C = 0.75*A + 0.25*B; Expected D = 0.25*A + 0.75*B\n")
cat(sprintf("All methods on same %d transcripts\n\n", length(ground_set)))

cat(sprintf("%-14s | %10s %10s | %10s %10s | %10s | %6s\n",
            "Method", "Pears(C)", "Pears(D)", "Spear(C)", "Spear(D)", "FC slope", "Det AB"))
cat(paste0(rep("-", 85), collapse = ""), "\n")

for (mname in names(methods)) {
  meth <- methods[[mname]]
  a <- get_on_ground(meth, "A")
  b <- get_on_ground(meth, "B")
  c_obs <- get_on_ground(meth, "C")
  d_obs <- get_on_ground(meth, "D")

  c_exp <- 0.75 * a + 0.25 * b
  d_exp <- 0.25 * a + 0.75 * b

  # Correlation on entire ground set (zeros included)
  pc <- cor(log2(c_obs + pseudo), log2(c_exp + pseudo), method = "pearson")
  pd <- cor(log2(d_obs + pseudo), log2(d_exp + pseudo), method = "pearson")
  sc <- cor(c_obs, c_exp, method = "spearman")
  sd <- cor(d_obs, d_exp, method = "spearman")

  # FC slope: observed C/D vs expected C/D
  obs_lfc <- log2(c_obs + pseudo) - log2(d_obs + pseudo)
  exp_lfc <- log2(c_exp + pseudo) - log2(d_exp + pseudo)
  big_fc <- abs(exp_lfc) > 0.5
  slope_str <- "N/A"
  if (sum(big_fc) > 100) {
    fit <- lm(obs_lfc[big_fc] ~ exp_lfc[big_fc])
    slope_str <- sprintf("%10.4f", coef(fit)[2])
  }

  det_ab <- sum(a > 0 & b > 0)

  cat(sprintf("%-14s | %10.4f %10.4f | %10.4f %10.4f | %s | %6d\n",
              mname, pc, pd, sc, sd, slope_str, det_ab))
}

# ============================================================
# Within-group concordance (no change needed — CV is per-method)
# ============================================================

compute_cv <- function(merged) {
  sample_cols <- grep("^tpm_", names(merged), value = TRUE)
  tpm_mat <- as.matrix(merged[, ..sample_cols])
  rm <- rowMeans(tpm_mat, na.rm = TRUE)
  rs <- apply(tpm_mat, 1, sd, na.rm = TRUE)
  data.table(target_name = merged$target_name, mean_tpm = rm, cv = rs / (rm + 1e-6))
}

cat("\n=== Within-Group Replicate Concordance (CV, mean TPM >= 1) ===\n\n")
cat(sprintf("%-14s | %8s %8s %8s %8s | %8s\n",
            "Method", "CV(A)", "CV(B)", "CV(C)", "CV(D)", "Mean"))
cat(paste0(rep("-", 14 + 4 * 11 + 11), collapse = ""), "\n")

for (mname in names(methods)) {
  meth <- methods[[mname]]
  cvs <- c()
  cat(sprintf("%-14s", mname))
  for (s in samples) {
    tpm_list <- list()
    for (r in reps) {
      df <- meth$reader(meth$path_fn(s, r))
      if (!is.null(df)) {
        setnames(df, "tpm", paste0("tpm_", r))
        tpm_list[[as.character(r)]] <- df
      }
    }
    if (length(tpm_list) == 0) { cat(sprintf(" | %8s", "N/A")); next }
    merged <- Reduce(function(a, b) merge(a, b, by = "target_name", all = TRUE), tpm_list)
    cv_df <- compute_cv(merged)
    expr <- cv_df[mean_tpm >= 1]
    med_cv <- median(expr$cv, na.rm = TRUE)
    cvs <- c(cvs, med_cv)
    cat(sprintf(" | %8.4f", med_cv))
  }
  cat(sprintf(" | %8.4f\n", mean(cvs)))
}

cat("\nDone.\n")
