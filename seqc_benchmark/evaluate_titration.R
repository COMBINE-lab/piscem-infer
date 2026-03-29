#!/usr/bin/env Rscript
#
# SEQC/MAQC-III Titration Consistency Evaluation
#
# Core test: Samples C (75%A+25%B) and D (25%A+75%B) have expected expression
# levels determined by the measured A and B concentrations. For each transcript:
#
#   expected_C = 0.75 * tpm_A + 0.25 * tpm_B
#   expected_D = 0.25 * tpm_A + 0.75 * tpm_B
#
# A good quantification method should produce C and D estimates that closely
# match these expected values. We evaluate:
#
# 1. Titration accuracy: Pearson/Spearman correlation of observed vs expected
# 2. Fold change recovery: slope of log2(C/D) observed vs expected
# 3. Within-group concordance: CV across replicates for A, B, C, D
# 4. Cross-method comparison: piscem-infer variants vs salmon vs kallisto

library(data.table)

benchdir <- "seqc_benchmark"
pseudo <- 0.01

# ============================================================
# Load quantification results
# ============================================================

# piscem-infer: TSV with columns target_name, len, eelen, tpm, ecount
read_piscem <- function(path) {
  if (!file.exists(path)) { warning("Not found: ", path); return(NULL) }
  df <- fread(path, sep = "\t")
  setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
  df[, .(target_name, tpm)]
}

# salmon: quant.sf with columns Name, Length, EffectiveLength, TPM, NumReads
read_salmon <- function(path) {
  if (!file.exists(path)) { warning("Not found: ", path); return(NULL) }
  df <- fread(path, sep = "\t")
  df[, .(target_name = Name, tpm = TPM)]
}

# kallisto: abundance.tsv with columns target_id, length, eff_length, est_counts, tpm
read_kallisto <- function(path) {
  if (!file.exists(path)) { warning("Not found: ", path); return(NULL) }
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
  "Adpt+Rescue" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_adapt_rescue", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Salmon" = list(
    reader = read_salmon,
    path_fn = function(s, r) file.path(benchdir, "quant_salmon", paste0(s, "_", r), "quant.sf")
  ),
  "Salmon def" = list(
    reader = read_salmon,
    path_fn = function(s, r) file.path(benchdir, "quant_salmon_nobias", paste0(s, "_", r), "quant.sf")
  ),
  "Soft Rescue" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_adapt_softrescue", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Salmon EM" = list(
    reader = read_salmon,
    path_fn = function(s, r) file.path(benchdir, "quant_salmon_em", paste0(s, "_", r), "quant.sf")
  ),
  "Kallisto" = list(
    reader = read_kallisto,
    path_fn = function(s, r) file.path(benchdir, "quant_kallisto", paste0(s, "_", r), "abundance.tsv")
  )
)

# ============================================================
# Load all data
# ============================================================

cat("Loading quantification results...\n")

load_sample_means <- function(meth, sample_type) {
  tpm_list <- list()
  for (r in reps) {
    path <- meth$path_fn(sample_type, r)
    df <- meth$reader(path)
    if (!is.null(df)) {
      tpm_list[[as.character(r)]] <- df
    }
  }
  if (length(tpm_list) == 0) return(NULL)

  # Merge all replicates
  merged <- Reduce(function(a, b) merge(a, b, by = "target_name", all = TRUE), tpm_list)
  tpm_cols <- setdiff(names(merged), "target_name")

  # Compute mean and CV
  tpm_mat <- as.matrix(merged[, ..tpm_cols])
  data.table(
    target_name = merged$target_name,
    mean_tpm = rowMeans(tpm_mat, na.rm = TRUE),
    sd_tpm = apply(tpm_mat, 1, sd, na.rm = TRUE),
    cv = apply(tpm_mat, 1, sd, na.rm = TRUE) / (rowMeans(tpm_mat, na.rm = TRUE) + 1e-6),
    n_reps = rowSums(!is.na(tpm_mat))
  )
}

# ============================================================
# Evaluation 1: Within-group replicate concordance
# ============================================================

cat("\n=== Within-Group Replicate Concordance (CV, mean TPM >= 1) ===\n\n")
cat(sprintf("%-14s", "Method"))
for (s in samples) cat(sprintf(" | %8s", paste0("CV(", s, ")")))
cat(sprintf(" | %8s\n", "Mean"))
cat(paste0(rep("-", 14 + 4 * 11 + 11), collapse = ""), "\n")

for (mname in names(methods)) {
  meth <- methods[[mname]]
  cvs <- c()
  cat(sprintf("%-14s", mname))
  for (s in samples) {
    sm <- load_sample_means(meth, s)
    if (is.null(sm)) {
      cat(sprintf(" | %8s", "N/A"))
      next
    }
    expr <- sm[mean_tpm >= 1]
    med_cv <- median(expr$cv, na.rm = TRUE)
    cvs <- c(cvs, med_cv)
    cat(sprintf(" | %8.4f", med_cv))
  }
  if (length(cvs) > 0) {
    cat(sprintf(" | %8.4f", mean(cvs)))
  } else {
    cat(sprintf(" | %8s", "N/A"))
  }
  cat("\n")
}

# ============================================================
# Evaluation 2: Titration accuracy (C and D vs expected)
# ============================================================

cat("\n=== Titration Accuracy: Observed vs Expected ===\n")
cat("Expected C = 0.75*A + 0.25*B; Expected D = 0.25*A + 0.75*B\n")
cat("(Transcripts with mean A or B TPM >= 1)\n\n")

cat(sprintf("%-14s | %10s %10s | %10s %10s | %10s\n",
            "Method", "Pears(C)", "Pears(D)", "Spear(C)", "Spear(D)", "FC slope"))
cat(paste0(rep("-", 80), collapse = ""), "\n")

for (mname in names(methods)) {
  meth <- methods[[mname]]
  A <- load_sample_means(meth, "A")
  B <- load_sample_means(meth, "B")
  C <- load_sample_means(meth, "C")
  D <- load_sample_means(meth, "D")

  if (is.null(A) || is.null(B) || is.null(C) || is.null(D)) {
    cat(sprintf("%-14s | %10s %10s | %10s %10s | %10s\n", mname, "N/A", "N/A", "N/A", "N/A", "N/A"))
    next
  }

  # Match target names across all four samples
  all_names <- Reduce(intersect, list(A$target_name, B$target_name, C$target_name, D$target_name))
  setkey(A, target_name); setkey(B, target_name); setkey(C, target_name); setkey(D, target_name)
  a <- A[all_names]$mean_tpm
  b <- B[all_names]$mean_tpm
  c_obs <- C[all_names]$mean_tpm
  d_obs <- D[all_names]$mean_tpm

  # Expected values from titration
  c_exp <- 0.75 * a + 0.25 * b
  d_exp <- 0.25 * a + 0.75 * b

  # Filter to transcripts expressed in at least one reference
  expr <- (a >= 1) | (b >= 1)

  # Correlations
  pc <- cor(log2(c_obs[expr] + pseudo), log2(c_exp[expr] + pseudo), method = "pearson")
  pd <- cor(log2(d_obs[expr] + pseudo), log2(d_exp[expr] + pseudo), method = "pearson")
  sc <- cor(c_obs[expr], c_exp[expr], method = "spearman")
  sd <- cor(d_obs[expr], d_exp[expr], method = "spearman")

  # Fold change C/D: observed vs expected
  obs_lfc <- log2(c_obs[expr] + pseudo) - log2(d_obs[expr] + pseudo)
  exp_lfc <- log2(c_exp[expr] + pseudo) - log2(d_exp[expr] + pseudo)

  # Only transcripts with non-trivial expected FC
  big_fc <- abs(exp_lfc) > 0.5
  if (sum(big_fc) > 100) {
    fit <- lm(obs_lfc[big_fc] ~ exp_lfc[big_fc])
    slope <- sprintf("%10.4f", coef(fit)[2])
  } else {
    slope <- sprintf("%10s", "N/A")
  }

  cat(sprintf("%-14s | %10.4f %10.4f | %10.4f %10.4f | %s\n",
              mname, pc, pd, sc, sd, slope))
}

# ============================================================
# Evaluation 3: Number of expressed transcripts
# ============================================================

cat("\n=== Number of Expressed Transcripts (TPM > 0) ===\n\n")
cat(sprintf("%-14s", "Method"))
for (s in samples) cat(sprintf(" | %8s", s))
cat("\n")
cat(paste0(rep("-", 14 + 4 * 11), collapse = ""), "\n")

for (mname in names(methods)) {
  meth <- methods[[mname]]
  cat(sprintf("%-14s", mname))
  for (s in samples) {
    sm <- load_sample_means(meth, s)
    if (is.null(sm)) {
      cat(sprintf(" | %8s", "N/A"))
    } else {
      cat(sprintf(" | %8d", sum(sm$mean_tpm > 0)))
    }
  }
  cat("\n")
}

cat("\nDone.\n")
