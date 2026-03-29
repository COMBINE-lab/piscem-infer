#!/usr/bin/env Rscript
# Evaluate within-condition replicate concordance across all 8 airway samples.

library(data.table)

benchdir <- "airway_benchmark"

# 4 untreated + 4 dexamethasone-treated
samples_untreated <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")
samples_treated   <- c("SRR1039509", "SRR1039513", "SRR1039517", "SRR1039521")
all_samples <- c(samples_untreated, samples_treated)

methods <- list(
  "Plain EM"       = list(dir = "quant/em",              subdir = FALSE),
  "NullSink .50e3" = list(dir = "quant/nullsink_s050_ec3", subdir = TRUE),
  "Hier. Cond"     = list(dir = "quant/hier_all_cond",   subdir = TRUE),
  "Hier. NoCond"   = list(dir = "quant/hier_all_nocond", subdir = TRUE)
)

load_method <- function(method_info) {
  dfs <- list()
  for (s in all_samples) {
    if (method_info$subdir) {
      path <- file.path(benchdir, method_info$dir, s, paste0(s, ".quant"))
    } else {
      path <- file.path(benchdir, method_info$dir, paste0(s, ".quant"))
    }
    if (!file.exists(path)) return(NULL)
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

cat("Loading...\n")
all_results <- list()
for (name in names(methods)) {
  all_results[[name]] <- load_method(methods[[name]])
}

# ---- Within-condition CV ----
for (cond_name in c("Untreated", "Treated")) {
  subs <- if (cond_name == "Untreated") samples_untreated else samples_treated
  cat(sprintf("\n=== %s Replicate Concordance (CV of TPM, mean TPM >= 1) ===\n\n", cond_name))
  cat(sprintf("%-14s | %8s %8s %8s | %8s | %8s\n",
              "Method", "Med CV", "Mean CV", "CV>1", "N expr", "N active"))
  cat(paste0(rep("-", 70), collapse = ""), "\n")

  for (name in names(all_results)) {
    m <- all_results[[name]]
    if (is.null(m)) { cat(sprintf("%-14s | MISSING\n", name)); next }
    cv_df <- compute_cv(m, subs)
    expr <- cv_df[mean_tpm >= 1]
    tpm_mat <- as.matrix(m[, ..subs])
    n_active <- mean(apply(tpm_mat, 2, function(x) sum(x > 0)))
    cat(sprintf("%-14s | %8.4f %8.4f %8d | %8d | %8.0f\n",
                name, median(expr$cv), mean(expr$cv), sum(expr$cv > 1), nrow(expr), n_active))
  }
}

# ---- Cross-condition: between-group differences ----
cat("\n=== Between-condition difference (mean |log2FC|, isoforms with mean TPM >= 1 in both) ===\n\n")
cat(sprintf("%-14s | %10s %10s %10s\n", "Method", "|log2FC|", "Cor(logFC)", "N both"))
cat(paste0(rep("-", 55), collapse = ""), "\n")

for (name in names(all_results)) {
  m <- all_results[[name]]
  if (is.null(m)) next
  pseudo <- 0.01
  untrt_mean <- rowMeans(as.matrix(m[, ..samples_untreated]))
  trt_mean   <- rowMeans(as.matrix(m[, ..samples_treated]))
  both_expr <- untrt_mean >= 1 & trt_mean >= 1
  if (sum(both_expr) > 0) {
    lfc <- log2(trt_mean[both_expr] + pseudo) - log2(untrt_mean[both_expr] + pseudo)
    cat(sprintf("%-14s | %10.4f %10s %10d\n",
                name, mean(abs(lfc)), "", sum(both_expr)))
  }
}

# ---- TPM distribution ----
cat("\n=== TPM distribution (first untreated sample) ===\n")
for (name in names(all_results)) {
  m <- all_results[[name]]
  if (is.null(m)) next
  tpm <- m[[samples_untreated[1]]]
  cat(sprintf("  %-14s: >0=%.0f, >0.1=%.0f, >1=%.0f, >10=%.0f\n",
              name, sum(tpm>0), sum(tpm>0.1), sum(tpm>1), sum(tpm>10)))
}

cat("\nDone.\n")
