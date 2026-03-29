#!/usr/bin/env Rscript
# Evaluate group LASSO vs plain EM on airway replicate concordance.

library(data.table)

benchdir <- "airway_benchmark"
samples <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")

methods <- list(
  "Plain EM"     = "quant/em",
  "Hier. Multi"  = "quant/hier",
  "Hier. NoCond" = "quant/hier_nocond",
  "Group LASSO"  = "quant/gl"
)

load_method <- function(method_dir) {
  dfs <- list()
  for (s in samples) {
    path <- file.path(benchdir, method_dir, paste0(s, "/", s, ".quant"))
    if (!file.exists(path)) {
      # Try without sample subdir
      path <- file.path(benchdir, method_dir, paste0(s, ".quant"))
    }
    if (!file.exists(path)) return(NULL)
    df <- fread(path, sep = "\t")
    setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
    dfs[[s]] <- df[, .(target_name, tpm)]
    setnames(dfs[[s]], "tpm", s)
  }
  Reduce(function(a, b) merge(a, b, by = "target_name", all = TRUE), dfs)
}

compute_cv <- function(merged) {
  tpm_mat <- as.matrix(merged[, ..samples])
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

# ---- Global CV ----
cat("\n=== Replicate Concordance (CV of TPM, isoforms with mean TPM >= 1) ===\n\n")
cat(sprintf("%-14s | %8s %8s %8s | %8s | %8s\n",
            "Method", "Med CV", "Mean CV", "CV>1", "N expr", "N active"))
cat(paste0(rep("-", 70), collapse = ""), "\n")

for (name in names(all_results)) {
  m <- all_results[[name]]
  if (is.null(m)) { cat(sprintf("%-14s | MISSING\n", name)); next }
  cv_df <- compute_cv(m)
  expr <- cv_df[mean_tpm >= 1]
  tpm_mat <- as.matrix(m[, ..samples])
  n_active <- mean(apply(tpm_mat, 2, function(x) sum(x > 0)))
  cat(sprintf("%-14s | %8.4f %8.4f %8d | %8d | %8.0f\n",
              name, median(expr$cv), mean(expr$cv), sum(expr$cv > 1), nrow(expr), n_active))
}

# ---- TPM distribution ----
cat("\n=== TPM distribution (sample 1) ===\n")
for (name in names(all_results)) {
  m <- all_results[[name]]
  if (is.null(m)) next
  tpm <- m[[samples[1]]]
  cat(sprintf("  %-14s: >0=%.0f, >0.1=%.0f, >1=%.0f, >10=%.0f, total_tpm=%.0f\n",
              name, sum(tpm>0), sum(tpm>0.1), sum(tpm>1), sum(tpm>10), sum(tpm)))
}

# ---- Correlation between methods ----
cat("\n=== Method agreement (log2 TPM correlation, sample 1) ===\n")
if (length(all_results) >= 2) {
  m1 <- all_results[[1]]; m2 <- all_results[[2]]
  if (!is.null(m1) && !is.null(m2)) {
    both <- merge(m1[, .(target_name, em=get(samples[1]))],
                  m2[, .(target_name, gl=get(samples[1]))],
                  by = "target_name")
    pseudo <- 0.01
    expr <- both[em > 1 | gl > 1]
    r <- cor(log2(expr$em + pseudo), log2(expr$gl + pseudo))
    cat(sprintf("  Pearson(log2 TPM, expressed): %.4f (n=%d)\n", r, nrow(expr)))
  }
}

cat("\nDone.\n")
