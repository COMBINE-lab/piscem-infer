#!/usr/bin/env Rscript
# Evaluate adaptive variance vs other methods on 8 airway samples.
# Same metrics as evaluate_all8.R but adds NoCond+AV method.

library(data.table)

benchdir <- "airway_benchmark"

samples_untreated <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")
samples_treated   <- c("SRR1039509", "SRR1039513", "SRR1039517", "SRR1039521")
all_samples <- c(samples_untreated, samples_treated)

methods <- list(
  "Plain EM"       = list(dir = "quant/em",                  subdir = FALSE),
  "Hier. Cond"     = list(dir = "quant/hier_all_cond",       subdir = TRUE),
  "Hier. NoCond"   = list(dir = "quant/hier_all_nocond",     subdir = TRUE),
  "NoCond+AV"      = list(dir = "quant/hier_all_nocond_av",  subdir = TRUE)
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
  cat(sprintf("%-14s | %8s %8s %8s | %8s\n",
              "Method", "Med CV", "Mean CV", "CV>1", "N expr"))
  cat(paste0(rep("-", 58), collapse = ""), "\n")

  for (name in names(all_results)) {
    m <- all_results[[name]]
    if (is.null(m)) { cat(sprintf("%-14s | MISSING\n", name)); next }
    cv_df <- compute_cv(m, subs)
    expr <- cv_df[mean_tpm >= 1]
    cat(sprintf("%-14s | %8.4f %8.4f %8d | %8d\n",
                name, median(expr$cv), mean(expr$cv), sum(expr$cv > 1), nrow(expr)))
  }
}

# ---- Cross-condition fold change ----
cat("\n=== Between-condition fold change (mean TPM >= 1 in both conditions) ===\n\n")
cat(sprintf("%-14s | %10s %10s %10s %10s\n",
            "Method", "|log2FC|", "SD(log2FC)", "Atten.Slope", "N both"))
cat(paste0(rep("-", 65), collapse = ""), "\n")

# Use plain EM as reference for fold changes
ref <- all_results[["Plain EM"]]
if (!is.null(ref)) {
  pseudo <- 0.01
  ref_untrt <- rowMeans(as.matrix(ref[, ..samples_untreated]))
  ref_trt   <- rowMeans(as.matrix(ref[, ..samples_treated]))
  ref_both  <- ref_untrt >= 1 & ref_trt >= 1
  ref_lfc   <- log2(ref_trt[ref_both] + pseudo) - log2(ref_untrt[ref_both] + pseudo)
}

for (name in names(all_results)) {
  m <- all_results[[name]]
  if (is.null(m)) next

  untrt_mean <- rowMeans(as.matrix(m[, ..samples_untreated]))
  trt_mean   <- rowMeans(as.matrix(m[, ..samples_treated]))
  both_expr <- untrt_mean >= 1 & trt_mean >= 1

  if (sum(both_expr) > 0) {
    lfc <- log2(trt_mean[both_expr] + pseudo) - log2(untrt_mean[both_expr] + pseudo)
    mean_abs_lfc <- mean(abs(lfc))
    sd_lfc <- sd(lfc)

    # Attenuation slope vs plain EM
    slope_str <- ""
    if (!is.null(ref) && name != "Plain EM") {
      # Match transcripts between methods
      m_names <- m$target_name[both_expr]
      ref_names <- ref$target_name[ref_both]
      common <- intersect(m_names, ref_names)
      if (length(common) > 100) {
        m_lfc <- lfc[match(common, m_names)]
        r_lfc <- ref_lfc[match(common, ref_names)]
        fit <- lm(m_lfc ~ r_lfc)
        slope_str <- sprintf("%.3f", coef(fit)[2])
      }
    }

    cat(sprintf("%-14s | %10.4f %10.4f %10s %10d\n",
                name, mean_abs_lfc, sd_lfc, slope_str, sum(both_expr)))
  }
}

# ---- Correlation with plain EM ----
cat("\n=== Spearman correlation of TPM with Plain EM (per sample) ===\n\n")
em_data <- all_results[["Plain EM"]]
if (!is.null(em_data)) {
  for (name in names(all_results)) {
    if (name == "Plain EM") next
    m <- all_results[[name]]
    if (is.null(m)) next
    cors <- c()
    for (s in all_samples) {
      # Match transcripts
      common <- intersect(em_data$target_name, m$target_name)
      em_tpm <- em_data[[s]][match(common, em_data$target_name)]
      m_tpm  <- m[[s]][match(common, m$target_name)]
      cors <- c(cors, cor(em_tpm, m_tpm, method = "spearman"))
    }
    cat(sprintf("  %-14s: mean Spearman = %.4f (range: %.4f - %.4f)\n",
                name, mean(cors), min(cors), max(cors)))
  }
}

cat("\nDone.\n")
