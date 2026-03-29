#!/usr/bin/env Rscript
# Investigate: does the no-condition model properly estimate high variance
# for DE transcripts and shrink them less?

library(data.table)

benchdir <- "airway_benchmark"
samples_untreated <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")
samples_treated   <- c("SRR1039509", "SRR1039513", "SRR1039517", "SRR1039521")
all_samples <- c(samples_untreated, samples_treated)

load_quant <- function(base_dir, subdir = TRUE) {
  dfs <- list()
  for (s in all_samples) {
    if (subdir) {
      path <- file.path(benchdir, base_dir, s, paste0(s, ".quant"))
    } else {
      path <- file.path(benchdir, base_dir, paste0(s, ".quant"))
    }
    if (!file.exists(path)) return(NULL)
    df <- fread(path, sep = "\t")
    setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
    dfs[[s]] <- df[, .(target_name, tpm)]
    setnames(dfs[[s]], "tpm", s)
  }
  Reduce(function(a, b) merge(a, b, by = "target_name", all = TRUE), dfs)
}

em   <- load_quant("quant/em", subdir = FALSE)
cond <- load_quant("quant/hier_all_cond", subdir = TRUE)
nocond <- load_quant("quant/hier_all_nocond", subdir = TRUE)

pseudo <- 0.01

# Compute fold changes from EM (ground truth proxy)
em_untrt <- rowMeans(as.matrix(em[, ..samples_untreated]))
em_trt   <- rowMeans(as.matrix(em[, ..samples_treated]))
em_lfc   <- log2(em_trt + pseudo) - log2(em_untrt + pseudo)

# Focus on transcripts with clear DE signal in EM (|log2FC| > 1 and mean > 1)
both_expr <- em_untrt >= 1 & em_trt >= 1
big_de <- both_expr & abs(em_lfc) > 1

cat(sprintf("Transcripts with |log2FC| > 1 and both conditions expressed: %d\n\n", sum(big_de)))

# For these DE transcripts, compare fold changes across methods
cond_untrt <- rowMeans(as.matrix(cond[, ..samples_untreated]))
cond_trt   <- rowMeans(as.matrix(cond[, ..samples_treated]))
cond_lfc   <- log2(cond_trt + pseudo) - log2(cond_untrt + pseudo)

nocond_untrt <- rowMeans(as.matrix(nocond[, ..samples_untreated]))
nocond_trt   <- rowMeans(as.matrix(nocond[, ..samples_treated]))
nocond_lfc   <- log2(nocond_trt + pseudo) - log2(nocond_untrt + pseudo)

cat("=== Fold change preservation for DE transcripts (|EM log2FC| > 1) ===\n\n")
cat(sprintf("%-14s | %10s %10s %10s %10s\n",
            "Method", "Mean|LFC|", "Med|LFC|", "Cor w/ EM", "Attenuation"))
cat(paste0(rep("-", 65), collapse = ""), "\n")

for (info in list(
  list(name = "Plain EM",     lfc = em_lfc),
  list(name = "Hier. Cond",   lfc = cond_lfc),
  list(name = "Hier. NoCond", lfc = nocond_lfc)
)) {
  de_lfc <- info$lfc[big_de]
  em_de_lfc <- em_lfc[big_de]
  r <- cor(de_lfc, em_de_lfc)
  # Attenuation: regression slope of method LFC on EM LFC
  slope <- coef(lm(de_lfc ~ em_de_lfc))[2]
  cat(sprintf("%-14s | %10.4f %10.4f %10.4f %10.4f\n",
              info$name, mean(abs(de_lfc)), median(abs(de_lfc)), r, slope))
}

# Look at variance of TPM across all 8 samples for DE vs non-DE
cat("\n=== Cross-sample variance (all 8 samples, TPM > 1 in at least one condition) ===\n\n")
em_all <- as.matrix(em[, ..all_samples])
em_var <- apply(em_all, 1, var)
cat(sprintf("  DE transcripts (n=%d): mean var = %.1f, median var = %.1f\n",
            sum(big_de), mean(em_var[big_de]), median(em_var[big_de])))
non_de <- both_expr & abs(em_lfc) < 0.5
cat(sprintf("  Non-DE (|LFC|<0.5, n=%d): mean var = %.1f, median var = %.1f\n",
            sum(non_de), mean(em_var[non_de]), median(em_var[non_de])))

# Look at specific examples
cat("\n=== Example DE transcripts ===\n")
cat(sprintf("%-30s | %8s %8s | %8s %8s | %8s %8s | %8s\n",
            "Transcript", "EM_ctrl", "EM_trt", "Cond_c", "Cond_t", "NoC_c", "NoC_t", "EM_LFC"))
cat(paste0(rep("-", 120), collapse = ""), "\n")

# Pick top 10 DE transcripts by |log2FC|
de_idx <- which(big_de)
de_ord <- de_idx[order(abs(em_lfc[de_idx]), decreasing = TRUE)]
for (i in head(de_ord, 15)) {
  cat(sprintf("%-30s | %8.1f %8.1f | %8.1f %8.1f | %8.1f %8.1f | %8.2f\n",
              substr(em$target_name[i], 1, 30),
              em_untrt[i], em_trt[i],
              cond_untrt[i], cond_trt[i],
              nocond_untrt[i], nocond_trt[i],
              em_lfc[i]))
}

# Scatter: nocond LFC vs EM LFC for DE transcripts
cat("\n=== Distribution of attenuation ratio (nocond LFC / EM LFC) for DE transcripts ===\n")
ratio <- nocond_lfc[big_de] / em_lfc[big_de]
ratio <- ratio[is.finite(ratio)]
cat(sprintf("  Quantiles: 5%%=%.3f 25%%=%.3f 50%%=%.3f 75%%=%.3f 95%%=%.3f\n",
            quantile(ratio, 0.05), quantile(ratio, 0.25), quantile(ratio, 0.50),
            quantile(ratio, 0.75), quantile(ratio, 0.95)))

cat("\nDone.\n")
