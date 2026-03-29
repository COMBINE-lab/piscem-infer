#!/usr/bin/env Rscript
# Deep dive on ENST00000449131 and sigma estimates

library(data.table)

benchdir <- "airway_benchmark"
samples_untreated <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")
samples_treated   <- c("SRR1039509", "SRR1039513", "SRR1039517", "SRR1039521")
all_samples <- c(samples_untreated, samples_treated)

read_one <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- fread(path, sep = "\t")
  setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
  df
}

target_pattern <- "ENST00000449131"
suspects <- c("ENST00000449131", "ENST00000639948", "ENST00000853539",
              "ENST00000547798", "ENST00000912395", "ENST00000282470")

cat("=== Per-sample TPM for suspect transcripts ===\n\n")

for (target in suspects) {
  cat(sprintf("--- %s ---\n", target))
  cat(sprintf("%-14s | %10s %10s %10s %10s | %10s %10s %10s %10s\n",
              "Method",
              "U:508", "U:512", "U:516", "U:520",
              "T:509", "T:513", "T:517", "T:521"))
  cat(paste0(rep("-", 110), collapse = ""), "\n")

  for (info in list(
    list(name = "Plain EM", dir = "quant/em", subdir = FALSE),
    list(name = "Hier. Cond", dir = "quant/hier_all_cond", subdir = TRUE),
    list(name = "Hier. NoCond", dir = "quant/hier_all_nocond", subdir = TRUE)
  )) {
    vals <- numeric(8)
    for (si in seq_along(all_samples)) {
      s <- all_samples[si]
      if (info$subdir) {
        path <- file.path(benchdir, info$dir, s, paste0(s, ".quant"))
      } else {
        path <- file.path(benchdir, info$dir, paste0(s, ".quant"))
      }
      df <- read_one(path)
      if (is.null(df)) { vals[si] <- NA; next }
      idx <- grep(target, df$target_name, fixed = TRUE)
      if (length(idx) == 0) { vals[si] <- NA; next }
      vals[si] <- df$tpm[idx[1]]
    }
    cat(sprintf("%-14s | %10.2f %10.2f %10.2f %10.2f | %10.2f %10.2f %10.2f %10.2f\n",
                info$name, vals[1], vals[2], vals[3], vals[4],
                vals[5], vals[6], vals[7], vals[8]))
  }
  cat("\n")
}

# Also show estimated counts
cat("\n=== Estimated counts for ENST00000449131 ===\n\n")
cat(sprintf("%-14s | %10s %10s %10s %10s | %10s %10s %10s %10s\n",
            "Method",
            "U:508", "U:512", "U:516", "U:520",
            "T:509", "T:513", "T:517", "T:521"))
cat(paste0(rep("-", 110), collapse = ""), "\n")

for (info in list(
  list(name = "Plain EM", dir = "quant/em", subdir = FALSE),
  list(name = "Hier. Cond", dir = "quant/hier_all_cond", subdir = TRUE),
  list(name = "Hier. NoCond", dir = "quant/hier_all_nocond", subdir = TRUE)
)) {
  vals <- numeric(8)
  for (si in seq_along(all_samples)) {
    s <- all_samples[si]
    if (info$subdir) {
      path <- file.path(benchdir, info$dir, s, paste0(s, ".quant"))
    } else {
      path <- file.path(benchdir, info$dir, paste0(s, ".quant"))
    }
    df <- read_one(path)
    if (is.null(df)) { vals[si] <- NA; next }
    idx <- grep("ENST00000449131", df$target_name, fixed = TRUE)
    if (length(idx) == 0) { vals[si] <- NA; next }
    vals[si] <- df$ecount[idx[1]]
  }
  cat(sprintf("%-14s | %10.1f %10.1f %10.1f %10.1f | %10.1f %10.1f %10.1f %10.1f\n",
              info$name, vals[1], vals[2], vals[3], vals[4],
              vals[5], vals[6], vals[7], vals[8]))
}

cat("\nDone.\n")
