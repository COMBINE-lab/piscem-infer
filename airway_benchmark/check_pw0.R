#!/usr/bin/env Rscript
library(data.table)
benchdir <- "airway_benchmark"
all_samples <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520",
                 "SRR1039509", "SRR1039513", "SRR1039517", "SRR1039521")

target <- "ENST00000449131"
cat(sprintf("=== %s (BEST1-202) per-sample TPM ===\n\n", target))
cat(sprintf("%-16s | %10s %10s %10s %10s | %10s %10s %10s %10s\n",
            "Method", "U:508", "U:512", "U:516", "U:520",
            "T:509", "T:513", "T:517", "T:521"))
cat(paste0(rep("-", 116), collapse = ""), "\n")

for (info in list(
  list(name = "Plain EM",       dir = "quant/em",                subdir = FALSE),
  list(name = "Hier pw=0.25",   dir = "quant/hier_all_cond",     subdir = TRUE),
  list(name = "Hier pw=0.00",   dir = "quant/hier_pw0", subdir = TRUE)
)) {
  vals <- numeric(8)
  for (si in seq_along(all_samples)) {
    s <- all_samples[si]
    if (info$subdir) {
      path <- file.path(benchdir, info$dir, s, paste0(s, ".quant"))
    } else {
      path <- file.path(benchdir, info$dir, paste0(s, ".quant"))
    }
    df <- fread(path, sep = "\t")
    setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
    idx <- grep(target, df$target_name, fixed = TRUE)
    vals[si] <- if (length(idx) > 0) df$tpm[idx[1]] else NA
  }
  cat(sprintf("%-16s | %10.2f %10.2f %10.2f %10.2f | %10.2f %10.2f %10.2f %10.2f\n",
              info$name, vals[1], vals[2], vals[3], vals[4],
              vals[5], vals[6], vals[7], vals[8]))
}
cat("\nDone.\n")
