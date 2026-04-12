#!/usr/bin/env Rscript

library(data.table)

simdir <- "sim_data_gencode"
mask_dir <- file.path(simdir, "oracle_masks")
dir.create(mask_dir, showWarnings = FALSE, recursive = TRUE)

gt <- fread(file.path(simdir, "ground_truth.csv"))
gt[, short_id := sub(" .*", "", transcript_id)]
sample_info <- fread(file.path(simdir, "sample_info.csv"))

for (i in seq_len(nrow(sample_info))) {
  sn <- sample_info$sample_name[i]
  cond <- sample_info$condition[i]
  tpm_col <- paste0("expected_tpm_", cond)
  keep <- gt$short_id[gt[[tpm_col]] > 0]
  writeLines(keep, file.path(mask_dir, paste0(sn, ".txt")))
  cat(sprintf("%s: wrote %d transcript IDs\n", sn, length(keep)))
}
