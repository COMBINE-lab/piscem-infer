#!/usr/bin/env Rscript
#
# Simulate multi-sample RNA-seq reads using polyester.
#
# This script generates paired-end FASTA reads for a multi-condition
# experiment with known ground truth abundances, suitable for validating
# piscem-infer's hierarchical multi-sample quantification.
#
# Requirements:
#   install.packages("BiocManager")
#   BiocManager::install("polyester")
#   BiocManager::install("Biostrings")
#
# Usage:
#   Rscript simulate_reads.R [output_dir] [transcriptome_fasta]
#
# If no transcriptome FASTA is provided, uses polyester's built-in chr22 data.

suppressPackageStartupMessages({
  library(polyester)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "sim_data"
fasta_file <- if (length(args) >= 2) args[2] else NULL

set.seed(42)

# ---- Configuration ----

num_reps <- 3          # replicates per condition
num_conditions <- 2     # control vs treatment
readlen <- 100
fraglen <- 250
fragsd <- 25
error_rate <- 0.005

# Number of transcripts to simulate (subset if using built-in data)
num_tx <- 200

# ---- Load transcriptome ----

if (is.null(fasta_file)) {
  # Use polyester's built-in chr22 transcripts
  fasta_file <- system.file("extdata", "chr22.fa", package = "polyester")
  cat("Using built-in chr22 transcriptome:", fasta_file, "\n")
}

fasta <- readDNAStringSet(fasta_file)
cat("Loaded", length(fasta), "transcripts\n")

# Subset to num_tx transcripts (pick the longest ones for more realistic mapping)
if (length(fasta) > num_tx) {
  tx_lengths <- width(fasta)
  keep_idx <- order(tx_lengths, decreasing = TRUE)[1:num_tx]
  fasta <- fasta[keep_idx]
  cat("Subsetted to", num_tx, "longest transcripts\n")
}

tx_names <- names(fasta)
tx_lengths <- width(fasta)
ntx <- length(fasta)

# Write the subsetted FASTA for use with piscem
ref_fasta <- file.path(outdir, "reference.fa")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
writeXStringSet(fasta, ref_fasta)
cat("Wrote reference FASTA:", ref_fasta, "\n")

# ---- Define ground truth abundances ----

# Base TPM-like abundances (log-normal distributed)
log_base_abundance <- rnorm(ntx, mean = 3, sd = 1.5)
base_abundance <- exp(log_base_abundance)
base_abundance <- base_abundance / sum(base_abundance)  # normalize to proportions

# reads_per_transcript: baseline read count (will be multiplied by fold changes)
total_reads <- 500000  # total reads per sample (approx)
reads_per_tx <- round(base_abundance * total_reads)
reads_per_tx[reads_per_tx < 1] <- 1  # ensure at least 1 read per transcript

# ---- Define differential expression ----

# Select DE transcripts: ~20% of transcripts will be DE
num_de <- round(ntx * 0.2)
de_indices <- sample(1:ntx, num_de)

# Split DE transcripts: half up in treatment, half down
num_de_up <- round(num_de / 2)
de_up <- de_indices[1:num_de_up]          # upregulated in treatment
de_down <- de_indices[(num_de_up+1):num_de]  # downregulated in treatment

# Fold changes: log2FC between 1 and 3 (so actual FC between 2x and 8x)
fc_up <- 2^runif(length(de_up), min = 1, max = 3)
fc_down <- 1 / (2^runif(length(de_down), min = 1, max = 3))

# Build fold change matrix: ntx rows x 2 columns (control, treatment)
fold_changes <- matrix(1, nrow = ntx, ncol = num_conditions)
fold_changes[de_up, 2] <- fc_up       # treatment column
fold_changes[de_down, 2] <- fc_down   # treatment column

cat("Differential expression setup:\n")
cat("  Total transcripts:", ntx, "\n")
cat("  DE transcripts:", num_de, "(", num_de_up, "up,", num_de - num_de_up, "down in treatment)\n")
cat("  Fold change range (up):", round(range(fc_up), 2), "\n")
cat("  Fold change range (down):", round(range(fc_down), 2), "\n")

# ---- Save ground truth ----

ground_truth <- data.frame(
  transcript_id = tx_names,
  length = tx_lengths,
  base_reads = reads_per_tx,
  fc_control = fold_changes[, 1],
  fc_treatment = fold_changes[, 2],
  is_de = 1:ntx %in% de_indices,
  de_direction = ifelse(1:ntx %in% de_up, "up",
                        ifelse(1:ntx %in% de_down, "down", "none")),
  stringsAsFactors = FALSE
)

# Compute expected TPM per condition
for (cond_idx in 1:num_conditions) {
  cond_name <- c("control", "treatment")[cond_idx]
  expected_reads <- reads_per_tx * fold_changes[, cond_idx]
  eff_len <- pmax(tx_lengths - fraglen + 1, 1)  # approximate effective length
  rpk <- expected_reads / (eff_len / 1000)
  tpm <- rpk / sum(rpk) * 1e6
  ground_truth[[paste0("expected_reads_", cond_name)]] <- round(expected_reads, 1)
  ground_truth[[paste0("expected_tpm_", cond_name)]] <- round(tpm, 2)
}

gt_file <- file.path(outdir, "ground_truth.csv")
write.csv(ground_truth, gt_file, row.names = FALSE)
cat("Wrote ground truth:", gt_file, "\n")

# ---- Simulate reads ----

# polyester simulates reads for all replicates of all groups at once.
# It creates files: sample_01_1.fasta, sample_01_2.fasta, ...
# Samples 1..num_reps are group 1 (control), num_reps+1..2*num_reps are group 2 (treatment)

sim_outdir <- file.path(outdir, "reads")
dir.create(sim_outdir, recursive = TRUE, showWarnings = FALSE)

cat("\nSimulating reads...\n")
cat("  Replicates per condition:", num_reps, "\n")
cat("  Total samples:", num_reps * num_conditions, "\n")
cat("  Read length:", readlen, "\n")
cat("  Fragment length:", fraglen, "+/-", fragsd, "\n")
cat("  Error rate:", error_rate, "\n")

simulate_experiment(
  fasta = ref_fasta,
  reads_per_transcript = reads_per_tx,
  num_reps = c(num_reps, num_reps),  # reps per group
  fold_changes = fold_changes,
  paired = TRUE,
  readlen = readlen,
  fraglen = fraglen,
  fragsd = fragsd,
  error_rate = error_rate,
  outdir = sim_outdir,
  write_info = TRUE,
  seed = 42
)

cat("\nSimulation complete.\n")

# ---- Create manifest for piscem-infer ----

# polyester names files: sample_01_1.fasta, sample_01_2.fasta, etc.
# Samples 1..num_reps = control, (num_reps+1)..2*num_reps = treatment

manifest_rows <- list()
total_samples <- num_reps * num_conditions

for (i in 1:total_samples) {
  sample_id <- sprintf("sample_%02d", i)
  condition <- if (i <= num_reps) "control" else "treatment"
  r1 <- file.path(sim_outdir, paste0(sample_id, "_1.fasta"))
  r2 <- file.path(sim_outdir, paste0(sample_id, "_2.fasta"))
  sample_out <- file.path(outdir, "quant", sample_id)

  if (!file.exists(r1)) {
    cat("WARNING: Expected file not found:", r1, "\n")
    next
  }

  # Count reads in the FASTA file
  n_reads <- length(readDNAStringSet(r1))
  cat(sprintf("  %s (%s): %d read pairs\n", sample_id, condition, n_reads))

  manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
    sample_name = sample_id,
    condition = condition,
    r1_path = normalizePath(r1, mustWork = FALSE),
    r2_path = normalizePath(r2, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

manifest_df <- do.call(rbind, manifest_rows)
manifest_file <- file.path(outdir, "sample_info.csv")
write.csv(manifest_df, manifest_file, row.names = FALSE)
cat("\nWrote sample info:", manifest_file, "\n")

# Also create the piscem-infer multi-quant manifest
# (rad_path and output_dir will be filled in by the pipeline script)
pq_manifest <- data.frame(
  sample_name = manifest_df$sample_name,
  condition = manifest_df$condition,
  rad_path = file.path(outdir, "mapped", manifest_df$sample_name, manifest_df$sample_name),
  output_dir = file.path(outdir, "quant", manifest_df$sample_name),
  stringsAsFactors = FALSE
)
pq_manifest_file <- file.path(outdir, "manifest.csv")
write.csv(pq_manifest, pq_manifest_file, row.names = FALSE)
cat("Wrote piscem-infer manifest:", pq_manifest_file, "\n")

cat("\n=== Done ===\n")
cat("Output directory:", normalizePath(outdir), "\n")
cat("Next steps:\n")
cat("  1. Build piscem index: piscem build -s", ref_fasta, "-k 31 -m 19 -t 8 -o", file.path(outdir, "index/ref"), "\n")
cat("  2. Run the pipeline: bash simulation/run_pipeline.sh", outdir, "\n")
