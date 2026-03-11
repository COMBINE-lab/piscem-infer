#!/usr/bin/env Rscript
#
# Simulate multi-sample RNA-seq reads using polyester.
#
# Generates paired-end FASTA reads for a multi-condition experiment with known
# ground truth abundances, suitable for validating piscem-infer's hierarchical
# multi-sample quantification.
#
# Polyester only simulates reads for transcripts with >0 abundance, but the
# full transcriptome is used as the reference for mapping and quantification.
# This tests the quantifier's ability to correctly assign zero abundance to
# unexpressed transcripts in the presence of a large background.
#
# Requirements:
#   install.packages("BiocManager")
#   BiocManager::install("polyester")
#   BiocManager::install("Biostrings")
#
# Usage:
#   Rscript simulate_reads.R [output_dir] [transcriptome_fasta]
#
# If no transcriptome FASTA is provided, uses polyester's built-in chr22 data
# with a small subset.

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
total_reads <- 5000000  # total reads per sample (approx)

# Number of expressed transcripts (subset of total transcriptome)
num_expressed <- 5000

# ---- Load transcriptome ----

small_mode <- is.null(fasta_file)

if (small_mode) {
  # Use polyester's built-in chr22 transcripts (small test mode)
  fasta_file <- system.file("extdata", "chr22.fa", package = "polyester")
  cat("Using built-in chr22 transcriptome:", fasta_file, "\n")
  num_expressed <- 200
  total_reads <- 500000
}

fasta <- readDNAStringSet(fasta_file)
ntx_total <- length(fasta)
cat("Loaded", ntx_total, "transcripts from", basename(fasta_file), "\n")

tx_names <- names(fasta)
tx_lengths <- width(fasta)

# ---- Select expressed transcripts ----

# Clamp num_expressed to available transcripts
num_expressed <- min(num_expressed, ntx_total)

# Select expressed transcripts with a bias toward longer ones (more realistic)
# Use a weighted sampling: weight = sqrt(length) to favor longer transcripts
# but not exclusively (shorter transcripts can be expressed too)
if (ntx_total > num_expressed) {
  weights <- sqrt(pmax(tx_lengths, 1))
  expressed_idx <- sort(sample(1:ntx_total, num_expressed, prob = weights))
  cat("Selected", num_expressed, "expressed transcripts out of", ntx_total, "total\n")
} else {
  expressed_idx <- 1:ntx_total
  cat("All", ntx_total, "transcripts will be expressed\n")
}

# Boolean mask for expressed transcripts
is_expressed <- rep(FALSE, ntx_total)
is_expressed[expressed_idx] <- TRUE

# ---- Write full reference FASTA ----

ref_fasta <- file.path(outdir, "reference.fa")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
writeXStringSet(fasta, ref_fasta)
cat("Wrote full reference FASTA:", ref_fasta, "(", ntx_total, "transcripts )\n")

# Write expressed-only FASTA for polyester (it only handles expressed transcripts)
expressed_fasta <- fasta[expressed_idx]
expressed_fasta_file <- file.path(outdir, "expressed_transcripts.fa")
writeXStringSet(expressed_fasta, expressed_fasta_file)
cat("Wrote expressed FASTA:", expressed_fasta_file, "(", num_expressed, "transcripts )\n")

# ---- Define ground truth abundances ----

# Log-normal abundance distribution for expressed transcripts
log_base_abundance <- rnorm(num_expressed, mean = 3, sd = 1.5)
base_abundance <- exp(log_base_abundance)
base_abundance <- base_abundance / sum(base_abundance)  # normalize to proportions

# reads_per_transcript for expressed genes
reads_per_tx_expressed <- round(base_abundance * total_reads)
reads_per_tx_expressed[reads_per_tx_expressed < 1] <- 1

# Full vector: 0 for unexpressed, positive for expressed
reads_per_tx <- rep(0, ntx_total)
reads_per_tx[expressed_idx] <- reads_per_tx_expressed

cat("\nAbundance summary (expressed transcripts):\n")
cat("  Min reads:", min(reads_per_tx_expressed), "\n")
cat("  Median reads:", median(reads_per_tx_expressed), "\n")
cat("  Max reads:", max(reads_per_tx_expressed), "\n")
cat("  Total reads per sample (approx):", sum(reads_per_tx_expressed), "\n")

# ---- Define differential expression ----

# DE only among expressed transcripts
num_de <- round(num_expressed * 0.2)
de_expressed_pos <- sample(1:num_expressed, num_de)  # positions within expressed set
de_indices <- expressed_idx[de_expressed_pos]          # positions within full transcriptome

# Split: half up, half down in treatment
num_de_up <- round(num_de / 2)
de_up <- de_indices[1:num_de_up]
de_down <- de_indices[(num_de_up + 1):num_de]

# Fold changes
fc_up <- 2^runif(length(de_up), min = 1, max = 3)
fc_down <- 1 / (2^runif(length(de_down), min = 1, max = 3))

# Full fold change matrix: ntx_total rows x 2 columns
fold_changes <- matrix(1, nrow = ntx_total, ncol = num_conditions)
fold_changes[de_up, 2] <- fc_up
fold_changes[de_down, 2] <- fc_down

cat("\nDifferential expression setup:\n")
cat("  Total transcripts:", ntx_total, "\n")
cat("  Expressed transcripts:", num_expressed, "\n")
cat("  DE transcripts:", num_de, "(", num_de_up, "up,", num_de - num_de_up, "down )\n")
cat("  Fold change range (up):", round(range(fc_up), 2), "\n")
cat("  Fold change range (down):", round(range(fc_down), 2), "\n")

# ---- Save ground truth ----

ground_truth <- data.frame(
  transcript_id = tx_names,
  length = tx_lengths,
  is_expressed = is_expressed,
  base_reads = reads_per_tx,
  fc_control = fold_changes[, 1],
  fc_treatment = fold_changes[, 2],
  is_de = 1:ntx_total %in% de_indices,
  de_direction = ifelse(1:ntx_total %in% de_up, "up",
                        ifelse(1:ntx_total %in% de_down, "down", "none")),
  stringsAsFactors = FALSE
)

# Compute expected TPM per condition (for ALL transcripts)
for (cond_idx in 1:num_conditions) {
  cond_name <- c("control", "treatment")[cond_idx]
  expected_reads <- reads_per_tx * fold_changes[, cond_idx]
  eff_len <- pmax(tx_lengths - fraglen + 1, 1)
  rpk <- expected_reads / (eff_len / 1000)
  tpm <- rpk / sum(rpk) * 1e6
  ground_truth[[paste0("expected_reads_", cond_name)]] <- round(expected_reads, 1)
  ground_truth[[paste0("expected_tpm_", cond_name)]] <- round(tpm, 2)
}

gt_file <- file.path(outdir, "ground_truth.csv")
write.csv(ground_truth, gt_file, row.names = FALSE)
cat("Wrote ground truth:", gt_file, "\n")
cat("  (includes", sum(!is_expressed), "unexpressed transcripts with TPM=0)\n")

# ---- Simulate reads ----

# polyester simulates reads for all replicates of all groups at once.
# IMPORTANT: polyester only accepts transcripts with reads_per_transcript > 0.
# We pass only the expressed subset and its fold changes.

# Extract fold changes for expressed transcripts only
fc_expressed <- fold_changes[expressed_idx, , drop = FALSE]

sim_outdir <- file.path(outdir, "reads")
dir.create(sim_outdir, recursive = TRUE, showWarnings = FALSE)

cat("\nSimulating reads...\n")
cat("  Expressed transcripts:", num_expressed, "\n")
cat("  Replicates per condition:", num_reps, "\n")
cat("  Total samples:", num_reps * num_conditions, "\n")
cat("  Read length:", readlen, "\n")
cat("  Fragment length:", fraglen, "+/-", fragsd, "\n")
cat("  Error rate:", error_rate, "\n")

simulate_experiment(
  fasta = expressed_fasta_file,  # only expressed transcripts
  reads_per_transcript = reads_per_tx_expressed,
  num_reps = c(num_reps, num_reps),
  fold_changes = fc_expressed,
  paired = TRUE,
  readlen = readlen,
  fraglen = fraglen,
  fragsd = fragsd,
  error_rate = error_rate,
  outdir = sim_outdir,
  seed = 42
)

cat("\nSimulation complete.\n")

# ---- Create manifest for piscem-infer ----

manifest_rows <- list()
total_samples <- num_reps * num_conditions

for (i in 1:total_samples) {
  sample_id <- sprintf("sample_%02d", i)
  condition <- if (i <= num_reps) "control" else "treatment"
  r1 <- file.path(sim_outdir, paste0(sample_id, "_1.fasta"))
  r2 <- file.path(sim_outdir, paste0(sample_id, "_2.fasta"))

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
