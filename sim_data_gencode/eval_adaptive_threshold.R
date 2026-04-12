#!/usr/bin/env Rscript
# Prototype: adaptive EC support threshold based on transcript ambiguity.
#
# For each transcript, compute the average EC size of its contributing ECs.
# Transcripts in highly ambiguous regions (large average EC) require more
# EC support to be called expressed.
#
# We test this by post-hoc filtering the Sel+Cons Support results with
# different threshold strategies.

library(data.table)

simdir <- "sim_data_gencode"
pseudo <- 0.01

gt <- fread(file.path(simdir, "ground_truth.csv"))
gt[, short_id := sub(" .*", "", transcript_id)]
sample_info <- fread(file.path(simdir, "sample_info.csv"))

read_quant <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- fread(path, sep = "\t"); setnames(df, c("target_name","len","eelen","tpm","ecount")); df
}

# For the adaptive threshold, we need per-transcript ambiguity info.
# Since we don't have direct access to the EC graph from R, we'll use a proxy:
# the number of transcripts with nonzero TPM that have "similar" expression patterns.
#
# Actually, a simpler proxy available from the quant files: transcripts at complex
# loci (many isoforms per gene) need higher thresholds. We can parse gene names
# from the GENCODE transcript IDs.

parse_gene <- function(target_name) {
  # GENCODE format: ENST...|ENSG...|...|...|name-NNN|gene_name|len|...
  sapply(strsplit(target_name, "\\|"), function(x) if (length(x) >= 6) x[6] else NA)
}

# Load one sample to get gene structure
qdf <- read_quant(file.path(simdir, "quant_consensus_sel_support/sample_01/sample_01.quant"))
qdf$gene <- parse_gene(qdf$target_name)
gene_iso_count <- qdf[, .N, by = gene]
setnames(gene_iso_count, "N", "n_isoforms")
qdf <- merge(qdf, gene_iso_count, by = "gene")

# Define adaptive threshold based on gene complexity
# Simple rule: threshold = max(2, floor(log2(n_isoforms)))
qdf$adaptive_thresh <- pmax(2, floor(log2(qdf$n_isoforms)))

cat("Adaptive threshold distribution:\n")
cat(sprintf("  n_isoforms  1: threshold = %d\n", max(2, floor(log2(1)))))
cat(sprintf("  n_isoforms  2-3: threshold = %d\n", max(2, floor(log2(2)))))
cat(sprintf("  n_isoforms  4-7: threshold = %d\n", max(2, floor(log2(4)))))
cat(sprintf("  n_isoforms  8-15: threshold = %d\n", max(2, floor(log2(8)))))
cat(sprintf("  n_isoforms 16-31: threshold = %d\n", max(2, floor(log2(16)))))
cat(sprintf("  n_isoforms 32-63: threshold = %d\n", max(2, floor(log2(32)))))
cat(sprintf("  n_isoforms 64+: threshold = %d\n", max(2, floor(log2(64)))))

# Now test: post-hoc re-filter with adaptive threshold.
# We need the per-transcript EC support from Phase 1 to apply adaptive thresholds.
# Since we don't have that directly, we'll simulate by applying threshold to
# per-sample expressed counts.
#
# Alternative approach: just test fixed thresholds at different levels to see
# what the optimal global threshold would be.

cat("\n=== Fixed threshold sweep ===\n")
cat("(Using Sel+Cons Support results with different fixed thresholds via post-hoc TPM filtering)\n\n")

# Actually, let's approach this differently: compare Sel+Cons Support at different
# min-ec-support values by running the tool. But we can approximate by looking at
# how many FP would be eliminated at higher thresholds by examining the
# intersection of called transcripts with ground truth.

# For now, let's just look at the gene-complexity stratified TP/FP rates
# to understand WHERE the FP come from.

ctrl_samples  <- sample_info$sample_name[sample_info$condition == "control"]
treat_samples <- sample_info$sample_name[sample_info$condition == "treatment"]

cat("\n=== Gene-complexity stratified TP/FP (Sel+Cons Support, sample_01) ===\n\n")

m <- match(qdf$target_name, gt$short_id)
valid <- !is.na(m)
true_tpm <- gt$expected_tpm_control[m[valid]]
est_tpm <- qdf$tpm[valid]
n_iso <- qdf$n_isoforms[valid]

tiers <- list(
  "1 iso"     = c(1, 1),
  "2-5 iso"   = c(2, 5),
  "6-10 iso"  = c(6, 10),
  "11-20 iso" = c(11, 20),
  "21-50 iso" = c(21, 50),
  "51+ iso"   = c(51, 9999)
)

cat(sprintf("%-12s | %6s %6s %6s | %8s %8s\n",
            "Tier", "TP", "FP", "FN", "Prec", "Recall"))
cat(paste0(rep("-", 55), collapse=""), "\n")

for (tier_name in names(tiers)) {
  bounds <- tiers[[tier_name]]
  sel <- n_iso >= bounds[1] & n_iso <= bounds[2]
  tp <- sum(est_tpm[sel] > 0 & true_tpm[sel] > 0)
  fp <- sum(est_tpm[sel] > 0 & true_tpm[sel] == 0)
  fn <- sum(est_tpm[sel] == 0 & true_tpm[sel] > 0)
  prec <- if (tp + fp > 0) tp / (tp + fp) else NA
  rec <- if (tp + fn > 0) tp / (tp + fn) else NA
  cat(sprintf("%-12s | %6d %6d %6d | %8.4f %8.4f\n",
              tier_name, tp, fp, fn,
              ifelse(is.na(prec), 0, prec),
              ifelse(is.na(rec), 0, rec)))
}

cat("\n=== Conclusion ===\n")
cat("If FP concentrate at complex loci, adaptive thresholds would help.\n")
cat("If FP are spread evenly, a global threshold change is sufficient.\n")

cat("\nDone.\n")
