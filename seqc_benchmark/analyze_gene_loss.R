#!/usr/bin/env Rscript
# Analyze which TaqMan genes are lost by consensus filtering and why.

library(data.table)

benchdir <- "seqc_benchmark"
pseudo <- 0.01

# Load TaqMan
annot <- fread(cmd = "grep -v '^[#!^]' seqc_benchmark/GPL4097.annot", sep = "\t", header = TRUE,
               select = c("ID", "Gene symbol"))
setnames(annot, c("taqman_id", "gene_symbol"))
annot <- annot[gene_symbol != "" & gene_symbol != "---"]

expr_lines <- readLines("seqc_benchmark/taqman_raw.txt")
data_start <- grep('^"ID_REF"', expr_lines)
expr_dt <- fread(text = expr_lines[data_start:length(expr_lines)], sep = "\t", header = TRUE)
setnames(expr_dt, 1, "taqman_id")
sample_labels <- c(paste0("A_", 1:4), paste0("B_", 1:4), paste0("C_", 1:4), paste0("D_", 1:4))
setnames(expr_dt, names(expr_dt)[-1], sample_labels)
expr_dt[, taq_A := rowMeans(.SD), .SDcols = paste0("A_", 1:4)]
expr_dt[, taq_B := rowMeans(.SD), .SDcols = paste0("B_", 1:4)]
taqman <- merge(expr_dt[, .(taqman_id, taq_A, taq_B)], annot, by = "taqman_id")
taqman <- taqman[!duplicated(gene_symbol)]

# Parse gene from GENCODE transcript names
parse_gene <- function(target_name) {
  sapply(strsplit(target_name, "\\|"), function(x) if (length(x) >= 6) x[6] else NA)
}

# Load gene-level TPM for different methods
gene_tpm <- function(dir, sample_key, reader_type = "piscem") {
  if (reader_type == "piscem") {
    path <- file.path(benchdir, dir, paste0(sample_key, ".quant"))
    if (!file.exists(path)) path <- file.path(benchdir, dir, sample_key, paste0(sample_key, ".quant"))
    if (!file.exists(path)) return(NULL)
    df <- fread(path, sep = "\t")
    setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
  } else {
    path <- file.path(benchdir, dir, sample_key, "quant.sf")
    if (!file.exists(path)) return(NULL)
    df <- fread(path, sep = "\t")
    df <- df[, .(target_name = Name, tpm = TPM)]
  }
  df[, gene := parse_gene(target_name)]
  df[!is.na(gene), .(gene_tpm = sum(tpm)), by = gene]
}

# Mean gene TPM across 4 replicates for sample A
mean_gene_tpm <- function(dir, sample_type, reader_type = "piscem") {
  dfs <- list()
  for (r in 1:4) {
    key <- paste0(sample_type, "_", r)
    df <- gene_tpm(dir, key, reader_type)
    if (!is.null(df)) {
      setnames(df, "gene_tpm", paste0("r", r))
      dfs[[r]] <- df
    }
  }
  if (length(dfs) == 0) return(NULL)
  merged <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), dfs)
  tpm_cols <- setdiff(names(merged), "gene")
  merged$mean_tpm <- rowMeans(as.matrix(merged[, ..tpm_cols]), na.rm = TRUE)
  merged[, .(gene, mean_tpm)]
}

# Load sample A for each method
em_a <- mean_gene_tpm("quant_em", "A")
sel_a <- mean_gene_tpm("quant_sel_support", "A")
adapt_a <- mean_gene_tpm("quant_sel_adapt", "A")
salmon_a <- mean_gene_tpm("quant_salmon", "A", "salmon")

# Match to TaqMan genes detected in A
taq_detected <- taqman[taq_A > 0.001]
cat(sprintf("TaqMan genes detected in A: %d\n", nrow(taq_detected)))

# For each method, check which TaqMan genes are detected (gene TPM > 0)
check_detection <- function(gene_df, method_name) {
  if (is.null(gene_df)) return(NULL)
  matched <- merge(taq_detected, gene_df, by.x = "gene_symbol", by.y = "gene", all.x = TRUE)
  matched[is.na(mean_tpm), mean_tpm := 0]
  detected <- matched[mean_tpm > 0]
  missed <- matched[mean_tpm == 0]
  cat(sprintf("\n%s: %d/%d TaqMan genes detected (%.1f%%), %d missed\n",
              method_name, nrow(detected), nrow(taq_detected),
              100 * nrow(detected) / nrow(taq_detected), nrow(missed)))

  # Characterize missed genes by TaqMan expression level
  if (nrow(missed) > 0) {
    cat("  Missed genes by TaqMan expression tier:\n")
    tiers <- list("taq > 0.1" = 0.1, "taq > 0.01" = 0.01, "taq > 0.001" = 0.001)
    for (tier_name in names(tiers)) {
      thresh <- tiers[[tier_name]]
      n_tier <- sum(missed$taq_A > thresh)
      cat(sprintf("    %s: %d missed\n", tier_name, n_tier))
    }
  }
  missed
}

em_missed <- check_detection(em_a, "Plain EM")
sel_missed <- check_detection(sel_a, "Sel+Support")
adapt_missed <- check_detection(adapt_a, "Sel+Adaptive")
salmon_missed <- check_detection(salmon_a, "Salmon")

# How many genes are missed by consensus but detected by plain EM?
cat("\n=== Genes lost by consensus filtering ===\n")
consensus_lost <- merge(
  sel_a[, .(gene, sel_tpm = mean_tpm)],
  em_a[, .(gene, em_tpm = mean_tpm)],
  by = "gene", all = TRUE
)
consensus_lost[is.na(sel_tpm), sel_tpm := 0]
consensus_lost[is.na(em_tpm), em_tpm := 0]

# Match to TaqMan
taq_match <- merge(consensus_lost, taq_detected, by.x = "gene", by.y = "gene_symbol")
lost_by_consensus <- taq_match[em_tpm > 0 & sel_tpm == 0]
cat(sprintf("\nGenes detected by EM but zeroed by Sel+Support: %d\n", nrow(lost_by_consensus)))
cat(sprintf("  Of these with TaqMan > 0.1: %d\n", sum(lost_by_consensus$taq_A > 0.1)))
cat(sprintf("  Of these with TaqMan > 0.01: %d\n", sum(lost_by_consensus$taq_A > 0.01)))

# What EM TPM did these have?
cat("\nDistribution of EM TPM for lost genes:\n")
cat(sprintf("  median: %.2f\n", median(lost_by_consensus$em_tpm)))
cat(sprintf("  mean: %.2f\n", mean(lost_by_consensus$em_tpm)))
cat(sprintf("  max: %.2f\n", max(lost_by_consensus$em_tpm)))
cat(sprintf("  EM TPM > 1: %d\n", sum(lost_by_consensus$em_tpm > 1)))
cat(sprintf("  EM TPM > 10: %d\n", sum(lost_by_consensus$em_tpm > 10)))

# What are the highest-TaqMan genes that we miss?
cat("\nTop 20 highest-TaqMan genes lost by Sel+Support:\n")
cat(sprintf("%-15s %10s %10s\n", "Gene", "TaqMan(A)", "EM TPM"))
cat(paste0(rep("-", 38), collapse = ""), "\n")
top_lost <- lost_by_consensus[order(-taq_A)][1:min(20, nrow(lost_by_consensus))]
for (i in 1:nrow(top_lost)) {
  cat(sprintf("%-15s %10.4f %10.2f\n", top_lost$gene[i], top_lost$taq_A[i], top_lost$em_tpm[i]))
}

# How many isoforms do the lost genes have?
cat("\n=== Isoform count for lost genes ===\n")
# Load one quant file to get transcript-gene mapping
qdf <- fread(file.path(benchdir, "quant_em/A_1.quant"), sep = "\t")
setnames(qdf, c("target_name", "len", "eelen", "tpm", "ecount"))
qdf[, gene := parse_gene(target_name)]
iso_count <- qdf[!is.na(gene), .N, by = gene]
setnames(iso_count, "N", "n_isoforms")

lost_with_iso <- merge(lost_by_consensus, iso_count, by = "gene")
cat(sprintf("Median isoforms for lost genes: %.0f\n", median(lost_with_iso$n_isoforms)))
cat(sprintf("Median isoforms for all TaqMan genes: %.0f\n",
            median(merge(taq_detected, iso_count, by.x = "gene_symbol", by.y = "gene")$n_isoforms)))

cat("\nLost genes by isoform tier:\n")
for (tier in c(1, 2, 5, 10, 20)) {
  n <- sum(lost_with_iso$n_isoforms <= tier)
  cat(sprintf("  n_isoforms <= %d: %d genes\n", tier, n))
}

cat("\nDone.\n")
