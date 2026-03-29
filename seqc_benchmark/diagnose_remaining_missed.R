#!/usr/bin/env Rscript
# What are the ~170 TaqMan genes that Cond+Gene10 misses but plain EM detects?

library(data.table)
benchdir <- "seqc_benchmark"

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

parse_gene <- function(tn) sapply(strsplit(tn, "\\|"), function(x) if (length(x) >= 6) x[6] else NA)

# Gene-level TPM from each method (mean of 4 reps, sample A)
gene_tpm_from <- function(dir, reader_type = "piscem") {
  dfs <- list()
  for (r in 1:4) {
    key <- paste0("A_", r)
    if (reader_type == "piscem") {
      path <- file.path(benchdir, dir, paste0(key, ".quant"))
      if (!file.exists(path)) path <- file.path(benchdir, dir, key, paste0(key, ".quant"))
      if (!file.exists(path)) next
      df <- fread(path, sep = "\t"); setnames(df, c("target_name","len","eelen","tpm","ecount"))
    } else {
      path <- file.path(benchdir, dir, key, "quant.sf")
      if (!file.exists(path)) next
      df <- fread(path, sep = "\t"); df <- df[, .(target_name = Name, tpm = TPM)]
    }
    df[, gene := parse_gene(target_name)]
    gene_df <- df[!is.na(gene), .(tpm = sum(tpm)), by = gene]
    setnames(gene_df, "tpm", paste0("r", r))
    dfs[[r]] <- gene_df
  }
  merged <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), dfs)
  tpm_cols <- grep("^r", names(merged), value = TRUE)
  merged$mean_tpm <- rowMeans(as.matrix(merged[, ..tpm_cols]), na.rm = TRUE)
  merged[, .(gene, mean_tpm)]
}

em <- gene_tpm_from("quant_em")
condgene <- gene_tpm_from("quant_sel_adapt_condgene")
adapt <- gene_tpm_from("quant_sel_adapt")

# TaqMan genes detected in A
taq_A <- taqman[taq_A > 0.001]

# Merge all
comp <- merge(taq_A, em, by.x = "gene_symbol", by.y = "gene", all.x = TRUE)
setnames(comp, "mean_tpm", "em_tpm")
comp <- merge(comp, condgene, by.x = "gene_symbol", by.y = "gene", all.x = TRUE)
setnames(comp, "mean_tpm", "cg_tpm")
comp <- merge(comp, adapt, by.x = "gene_symbol", by.y = "gene", all.x = TRUE)
setnames(comp, "mean_tpm", "adapt_tpm")
comp[is.na(em_tpm), em_tpm := 0]
comp[is.na(cg_tpm), cg_tpm := 0]
comp[is.na(adapt_tpm), adapt_tpm := 0]

# Genes detected by EM but missed by Cond+Gene10
missed <- comp[em_tpm > 0 & cg_tpm == 0]
cat(sprintf("Genes detected by EM but missed by Cond+Gene10: %d\n", nrow(missed)))
cat(sprintf("  Of these also missed by Sel+Adapt: %d\n", sum(missed$adapt_tpm == 0)))
cat(sprintf("  Of these detected by Sel+Adapt: %d\n", sum(missed$adapt_tpm > 0)))

# Expression distribution
cat("\n=== Missed genes by EM TPM ===\n")
for (tier in c(0, 0.1, 1, 5, 10, 50)) {
  n <- sum(missed$em_tpm > tier)
  cat(sprintf("  EM TPM > %5.1f: %d genes\n", tier, n))
}

cat("\n=== Missed genes by TaqMan level ===\n")
for (tier in c(0.001, 0.01, 0.05, 0.1, 0.5, 1.0)) {
  n <- sum(missed$taq_A > tier)
  cat(sprintf("  TaqMan > %5.3f: %d genes\n", tier, n))
}

cat("\n=== Top 30 missed genes (by EM TPM) ===\n")
cat(sprintf("%-12s %8s %8s %8s %8s\n", "Gene", "EM_TPM", "CG_TPM", "Adapt_T", "TaqMan"))
cat(paste0(rep("-", 50), collapse = ""), "\n")
top <- missed[order(-em_tpm)][1:min(30, nrow(missed))]
for (i in 1:nrow(top)) {
  r <- top[i]
  cat(sprintf("%-12s %8.2f %8.2f %8.2f %8.4f\n",
              r$gene_symbol, r$em_tpm, r$cg_tpm, r$adapt_tpm, r$taq_A))
}

# What fraction of these missed genes have very low TaqMan?
cat(sprintf("\n=== Summary ===\n"))
cat(sprintf("Missed with TaqMan < 0.01: %d / %d (%.0f%%) — near TaqMan detection limit\n",
            sum(missed$taq_A < 0.01), nrow(missed), 100 * sum(missed$taq_A < 0.01) / nrow(missed)))
cat(sprintf("Missed with TaqMan < 0.05: %d / %d (%.0f%%)\n",
            sum(missed$taq_A < 0.05), nrow(missed), 100 * sum(missed$taq_A < 0.05) / nrow(missed)))
cat(sprintf("Missed with EM TPM < 1: %d / %d (%.0f%%) — very low RNA-seq expression\n",
            sum(missed$em_tpm < 1), nrow(missed), 100 * sum(missed$em_tpm < 1) / nrow(missed)))

cat("\nDone.\n")
