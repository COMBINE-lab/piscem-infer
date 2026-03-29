#!/usr/bin/env Rscript
library(data.table)
benchdir <- "seqc_benchmark"; pseudo <- 0.01
parse_gene <- function(tn) sapply(strsplit(tn, "\\|"), function(x) if (length(x) >= 6) x[6] else NA)

annot <- fread(cmd = "grep -v '^[#!^]' seqc_benchmark/GPL4097.annot", sep = "\t", header = TRUE, select = c("ID", "Gene symbol"))
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
taq_both <- taqman[taq_A > 0.001 & taq_B > 0.001]

gene_tpm <- function(dir, st) {
  dfs <- list()
  for (r in 1:4) {
    key <- paste0(st, "_", r)
    path <- file.path(benchdir, dir, key, paste0(key, ".quant"))
    if (!file.exists(path)) path <- file.path(benchdir, dir, paste0(key, ".quant"))
    if (!file.exists(path)) next
    df <- fread(path, sep = "\t"); setnames(df, c("target_name","len","eelen","tpm","ecount"))
    df[, gene := parse_gene(target_name)]
    g <- df[!is.na(gene), .(tpm = sum(tpm)), by = gene]
    setnames(g, "tpm", paste0("r", r))
    dfs[[r]] <- g
  }
  merged <- Reduce(function(a,b) merge(a,b,by="gene",all=TRUE), dfs)
  tpm_cols <- grep("^r", names(merged), value = TRUE)
  merged$mean_tpm <- rowMeans(as.matrix(merged[, ..tpm_cols]), na.rm = TRUE)
  merged[, .(gene, mean_tpm)]
}

cat("=== FC slope: ALL TaqMan genes in both (uniform, 802 genes) ===\n")
cat("=== vs. RESTRICTED to genes detected by the method in both conditions ===\n\n")

methods <- list(
  c("Plain EM", "quant_em"),
  c("Sel+Adapt", "quant_sel_adapt"),
  c("Cond+Gene", "quant_sel_adapt_condgene"),
  c("Adpt+Rescue", "quant_sel_adapt_rescue")
)

cat(sprintf("%-14s | %10s %6s | %10s %6s\n", "Method", "Uniform", "N", "Both-det", "N"))
cat(paste0(rep("-", 60), collapse = ""), "\n")

for (m in methods) {
  mname <- m[1]; dir <- m[2]
  ga <- gene_tpm(dir, "A"); gb <- gene_tpm(dir, "B")
  merged <- merge(ga, gb, by = "gene", suffixes = c("_A", "_B"))

  # Uniform: all 802 TaqMan genes, using 0 for undetected
  m_all <- match(taq_both$gene_symbol, merged$gene)
  tpm_a_all <- rep(0.0, nrow(taq_both)); tpm_b_all <- rep(0.0, nrow(taq_both))
  found <- !is.na(m_all)
  tpm_a_all[found] <- merged$mean_tpm_A[m_all[found]]
  tpm_b_all[found] <- merged$mean_tpm_B[m_all[found]]
  est_lfc_all <- log2(tpm_a_all + pseudo) - log2(tpm_b_all + pseudo)
  taq_lfc <- log2(taq_both$taq_A + pseudo) - log2(taq_both$taq_B + pseudo)
  fit_all <- lm(est_lfc_all ~ taq_lfc)

  # Both-detected: only genes where method has TPM > 0 in both A and B
  matched <- merge(taq_both, merged[mean_tpm_A > 0 & mean_tpm_B > 0], by.x = "gene_symbol", by.y = "gene")
  est_lfc_bd <- log2(matched$mean_tpm_A + pseudo) - log2(matched$mean_tpm_B + pseudo)
  taq_lfc_bd <- log2(matched$taq_A + pseudo) - log2(matched$taq_B + pseudo)
  fit_bd <- lm(est_lfc_bd ~ taq_lfc_bd)

  cat(sprintf("%-14s | %10.4f %6d | %10.4f %6d\n",
              mname, coef(fit_all)[2], nrow(taq_both), coef(fit_bd)[2], nrow(matched)))
}
