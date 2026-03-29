#!/usr/bin/env Rscript
# Diagnose: what's different about TaqMan genes missed by consensus filtering?
# Examine Phase 1 EM estimates, EC structure, and per-isoform evidence
# to understand what signal exists for recovery.

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
taqman <- merge(expr_dt[, .(taqman_id, taq_A)], annot, by = "taqman_id")
taqman <- taqman[!duplicated(gene_symbol)]

# Parse gene from GENCODE names
parse_gene <- function(target_name) {
  sapply(strsplit(target_name, "\\|"), function(x) if (length(x) >= 6) x[6] else NA)
}

# Load Phase 1 plain EM results (sample A, all 4 replicates)
cat("Loading Phase 1 EM estimates for sample A...\n")
em_txps <- list()
for (r in 1:4) {
  path <- file.path(benchdir, "quant_em", paste0("A_", r, ".quant"))
  df <- fread(path, sep = "\t")
  setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
  df$gene <- parse_gene(df$target_name)
  setnames(df, "tpm", paste0("tpm_", r))
  em_txps[[r]] <- df[, c("target_name", "gene", paste0("tpm_", r)), with = FALSE]
}
em <- Reduce(function(a, b) merge(a, b, by = c("target_name", "gene")), em_txps)
em$mean_tpm <- rowMeans(em[, .(tpm_1, tpm_2, tpm_3, tpm_4)])

# Load consensus Sel+Adapt results (sample A)
cons_txps <- list()
for (r in 1:4) {
  path <- file.path(benchdir, "quant_sel_adapt", paste0("A_", r), paste0("A_", r, ".quant"))
  df <- fread(path, sep = "\t")
  setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
  setnames(df, "tpm", paste0("cons_tpm_", r))
  cons_txps[[r]] <- df[, c("target_name", paste0("cons_tpm_", r)), with = FALSE]
}
cons <- Reduce(function(a, b) merge(a, b, by = "target_name"), cons_txps)
cons$cons_mean_tpm <- rowMeans(cons[, .(cons_tpm_1, cons_tpm_2, cons_tpm_3, cons_tpm_4)])

# Merge EM and consensus at transcript level
txp <- merge(em[, .(target_name, gene, mean_tpm)], cons[, .(target_name, cons_mean_tpm)], by = "target_name")

# Gene-level aggregation
gene_em <- txp[!is.na(gene), .(
  gene_tpm_em = sum(mean_tpm),
  n_isoforms = .N,
  n_isoforms_expressed = sum(mean_tpm > 0),
  max_isoform_tpm = max(mean_tpm),
  top_isoform_frac = max(mean_tpm) / (sum(mean_tpm) + 1e-10)
), by = gene]

gene_cons <- txp[!is.na(gene), .(
  gene_tpm_cons = sum(cons_mean_tpm),
  n_isoforms_consensus = sum(cons_mean_tpm > 0)
), by = gene]

gene <- merge(gene_em, gene_cons, by = "gene")

# Match to TaqMan
taq_a_detected <- taqman[taq_A > 0.001]
gene_taq <- merge(gene, taq_a_detected, by.x = "gene", by.y = "gene_symbol")

# Classify: detected by consensus vs missed
gene_taq$detected <- gene_taq$gene_tpm_cons > 0
gene_taq$missed <- !gene_taq$detected

cat(sprintf("\n=== TaqMan genes in sample A: %d ===\n", nrow(gene_taq)))
cat(sprintf("Detected by Sel+Adapt: %d\n", sum(gene_taq$detected)))
cat(sprintf("Missed by Sel+Adapt: %d\n", sum(gene_taq$missed)))

# Compare properties of detected vs missed genes
cat("\n=== Properties of detected vs missed genes ===\n\n")
cat(sprintf("%-25s | %12s %12s\n", "Property", "Detected", "Missed"))
cat(paste0(rep("-", 55), collapse = ""), "\n")

props <- list(
  "N genes" = function(d) nrow(d),
  "Median gene TPM (EM)" = function(d) median(d$gene_tpm_em),
  "Median n_isoforms" = function(d) median(d$n_isoforms),
  "Median n_iso expressed" = function(d) median(d$n_isoforms_expressed),
  "Median max_iso TPM" = function(d) median(d$max_isoform_tpm),
  "Median top_iso fraction" = function(d) median(d$top_isoform_frac),
  "Median TaqMan" = function(d) median(d$taq_A),
  "Gene TPM > 10 (EM)" = function(d) sum(d$gene_tpm_em > 10),
  "Gene TPM > 100 (EM)" = function(d) sum(d$gene_tpm_em > 100),
  "top_iso_frac < 0.5" = function(d) sum(d$top_isoform_frac < 0.5),
  "top_iso_frac < 0.25" = function(d) sum(d$top_isoform_frac < 0.25),
  "n_isoforms > 10" = function(d) sum(d$n_isoforms > 10),
  "n_isoforms > 20" = function(d) sum(d$n_isoforms > 20)
)

det <- gene_taq[detected == TRUE]
mis <- gene_taq[missed == TRUE]
for (pname in names(props)) {
  fn <- props[[pname]]
  cat(sprintf("%-25s | %12s %12s\n", pname, format(fn(det), digits = 4), format(fn(mis), digits = 4)))
}

# Distribution of missed genes by expression level
cat("\n=== Missed genes by EM gene-level TPM tier ===\n\n")
tiers <- c(0, 1, 10, 50, 100, 500, Inf)
tier_labels <- c("<1", "1-10", "10-50", "50-100", "100-500", ">500")
for (i in 1:(length(tiers)-1)) {
  in_tier <- mis[gene_tpm_em >= tiers[i] & gene_tpm_em < tiers[i+1]]
  cat(sprintf("  EM gene TPM %8s: %3d missed genes (median %d isoforms, top_iso_frac %.2f)\n",
              tier_labels[i], nrow(in_tier),
              ifelse(nrow(in_tier) > 0, as.integer(median(in_tier$n_isoforms)), 0L),
              ifelse(nrow(in_tier) > 0, median(in_tier$top_isoform_frac), 0)))
}

# Top 20 most-expressed missed genes
cat("\n=== Top 20 most-expressed missed genes (by EM gene TPM) ===\n\n")
cat(sprintf("%-12s %8s %8s %6s %6s %8s\n",
            "Gene", "EM_TPM", "TaqMan", "#Iso", "#Expr", "TopFrac"))
cat(paste0(rep("-", 55), collapse = ""), "\n")
top_missed <- mis[order(-gene_tpm_em)][1:min(20, nrow(mis))]
for (i in 1:nrow(top_missed)) {
  r <- top_missed[i]
  cat(sprintf("%-12s %8.1f %8.4f %6d %6d %8.3f\n",
              r$gene, r$gene_tpm_em, r$taq_A, r$n_isoforms,
              r$n_isoforms_expressed, r$top_isoform_frac))
}

cat("\nDone.\n")
