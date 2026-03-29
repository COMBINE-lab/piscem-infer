#!/usr/bin/env Rscript
#
# Gene-level evaluation against TaqMan qRT-PCR ground truth.
# Aggregates transcript-level TPMs to gene level and correlates
# with TaqMan expression measurements for ~1000 genes.

library(data.table)

benchdir <- "seqc_benchmark"
pseudo <- 0.01

# ============================================================
# Load TaqMan data
# ============================================================

cat("Loading TaqMan data...\n")

# Platform annotation: ID -> Gene Symbol
annot <- fread(cmd = "grep -v '^[#!^]' seqc_benchmark/GPL4097.annot", sep = "\t", header = TRUE,
               select = c("ID", "Gene symbol"))
setnames(annot, c("taqman_id", "gene_symbol"))
annot <- annot[gene_symbol != "" & gene_symbol != "---"]

# Expression matrix
expr_lines <- readLines("seqc_benchmark/taqman_raw.txt")
data_start <- grep('^"ID_REF"', expr_lines)
expr_dt <- fread(text = expr_lines[data_start:length(expr_lines)], sep = "\t", header = TRUE)
setnames(expr_dt, 1, "taqman_id")

# Map column names to samples: A1-A4, B1-B4, C1-C4, D1-D4
# Columns are GSM129638-653 in order: A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,D4
sample_cols <- names(expr_dt)[-1]
sample_labels <- c(paste0("A_", 1:4), paste0("B_", 1:4), paste0("C_", 1:4), paste0("D_", 1:4))
setnames(expr_dt, sample_cols, sample_labels)

# Compute mean per sample type
expr_dt[, taq_A := rowMeans(.SD), .SDcols = paste0("A_", 1:4)]
expr_dt[, taq_B := rowMeans(.SD), .SDcols = paste0("B_", 1:4)]
expr_dt[, taq_C := rowMeans(.SD), .SDcols = paste0("C_", 1:4)]
expr_dt[, taq_D := rowMeans(.SD), .SDcols = paste0("D_", 1:4)]

# Merge with gene symbols
taqman <- merge(expr_dt[, .(taqman_id, taq_A, taq_B, taq_C, taq_D)], annot, by = "taqman_id")
# Remove duplicates (keep first)
taqman <- taqman[!duplicated(gene_symbol)]
cat(sprintf("Loaded %d TaqMan genes with expression data\n", nrow(taqman)))

# ============================================================
# RNA-seq: aggregate transcript TPMs to gene level
# ============================================================

# Parse gene name from GENCODE transcript IDs
# Format: ENST...|ENSG...|...|...|name-NNN|gene_name|len|...
parse_gene <- function(target_name) {
  sapply(strsplit(target_name, "\\|"), function(x) if (length(x) >= 6) x[6] else NA)
}

read_piscem <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- fread(path, sep = "\t")
  setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
  df[, gene := parse_gene(target_name)]
  df[!is.na(gene), .(tpm = sum(tpm)), by = gene]
}

read_salmon <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- fread(path, sep = "\t")
  df[, gene := parse_gene(Name)]
  df[!is.na(gene), .(tpm = sum(TPM)), by = gene]
}

read_kallisto <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- fread(path, sep = "\t")
  df[, gene := parse_gene(target_id)]
  df[!is.na(gene), .(tpm = sum(tpm)), by = gene]
}

samples <- c("A", "B", "C", "D")
reps <- 1:4

methods <- list(
  "Plain EM" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_em", paste0(s, "_", r, ".quant"))
  ),
  "Sel+Supp" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_support", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Sel+Adapt" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_adapt", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Adpt+Rescue" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_adapt_rescue", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Salmon" = list(
    reader = read_salmon,
    path_fn = function(s, r) file.path(benchdir, "quant_salmon", paste0(s, "_", r), "quant.sf")
  ),
  "Salmon noBias" = list(
    reader = read_salmon,
    path_fn = function(s, r) file.path(benchdir, "quant_salmon_nobias", paste0(s, "_", r), "quant.sf")
  ),
  "Soft Rescue" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_adapt_softrescue", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Salmon EM" = list(
    reader = read_salmon,
    path_fn = function(s, r) file.path(benchdir, "quant_salmon_em", paste0(s, "_", r), "quant.sf")
  ),
  "Kallisto" = list(
    reader = read_kallisto,
    path_fn = function(s, r) file.path(benchdir, "quant_kallisto", paste0(s, "_", r), "abundance.tsv")
  )
)

# Compute gene-level mean TPM per sample type
gene_means <- function(meth, sample_type) {
  gene_list <- list()
  for (r in reps) {
    path <- meth$path_fn(sample_type, r)
    df <- meth$reader(path)
    if (!is.null(df)) {
      setnames(df, "tpm", paste0("tpm_", r))
      gene_list[[as.character(r)]] <- df
    }
  }
  if (length(gene_list) == 0) return(NULL)
  merged <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), gene_list)
  tpm_cols <- setdiff(names(merged), "gene")
  merged$mean_tpm <- rowMeans(as.matrix(merged[, ..tpm_cols]), na.rm = TRUE)
  merged[, .(gene, mean_tpm)]
}

# ============================================================
# Correlate with TaqMan
# ============================================================

cat("\n=== Gene-Level Correlation with TaqMan qRT-PCR ===\n")
cat("(Pearson on log2 scale, genes with TaqMan > 0.001 and RNA-seq TPM > 0)\n\n")

cat(sprintf("%-14s | %10s %10s %10s %10s | %10s %10s | %6s\n",
            "Method", "Pears(A)", "Pears(B)", "Pears(C)", "Pears(D)",
            "Spear(A)", "Spear(B)", "N genes"))
cat(paste0(rep("-", 95), collapse = ""), "\n")

for (mname in names(methods)) {
  meth <- methods[[mname]]
  pearson_vals <- c()
  spearman_vals <- c()
  n_genes <- NA

  for (s in c("A", "B", "C", "D")) {
    gm <- gene_means(meth, s)
    if (is.null(gm)) { pearson_vals <- c(pearson_vals, NA); spearman_vals <- c(spearman_vals, NA); next }

    # Match to TaqMan
    matched <- merge(gm, taqman[, .(gene_symbol, taq = get(paste0("taq_", s)))],
                     by.x = "gene", by.y = "gene_symbol")
    # Filter to detected in both
    valid <- matched$taq > 0.001 & matched$mean_tpm > 0
    if (sum(valid) < 10) { pearson_vals <- c(pearson_vals, NA); spearman_vals <- c(spearman_vals, NA); next }

    pc <- cor(log2(matched$mean_tpm[valid] + pseudo), log2(matched$taq[valid] + pseudo), method = "pearson")
    sc <- cor(matched$mean_tpm[valid], matched$taq[valid], method = "spearman")
    pearson_vals <- c(pearson_vals, pc)
    spearman_vals <- c(spearman_vals, sc)
    if (s == "A") n_genes <- sum(valid)
  }

  cat(sprintf("%-14s | %10.4f %10.4f %10.4f %10.4f | %10.4f %10.4f | %6d\n",
              mname,
              pearson_vals[1], pearson_vals[2], pearson_vals[3], pearson_vals[4],
              spearman_vals[1], spearman_vals[2],
              ifelse(is.na(n_genes), 0, n_genes)))
}

# ============================================================
# TaqMan FC comparison (A vs B)
# ============================================================

cat("\n=== Gene-Level Fold Change: A vs B (RNA-seq vs TaqMan) ===\n")
cat("(Genes detected by TaqMan in both A and B, TaqMan > 0.001)\n\n")

cat(sprintf("%-14s | %10s %10s | %6s\n",
            "Method", "Pears(FC)", "FC slope", "N genes"))
cat(paste0(rep("-", 50), collapse = ""), "\n")

for (mname in names(methods)) {
  meth <- methods[[mname]]
  gm_a <- gene_means(meth, "A")
  gm_b <- gene_means(meth, "B")
  if (is.null(gm_a) || is.null(gm_b)) {
    cat(sprintf("%-14s | %10s %10s | %6s\n", mname, "N/A", "N/A", "N/A"))
    next
  }

  # Merge RNA-seq A and B
  rnaseq <- merge(gm_a, gm_b, by = "gene", suffixes = c("_A", "_B"))

  # Merge with TaqMan
  matched <- merge(rnaseq, taqman[, .(gene_symbol, taq_A, taq_B)],
                   by.x = "gene", by.y = "gene_symbol")

  # Filter to detected in both by TaqMan
  valid <- matched$taq_A > 0.001 & matched$taq_B > 0.001 &
           matched$mean_tpm_A > 0 & matched$mean_tpm_B > 0

  if (sum(valid) < 50) {
    cat(sprintf("%-14s | %10s %10s | %6s\n", mname, "N/A", "N/A", "N/A"))
    next
  }

  rnaseq_lfc <- log2(matched$mean_tpm_A[valid] + pseudo) - log2(matched$mean_tpm_B[valid] + pseudo)
  taqman_lfc <- log2(matched$taq_A[valid] + pseudo) - log2(matched$taq_B[valid] + pseudo)

  pc <- cor(rnaseq_lfc, taqman_lfc, method = "pearson")
  fit <- lm(rnaseq_lfc ~ taqman_lfc)
  slope <- coef(fit)[2]

  cat(sprintf("%-14s | %10.4f %10.4f | %6d\n",
              mname, pc, slope, sum(valid)))
}

cat("\nDone.\n")
