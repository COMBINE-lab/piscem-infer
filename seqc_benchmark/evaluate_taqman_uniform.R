#!/usr/bin/env Rscript
#
# Gene-level evaluation against TaqMan qRT-PCR ground truth.
# UNIFORM GROUND SET: all methods evaluated on the same set of TaqMan genes.
# Genes not detected by a method get TPM=0 (counted as FN in correlation).

library(data.table)

benchdir <- "seqc_benchmark"
pseudo <- 0.01

# ============================================================
# Load TaqMan data
# ============================================================

cat("Loading TaqMan data...\n")

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
expr_dt[, taq_C := rowMeans(.SD), .SDcols = paste0("C_", 1:4)]
expr_dt[, taq_D := rowMeans(.SD), .SDcols = paste0("D_", 1:4)]
taqman <- merge(expr_dt[, .(taqman_id, taq_A, taq_B, taq_C, taq_D)], annot, by = "taqman_id")
taqman <- taqman[!duplicated(gene_symbol)]

# ============================================================
# Define uniform ground set
# ============================================================

# Ground set: TaqMan genes detected in at least one sample (TaqMan > 0.001)
ground_set <- taqman[taq_A > 0.001 | taq_B > 0.001 | taq_C > 0.001 | taq_D > 0.001]
cat(sprintf("Uniform ground set: %d TaqMan genes\n", nrow(ground_set)))

# For sample-specific evaluation, use genes detected by TaqMan in that sample
ground_A <- taqman[taq_A > 0.001]$gene_symbol
ground_B <- taqman[taq_B > 0.001]$gene_symbol
ground_AB <- intersect(ground_A, ground_B)  # for FC evaluation
cat(sprintf("  Sample A: %d genes, Sample B: %d genes, both: %d\n",
            length(ground_A), length(ground_B), length(ground_AB)))

# ============================================================
# RNA-seq: aggregate transcript TPMs to gene level
# ============================================================

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
  "Soft Rescue" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_adapt_softrescue", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Salmon" = list(
    reader = read_salmon,
    path_fn = function(s, r) file.path(benchdir, "quant_salmon", paste0(s, "_", r), "quant.sf")
  ),
  "Cond+Gene" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_adapt_condgene", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
  ),
  "Adpt+GeneRsc" = list(
    reader = read_piscem,
    path_fn = function(s, r) file.path(benchdir, "quant_sel_adapt_generescue", paste0(s, "_", r), paste0(s, "_", r, ".quant"))
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

# Get gene-level mean TPM for a method+sample, evaluated on a fixed gene set.
# Returns TPM=0 for genes in the ground set that the method doesn't detect.
gene_means_on_ground <- function(meth, sample_type, gene_set) {
  gene_list <- list()
  for (r in reps) {
    path <- meth$path_fn(sample_type, r)
    df <- meth$reader(path)
    if (!is.null(df)) {
      setnames(df, "tpm", paste0("tpm_", r))
      gene_list[[as.character(r)]] <- df
    }
  }
  if (length(gene_list) == 0) return(rep(0, length(gene_set)))
  merged <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), gene_list)
  tpm_cols <- setdiff(names(merged), "gene")
  merged$mean_tpm <- rowMeans(as.matrix(merged[, ..tpm_cols]), na.rm = TRUE)

  # Map to ground set (0 for missing genes)
  m <- match(gene_set, merged$gene)
  result <- rep(0.0, length(gene_set))
  found <- !is.na(m)
  result[found] <- merged$mean_tpm[m[found]]
  result
}

# ============================================================
# Evaluation 1: Correlation with TaqMan on uniform ground set
# ============================================================

cat("\n=== Gene-Level Correlation with TaqMan (UNIFORM ground set) ===\n")
cat("All methods evaluated on the same gene set; undetected genes = TPM 0\n\n")

cat(sprintf("%-14s | %10s %10s %10s %10s | %6s %6s\n",
            "Method", "Pears(A)", "Pears(B)", "Spear(A)", "Spear(B)", "Det A", "Det B"))
cat(paste0(rep("-", 80), collapse = ""), "\n")

for (mname in names(methods)) {
  meth <- methods[[mname]]
  tpm_a <- gene_means_on_ground(meth, "A", ground_A)
  tpm_b <- gene_means_on_ground(meth, "B", ground_B)
  taq_a <- taqman[match(ground_A, gene_symbol)]$taq_A
  taq_b <- taqman[match(ground_B, gene_symbol)]$taq_B

  pa <- cor(log2(tpm_a + pseudo), log2(taq_a + pseudo), method = "pearson")
  pb <- cor(log2(tpm_b + pseudo), log2(taq_b + pseudo), method = "pearson")
  sa <- cor(tpm_a, taq_a, method = "spearman")
  sb <- cor(tpm_b, taq_b, method = "spearman")
  det_a <- sum(tpm_a > 0)
  det_b <- sum(tpm_b > 0)

  cat(sprintf("%-14s | %10.4f %10.4f %10.4f %10.4f | %6d %6d\n",
              mname, pa, pb, sa, sb, det_a, det_b))
}

# ============================================================
# Evaluation 2: FC slope on uniform ground set
# ============================================================

cat("\n=== Gene-Level FC: A vs B (UNIFORM ground set) ===\n")
cat("All methods on same %d genes detected by TaqMan in both A and B\n\n", length(ground_AB))

cat(sprintf("%-14s | %10s %10s | %6s\n", "Method", "Pears(FC)", "FC slope", "Det AB"))
cat(paste0(rep("-", 45), collapse = ""), "\n")

for (mname in names(methods)) {
  meth <- methods[[mname]]
  tpm_a <- gene_means_on_ground(meth, "A", ground_AB)
  tpm_b <- gene_means_on_ground(meth, "B", ground_AB)
  taq_a <- taqman[match(ground_AB, gene_symbol)]$taq_A
  taq_b <- taqman[match(ground_AB, gene_symbol)]$taq_B

  rnaseq_lfc <- log2(tpm_a + pseudo) - log2(tpm_b + pseudo)
  taqman_lfc <- log2(taq_a + pseudo) - log2(taq_b + pseudo)

  pc <- cor(rnaseq_lfc, taqman_lfc, method = "pearson")
  fit <- lm(rnaseq_lfc ~ taqman_lfc)
  slope <- coef(fit)[2]
  det_ab <- sum(tpm_a > 0 & tpm_b > 0)

  cat(sprintf("%-14s | %10.4f %10.4f | %6d\n", mname, pc, slope, det_ab))
}

# ============================================================
# Evaluation 3: Detection rate
# ============================================================

cat("\n=== Detection Rate on TaqMan Ground Set ===\n\n")
cat(sprintf("%-14s | %10s %10s %10s\n", "Method", "Det A", "Det B", "Det A&B"))
cat(paste0(rep("-", 50), collapse = ""), "\n")

for (mname in names(methods)) {
  meth <- methods[[mname]]
  tpm_a <- gene_means_on_ground(meth, "A", ground_A)
  tpm_b <- gene_means_on_ground(meth, "B", ground_B)
  tpm_ab_a <- gene_means_on_ground(meth, "A", ground_AB)
  tpm_ab_b <- gene_means_on_ground(meth, "B", ground_AB)

  cat(sprintf("%-14s | %6d/%3d %6d/%3d %6d/%3d\n",
              mname,
              sum(tpm_a > 0), length(ground_A),
              sum(tpm_b > 0), length(ground_B),
              sum(tpm_ab_a > 0 & tpm_ab_b > 0), length(ground_AB)))
}

cat("\nDone.\n")
