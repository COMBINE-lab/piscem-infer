#!/usr/bin/env Rscript
# Diagnose FC overcorrection in condition rescue.
# Compare gene-level FC from Phase 1 (pre-masking) vs Phase 2 (post-masking)
# to understand where the inflation comes from.

library(data.table)
benchdir <- "seqc_benchmark"
pseudo <- 0.01

parse_gene <- function(tn) sapply(strsplit(tn, "\\|"), function(x) if (length(x) >= 6) x[6] else NA)

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

# Gene-level TPM helper
gene_tpm <- function(dir, sample_type, reader_type = "piscem") {
  dfs <- list()
  for (r in 1:4) {
    key <- paste0(sample_type, "_", r)
    if (reader_type == "piscem") {
      path <- file.path(benchdir, dir, paste0(key, ".quant"))
      if (!file.exists(path)) path <- file.path(benchdir, dir, key, paste0(key, ".quant"))
    } else {
      path <- file.path(benchdir, dir, key, "quant.sf")
    }
    if (!file.exists(path)) next
    if (reader_type == "piscem") {
      df <- fread(path, sep = "\t"); setnames(df, c("target_name","len","eelen","tpm","ecount"))
    } else {
      df <- fread(path, sep = "\t"); df <- df[, .(target_name = Name, tpm = TPM)]
    }
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

# Load several methods
methods <- list(
  "Plain EM" = list(dir = "quant_em"),
  "Sel+Adapt" = list(dir = "quant_sel_adapt"),
  "Cond+Gene" = list(dir = "quant_sel_adapt_condgene")
)

# Compute gene-level FC (A vs B) for each method
taq_both <- taqman[taq_A > 0.001 & taq_B > 0.001]
taq_lfc <- log2(taq_both$taq_A + pseudo) - log2(taq_both$taq_B + pseudo)

cat("=== Gene-level FC (A vs B) stratified by TaqMan FC magnitude ===\n\n")

for (mname in names(methods)) {
  m <- methods[[mname]]
  ga <- gene_tpm(m$dir, "A")
  gb <- gene_tpm(m$dir, "B")
  merged <- merge(ga, gb, by = "gene", suffixes = c("_A", "_B"))
  matched <- merge(taq_both, merged, by.x = "gene_symbol", by.y = "gene")

  est_lfc <- log2(matched$mean_tpm_A + pseudo) - log2(matched$mean_tpm_B + pseudo)
  taq_l <- log2(matched$taq_A + pseudo) - log2(matched$taq_B + pseudo)

  cat(sprintf("=== %s ===\n", mname))

  # Stratify by TaqMan FC magnitude
  for (tier_name in c("|FC| < 0.5", "0.5 <= |FC| < 2", "2 <= |FC| < 5", "|FC| >= 5")) {
    sel <- switch(tier_name,
      "|FC| < 0.5" = abs(taq_l) < 0.5,
      "0.5 <= |FC| < 2" = abs(taq_l) >= 0.5 & abs(taq_l) < 2,
      "2 <= |FC| < 5" = abs(taq_l) >= 2 & abs(taq_l) < 5,
      "|FC| >= 5" = abs(taq_l) >= 5
    )
    if (sum(sel) < 10) {
      cat(sprintf("  %18s: n=%d (too few)\n", tier_name, sum(sel)))
      next
    }
    fit <- lm(est_lfc[sel] ~ taq_l[sel])
    slope <- coef(fit)[2]
    mean_abs_est <- mean(abs(est_lfc[sel]))
    mean_abs_taq <- mean(abs(taq_l[sel]))
    cat(sprintf("  %18s: n=%4d, slope=%.3f, mean|est|=%.2f, mean|taq|=%.2f\n",
                tier_name, sum(sel), slope, mean_abs_est, mean_abs_taq))
  }

  # Overall
  fit_all <- lm(est_lfc ~ taq_l)
  cat(sprintf("  %18s: n=%4d, slope=%.3f\n\n", "Overall", length(est_lfc), coef(fit_all)[2]))
}

# Direct comparison: which genes have the worst FC overcorrection?
cat("=== Genes with largest FC overcorrection (Cond+Gene vs TaqMan) ===\n\n")
ga <- gene_tpm("quant_sel_adapt_condgene", "A")
gb <- gene_tpm("quant_sel_adapt_condgene", "B")
merged <- merge(ga, gb, by = "gene", suffixes = c("_A", "_B"))
matched <- merge(taq_both, merged, by.x = "gene_symbol", by.y = "gene")
matched$est_lfc <- log2(matched$mean_tpm_A + pseudo) - log2(matched$mean_tpm_B + pseudo)
matched$taq_lfc <- log2(matched$taq_A + pseudo) - log2(matched$taq_B + pseudo)
matched$fc_error <- matched$est_lfc - matched$taq_lfc

# Top overcorrected (est FC too large)
cat("Top 15 most overcorrected (est FC too high relative to TaqMan):\n")
cat(sprintf("%-12s %8s %8s %8s %8s %8s\n", "Gene", "est_LFC", "taq_LFC", "error", "TPM_A", "TPM_B"))
cat(paste0(rep("-", 58), collapse = ""), "\n")
top_over <- matched[order(-fc_error)][1:15]
for (i in 1:nrow(top_over)) {
  r <- top_over[i]
  cat(sprintf("%-12s %8.2f %8.2f %8.2f %8.1f %8.1f\n",
              r$gene_symbol, r$est_lfc, r$taq_lfc, r$fc_error, r$mean_tpm_A, r$mean_tpm_B))
}

cat("\nDone.\n")
