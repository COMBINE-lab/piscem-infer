#!/usr/bin/env Rscript
library(data.table)
benchdir <- "airway_benchmark"
samples_u <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")

methods <- list(
  "Plain EM"       = "quant/em",
  "Sel+Supp(8)"    = "quant/sel_support_all8_fast",
  "Sel+Adapt(8)"   = "quant/sel_support_all8_adaptive",
  "Sel+Pos5(8)"    = "quant/sel_support_all8_pos5",
  "Pos5+Cov(8)"    = "quant/sel_support_all8_pos5_cov",
  "Pos5+CG(8)"     = "quant/sel_pos5_condgene"
)

load_method <- function(dir) {
  dfs <- list()
  for (s in samples_u) {
    path <- file.path(benchdir, dir, s, paste0(s, ".quant"))
    if (!file.exists(path)) path <- file.path(benchdir, dir, paste0(s, ".quant"))
    if (!file.exists(path)) return(NULL)
    df <- fread(path, sep="\t"); setnames(df, c("target_name","len","eelen","tpm","ecount"))
    dfs[[s]] <- df[, .(target_name, tpm)]; setnames(dfs[[s]], "tpm", s)
  }
  Reduce(function(a,b) merge(a,b,by="target_name",all=TRUE), dfs)
}

compute_cv <- function(merged) {
  mat <- as.matrix(merged[, ..samples_u])
  rm <- rowMeans(mat); rs <- apply(mat,1,sd)
  data.table(target_name=merged$target_name, mean_tpm=rm, cv=rs/(rm+1e-6))
}

parse_gene <- function(tn) sapply(strsplit(tn, "\\|"), function(x) if (length(x) >= 6) x[6] else NA)

cat("\n=== Untreated Replicate Concordance ===\n\n")
cat(sprintf("%-16s | %8s %8s %8s | %8s\n", "Method", "Med CV", "Mean CV", "CV>1", "N expr"))
cat(paste0(rep("-",60), collapse=""), "\n")
for (name in names(methods)) {
  m <- load_method(methods[[name]])
  if (is.null(m)) { cat(sprintf("%-16s | MISSING\n", name)); next }
  cv <- compute_cv(m); e <- cv[mean_tpm >= 1]
  cat(sprintf("%-16s | %8.4f %8.4f %8d | %8d\n", name,
              median(e$cv), mean(e$cv), sum(e$cv>1), nrow(e)))
}

# Gene complexity stratification
ref <- load_method(methods[[1]])
ref$gene_name <- parse_gene(ref$target_name)
gene_iso_count <- ref[, .N, by = gene_name]; setnames(gene_iso_count, "N", "n_isoforms")

tiers <- list("1 iso"=c(1,1), "2-5"=c(2,5), "6-10"=c(6,10), "11-20"=c(11,20), "21-50"=c(21,50), "51+"=c(51,9999))
cat("\n=== Stratified by gene complexity (median CV, TPM >= 1) ===\n\n")
cat(sprintf("%-16s", "Method"))
for (t in names(tiers)) cat(sprintf(" | %8s", t))
cat("\n"); cat(paste0(rep("-", 16 + 6*11), collapse=""), "\n")

for (name in names(methods)) {
  m <- load_method(methods[[name]])
  if (is.null(m)) next
  cv <- compute_cv(m)
  cv$gene_name <- parse_gene(cv$target_name)
  cv <- merge(cv, gene_iso_count, by="gene_name")
  cat(sprintf("%-16s", name))
  for (tn in names(tiers)) {
    bounds <- tiers[[tn]]
    tier <- cv[n_isoforms >= bounds[1] & n_isoforms <= bounds[2] & mean_tpm >= 1]
    cat(sprintf(" | %8.4f", median(tier$cv)))
  }
  cat("\n")
}
cat("\nDone.\n")
