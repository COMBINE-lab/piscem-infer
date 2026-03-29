#!/usr/bin/env Rscript
# Compare SQUAREM vs plain EM for Sel+Support on airway 8-sample.

library(data.table)
benchdir <- "airway_benchmark"
samples_u <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")
samples_t <- c("SRR1039509", "SRR1039513", "SRR1039517", "SRR1039521")

methods <- list(
  "Plain EM"       = "quant/em",
  "Sel+Supp(8)"    = "quant/sel_support_all8",
  "Sel+Supp(8)Sq"  = "quant/sel_support_all8_sq",
  "Sel+Supp(8)Sq3" = "quant/sel_support_all8_sq3",
  "Sel+Supp(8)Fast" = "quant/sel_support_all8_fast"
)

load_method <- function(dir, sample_list) {
  dfs <- list()
  for (s in sample_list) {
    path <- file.path(benchdir, dir, s, paste0(s, ".quant"))
    if (!file.exists(path)) path <- file.path(benchdir, dir, paste0(s, ".quant"))
    if (!file.exists(path)) return(NULL)
    df <- fread(path, sep="\t"); setnames(df, c("target_name","len","eelen","tpm","ecount"))
    dfs[[s]] <- df[, .(target_name, tpm)]; setnames(dfs[[s]], "tpm", s)
  }
  Reduce(function(a,b) merge(a,b,by="target_name",all=TRUE), dfs)
}

compute_cv <- function(merged, subs) {
  mat <- as.matrix(merged[, ..subs])
  rm <- rowMeans(mat); rs <- apply(mat,1,sd)
  data.table(target_name=merged$target_name, mean_tpm=rm, cv=rs/(rm+1e-6))
}

for (cond in c("Untreated", "Treated")) {
  subs <- if (cond == "Untreated") samples_u else samples_t
  cat(sprintf("\n=== %s Replicate Concordance ===\n\n", cond))
  cat(sprintf("%-18s | %8s %8s %8s | %8s\n", "Method", "Med CV", "Mean CV", "CV>1", "N expr"))
  cat(paste0(rep("-",62), collapse=""), "\n")
  for (name in names(methods)) {
    m <- load_method(methods[[name]], c(samples_u, samples_t))
    if (is.null(m)) { m <- load_method(methods[[name]], subs) }
    if (is.null(m)) { cat(sprintf("%-18s | MISSING\n", name)); next }
    cv <- compute_cv(m, subs); e <- cv[mean_tpm >= 1]
    cat(sprintf("%-18s | %8.4f %8.4f %8d | %8d\n", name,
                median(e$cv), mean(e$cv), sum(e$cv>1), nrow(e)))
  }
}
cat("\nDone.\n")
