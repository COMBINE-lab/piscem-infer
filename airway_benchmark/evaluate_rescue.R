#!/usr/bin/env Rscript
# Evaluate condition rescue vs baselines on airway 8-sample.

library(data.table)
benchdir <- "airway_benchmark"
samples_u <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")
samples_t <- c("SRR1039509", "SRR1039513", "SRR1039517", "SRR1039521")
all_samples <- c(samples_u, samples_t)
pseudo <- 0.01

methods <- list(
  "Plain EM"         = "quant/em",
  "Sel+Supp(8)"      = "quant/sel_support_all8_fast",
  "Sel+Supp(8)Resc"  = "quant/sel_support_all8_rescue",
  "Sel+Supp(8)Adpt"  = "quant/sel_support_all8_adaptive",
  "Adpt+Rescue"      = "quant/sel_support_all8_adapt_rescue"
)

load_method <- function(dir) {
  dfs <- list()
  for (s in all_samples) {
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

cat("Loading...\n")
all_results <- list()
for (name in names(methods)) {
  cat("  ", name, "\n")
  all_results[[name]] <- load_method(methods[[name]])
}

# Within-condition CV
for (cond in c("Untreated", "Treated")) {
  subs <- if (cond == "Untreated") samples_u else samples_t
  cat(sprintf("\n=== %s Replicate Concordance ===\n\n", cond))
  cat(sprintf("%-18s | %8s %8s %8s | %8s\n", "Method", "Med CV", "Mean CV", "CV>1", "N expr"))
  cat(paste0(rep("-",62), collapse=""), "\n")
  for (name in names(all_results)) {
    m <- all_results[[name]]
    if (is.null(m)) { cat(sprintf("%-18s | MISSING\n", name)); next }
    cv <- compute_cv(m, subs); e <- cv[mean_tpm >= 1]
    cat(sprintf("%-18s | %8.4f %8.4f %8d | %8d\n", name,
                median(e$cv), mean(e$cv), sum(e$cv>1), nrow(e)))
  }
}

# Between-condition fold change
cat("\n=== Between-condition fold change (dex vs untreated) ===\n\n")
cat(sprintf("%-18s | %10s %10s %10s | %8s\n", "Method", "Mean|LFC|", "Med|LFC|", "SD(LFC)", "N both"))
cat(paste0(rep("-",75), collapse=""), "\n")
for (name in names(all_results)) {
  m <- all_results[[name]]
  if (is.null(m)) next
  untrt_mean <- rowMeans(as.matrix(m[, ..samples_u]))
  trt_mean   <- rowMeans(as.matrix(m[, ..samples_t]))
  both_expr <- untrt_mean >= 1 & trt_mean >= 1
  if (sum(both_expr) > 0) {
    lfc <- log2(trt_mean[both_expr] + pseudo) - log2(untrt_mean[both_expr] + pseudo)
    cat(sprintf("%-18s | %10.4f %10.4f %10.4f | %8d\n",
                name, mean(abs(lfc)), median(abs(lfc)), sd(lfc), sum(both_expr)))
  }
}

cat("\nDone.\n")
