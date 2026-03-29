#!/usr/bin/env Rscript
#
# Evaluate structural selection + consensus methods on airway data.
# Compares replicate concordance (CV) across 4 untreated replicates.

library(data.table)

benchdir <- "airway_benchmark"
samples <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")

methods <- list(
  "Plain EM"       = "quant/em",
  "EM+Sel"         = "quant/sel",
  "NullSink .50"   = "quant/nullsink_u4",
  "Cons TPM"       = "quant/consensus_tpm",
  "Cons Support"   = "quant/consensus_support",
  "Sel+Cons TPM"   = "quant/sel_tpm",
  "Sel+Cons Supp"  = "quant/sel_support"
)

# ---- Load quant files ----
cat("Loading quantification results...\n")

load_method <- function(method_dir) {
  dfs <- list()
  for (s in samples) {
    # Try both naming conventions
    path <- file.path(benchdir, method_dir, s, paste0(s, ".quant"))
    if (!file.exists(path)) {
      path <- file.path(benchdir, method_dir, paste0(s, ".quant"))
    }
    if (!file.exists(path)) {
      cat("  MISSING:", path, "\n")
      return(NULL)
    }
    df <- fread(path, sep = "\t")
    setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
    dfs[[s]] <- df[, .(target_name, tpm)]
    setnames(dfs[[s]], "tpm", s)
  }
  Reduce(function(a, b) merge(a, b, by = "target_name", all = TRUE), dfs)
}

all_results <- list()
for (name in names(methods)) {
  cat("  Loading", name, "...\n")
  all_results[[name]] <- load_method(methods[[name]])
}

# ---- Gene complexity ----
parse_gene <- function(target_name) {
  sapply(strsplit(target_name, "\\|"), function(x) x[6])
}

ref <- all_results[[1]]
if (is.null(ref)) stop("Plain EM baseline not found")
ref$gene_name <- parse_gene(ref$target_name)
gene_iso_count <- ref[, .N, by = gene_name]
setnames(gene_iso_count, "N", "n_isoforms")

# ---- Concordance ----
cat("\n=== Replicate Concordance (CV of TPM across 4 untreated replicates) ===\n\n")

compute_cv <- function(merged, min_mean_tpm = 1.0) {
  tpm_mat <- as.matrix(merged[, ..samples])
  row_mean <- rowMeans(tpm_mat)
  row_sd <- apply(tpm_mat, 1, sd)
  cv <- row_sd / (row_mean + 1e-6)
  data.table(
    target_name = merged$target_name,
    mean_tpm = row_mean,
    sd_tpm = row_sd,
    cv = cv
  )
}

cat(sprintf("%-16s | %8s %8s %8s | %8s %8s | %8s\n",
            "Method", "Med CV", "Mean CV", "CV>1",
            "Med CV>5", "Med CV>50", "N expr"))
cat(paste0(rep("-", 90), collapse = ""), "\n")

for (name in names(all_results)) {
  m <- all_results[[name]]
  if (is.null(m)) { cat(sprintf("%-16s | MISSING\n", name)); next }

  cv_df <- compute_cv(m)
  expr <- cv_df[mean_tpm >= 1]
  expr5 <- cv_df[mean_tpm >= 5]
  expr50 <- cv_df[mean_tpm >= 50]

  cat(sprintf("%-16s | %8.4f %8.4f %8d | %8.4f %8.4f | %8d\n",
              name,
              median(expr$cv), mean(expr$cv), sum(expr$cv > 1),
              median(expr5$cv), median(expr50$cv),
              nrow(expr)))
}

# ---- Stratified by gene complexity ----
cat("\n=== Stratified by gene complexity (isoforms per gene) ===\n")
cat("(Median CV of expressed isoforms, TPM >= 1)\n\n")

tiers <- list(
  "1 iso"     = c(1, 1),
  "2-5 iso"   = c(2, 5),
  "6-10 iso"  = c(6, 10),
  "11-20 iso" = c(11, 20),
  "21-50 iso" = c(21, 50),
  "51+ iso"   = c(51, 9999)
)

cat(sprintf("%-16s", "Method"))
for (t in names(tiers)) cat(sprintf(" | %10s", t))
cat("\n")
cat(paste0(rep("-", 18 + length(tiers) * 13), collapse = ""), "\n")

for (name in names(all_results)) {
  m <- all_results[[name]]
  if (is.null(m)) next

  cv_df <- compute_cv(m)
  cv_df$gene_name <- parse_gene(cv_df$target_name)
  cv_df <- merge(cv_df, gene_iso_count, by = "gene_name")

  cat(sprintf("%-16s", name))
  for (tier_name in names(tiers)) {
    bounds <- tiers[[tier_name]]
    tier <- cv_df[n_isoforms >= bounds[1] & n_isoforms <= bounds[2] & mean_tpm >= 1]
    if (nrow(tier) > 0) {
      cat(sprintf(" | %10.4f", median(tier$cv)))
    } else {
      cat(sprintf(" | %10s", "n/a"))
    }
  }
  cat("\n")
}

# ---- Paired comparison ----
cat("\n=== Paired Wilcoxon: is method CV significantly different from Plain EM? ===\n\n")

em <- all_results[["Plain EM"]]
if (!is.null(em)) {
  for (comp_name in names(all_results)) {
    if (comp_name == "Plain EM") next
    comp <- all_results[[comp_name]]
    if (is.null(comp)) next

    cv_em <- compute_cv(em)
    cv_comp <- compute_cv(comp)

    both <- merge(cv_em, cv_comp, by = "target_name", suffixes = c(".em", ".comp"))
    both <- both[mean_tpm.em >= 1 & mean_tpm.comp >= 1]

    wt <- wilcox.test(both$cv.em, both$cv.comp, paired = TRUE, alternative = "greater")
    n_better <- sum(both$cv.comp < both$cv.em)
    n_worse <- sum(both$cv.comp > both$cv.em)

    cat(sprintf("  %-16s: p=%.2e, %d improved / %d worsened / %d tied (n=%d)\n",
                comp_name, wt$p.value, n_better, n_worse,
                nrow(both) - n_better - n_worse, nrow(both)))
  }
}

cat("\n=== Done ===\n")
