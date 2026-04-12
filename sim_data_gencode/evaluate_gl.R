#!/usr/bin/env Rscript
#
# Evaluate group LASSO vs single-sample EM vs hierarchical multi-sample.

simdir <- "sim_data_gencode"

# ---- Load ground truth ----
gt <- read.csv(file.path(simdir, "ground_truth.csv"), stringsAsFactors = FALSE)
gt$short_id <- sub(" .*", "", gt$transcript_id)
if (!"is_expressed" %in% colnames(gt)) gt$is_expressed <- gt$base_reads > 0

num_expressed <- sum(gt$is_expressed)
num_total <- nrow(gt)
cat(sprintf("Ground truth: %d transcripts (%d expressed, %d unexpressed)\n\n",
            num_total, num_expressed, num_total - num_expressed))

# ---- Load sample info ----
sample_info <- read.csv(file.path(simdir, "sample_info.csv"), stringsAsFactors = FALSE)

read_quant <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- read.delim(path, stringsAsFactors = FALSE)
  colnames(df) <- c("target_name", "len", "eelen", "tpm", "ecount")
  df
}

# ---- Methods to compare ----
methods <- list(
  "Single EM"      = file.path(simdir, "quant_single"),
  "Hier. Multi"    = file.path(simdir, "quant_multi"),
  "Hier. NoCond"   = file.path(simdir, "quant_hier_nocond"),
  "Group LASSO"    = file.path(simdir, "quant_gl")
)

# ---- Evaluate ----
pseudo <- 0.01
fp_thresh <- 1.0

cat(sprintf("%-14s | %8s %8s %8s | %8s %8s | %6s %6s | %6s\n",
            "Method", "P(log2)", "Spearman", "MARD", "Expr P", "Expr S",
            "FP", "FP%", "Active"))
cat(paste0(rep("-", 95), collapse = ""), "\n")

for (mname in names(methods)) {
  mdir <- methods[[mname]]

  all_pearson <- c(); all_spearman <- c(); all_mard <- c()
  expr_pearson <- c(); expr_spearman <- c()
  all_fp <- 0; all_unexpr <- 0
  n_active <- c()

  for (i in 1:nrow(sample_info)) {
    sname <- sample_info$sample_name[i]
    condition <- sample_info$condition[i]
    qpath <- file.path(mdir, sname, paste0(sname, ".quant"))
    q <- read_quant(qpath)
    if (is.null(q)) next

    tpm_col <- paste0("expected_tpm_", condition)
    m <- match(q$target_name, gt$short_id)
    valid <- !is.na(m)
    est <- q$tpm[valid]
    tru <- gt[[tpm_col]][m[valid]]
    is_expr <- gt$is_expressed[m[valid]]

    all_pearson <- c(all_pearson, cor(log2(est+pseudo), log2(tru+pseudo)))
    all_spearman <- c(all_spearman, cor(est, tru, method="spearman"))
    all_mard <- c(all_mard, median(abs(est-tru)/(tru+pseudo)))

    if (sum(is_expr) > 0) {
      expr_pearson <- c(expr_pearson, cor(log2(est[is_expr]+pseudo), log2(tru[is_expr]+pseudo)))
      expr_spearman <- c(expr_spearman, cor(est[is_expr], tru[is_expr], method="spearman"))
    }

    all_fp <- all_fp + sum(est[!is_expr] > fp_thresh)
    all_unexpr <- all_unexpr + sum(!is_expr)
    n_active <- c(n_active, sum(est > 0))
  }

  cat(sprintf("%-14s | %8.4f %8.4f %8.4f | %8.4f %8.4f | %6d %5.1f%% | %6.0f\n",
              mname,
              mean(all_pearson), mean(all_spearman), mean(all_mard),
              mean(expr_pearson), mean(expr_spearman),
              all_fp, 100*all_fp/all_unexpr,
              mean(n_active)))
}

cat("\n")
