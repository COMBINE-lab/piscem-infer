#!/usr/bin/env Rscript
library(data.table)
simdir <- "sim_data_gencode"; pseudo <- 0.01
gt <- fread(file.path(simdir, "ground_truth.csv"))
gt[, short_id := sub(" .*", "", transcript_id)]
sample_info <- fread(file.path(simdir, "sample_info.csv"))
ctrl_samples <- sample_info$sample_name[sample_info$condition == "control"]
treat_samples <- sample_info$sample_name[sample_info$condition == "treatment"]
read_quant <- function(path) { if (!file.exists(path)) return(NULL); df <- fread(path, sep="\t"); setnames(df, c("target_name","len","eelen","tpm","ecount")); df }

methods <- list(
  list(name="Single EM",      dir="quant_em"),
  list(name="Sel+Supp",       dir="quant_consensus_sel_support"),
  list(name="Sel+Supp Resc",  dir="quant_consensus_sel_support_rescue"),
  list(name="Sel+Supp Adapt", dir="quant_consensus_sel_support_adaptive"),
  list(name="Adpt+Rescue",    dir="quant_consensus_sel_support_adapt_rescue"),
  list(name="Soft Rescue",    dir="quant_consensus_sel_support_soft_rescue"),
  list(name="Adpt+GeneResc", dir="quant_consensus_sel_adapt_generescue"),
  list(name="Adpt+Pos5",    dir="quant_consensus_sel_adapt_pos5"),
  list(name="Pos5+CG",     dir="quant_consensus_sel_pos5_condgene")
)

cat(sprintf("%-16s | %6s %6s %6s | %8s %8s %8s | %8s %8s %8s %8s\n",
            "Method","TP","FP","FN","Prec","Recall","F1","Pearson","Spearman","RMSE","FC slope"))
cat(paste0(rep("-",120),collapse=""),"\n")

for (meth in methods) {
  all_p <- c(); all_s <- c(); all_r <- c(); all_tp <- c(); all_fp <- c(); all_fn <- c()
  for (i in 1:nrow(sample_info)) {
    sn <- sample_info$sample_name[i]; cond <- sample_info$condition[i]
    qdf <- read_quant(file.path(simdir,meth$dir,sn,paste0(sn,".quant")))
    if (is.null(qdf)) next
    m <- match(qdf$target_name, gt$short_id); valid <- !is.na(m)
    true_t <- gt[[paste0("expected_tpm_",cond)]][m[valid]]; est <- qdf$tpm[valid]
    all_p <- c(all_p, cor(log2(est+pseudo),log2(true_t+pseudo),method="pearson"))
    all_s <- c(all_s, cor(est,true_t,method="spearman"))
    all_r <- c(all_r, sqrt(mean((log2(est+1)-log2(true_t+1))^2)))
    all_tp <- c(all_tp, sum(est>0 & true_t>0))
    all_fp <- c(all_fp, sum(est>0 & true_t==0))
    all_fn <- c(all_fn, sum(est==0 & true_t>0))
  }
  tp <- mean(all_tp); fp <- mean(all_fp); fn_ <- mean(all_fn)
  prec <- tp/(tp+fp); rec <- tp/(tp+fn_); f1 <- 2*prec*rec/(prec+rec)
  ctrl_t <- list(); treat_t <- list()
  for (sn in ctrl_samples) { qdf <- read_quant(file.path(simdir,meth$dir,sn,paste0(sn,".quant"))); if (!is.null(qdf)) ctrl_t[[sn]] <- qdf$tpm }
  for (sn in treat_samples) { qdf <- read_quant(file.path(simdir,meth$dir,sn,paste0(sn,".quant"))); if (!is.null(qdf)) treat_t[[sn]] <- qdf$tpm }
  slope_str <- "   N/A"
  if (length(ctrl_t)>0 && length(treat_t)>0) {
    mc <- Reduce("+",ctrl_t)/length(ctrl_t); mt <- Reduce("+",treat_t)/length(treat_t)
    elfc <- log2(mt+pseudo)-log2(mc+pseudo)
    rn <- read_quant(file.path(simdir,meth$dir,ctrl_samples[1],paste0(ctrl_samples[1],".quant")))$target_name
    m <- match(rn,gt$short_id); v <- !is.na(m)
    tlfc <- log2(gt$expected_tpm_treatment[m[v]]+pseudo)-log2(gt$expected_tpm_control[m[v]]+pseudo)
    de <- rn[v] %in% gt$short_id[gt$is_de]; big <- abs(tlfc)>1
    if (sum(de&big)>10) { fit <- lm(elfc[v][de&big]~tlfc[de&big]); slope_str <- sprintf("%6.3f",coef(fit)[2]) }
  }
  cat(sprintf("%-16s | %6d %6d %6d | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f %8s\n",
              meth$name, round(tp), round(fp), round(fn_), prec, rec, f1,
              mean(all_p), mean(all_s), mean(all_r), slope_str))
}
