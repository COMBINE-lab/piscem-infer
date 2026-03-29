#!/usr/bin/env Rscript
# Comprehensive GENCODE evaluation: all methods including salmon and kallisto.
library(data.table)
simdir <- "sim_data_gencode"; pseudo <- 0.01
gt <- fread(file.path(simdir, "ground_truth.csv"))
gt[, short_id := sub(" .*", "", transcript_id)]
sample_info <- fread(file.path(simdir, "sample_info.csv"))
ctrl_samples <- sample_info$sample_name[sample_info$condition == "control"]
treat_samples <- sample_info$sample_name[sample_info$condition == "treatment"]

read_piscem <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- fread(path, sep="\t"); setnames(df, c("target_name","len","eelen","tpm","ecount")); df
}
read_salmon <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- fread(path, sep="\t"); df[, .(target_name=Name, len=Length, eelen=EffectiveLength, tpm=TPM, ecount=NumReads)]
}
read_kallisto <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- fread(path, sep="\t"); df[, .(target_name=target_id, len=length, eelen=eff_length, tpm, ecount=est_counts)]
}

methods <- list(
  list(name="Single EM",  reader=read_piscem,  fn=function(sn) file.path(simdir,"quant_em",sn,paste0(sn,".quant"))),
  list(name="Salmon",     reader=read_salmon,  fn=function(sn) file.path(simdir,"quant_salmon",sn,"quant.sf")),
  list(name="Salmon EM",  reader=read_salmon,  fn=function(sn) file.path(simdir,"quant_salmon_em",sn,"quant.sf")),
  list(name="Kallisto",   reader=read_kallisto, fn=function(sn) file.path(simdir,"quant_kallisto",sn,"abundance.tsv")),
  list(name="Sel+Adapt",  reader=read_piscem,  fn=function(sn) file.path(simdir,"quant_consensus_sel_support_adaptive",sn,paste0(sn,".quant"))),
  list(name="Sel+Pos5",   reader=read_piscem,  fn=function(sn) file.path(simdir,"quant_consensus_sel_adapt_pos5",sn,paste0(sn,".quant"))),
  list(name="Pos5+CG",    reader=read_piscem,  fn=function(sn) file.path(simdir,"quant_consensus_sel_pos5_condgene",sn,paste0(sn,".quant")))
)

cat(sprintf("%-14s | %6s %6s %6s | %8s %8s %8s | %8s %8s %8s %8s\n",
            "Method","TP","FP","FN","Prec","Recall","F1","Pearson","Spearman","RMSE","FC slope"))
cat(paste0(rep("-",120),collapse=""),"\n")

for (meth in methods) {
  all_p <- c(); all_s <- c(); all_r <- c(); all_tp <- c(); all_fp <- c(); all_fn <- c()
  for (i in 1:nrow(sample_info)) {
    sn <- sample_info$sample_name[i]; cond <- sample_info$condition[i]
    qdf <- meth$reader(meth$fn(sn))
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
  for (sn in ctrl_samples) { qdf <- meth$reader(meth$fn(sn)); if (!is.null(qdf)) ctrl_t[[sn]] <- qdf$tpm }
  for (sn in treat_samples) { qdf <- meth$reader(meth$fn(sn)); if (!is.null(qdf)) treat_t[[sn]] <- qdf$tpm }
  slope_str <- "   N/A"
  if (length(ctrl_t)>0 && length(treat_t)>0) {
    mc <- Reduce("+",ctrl_t)/length(ctrl_t); mt <- Reduce("+",treat_t)/length(treat_t)
    elfc <- log2(mt+pseudo)-log2(mc+pseudo)
    rn <- meth$reader(meth$fn(ctrl_samples[1]))$target_name
    m <- match(rn,gt$short_id); v <- !is.na(m)
    tlfc <- log2(gt$expected_tpm_treatment[m[v]]+pseudo)-log2(gt$expected_tpm_control[m[v]]+pseudo)
    de <- rn[v] %in% gt$short_id[gt$is_de]; big <- abs(tlfc)>1
    if (sum(de&big)>10) { fit <- lm(elfc[v][de&big]~tlfc[de&big]); slope_str <- sprintf("%6.3f",coef(fit)[2]) }
  }
  cat(sprintf("%-14s | %6d %6d %6d | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f %8s\n",
              meth$name, round(tp), round(fp), round(fn_), prec, rec, f1,
              mean(all_p), mean(all_s), mean(all_r), slope_str))
}
