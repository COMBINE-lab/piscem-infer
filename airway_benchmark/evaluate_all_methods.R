#!/usr/bin/env Rscript
library(data.table)
benchdir <- "airway_benchmark"
samples_u <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")
samples_t <- c("SRR1039509", "SRR1039513", "SRR1039517", "SRR1039521")

read_piscem <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- fread(path, sep="\t"); setnames(df, c("target_name","len","eelen","tpm","ecount"))
  df[, .(target_name, tpm)]
}
read_salmon <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- fread(path, sep="\t"); df[, .(target_name = Name, tpm = TPM)]
}
read_kallisto <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- fread(path, sep="\t"); df[, .(target_name = target_id, tpm)]
}

methods <- list(
  "Plain EM"     = list(reader=read_piscem, fn=function(s) file.path(benchdir,"quant/em",paste0(s,".quant"))),
  "Salmon"       = list(reader=read_salmon, fn=function(s) file.path(benchdir,"quant/salmon",s,"quant.sf")),
  "Salmon EM"    = list(reader=read_salmon, fn=function(s) file.path(benchdir,"quant/salmon_em",s,"quant.sf")),
  "Kallisto"     = list(reader=read_kallisto, fn=function(s) file.path(benchdir,"quant/kallisto",s,"abundance.tsv")),
  "Sel+Adapt"    = list(reader=read_piscem, fn=function(s) file.path(benchdir,"quant/sel_support_all8_adaptive",s,paste0(s,".quant"))),
  "Sel+Pos5"     = list(reader=read_piscem, fn=function(s) file.path(benchdir,"quant/sel_support_all8_pos5",s,paste0(s,".quant"))),
  "Pos5+CG"      = list(reader=read_piscem, fn=function(s) file.path(benchdir,"quant/sel_pos5_condgene",s,paste0(s,".quant")))
)

load_method <- function(meth, sample_list) {
  dfs <- list()
  for (s in sample_list) {
    df <- meth$reader(meth$fn(s))
    if (!is.null(df)) { dfs[[s]] <- df; setnames(dfs[[s]], "tpm", s) }
  }
  if (length(dfs) == 0) return(NULL)
  Reduce(function(a,b) merge(a,b,by="target_name",all=TRUE), dfs)
}

compute_cv <- function(merged, subs) {
  mat <- as.matrix(merged[, ..subs])
  rm <- rowMeans(mat, na.rm=TRUE); rs <- apply(mat,1,sd,na.rm=TRUE)
  data.table(target_name=merged$target_name, mean_tpm=rm, cv=rs/(rm+1e-6))
}

parse_gene <- function(tn) sapply(strsplit(tn,"\\|"), function(x) if(length(x)>=6) x[6] else NA)

# Concordance
for (cond in c("Untreated","Treated")) {
  subs <- if(cond=="Untreated") samples_u else samples_t
  cat(sprintf("\n=== %s Replicate Concordance (CV, mean TPM >= 1) ===\n\n", cond))
  cat(sprintf("%-14s | %8s %8s %8s | %8s\n","Method","Med CV","Mean CV","CV>1","N expr"))
  cat(paste0(rep("-",58),collapse=""),"\n")
  for (name in names(methods)) {
    m <- load_method(methods[[name]], c(samples_u, samples_t))
    if (is.null(m)) { cat(sprintf("%-14s | MISSING\n",name)); next }
    cv <- compute_cv(m, subs); e <- cv[mean_tpm>=1]
    cat(sprintf("%-14s | %8.4f %8.4f %8d | %8d\n",name,
                median(e$cv),mean(e$cv),sum(e$cv>1),nrow(e)))
  }
}

# Gene complexity stratification (untreated)
ref <- load_method(methods[[1]], c(samples_u,samples_t))
ref$gene_name <- parse_gene(ref$target_name)
gene_iso_count <- ref[, .N, by=gene_name]; setnames(gene_iso_count,"N","n_isoforms")
tiers <- list("1 iso"=c(1,1),"2-5"=c(2,5),"6-10"=c(6,10),"11-20"=c(11,20),"21-50"=c(21,50),"51+"=c(51,9999))

cat("\n=== Untreated: Stratified by gene complexity (median CV, TPM >= 1) ===\n\n")
cat(sprintf("%-14s","")); for(t in names(tiers)) cat(sprintf(" | %8s",t)); cat("\n")
cat(paste0(rep("-",14+6*11),collapse=""),"\n")
for (name in names(methods)) {
  m <- load_method(methods[[name]], c(samples_u,samples_t))
  if (is.null(m)) next
  cv <- compute_cv(m, samples_u)
  cv$gene_name <- parse_gene(cv$target_name)
  cv <- merge(cv, gene_iso_count, by="gene_name")
  cat(sprintf("%-14s",name))
  for (tn in names(tiers)) {
    bounds <- tiers[[tn]]
    tier <- cv[n_isoforms>=bounds[1] & n_isoforms<=bounds[2] & mean_tpm>=1]
    cat(sprintf(" | %8.4f",median(tier$cv)))
  }
  cat("\n")
}

# Between-condition FC
cat("\n=== Between-condition fold change (dex vs untreated) ===\n\n")
pseudo <- 0.01
cat(sprintf("%-14s | %10s %10s | %8s\n","Method","Mean|LFC|","SD(LFC)","N both"))
cat(paste0(rep("-",50),collapse=""),"\n")
for (name in names(methods)) {
  m <- load_method(methods[[name]], c(samples_u,samples_t))
  if (is.null(m)) next
  u_mean <- rowMeans(as.matrix(m[,..samples_u]),na.rm=TRUE)
  t_mean <- rowMeans(as.matrix(m[,..samples_t]),na.rm=TRUE)
  both <- u_mean>=1 & t_mean>=1
  if (sum(both)>0) {
    lfc <- log2(t_mean[both]+pseudo) - log2(u_mean[both]+pseudo)
    cat(sprintf("%-14s | %10.4f %10.4f | %8d\n",name,mean(abs(lfc)),sd(lfc),sum(both)))
  }
}
cat("\nDone.\n")
