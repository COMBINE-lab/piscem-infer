#!/usr/bin/env Rscript
# Check if ENST00000449131 shares EQ classes with other transcripts
# by looking at correlation of counts across methods

library(data.table)

benchdir <- "airway_benchmark"
all_samples <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520",
                 "SRR1039509", "SRR1039513", "SRR1039517", "SRR1039521")

# Load plain EM and hier cond for sample SRR1039508 (untreated)
em <- fread(file.path(benchdir, "quant/em/SRR1039508.quant"), sep = "\t")
setnames(em, c("target_name", "len", "eelen", "tpm", "ecount"))

hier <- fread(file.path(benchdir, "quant/hier_all_cond/SRR1039508/SRR1039508.quant"), sep = "\t")
setnames(hier, c("target_name", "len", "eelen", "tpm", "ecount"))

# ENST00000449131 is BEST1-202 gene
# Find all BEST1 transcripts
best1_em <- em[grepl("BEST1", target_name)]
best1_hier <- hier[grepl("BEST1", target_name)]

cat("=== All BEST1 transcripts in SRR1039508 (untreated) ===\n\n")
cat(sprintf("%-80s | %10s %10s | %10s %10s\n",
            "Transcript", "EM_TPM", "EM_ecount", "Hier_TPM", "Hier_ecount"))
cat(paste0(rep("-", 130), collapse = ""), "\n")

both <- merge(em[grepl("BEST1", target_name), .(target_name, em_tpm=tpm, em_ec=ecount)],
              hier[grepl("BEST1", target_name), .(target_name, h_tpm=tpm, h_ec=ecount)],
              by = "target_name")
for (i in 1:nrow(both)) {
  cat(sprintf("%-80s | %10.2f %10.1f | %10.2f %10.1f\n",
              substr(both$target_name[i], 1, 80),
              both$em_tpm[i], both$em_ec[i],
              both$h_tpm[i], both$h_ec[i]))
}

# Also look at the gene: ENSG00000167995
cat("\n=== All transcripts of gene ENSG00000167995 ===\n\n")
gene_em <- em[grepl("ENSG00000167995", target_name)]
gene_hier <- hier[grepl("ENSG00000167995", target_name)]

both2 <- merge(gene_em[, .(target_name, em_tpm=tpm, em_ec=ecount)],
               gene_hier[, .(target_name, h_tpm=tpm, h_ec=ecount)],
               by = "target_name")
for (i in 1:nrow(both2)) {
  cat(sprintf("%-80s | %10.2f %10.1f | %10.2f %10.1f\n",
              substr(both2$target_name[i], 1, 80),
              both2$em_tpm[i], both2$em_ec[i],
              both2$h_tpm[i], both2$h_ec[i]))
}

# Now compare the BIGGEST changes between EM and Hier for SRR1039508
cat("\n=== Top 20 transcripts with largest absolute TPM increase (EM -> Hier, SRR1039508) ===\n")
both_all <- merge(em[, .(target_name, em_tpm=tpm, em_ec=ecount)],
                  hier[, .(target_name, h_tpm=tpm, h_ec=ecount)],
                  by = "target_name")
both_all[, diff := h_tpm - em_tpm]
setorder(both_all, -diff)

cat(sprintf("\n%-80s | %10s %10s | %10s\n",
            "Transcript", "EM_TPM", "Hier_TPM", "Diff"))
cat(paste0(rep("-", 115), collapse = ""), "\n")
for (i in 1:20) {
  cat(sprintf("%-80s | %10.2f %10.2f | %10.2f\n",
              substr(both_all$target_name[i], 1, 80),
              both_all$em_tpm[i], both_all$h_tpm[i], both_all$diff[i]))
}

cat("\nDone.\n")
