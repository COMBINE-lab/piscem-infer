#!/usr/bin/env Rscript
# Scatter plots of estimated vs true log2 fold change (treatment/control)
# per transcript, for each quantification method.
# Includes "EM+Selection" = single EM filtered to NoCond+AV expressed set.

library(ggplot2)
library(data.table)
library(ggnewscale)

simdir <- "sim_data_gencode"
pseudo <- 0.01

# ---- Load ground truth ----
gt <- fread(file.path(simdir, "ground_truth.csv"))
gt[, short_id := sub(" .*", "", transcript_id)]
if (!"is_expressed" %in% colnames(gt)) gt[, is_expressed := base_reads > 0]

sample_info <- fread(file.path(simdir, "sample_info.csv"))
ctrl_samples  <- sample_info$sample_name[sample_info$condition == "control"]
treat_samples <- sample_info$sample_name[sample_info$condition == "treatment"]

read_quant <- function(path) {
  if (!file.exists(path)) { warning("Not found: ", path); return(NULL) }
  df <- fread(path, sep = "\t")
  setnames(df, c("target_name", "len", "eelen", "tpm", "ecount"))
  df
}

# ---- Get NoCond+AV expressed set (union across all samples) ----
av_expressed <- character(0)
for (i in 1:nrow(sample_info)) {
  sn <- sample_info$sample_name[i]
  qpath <- file.path(simdir, "quant_hier_nocond_av", sn, paste0(sn, ".quant"))
  qdf <- read_quant(qpath)
  if (!is.null(qdf)) av_expressed <- union(av_expressed, qdf$target_name[qdf$tpm > 0])
}

methods <- list(
  list(name = "Single EM",     dir = "quant_em",                filter_set = NULL),
  list(name = "Cons TPM",      dir = "quant_consensus",         filter_set = NULL),
  list(name = "Cons UES",      dir = "quant_consensus_ues",     filter_set = NULL),
  list(name = "Cons Support",  dir = "quant_consensus_support", filter_set = NULL),
  list(name = "NoCond+AV",     dir = "quant_hier_nocond_av",    filter_set = NULL)
)
method_levels <- c("Single EM", "Cons TPM", "Cons UES", "Cons Support", "NoCond+AV")

# ---- True fold changes ----
true_lfc <- log2(gt$expected_tpm_treatment + pseudo) - log2(gt$expected_tpm_control + pseudo)
expressed <- gt$expected_tpm_control > 0 | gt$expected_tpm_treatment > 0

# ---- Helper: read and optionally filter TPM for a sample ----
get_tpm <- function(meth, sname) {
  qpath <- file.path(simdir, meth$dir, sname, paste0(sname, ".quant"))
  qdf <- read_quant(qpath)
  if (is.null(qdf)) return(NULL)

  est_tpm <- qdf$tpm
  if (!is.null(meth$filter_set)) {
    keep <- qdf$target_name %in% meth$filter_set
    est_tpm <- est_tpm * keep
    total <- sum(est_tpm)
    if (total > 0) est_tpm <- est_tpm * 1e6 / total
  }

  list(tpm = est_tpm, names = qdf$target_name)
}

# ---- Compute estimated fold changes per method ----
all_data <- list()

for (meth in methods) {
  ctrl_tpms <- list()
  ref_names <- NULL
  for (sn in ctrl_samples) {
    res <- get_tpm(meth, sn)
    if (!is.null(res)) {
      ctrl_tpms[[sn]] <- res$tpm
      if (is.null(ref_names)) ref_names <- res$names
    }
  }

  treat_tpms <- list()
  for (sn in treat_samples) {
    res <- get_tpm(meth, sn)
    if (!is.null(res)) treat_tpms[[sn]] <- res$tpm
  }

  if (length(ctrl_tpms) == 0 || length(treat_tpms) == 0) next

  mean_ctrl  <- Reduce("+", ctrl_tpms)  / length(ctrl_tpms)
  mean_treat <- Reduce("+", treat_tpms) / length(treat_tpms)
  est_lfc    <- log2(mean_treat + pseudo) - log2(mean_ctrl + pseudo)

  m <- match(ref_names, gt$short_id)
  valid <- !is.na(m) & expressed[m]

  all_data[[meth$name]] <- data.table(
    method    = meth$name,
    true_lfc  = true_lfc[m[valid]],
    est_lfc   = est_lfc[valid],
    is_de     = gt$is_de[m[valid]]
  )
}

df <- rbindlist(all_data)
df[, method := factor(method, levels = method_levels)]

# ---- Stats per method ----
stats <- df[, {
  p_all <- cor(true_lfc, est_lfc, method = "pearson")
  de_big <- is_de & abs(true_lfc) > 1
  if (sum(de_big) > 10) {
    fit <- lm(est_lfc[de_big] ~ true_lfc[de_big])
    slope <- coef(fit)[2]
  } else {
    slope <- NA_real_
  }
  list(pearson = p_all, slope = slope, n_de = sum(de_big))
}, by = method]

stats[, label := sprintf("slope = %.3f\nr = %.3f", slope, pearson)]

# ---- Plot ----
lim <- max(abs(c(df$true_lfc, df$est_lfc))) * 1.05

df_nonde <- df[is_de == FALSE]
df_de    <- df[is_de == TRUE]

p <- ggplot(mapping = aes(x = true_lfc, y = est_lfc)) +
  geom_hex(data = df_nonde, bins = 60, alpha = 0.85) +
  scale_fill_gradientn(
    colours = c("#f7f7f7", "#d9d9d9", "#bdbdbd", "#969696", "#636363", "#252525"),
    trans = "log10",
    name = "non-DE",
    guide = guide_colorbar(barwidth = 5, barheight = 0.5, order = 2)
  ) +
  new_scale_fill() +
  geom_point(data = df_de, colour = "#e31a1c", size = 0.5, alpha = 0.5, shape = 16) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              colour = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_text(
    data = stats,
    aes(x = -lim * 0.95, y = lim * 0.90, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.2,
    colour = "grey20", fontface = "italic", lineheight = 0.9
  ) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  facet_wrap(~ method, nrow = 1) +
  labs(
    x = expression("True" ~ log[2] ~ "fold change (treatment / control)"),
    y = expression("Estimated" ~ log[2] ~ "fold change"),
    title = "Gencode simulation: fold change recovery by method"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 11),
    legend.position  = "bottom",
    legend.box       = "horizontal",
    plot.title       = element_text(size = 13, face = "bold"),
    aspect.ratio     = 1
  )

ggsave(file.path(simdir, "scatter_foldchange.pdf"), p, width = 17, height = 4.8, device = pdf)
cat("Saved:", file.path(simdir, "scatter_foldchange.pdf"), "\n")
