#!/usr/bin/env Rscript
# Log scatter plots of predicted vs true TPM per transcript
# using log2(1+x) transform to include zeros.
# Top row: scatter; Bottom row: MA (Bland-Altman) plot showing residuals vs mean.
# Includes "EM+Selection" = single EM filtered to NoCond+AV expressed set.

library(ggplot2)
library(data.table)
library(hexbin)
library(ggnewscale)
library(patchwork)

simdir <- "sim_data_gencode"

# ---- Load ground truth ----
gt <- fread(file.path(simdir, "ground_truth.csv"))
gt[, short_id := sub(" .*", "", transcript_id)]

sample_info <- fread(file.path(simdir, "sample_info.csv"))

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
cat(sprintf("NoCond+AV expressed set: %d transcripts\n", length(av_expressed)))

methods <- list(
  list(name = "Single EM",     dir = "quant_em",                subdir = TRUE),
  list(name = "Cons TPM",      dir = "quant_consensus",         subdir = TRUE),
  list(name = "Cons UES",      dir = "quant_consensus_ues",     subdir = TRUE),
  list(name = "Cons Support",  dir = "quant_consensus_support", subdir = TRUE),
  list(name = "NoCond+AV",     dir = "quant_hier_nocond_av",    subdir = TRUE)
)

method_levels <- c("Single EM", "Cons TPM", "Cons UES", "Cons Support", "NoCond+AV")

# ---- Helper: build data for one sample ----
build_sample_data <- function(meth, sname, condition) {
  tpm_col <- paste0("expected_tpm_", condition)

  qpath <- file.path(simdir, meth$dir, sname, paste0(sname, ".quant"))
  qdf <- read_quant(qpath)
  if (is.null(qdf)) return(NULL)

  est_tpm <- qdf$tpm

  # Apply selection filter if specified
  if (!is.null(meth$filter_set)) {
    keep <- qdf$target_name %in% meth$filter_set
    est_tpm <- est_tpm * keep
    total <- sum(est_tpm)
    if (total > 0) est_tpm <- est_tpm * 1e6 / total
  }

  m <- match(qdf$target_name, gt$short_id)
  valid <- !is.na(m)

  data.table(
    method   = meth$name,
    true_tpm = gt[[tpm_col]][m[valid]],
    est_tpm  = est_tpm[valid],
    is_de    = gt$is_de[m[valid]]
  )
}

# log2(1+x) transform for ggplot
log2p1_trans <- scales::trans_new(
  name = "log2p1",
  transform = function(x) log2(x + 1),
  inverse = function(x) 2^x - 1,
  breaks = function(x) c(0, 1, 10, 100, 1000, 10000, 1e5, 1e6)
)

# ---- Composite plot function ----
make_composite <- function(df, title_str, outfile, width = 20, height = 8) {
  df[, method := factor(method, levels = method_levels)]

  # Compute both correlations
  pseudo <- 0.01
  stats <- df[, {
    p <- cor(log2(true_tpm + pseudo), log2(est_tpm + pseudo), method = "pearson")
    s <- cor(true_tpm, est_tpm, method = "spearman")
    list(pearson = p, spearman = s)
  }, by = method]
  stats[, label := paste0("r = ", sprintf("%.3f", pearson), "\nrho = ", sprintf("%.3f", spearman))]

  # Split non-DE and DE
  df_nonde <- df[is_de == FALSE]
  df_de    <- df[is_de == TRUE]

  # Axis limits
  max_val <- max(c(df$true_tpm, df$est_tpm))
  ax_max <- max_val * 1.2
  ax_breaks <- c(0, 1, 10, 100, 1000, 10000, 1e5)
  ax_breaks <- ax_breaks[ax_breaks <= ax_max]
  ax_labels <- function(x) {
    ifelse(x == 0, "0",
      ifelse(x < 1000, as.character(round(x)),
        scales::label_comma()(x)))
  }

  ann_x <- 2^(0.3) - 1
  ann_y <- 2^(log2(ax_max + 1) - 0.3) - 1

  # ======== Top row: scatter plots ========
  p_scatter <- ggplot(mapping = aes(x = true_tpm, y = est_tpm)) +
    geom_hex(data = df_nonde, bins = 70, alpha = 0.9) +
    scale_fill_gradientn(
      colours = c("#f7f7f7", "#d9d9d9", "#bdbdbd", "#969696", "#636363", "#252525"),
      trans = "log10",
      name = "non-DE",
      guide = guide_colorbar(barwidth = 5, barheight = 0.5)
    ) +
    new_scale_fill() +
    geom_point(data = df_de, colour = "#e31a1c", size = 0.4, alpha = 0.5, shape = 16) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                colour = "grey40", linewidth = 0.4) +
    scale_x_continuous(trans = log2p1_trans, breaks = ax_breaks, labels = ax_labels,
                       limits = c(0, ax_max)) +
    scale_y_continuous(trans = log2p1_trans, breaks = ax_breaks, labels = ax_labels,
                       limits = c(0, ax_max)) +
    geom_text(
      data = stats,
      aes(x = ann_x, y = ann_y, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.0,
      colour = "grey20", fontface = "italic", lineheight = 0.9
    ) +
    facet_wrap(~ method, nrow = 1) +
    labs(x = "True TPM", y = "Estimated TPM") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold", size = 11),
      legend.position  = "none",
      aspect.ratio     = 1
    )

  # ======== Bottom row: MA plot ========
  df[, M := log2(est_tpm + 1) - log2(true_tpm + 1)]
  df[, A := (log2(est_tpm + 1) + log2(true_tpm + 1)) / 2]

  df_nonde_ma <- df[is_de == FALSE]
  df_de_ma    <- df[is_de == TRUE]

  ma_stats <- df[, {
    rmse <- sqrt(mean(M^2))
    list(rmse = rmse)
  }, by = method]
  ma_stats[, label := sprintf("RMSE = %.2f", rmse)]

  a_max <- max(df$A) + 0.5
  m_lim <- max(abs(df$M[is.finite(df$M)])) * 0.6
  m_lim <- min(m_lim, 15)

  p_ma <- ggplot(mapping = aes(x = A, y = M)) +
    geom_hex(data = df_nonde_ma, bins = 70, alpha = 0.9) +
    scale_fill_gradientn(
      colours = c("#f7f7f7", "#d9d9d9", "#bdbdbd", "#969696", "#636363", "#252525"),
      trans = "log10",
      name = "non-DE",
      guide = guide_colorbar(barwidth = 5, barheight = 0.5)
    ) +
    new_scale_fill() +
    geom_point(data = df_de_ma, colour = "#e31a1c", size = 0.4, alpha = 0.5, shape = 16) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    scale_x_continuous(limits = c(0, a_max)) +
    scale_y_continuous(limits = c(-m_lim, m_lim)) +
    geom_text(
      data = ma_stats,
      aes(x = 0.3, y = m_lim * 0.9, label = label),
      inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.0,
      colour = "grey20", fontface = "italic"
    ) +
    facet_wrap(~ method, nrow = 1) +
    labs(
      x = expression("A = mean " * log[2] * "(TPM + 1)"),
      y = expression("M = " * log[2] * "(est + 1) - " * log[2] * "(true + 1)")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text       = element_blank(),
      legend.position  = "none",
      aspect.ratio     = 1
    )

  combined <- p_scatter / p_ma +
    plot_annotation(title = title_str,
                    theme = theme(plot.title = element_text(size = 13, face = "bold")))

  ggsave(outfile, combined, width = width, height = height, device = pdf)
  cat("Saved:", outfile, "\n")
}

# ======================================================================
# Plot 1: single sample (sample_01, control)
# ======================================================================
cat("Building single-sample data...\n")
all_data <- list()
for (meth in methods) {
  d <- build_sample_data(meth, "sample_01", "control")
  if (!is.null(d)) all_data[[meth$name]] <- d
}
df1 <- rbindlist(all_data)

make_composite(
  df1,
  "Gencode simulation: estimated vs. true TPM (sample_01, control)",
  file.path(simdir, "scatter_pred_vs_true.pdf")
)

# ======================================================================
# Plot 2: all samples pooled
# ======================================================================
cat("Building all-sample data...\n")
all_pooled <- list()
for (meth in methods) {
  parts <- list()
  for (i in 1:nrow(sample_info)) {
    d <- build_sample_data(meth, sample_info$sample_name[i], sample_info$condition[i])
    if (!is.null(d)) parts[[length(parts) + 1]] <- d
  }
  if (length(parts) > 0) all_pooled[[meth$name]] <- rbindlist(parts)
}
df2 <- rbindlist(all_pooled)

make_composite(
  df2,
  "Gencode simulation: estimated vs. true TPM (all 6 samples pooled)",
  file.path(simdir, "scatter_pred_vs_true_allsamples.pdf")
)
