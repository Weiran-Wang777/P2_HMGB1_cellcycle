# -*- coding: utf-8 -*-
# P2 Figure 5 — Concomitant curve: HMGB1 downregulation vs G2M along viral load
# Fixed: simple numeric binning instead of cut() label parsing

set.seed(123)
suppressMessages({library(Seurat); library(ggplot2); library(patchwork)})

base <- "data"
if (!file.exists(file.path(base, "denv_cc.rds")))
  stop("Error: data/denv_cc.rds not found. See README for data download instructions.")
if (!file.exists(file.path(base, "zika_clustered.rds")))
  stop("Error: data/zika_clustered.rds not found. See README for data download instructions.")
denv <- readRDS(file.path(base, "denv_cc.rds"))
zika <- readRDS(file.path(base, "zika_clustered.rds"))
outdir <- "P2_HMGB1/figures/"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

theme_pub <- function(base=12) {
  theme_classic(base_size=base) +
  theme(text=element_text(family="sans"),
        axis.text=element_text(colour="black"),
        axis.ticks=element_line(colour="black"),
        legend.key=element_blank(),
        legend.title=element_text(size=base-2),
        plot.title=element_text(face="bold", size=base+1, hjust=0.5),
        plot.margin=margin(6,6,6,6))
}

# ---- Binned analysis (simple numeric approach) ----
make_binned_df <- function(seu, bin_width=2) {
  df <- FetchData(seu, vars=c("HMGB1","viral_log","G2M.Score"))
  colnames(df) <- c("HMGB1","viral_log","G2M.Score")

  # Assign numeric bin: floor(viral_log / bin_width) * bin_width + bin_width/2
  df$bin_mid <- floor(df$viral_log / bin_width) * bin_width + bin_width/2

  # Aggregate by bin
  agg <- aggregate(cbind(HMGB1, G2M.Score) ~ bin_mid, data=df, FUN=mean)
  agg$n <- as.numeric(table(df$bin_mid)[as.character(agg$bin_mid)])

  # Keep bins with >=20 cells
  agg <- agg[agg$n >= 20, ]
  agg <- agg[order(agg$bin_mid), ]
  rownames(agg) <- NULL
  return(agg)
}

binned_d <- make_binned_df(denv)
binned_z <- make_binned_df(zika)

cat("\n=== Binned summary ===\n")
cat("DENV:\n"); print(binned_d)
cat("\nZIKV:\n"); print(binned_z)

# ---- Dual-axis plot per virus ----
plot_concomitant <- function(binned, virus_label, color_hmg, color_g2m) {
  # Scale G2M to overlay on HMGB1 axis
  hmg_range <- range(binned$HMGB1)
  g2m_range <- range(binned$G2M.Score)
  binned$G2M_scaled <- (binned$G2M.Score - g2m_range[1]) / diff(g2m_range) * diff(hmg_range) + hmg_range[1]

  ggplot(binned, aes(x=bin_mid)) +
    geom_line(aes(y=HMGB1, color="HMGB1"), linewidth=0.9) +
    geom_point(aes(y=HMGB1, color="HMGB1"), size=2.5) +
    geom_line(aes(y=G2M_scaled, color="G2M.Score"), linewidth=0.9, linetype="dashed") +
    geom_point(aes(y=G2M_scaled, color="G2M.Score"), size=2.5, shape=17) +
    scale_color_manual(values=c("HMGB1"=color_hmg, "G2M.Score"=color_g2m),
                       name="Expression / Score") +
    labs(title=virus_label, x="Viral load (log2, binned)", y="Mean normalized expression") +
    theme_pub(12) +
    theme(legend.position="bottom")
}

p5_d <- plot_concomitant(binned_d, "DENV", "#C0392B", "#55A868")
p5_z <- plot_concomitant(binned_z, "ZIKV", "#1A6FA0", "#55A868")

combined <- p5_d | p5_z + plot_annotation(tag_levels='A')
ggsave(paste0(outdir, "Fig5_HMGB1_G2M_concomitant.png"), combined, width=7.5, height=3.75, dpi=600, bg="white")
ggsave(paste0(outdir, "Fig5_HMGB1_G2M_concomitant.pdf"), combined, width=7.5, height=3.75, device=cairo_pdf, bg="white")
cat("\nSaved Fig5 (concomitant HMGB1 down vs G2M along viral load)\n")
