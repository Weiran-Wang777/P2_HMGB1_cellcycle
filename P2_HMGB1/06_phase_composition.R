# -*- coding: utf-8 -*-
# P2 Figure 4 — Combined: A=HMGB1 by cell cycle phase (boxplots), B=phase composition by infection status
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
cols_phase <- c("G1"="#E8A33D", "S"="#4C72B0", "G2M"="#55A868")

# ========== Panel A: Boxplots — HMGB1 by cell cycle phase ==========
expr_d <- FetchData(denv, vars=c("HMGB1","Phase"))
expr_z <- FetchData(zika, vars=c("HMGB1","Phase"))
expr_d$Phase <- factor(expr_d$Phase, levels=c("G1","S","G2M"))
expr_z$Phase <- factor(expr_z$Phase, levels=c("G1","S","G2M"))
pA_d <- ggplot(expr_d, aes(Phase, HMGB1, fill=Phase)) +
  geom_boxplot(outlier.size=0.6, outlier.alpha=0.4, width=0.6) +
  scale_fill_manual(values=cols_phase, guide="none") +
  labs(title="DENV", x=NULL, y="HMGB1 expression") +
  theme_pub(13)
pA_z <- ggplot(expr_z, aes(Phase, HMGB1, fill=Phase)) +
  geom_boxplot(outlier.size=0.6, outlier.alpha=0.4, width=0.6) +
  scale_fill_manual(values=cols_phase, guide="none") +
  labs(title="ZIKV", x=NULL, y="HMGB1 expression") +
  theme_pub(13)
panel_A <- pA_d | pA_z

# ========== Panel B: Stacked bar — phase composition by infection status ==========
get_comp <- function(seu, virus_name) {
  meta <- seu@meta.data
  meta$group <- ifelse(meta$viral_log > 0.5, "Infected", "Uninfected")
  tab <- as.data.frame(prop.table(table(meta$group, meta$Phase), margin=1) * 100)
  colnames(tab) <- c("Group", "Phase", "Pct")
  tab$Virus <- virus_name
  tab
}
denv_comp <- get_comp(denv, "DENV")
zika_comp <- get_comp(zika, "ZIKV")
comp <- rbind(denv_comp, zika_comp)
comp$Phase <- factor(comp$Phase, levels=c("G1", "S", "G2M"))
comp$Group <- factor(comp$Group, levels=c("Uninfected", "Infected"))
comp$Virus <- factor(comp$Virus, levels=c("DENV", "ZIKV"))
panel_B <- ggplot(comp, aes(x=Group, y=Pct, fill=Phase)) +
  geom_bar(stat="identity", width=0.7) +
  facet_grid(~Virus) +
  scale_fill_manual(values=cols_phase) +
  geom_text(aes(label=sprintf("%.1f%%", Pct)), position=position_stack(vjust=0.5), size=2.5, color="white") +
  labs(x="", y="Proportion of cells (%)", fill="Cell cycle phase") +
  theme_classic(base_size=12) +
  theme(
    axis.text.x=element_text(angle=20, hjust=1),
    strip.text=element_text(face="bold", size=13),
    legend.position="bottom",
    panel.grid.major.x=element_blank(),
    plot.margin=margin(6,6,6,6)
  )

# ========== Combine: A on top, B below, with tags ==========
combined <- (panel_A / panel_B) + plot_annotation(tag_levels='A')
ggsave(paste0(outdir, "Fig4_HMGB1_cellcycle.png"), combined, width=7.5, height=6.7, dpi=600, bg="white")
ggsave(paste0(outdir, "Fig4_HMGB1_cellcycle.pdf"), combined, width=7.5, height=6.7, device=cairo_pdf, bg="white")
cat("Saved combined Fig4 (A=boxplots, B=phase composition by infection)\n")
