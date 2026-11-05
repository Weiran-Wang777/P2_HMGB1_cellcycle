# -*- coding: utf-8 -*-
# P2 figures — New main line: HMGB1 & cell-cycle decomposition
# Fig1: DAMP family (HMGB1/HMGB2/HMGN1/HMGN2) vs viral load gradient
# Fig2: DAMP genes by cell cycle phase (violin plots)
# Fig3: HMGB1 correlation with cell-cycle / proliferation genes

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

# ---- Auto-detect high-viral cluster (for reference, not used as grouping) ----
find_high_viral <- function(seu) {
  viral_by_cluster <- aggregate(seu$viral_log, by = list(cluster = seu$seurat_clusters),
                                FUN = median, na.rm = TRUE)
  hc <- viral_by_cluster$cluster[which.max(viral_by_cluster$x)]
  message("High-viral cluster: ", hc, " (median viral_log = ", round(max(viral_by_cluster$x), 2), ")")
  return(hc)
}
hc_d <- find_high_viral(denv)
hc_z <- find_high_viral(zika)

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
save_fig <- function(p, name, w, h) {
  # Cap width at 7.5 inches (190mm = journal double-column max)
  if (w > 7.5) { h <- h * (7.5 / w); w <- 7.5 }
  ggsave(paste0(outdir, name, ".png"), p, width=w, height=h, dpi=600, bg="white")
  ggsave(paste0(outdir, name, ".pdf"), p, width=w, height=h, device=cairo_pdf, bg="white")
  cat("Saved", name, " (", round(w*25.4), "x", round(h*25.4), "mm)\n", sep="")
}
cols_virus <- c("DENV"="#C0392B", "ZIKV"="#1A6FA0")
cols_phase <- c("G1"="#E8A33D", "S"="#4C72B0", "G2M"="#55A868")

# ---- Gene sets ----
damp_genes <- c("HMGB1","HMGB2","HMGN1","HMGN2")
cellcycle_genes <- c("MKI67","PCNA","TOP2A","CCNB1","CDK1","CCNA2",
                     "CDKN1A","GADD45A","WEE1")

# ================= Fig1: DAMP family vs viral load gradient =================
# 2 rows (DENV top, ZIKV bottom) × 4 columns (HMGB1/HMGB2/HMGN1/HMGN2)
plot_damp_gene <- function(seu, gene, virus_name, color) {
  df <- FetchData(seu, vars=c(gene, "viral_log"))
  colnames(df) <- c("expr", "viral_log")
  ggplot(df, aes(viral_log, expr)) +
    geom_point(alpha=0.25, size=0.8, colour=color) +
    geom_smooth(method="lm", colour="black", linewidth=0.6, se=TRUE, fill="grey70", alpha=0.3) +
    labs(title=gene, x="log2(viral RNA)", y=NULL) +
    theme_pub(11)
}

fig1_list <- list()
for (i in seq_along(damp_genes)) {
  g <- damp_genes[i]
  fig1_list[[paste0("d",i)]] <- plot_damp_gene(denv, g, "DENV", "#C0392B")
  fig1_list[[paste0("z",i)]] <- plot_damp_gene(zika, g, "ZIKV", "#1A6FA0")
}

# Assemble: row 1 = DENV, row 2 = ZIKV
fig1_top <- (fig1_list$d1 | fig1_list$d2 | fig1_list$d3 | fig1_list$d4)
fig1_bot <- (fig1_list$z1 | fig1_list$z2 | fig1_list$z3 | fig1_list$z4)
save_fig(fig1_top / fig1_bot + plot_annotation(tag_levels='A'), "Fig1_DAMP_viral_gradient", 12, 7)

# ================= Fig2: DAMP genes by cell cycle phase (violin) =================
# Each row = one DAMP gene, columns = DENV / ZIKV
plot_violin <- function(seu, gene, virus_label, color) {
  df <- FetchData(seu, vars=c(gene, "Phase"))
  colnames(df) <- c("expr", "Phase")
  df$Phase <- factor(df$Phase, levels=c("G1","S","G2M"))
  ggplot(df, aes(Phase, expr, fill=Phase)) +
    geom_violin(scale="width", width=0.8, alpha=0.8) +
    geom_boxplot(width=0.15, outlier.size=0.3, fill="white", alpha=0.7) +
    scale_fill_manual(values=cols_phase, guide="none") +
    labs(title=virus_label, x=NULL, y=gene) +
    theme_pub(11)
}

fig2_list <- list()
for (g in damp_genes) {
  fig2_list[[paste0(g,"_d")]] <- plot_violin(denv, g, "DENV", "#C0392B")
  fig2_list[[paste0(g,"_z")]] <- plot_violin(zika, g, "ZIKV", "#1A6FA0")
}

fig2_panel <- (fig2_list$HMGB1_d | fig2_list$HMGB1_z) /
              (fig2_list$HMGB2_d | fig2_list$HMGB2_z) /
              (fig2_list$HMGN1_d | fig2_list$HMGN1_z) /
              (fig2_list$HMGN2_d | fig2_list$HMGN2_z)
save_fig(fig2_panel + plot_annotation(tag_levels='A'), "Fig2_DAMP_cellcycle_violin", 8, 12)

# ================= Fig3: HMGB1 vs cell-cycle gene correlation =================
# ZIKV cells only — Spearman correlation of HMGB1 with proliferation / arrest genes
cc_genes_present <- intersect(cellcycle_genes, rownames(zika))
expr_cc <- FetchData(zika, vars=c("HMGB1", cc_genes_present))
cor_mat <- cor(expr_cc, method="spearman")
cor_vec <- cor_mat["HMGB1", cc_genes_present]
cor_df <- data.frame(gene = names(cor_vec), r = unname(cor_vec))
cor_df <- cor_df[order(cor_df$r), ]
cor_df$gene <- factor(cor_df$gene, levels=cor_df$gene)
cor_df$group <- ifelse(cor_df$r > 0, "Positive", "Negative")

p3 <- ggplot(cor_df, aes(gene, r, fill=group)) +
  geom_col(width=0.7) +
  geom_hline(yintercept=0, colour="black", linewidth=0.4) +
  scale_fill_manual(values=c("Positive"="#55A868", "Negative"="#C0392B"), guide="none") +
  coord_flip() +
  labs(x=NULL, y="Spearman r with HMGB1",
       title="ZIKV: HMGB1 vs cell-cycle / proliferation genes") +
  theme_pub(12) +
  theme(axis.text.y=element_text(size=10))
save_fig(p3, "Fig3_HMGB1_cellcycle_correlation", 7, 5.5)

cat("ALL P2 FIGURES REGENERATED (new main line)\n")
