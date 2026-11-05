# ============================================================
# P2 Supplementary Tables export
# Run from P2_HMGB1_cellcycle_code/ root
# ============================================================
set.seed(123)
suppressMessages(library(Seurat))

outdir <- "P2_HMGB1/supplementary"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

denv <- readRDS("data/denv_cc.rds")
zika <- readRDS("data/zika_clustered.rds")

find_high_viral <- function(seu) {
  vbc <- aggregate(seu$viral_log, by = list(cluster = seu$seurat_clusters), FUN = median, na.rm = TRUE)
  return(list(cluster=vbc$cluster[which.max(vbc$x)], median=max(vbc$x)))
}
hc_d <- find_high_viral(denv)
hc_z <- find_high_viral(zika)
cat(sprintf("DENV high-viral: cluster %s (median=%.2f)\n", hc_d$cluster, hc_d$median))
cat(sprintf("ZIKV high-viral: cluster %s (median=%.2f)\n", hc_z$cluster, hc_z$median))

# ================= S1: HMGB1 correlation with all tested genes =================
damp_genes <- c("HMGB1","HMGB2","HMGN1","HMGN2")
cc_genes <- c("MKI67","PCNA","TOP2A","CCNB1","CDK1","CCNA2",
               "CDKN1A","GADD45A","WEE1")
test_genes <- c(damp_genes, cc_genes)

corr_table <- function(seu, virus) {
  present <- intersect(test_genes, rownames(seu))
  expr <- FetchData(seu, vars=present)
  cor_mat <- cor(expr, method="spearman")
  hmgb1_cor <- cor_mat["HMGB1", present]
  df <- data.frame(
    Virus = virus,
    Gene = names(hmgb1_cor),
    Spearman_r = round(unname(hmgb1_cor), 4)
  )
  df[df$Gene != "HMGB1", ]
}

s1_d <- corr_table(denv, "DENV")
s1_z <- corr_table(zika, "ZIKV")
s1 <- rbind(s1_d, s1_z)
write.csv(s1, file.path(outdir, "TableS1_HMGB1_correlations.csv"), row.names=FALSE)
cat("S1 exported:", nrow(s1), "rows\n")
print(s1)

# ================= S2: Cell-cycle decomposition details =================
decomp_table <- function(seu, virus_name, high_cluster) {
  meta <- seu@meta.data
  expr <- FetchData(seu, vars="HMGB1")
  df <- data.frame(
    HMGB1 = expr$HMGB1,
    viral_log = meta$viral_log,
    S.Score = meta$S.Score,
    G2M.Score = meta$G2M.Score,
    cluster = meta$seurat_clusters,
    Phase = meta$Phase
  )
  infected <- df[df$viral_log > 0.5, ]

  m1 <- lm(HMGB1 ~ viral_log, data=infected)
  m2 <- lm(HMGB1 ~ viral_log + S.Score + G2M.Score, data=infected)

  # Raw vs adjusted diff
  m_cc <- lm(HMGB1 ~ S.Score + G2M.Score, data=df)
  df$HMGB1_resid <- resid(m_cc)
  high <- df[df$cluster == high_cluster, ]
  rest <- df[df$cluster != high_cluster, ]

  data.frame(
    Virus = virus_name,
    n_infected = nrow(infected),
    n_high_viral = nrow(high),
    high_viral_cluster = high_cluster,
    R2_viral_only = round(summary(m1)$r.squared, 4),
    R2_full = round(summary(m2)$r.squared, 4),
    R2_incremental = round(summary(m2)$r.squared - summary(m1)$r.squared, 4),
    beta_viral_full = round(coef(m2)["viral_log"], 4),
    beta_S_full = round(coef(m2)["S.Score"], 4),
    beta_G2M_full = round(coef(m2)["G2M.Score"], 4),
    raw_diff = round(mean(high$HMGB1) - mean(rest$HMGB1), 4),
    adjusted_diff = round(mean(high$HMGB1_resid) - mean(rest$HMGB1_resid), 4),
    pct_remaining = round(abs(mean(high$HMGB1_resid) - mean(rest$HMGB1_resid)) /
                         abs(mean(high$HMGB1) - mean(rest$HMGB1)) * 100, 1)
  )
}

s2_d <- decomp_table(denv, "DENV", hc_d$cluster)
s2_z <- decomp_table(zika, "ZIKV", hc_z$cluster)
s2 <- rbind(s2_d, s2_z)
write.csv(s2, file.path(outdir, "TableS2_cellcycle_decomposition.csv"), row.names=FALSE)
cat("\nS2 exported:\n")
print(s2)

# ================= S3: Phase composition by cluster =================
phase_comp <- function(seu, virus) {
  meta <- seu@meta.data
  agg <- aggregate(row.names(meta) ~ seurat_clusters + Phase, data=meta, FUN=length)
  colnames(agg) <- c("Cluster", "Phase", "n_cells")
  agg$Virus <- virus
  # Calculate percentage
  totals <- aggregate(n_cells ~ Cluster, data=agg, FUN=sum)
  colnames(totals) <- c("Cluster", "total")
  agg <- merge(agg, totals, by="Cluster")
  agg$pct <- round(agg$n_cells / agg$total * 100, 1)
  agg[, c("Virus", "Cluster", "Phase", "n_cells", "total", "pct")]
}

s3_d <- phase_comp(denv, "DENV")
s3_z <- phase_comp(zika, "ZIKV")
s3 <- rbind(s3_d, s3_z)
write.csv(s3, file.path(outdir, "TableS3_phase_composition.csv"), row.names=FALSE)
cat("\nS3 phase composition exported:", nrow(s3), "rows\n")

cat("\nALL P2 SUPPLEMENTARY TABLES DONE\n")
