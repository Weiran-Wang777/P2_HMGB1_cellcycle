# P2 pseudobulk: descriptive aggregation for HMGB1 and key genes
# Note: No biological replicates available. Results are directional only.
set.seed(123)
suppressPackageStartupMessages({library(Seurat)})

base <- "data"
if (!file.exists(file.path(base, "denv_cc.rds")))
  stop("Error: data/denv_cc.rds not found. See README for data download instructions.")
if (!file.exists(file.path(base, "zika_clustered.rds")))
  stop("Error: data/zika_clustered.rds not found. See README for data download instructions.")
denv <- readRDS(file.path(base, "denv_cc.rds"))
zika <- readRDS(file.path(base, "zika_clustered.rds"))

genes <- c("HMGB1","HMGB2","HMGN2","CXCL10","VCAM1","NFKB2","RELB","SOCS2")

pb_analysis <- function(seu, virus_name) {
  meta <- seu@meta.data
  counts <- GetAssayData(seu, layer="counts")
  meta$group <- paste0(meta$time, "_", ifelse(meta$viral_log > 0.5, "inf", "uninf"))
  results <- data.frame()
  for (g in genes) {
    if (!(g %in% rownames(counts))) next
    for (grp in unique(meta$group)) {
      cells <- rownames(meta)[meta$group == grp]
      total_counts <- sum(counts[g, cells])
      n_cells <- length(cells)
      cpc <- total_counts / n_cells
      results <- rbind(results, data.frame(
        virus = virus_name, gene = g, group = grp,
        n_cells = n_cells, total_counts = total_counts,
        counts_per_cell = round(cpc, 3)
      ))
    }
  }
  return(results)
}

pb_d <- pb_analysis(denv, "DENV")
pb_z <- pb_analysis(zika, "ZIKV")
pb <- rbind(pb_d, pb_z)

cat("\n=== PSEUDOBULK: 48h infected vs uninfected (counts/cell) ===\n")
pb_48 <- pb[grepl("48", pb$group), ]
for (g in genes) {
  if (!(g %in% pb_48$gene)) next
  d_inf <- pb_48[pb_48$virus=="DENV" & pb_48$gene==g & grepl("48_inf", pb_48$group), "counts_per_cell"]
  d_uninf <- pb_48[pb_48$virus=="DENV" & pb_48$gene==g & grepl("48_uninf", pb_48$group), "counts_per_cell"]
  z_inf <- pb_48[pb_48$virus=="ZIKV" & pb_48$gene==g & grepl("48_inf", pb_48$group), "counts_per_cell"]
  z_uninf <- pb_48[pb_48$virus=="ZIKV" & pb_48$gene==g & grepl("48_uninf", pb_48$group), "counts_per_cell"]
  cat(sprintf("  %-10s DENV: inf=%.3f uninf=%.3f  | ZIKV: inf=%.3f uninf=%.3f\n",
              g, ifelse(length(d_inf),d_inf,0), ifelse(length(d_uninf),d_uninf,0),
              ifelse(length(z_inf),z_inf,0), ifelse(length(z_uninf),z_uninf,0)))
}

write.csv(pb, "P2_HMGB1/p2_pseudobulk_results.csv", row.names=FALSE)
cat("\nSaved p2_pseudobulk_results.csv\n")
