# P2: Cell-cycle covariate analysis for HMGB1 downregulation
# Goal: quantify how much of HMGB1 downregulation is explained by cell cycle vs viral load
set.seed(123)
suppressPackageStartupMessages({library(Seurat)})

base <- "data"
if (!file.exists(file.path(base, "denv_cc.rds")))
  stop("Error: data/denv_cc.rds not found. See README for data download instructions.")
if (!file.exists(file.path(base, "zika_clustered.rds")))
  stop("Error: data/zika_clustered.rds not found. See README for data download instructions.")
denv <- readRDS(file.path(base, "denv_cc.rds"))
zika <- readRDS(file.path(base, "zika_clustered.rds"))

# ---- Auto-detect high-viral cluster ----
find_high_viral <- function(seu) {
  viral_by_cluster <- aggregate(seu$viral_log, by = list(cluster = seu$seurat_clusters),
                                FUN = median, na.rm = TRUE)
  hc <- viral_by_cluster$cluster[which.max(viral_by_cluster$x)]
  message("High-viral cluster: ", hc, " (median viral_log = ", round(max(viral_by_cluster$x), 2), ")")
  return(hc)
}

analyze_virus <- function(seu, virus_name, high_cluster) {
  cat("\n=====", virus_name, "=====\n")
  meta <- seu@meta.data
  expr <- FetchData(seu, vars = "HMGB1")
  df <- data.frame(
    HMGB1 = expr$HMGB1,
    viral_log = meta$viral_log,
    S.Score = meta$S.Score,
    G2M.Score = meta$G2M.Score,
    cluster = meta$seurat_clusters
  )
  # Infected cells (viral_log > 0.5, consistent with P2 manuscript)
  infected <- df[df$viral_log > 0.5, ]
  cat("Infected cells (viral_log > 0.5):", nrow(infected), "\n")

  # Model 1: HMGB1 ~ viral_load only
  m1 <- lm(HMGB1 ~ viral_log, data = infected)
  r2_viral <- summary(m1)$r.squared
  cat("Model 1 (viral_log only): R2 =", round(r2_viral, 4), "\n")

  # Model 2: HMGB1 ~ viral_load + cell cycle
  m2 <- lm(HMGB1 ~ viral_log + S.Score + G2M.Score, data = infected)
  r2_full <- summary(m2)$r.squared
  cat("Model 2 (viral_log + S.Score + G2M.Score): R2 =", round(r2_full, 4), "\n")

  r2_incremental <- r2_full - r2_viral
  cat("Incremental R2 from cell cycle:", round(r2_incremental, 4),
      "(", round(r2_incremental/r2_full*100, 1), "% of full model)\n")

  # Residual approach: regress HMGB1 on cell cycle, then compare
  m_cc <- lm(HMGB1 ~ S.Score + G2M.Score, data = df)
  df$HMGB1_resid <- resid(m_cc)
  high <- df[df$cluster == high_cluster, ]
  rest <- df[df$cluster != high_cluster, ]

  cat("\nHigh-viral cluster:", high_cluster, "(", nrow(high), "cells)\n")
  raw_diff <- mean(high$HMGB1) - mean(rest$HMGB1)
  adj_diff <- mean(high$HMGB1_resid) - mean(rest$HMGB1_resid)
  cat("  Raw HMGB1 diff:", round(raw_diff, 4), "\n")
  cat("  Cell-cycle-adjusted diff:", round(adj_diff, 4), "\n")
  if (abs(raw_diff) > 0.001) {
    pct_remaining <- abs(adj_diff) / abs(raw_diff) * 100
    cat("  Adjusted retains", round(pct_remaining, 1), "% of raw effect\n")
  }
  return(list(r2_viral = r2_viral, r2_full = r2_full,
              r2_incremental = r2_incremental,
              raw_diff = raw_diff, adj_diff = adj_diff))
}

denv_res <- analyze_virus(denv, "DENV", find_high_viral(denv))
zika_res <- analyze_virus(zika, "ZIKV", find_high_viral(zika))

# Save results
sink("P2_HMGB1/p2_cellcycle_results.txt")
cat("P2 Cell-cycle covariate analysis\n")
cat("================================\n\n")
cat("DENV: viral_log R2 =", round(denv_res$r2_viral, 4),
    " full R2 =", round(denv_res$r2_full, 4),
    " cell-cycle incremental R2 =", round(denv_res$r2_incremental, 4), "\n")
cat("ZIKV: viral_log R2 =", round(zika_res$r2_viral, 4),
    " full R2 =", round(zika_res$r2_full, 4),
    " cell-cycle incremental R2 =", round(zika_res$r2_incremental, 4), "\n\n")
cat("DENV high-viral vs rest: raw diff =", round(denv_res$raw_diff, 4),
    " cell-cycle-adjusted diff =", round(denv_res$adj_diff, 4), "\n")
cat("ZIKV high-viral vs rest: raw diff =", round(zika_res$raw_diff, 4),
    " cell-cycle-adjusted diff =", round(zika_res$adj_diff, 4), "\n")
sink()
cat("\nResults saved to p2_cellcycle_results.txt\n")
