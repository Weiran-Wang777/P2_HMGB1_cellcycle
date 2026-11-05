# ============================================================
# run_all.R — One-click reproduction of P2 (HMGB1 & cell-cycle)
# Run this from the repository root:  Rscript run_all.R
# ============================================================
set.seed(123)
cat("========================================\n")
cat("P2: HMGB1 downregulation and cell-cycle perturbation\n")
cat("========================================\n\n")

cat(">>> Step 1: Fig1-3 (DAMP gradient, violin, cell-cycle correlation) ...\n")
source("P2_HMGB1/01_qc_cluster.R")

cat("\n>>> Step 2: Cell-cycle decomposition (nested linear models) ...\n")
source("P2_HMGB1/05_cellcycle_decomp.R")

cat("\n>>> Step 3: Fig4 (boxplots + phase composition) ...\n")
source("P2_HMGB1/06_phase_composition.R")

cat("\n>>> Step 4: Fig5 (concomitant HMGB1 down vs G2M up) ...\n")
source("P2_HMGB1/08_fig5_concomitant.R")

cat("\n>>> Step 5: Pseudobulk descriptive analysis ...\n")
source("P2_HMGB1/07_pseudobulk.R")

cat("\n>>> Step 6: Recording session info (after all packages loaded) ...\n")
source("P2_HMGB1/00_sessionInfo.R")

cat("\n========================================\n")
cat("ALL P2 ANALYSES COMPLETE.\n")
cat("Figures saved to P2_HMGB1/figures/\n")
cat("========================================\n")
