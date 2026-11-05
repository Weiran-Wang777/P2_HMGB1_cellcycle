# 00_sessionInfo.R — Run this to record package versions for reproducibility
# Run at the end of the analysis and save output to session_info.txt

set.seed(123)
sink("P2_HMGB1/session_info.txt")
cat("Session info for P2 HMGB1 reanalysis\n")
cat("======================================\n\n")
print(sessionInfo())
cat("\n\nKey package versions:\n")
for (pkg in c("Seurat", "ggplot2", "patchwork")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %s: %s\n", pkg, as.character(packageVersion(pkg))))
  }
}
sink()
cat("Saved session_info.txt\n")
