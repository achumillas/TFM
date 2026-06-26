# ==============================================================================
# zzz.R
# ==============================================================================
# Package-level declarations loaded last (zzz prefix ensures load order).
#
# utils::globalVariables() suppresses R CMD check NOTEs about "no visible
# binding for global variable" that arise from ggplot2 NSE (non-standard
# evaluation), dplyr column references, and base functions used inside
# lapply/vapply closures that the static analyser cannot resolve.
# ==============================================================================

utils::globalVariables(c(
  # ggplot2 aes() column names used in plot_ROC_curves and plot_bootstrap_results
  "FPR", "TPR", "Variable",
  "Estimate", "BootMean", "lwr", "upr", "t0",
  ".data",
  # dplyr column references in allIndices
  "pAUC_FPR", "pAUC_TPR", "MCpAUC", "TpAUC", "NpAUC", "FpAUC", "pAUC",
  # base functions flagged by R CMD check inside closures
  "setNames"
))