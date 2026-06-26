# ==============================================================================
# 02_integration.R
# ==============================================================================
# Trapezoidal integration of the ROC curve partial area.
#
# Functions in this file:
#   pA_FPR - Integrate TPR over a FPR slice (high-specificity domain).
#   pA_TPR - Integrate Specificity over a TPR slice (high-sensitivity domain).
#
# Both functions receive the output of portion_ROC_FPR / portion_ROC_TPR
# (a 2-column matrix with interpolated boundary points) and return the raw,
# unscaled partial area. The scaled indexes (MCpA, TpA, NpA, FpA) in
# 04_index_engines.R build upon these raw values.
# ==============================================================================


# ------------------------------------------------------------------------------
# pA_FPR  (Trapezoidal integration — FPR axis)
# ------------------------------------------------------------------------------
# Description:
#   Computes the raw partial area under the ROC curve over a given FPR interval
#   using the trapezoidal rule. Integrates vertically (TPR over FPR).
#   Used internally by MCpAUC, TpAUC, and pAUC(axis = "x").
#
# Parameters:
#   fpr.proc - Numeric vector of FPR values (the sub-interval slice)
#   sen.proc - Numeric vector of corresponding TPR values
#
# Returns:
#   A single numeric value: the trapezoidal partial area.
#   Returns NA_real_ for degenerate input (< 2 points or NAs present).
# ------------------------------------------------------------------------------

pA_FPR = function(fpr.proc, sen.proc) {

  # Sort points by increasing FPR to ensure correct integration direction
  ord = order(fpr.proc)
  fpr = fpr.proc[ord]
  sen = sen.proc[ord]

  # Return NA for degenerate input (fewer than 2 points or any missing values)
  if (length(fpr) < 2 || anyNA(fpr) || anyNA(sen)) return(NA_real_)

  # Trapezoidal rule: sum of trapezoid areas over each FPR step
  dx = diff(fpr)
  return(as.numeric(sum(dx * (sen[-1] + sen[-length(sen)]) / 2)))
}


# ------------------------------------------------------------------------------
# pA_TPR  (Trapezoidal integration — TPR axis)
# ------------------------------------------------------------------------------
# Description:
#   Computes the raw partial area under the ROC curve over a given TPR interval
#   using the trapezoidal rule. Integrates horizontally (Specificity over TPR).
#   Used internally by NpAUC, FpAUC, and pAUC(axis = "y").
#
# Parameters:
#   tpr.proc - Numeric vector of TPR values (the sub-interval slice)
#   fpr.proc - Numeric vector of corresponding FPR values
#
# Returns:
#   A single numeric value: the trapezoidal partial area (in specificity units).
#   Returns NA_real_ for degenerate input.
# ------------------------------------------------------------------------------

pA_TPR = function(tpr.proc, fpr.proc) {

  # Sort points by increasing TPR to ensure correct integration direction
  ord  = order(tpr.proc)
  tpr  = tpr.proc[ord]
  fpr  = fpr.proc[ord]

  # Return NA for degenerate input
  if (length(tpr) < 2 || anyNA(tpr) || anyNA(fpr)) return(NA_real_)

  # Convert FPR to specificity (1 - FPR) for horizontal integration
  dy   = diff(tpr)
  spec = 1 - fpr

  # Trapezoidal rule: sum of trapezoid areas over each TPR step
  return(as.numeric(sum(dy * (spec[-1] + spec[-length(spec)]) / 2)))
}
