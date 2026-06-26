# ==============================================================================
# 03_classification.R
# ==============================================================================
# ROC curve shape classification.
#
# Functions in this file:
#   classification_Tp - Classify the ROC slice in the FPR domain based on the
#                       Positive Likelihood Ratio (PLR). Used by MCpA and TpA.
#   classification_Fp - Classify the ROC slice in the TPR domain based on the
#                       Negative Likelihood Ratio (NLR). Used by NpA and FpA.
#
# Classification is necessary because the correct normalisation formula for
# each scaled index depends on the shape of the ROC curve in the region of
# interest. A misclassified curve would produce a mathematically invalid index.
# ==============================================================================


# ------------------------------------------------------------------------------
# classification_Tp  (ROC shape classifier — FPR domain)
# ------------------------------------------------------------------------------
# Description:
#   Determines the shape of the ROC curve over the selected FPR interval.
#   Used by MCpA (proper check) and TpA (PLR-bounded vs proper vs other).
#
#   Returns a length-2 logical vector:
#     [1] PLR.bounded : partial PLR attains its minimum at FPR2 (the upper extreme)
#     [2] proper      : the curve lies entirely above the diagonal (TPR >= FPR)
#
# Parameters:
#   fpr.proc - FPR values of the ROC sub-interval
#   sen.proc - TPR values of the ROC sub-interval
#
# Returns:
#   logical[2]: c(PLR.bounded, proper)
# ------------------------------------------------------------------------------

classification_Tp = function(fpr.proc, sen.proc) {

  # A curve is proper if it lies entirely above the diagonal (TPR >= FPR everywhere)
  proper = all(sen.proc >= fpr.proc)

  # Compute increments relative to the first point to build the partial PLR
  d_fpr = fpr.proc - fpr.proc[1]
  d_sen = sen.proc - sen.proc[1]

  # Exclude points where the FPR denominator is zero (vertical segment at start)
  # to avoid Inf / NaN in the PLR ratio
  finite_idx = which(d_fpr > 0)

  if (length(finite_idx) < 1L) {
    # All FPR values are identical — PLR is undefined; treat as not PLR-bounded
    plr.bd = FALSE
  } else {
    plr = d_sen[finite_idx] / d_fpr[finite_idx]
    # PLR-bounded: the partial PLR must be non-decreasing toward the upper boundary FPR2
    plr.bd = all(plr >= plr[length(plr)])
  }

  # Return the two classification flags as a logical vector
  return(c(plr.bd, proper))
}


# ------------------------------------------------------------------------------
# classification_Fp  (ROC shape classifier — TPR domain)
# ------------------------------------------------------------------------------
# Description:
#   Determines the shape of the ROC curve over the selected TPR interval based
#   on the Negative Likelihood Ratio (NLR = (1 - TPR) / (1 - FPR)).
#   Used by NpA (bounded check) and FpA (case 1 / 2 / 3 selection).
#
#   Returns a length-3 logical vector (exactly one element is TRUE):
#     [1] bounded  : NLR(x) <= NLR0 for all x >= FPR0
#     [2] partial  : some NLR(x) > NLR0, but NLR(x) <= 1 everywhere
#     [3] improper : some NLR(x) > 1  (hook present)
#
# Parameters:
#   tpr.proc - TPR values of the ROC sub-interval (ascending, starting at TPR0)
#   fpr.proc - FPR values of the ROC sub-interval
#
# Returns:
#   logical[3]: c(bounded, partial, improper)
# ------------------------------------------------------------------------------

classification_Fp = function(tpr.proc, fpr.proc) {

  # Reference point: NLR at the lower boundary of the sensitivity interval (TPR0)
  TPR0 = tpr.proc[1]
  FPR0 = fpr.proc[1]
  NLR0 = (1 - TPR0) / (1 - FPR0)

  # Compute NLR at every point in the slice and remove non-finite values
  NLR_v = (1 - tpr.proc) / (1 - fpr.proc)
  NLR_v = NLR_v[is.finite(NLR_v)]

  # Classify the curve shape based on NLR behaviour
  bounded  = all(NLR_v <= NLR0)      # NLR stays below the reference throughout
  improper = any(NLR_v > 1)           # NLR exceeds 1 somewhere (hook present)
  partial  = !bounded && !improper    # intermediate case: partially proper

  # Return the three mutually exclusive classification flags
  return(c(bounded, partial, improper))
}
