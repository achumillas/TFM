# ==============================================================================
# 04_index_engines.R
# ==============================================================================
# Scalar index computation engines. Each function receives a sub-interval
# of the ROC curve (output of portion_ROC_FPR or portion_ROC_TPR) and returns
# the corresponding scaled partial-area index as a single numeric value.
#
# Functions in this file:
#   MCpA - McClish standardised pAUC (FPR domain, proper curves only).
#   TpA  - Tighter partial area index (FPR domain, any curve shape).
#   NpA  - Normalised pAUC (TPR domain, NLR-bounded curves only).
#   FpA  - Fitted partial area index (TPR domain, any curve shape).
#
# These functions are internal engines called by the exported *AUC functions
# in 08_exported_functions.R, and also by the bootstrap statistic functions
# in 05_bootstrap_stats.R.
# ==============================================================================


# ------------------------------------------------------------------------------
# MCpA  (McClish, 1989)
# ------------------------------------------------------------------------------
# Description:
#   Computes the McClish standardised partial AUC (MCpAUC) for a sub-interval
#   of the ROC curve on the FPR axis. The index rescales the raw pAUC to
#   [0.5, 1] using the area under the diagonal as the minimum reference and
#   the rectangle of height 1 as the maximum reference:
#
#       MCpAUC = (1 + (pAUC - minA) / (maxA - minA)) / 2
#
#   where  minA = (FPR_max^2 - FPR_min^2) / 2   (area under the diagonal)
#          maxA = FPR_max - FPR_min              (rectangle of height 1)
#
#   Returns NA for improper ROC curves (below the diagonal), where the index
#   is not mathematically defined.
#
# Parameters:
#   sen.proc - TPR values of the ROC sub-interval
#   fpr.proc - FPR values of the ROC sub-interval
#
# Returns:
#   Numeric scalar in [0.5, 1], or NA_real_ for improper curves or degenerate input.
#
# References:
#   McClish, D.K. (1989). Analyzing a portion of the ROC curve.
#   Medical Decision Making, 9(3), 190-195.
# ------------------------------------------------------------------------------

MCpA = function(sen.proc, fpr.proc) {

  # Return NA immediately for degenerate input
  if (length(fpr.proc) < 2 || anyNA(fpr.proc) || anyNA(sen.proc))
    return(NA_real_)

  # Compute raw pAUC and classify the curve shape
  pAUC = pA_FPR(fpr.proc, sen.proc)
  type = classification_Tp(fpr.proc, sen.proc)

  # MCpAUC is only defined for proper ROC curves (above the diagonal)
  if (!type[2]) return(NA_real_)

  fpr_min   = min(fpr.proc)
  fpr_max   = max(fpr.proc)
  fpr_range = fpr_max - fpr_min

  # Return NA for a degenerate single-point interval
  if (fpr_range == 0) return(NA_real_)

  # Upper reference: area of the rectangle of height 1 over the FPR interval
  maxA = fpr_range

  # Lower reference: area under the diagonal (integral of x from fpr_min to fpr_max)
  minA = (fpr_max^2 - fpr_min^2) / 2

  # Rescale pAUC to [0.5, 1] using the diagonal-to-rectangle range
  return((1 + (pAUC - minA) / (maxA - minA)) / 2)
}


# ------------------------------------------------------------------------------
# TpA  (Vivo et al., 2018)
# ------------------------------------------------------------------------------
# Description:
#   Computes the tighter partial area index (TpAUC) for a sub-interval on the
#   FPR axis. Unlike MCpAUC, TpAUC uses bounds derived from the actual TPR
#   values at the interval endpoints, making it valid for both proper and
#   improper curves and sensitive to crossing ROC curves with equal raw pAUC.
#
#   The lower bound (TpAUC.min) depends on the curve shape:
#     PLR-bounded → min = (1/2)(TPR1 + TPR2) * dFPR          (eq. 5, paper)
#     Proper      → min = max(TPR1*dFPR, sum(d(FPR^2))/2)     (eq. 4, paper)
#     Other       → min = TPR1 * dFPR                         (eq. 3, paper)
#
#   The upper bound is always: TpAUC.max = TPR2 * dFPR
#
#       TpAUC = (1 + (pAUC - TpAUC.min) / (TpAUC.max - TpAUC.min)) / 2
#
# Parameters:
#   fpr.proc - FPR values of the ROC sub-interval
#   sen.proc - TPR values of the ROC sub-interval
#
# Returns:
#   Numeric scalar in [0.5, 1] for any ROC curve shape.
#
# References:
#   Vivo, J.-M., Franco, M. & Vicari, D. (2018). Rethinking an ROC partial
#   area index for evaluating the classification performance at a high
#   specificity range. Advances in Data Analysis and Classification,
#   12(3), 683-704.
# ------------------------------------------------------------------------------

TpA = function(fpr.proc, sen.proc) {

  # Compute raw pAUC and classify the curve shape
  pAUC   = pA_FPR(fpr.proc, sen.proc)
  type   = classification_Tp(fpr.proc, sen.proc)

  # Total width of the FPR interval
  dFPR   = sum(diff(fpr.proc))

  # TPR at the upper and lower boundaries of the interval
  maxTPR = max(sen.proc)
  minTPR = min(sen.proc)

  # Upper bound: rectangle of width dFPR and height TPR2 (the maximum TPR)
  # Special case: if the curve is constant (TPR1 == TPR2), a warning is issued
  # and TpAUC.max is set to dFPR (consistent with the original paper supplementary code)
  if (min(fpr.proc) == max(fpr.proc)) {
    warning("Constant ROC curve over the prefixed FPR range")
    TpAUC.max = dFPR
  } else {
    TpAUC.max = dFPR * maxTPR
  }

  # Lower bound depends on the curve shape
  TpAUC.min =
    if   (minTPR == maxTPR) 0
    else if (type[1]) dFPR * mean(c(minTPR, maxTPR))                 # PLR-bounded
    else if (type[2]) max(dFPR * minTPR, sum(diff(fpr.proc^2)) / 2)  # proper
    else              dFPR * minTPR                                   # other

  # Guard against degenerate case where the curve collapses to zero
  if (maxTPR == 0) return(0)

  # Rescale pAUC to [0.5, 1] using the tighter bounds
  return((1 + (pAUC - TpAUC.min) / (TpAUC.max - TpAUC.min)) / 2)
}


# ------------------------------------------------------------------------------
# NpA  (Jiang et al., 1996)
# ------------------------------------------------------------------------------
# Description:
#   Computes the normalised partial AUC (NpAUC) for a sub-interval on the
#   TPR axis (high-sensitivity domain). The raw pAUC is divided by (1 - TPR0),
#   which is the width of the horizontal band above TPR0:
#
#       NpAUC = A_TPR0 / (1 - TPR0)
#
#   Valid only when the NLR is bounded throughout the region of interest.
#   Returns NA otherwise — use FpAUC for those cases.
#
# Parameters:
#   tpr.proc - TPR values of the ROC sub-interval (ascending, starting at TPR0)
#   fpr.proc - FPR values of the ROC sub-interval
#
# Returns:
#   Numeric scalar in [0, 1] when the NLR condition holds, NA_real_ otherwise.
#
# References:
#   Jiang, Y., Metz, C.E. & Nishikawa, R.M. (1996). A receiver operating
#   characteristic partial area index for highly sensitive diagnostic tests.
#   Radiology, 201(3), 745-750.
# ------------------------------------------------------------------------------

NpA = function(tpr.proc, fpr.proc) {

  # Compute raw horizontal-band pAUC and classify the curve shape
  pAUC = pA_TPR(tpr.proc, fpr.proc)
  type = classification_Fp(tpr.proc, fpr.proc)

  # TPR0 is the lower boundary of the sensitivity interval
  TPR0 = tpr.proc[1]

  # NpAUC is only defined when the NLR is bounded; return NA otherwise
  if (type[1]) return(pAUC / (1 - TPR0)) else return(NA_real_)
}


# ------------------------------------------------------------------------------
# FpA  (Franco & Vivo, 2021)
# ------------------------------------------------------------------------------
# Description:
#   Computes the fitted partial area index (FpAUC) for a sub-interval on the
#   TPR axis. FpAUC extends NpAUC by applying tighter NLR-based bounds derived
#   from the actual curve shape, making it valid for ANY ROC curve and always
#   returning a value in [0.5, 1].
#
#   Algorithm 1 from Franco & Vivo (2021) selects the formula based on
#   classification_Fp output:
#
#     Case 1 — Bounded   (type[1]): NLR(x) <= NLR0 for all x >= FPR0
#       lower bound = (1/2)(1-FPR0)(1-TPR0)
#       FpAUC = pAUC / [(1-FPR0)(1-TPR0)]                               (eq. 9)
#
#     Case 2 — Partially proper (type[2]): NLR(x) <= 1 but some > NLR0
#       lower bound = (1/2)(1-TPR0)^2
#       FpAUC = [pAUC + (1-TPR0)(TPR0-FPR0)] /
#               [(1+TPR0-2*FPR0)(1-TPR0)]                               (eq. 10)
#
#     Case 3 — Improper  (type[3]): some NLR(x) > 1  (hook present)
#       lower bound = 0
#       FpAUC = [pAUC + (1-TPR0)(1-FPR0)] /
#               [2(1-FPR0)(1-TPR0)]                                     (eq. 11)
#
# Parameters:
#   tpr.proc - TPR values of the ROC sub-interval (ascending, starting at TPR0)
#   fpr.proc - FPR values of the ROC sub-interval
#
# Returns:
#   Numeric scalar in [0.5, 1] for any ROC curve shape.
#
# References:
#   Franco, M. & Vivo, J.-M. (2021). Evaluating the performances of biomarkers
#   over a restricted domain of high sensitivity. Mathematics, 9(21), 2826.
# ------------------------------------------------------------------------------

# ── shapepROC ─────────────────────────────────────────────────────────────────
# Classify the behaviour of the RVN over the sensitivity sub-interval.
# The last point (typically (1,1)) is excluded to avoid calculating RVN = 0/0 = NaN.
# Returns one of: "BpNLR", "pProp", "other""
shapepROC = function(fpr.p, sen.p) {
  # We exclude the last point to avoid calculating RVN = (1-1)/(1-1) = NaN
  fpr.p = fpr.p[-length(fpr.p)]
  sen.p = sen.p[-length(sen.p)]
  
  # We calculate the RVN at each remaining point and the reference RVN at TVP0 (the first point)
  nlr.p = (1 - sen.p) / (1 - fpr.p)
  nlr0  = nlr.p[1]
  
  # We classify based on the behaviour of the RVN over the interval
  if (all(nlr.p <= nlr0) && is.finite(nlr0)) return("BpNLR")   # RVN bounded by RVN0
  else if (all(nlr.p <= 1))                  return("pProp")   # RVN ≤ 1 but exceeds RVN0
  else                                       return("other")   # RVN > 1 at some point (hook)
}

FpA = function(tpr.proc, fpr.proc) {

  # Reference points at the lower boundary of the sensitivity interval
  TPR0 = tpr.proc[1]
  FPR0 = fpr.proc[1]

  # Compute raw horizontal-band pAUC
  pAUC = pA_TPR(tpr.proc, fpr.proc)

  # Classify the curve shape based on NLR behaviour (excluding the last point)
  sproc = shapepROC(fpr.proc, tpr.proc)

  switch(sproc,

    BpNLR = {
      # Case 1 — NLR bounded: NLR(x) <= NLR0 throughout the interval
      # Lower bound = (1/2)(1-FPR0)(1-TPR0); upper bound = (1-FPR0)(1-TPR0)
      # FpAUC = pAUC / [(1-FPR0)(1-TPR0)]                          (eq. 9)
      pAUC / ((1 - FPR0) * (1 - TPR0))
    },

    pProp = {
      # Case 2 — Partially proper: NLR(x) <= 1 everywhere but exceeds NLR0 somewhere
      # Lower bound = (1/2)(1-TPR0)^2; upper bound = (1-FPR0)(1-TPR0)
      # FpAUC = [pAUC + (1-TPR0)(TPR0-FPR0)] / [(1+TPR0-2*FPR0)(1-TPR0)]  (eq. 10)
      (pAUC + (1 - TPR0) * (TPR0 - FPR0)) /
        ((1 + TPR0 - 2 * FPR0) * (1 - TPR0))
    },

    other = {
      # Case 3 — Improper: NLR(x) > 1 for some x (hook present in the curve)
      # Lower bound = 0; upper bound = (1-FPR0)(1-TPR0)
      # FpAUC = [pAUC + (1-TPR0)(1-FPR0)] / [2(1-FPR0)(1-TPR0)]   (eq. 11)
      (pAUC + (1 - TPR0) * (1 - FPR0)) /
        (2 * (1 - FPR0) * (1 - TPR0))
    }
  )
}
