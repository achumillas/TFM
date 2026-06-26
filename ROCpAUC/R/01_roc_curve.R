# ==============================================================================
# 01_roc_curve.R
# ==============================================================================
# Dataset preparation, empirical ROC curve construction, and sub-interval
# extraction. These are the entry-point functions called by every exported
# index function before any computation takes place.
#
# Functions in this file:
#   prepare_ROC_dataset - Validate and standardise input (data.frame or SE).
#   pointsCurve         - Compute (FPR, TPR) pairs from outcome + predictor.
#   portion_ROC_FPR     - Extract the ROC slice in a given FPR interval.
#   portion_ROC_TPR     - Extract the ROC slice in a given TPR interval.
# ==============================================================================


# ------------------------------------------------------------------------------
# prepare_ROC_dataset
# ------------------------------------------------------------------------------
# Description:
#   Validates and standardises the input dataset for all ROC functions.
#   Accepts either a plain data.frame or a Bioconductor SummarizedExperiment.
#
#   For SummarizedExperiment inputs the function:
#     - Extracts the binary outcome vector from colData[[variable]].
#     - Subsets the assay to the features named in 'selection'.
#     - Transposes (samples as rows) and applies mean-centring and
#       unit-variance scaling so predictors are on a comparable scale.
#   For data.frame inputs no transformation is applied; the data are used as-is.
#
#   The returned data.frame always has:
#     - Column 1  : binary outcome (0/1 or 1/2 coding accepted downstream)
#     - Columns 2+: predictor (biomarker) values
#
# Parameters:
#   dataset   - Input data.frame OR SummarizedExperiment object
#   selection - Character vector of feature/gene names  (SE input only)
#   variable  - Name of the outcome column in colData   (SE input only)
#
# Returns:
#   A data.frame ready for ROC analysis with at least 2 columns.
# ------------------------------------------------------------------------------

prepare_ROC_dataset = function(dataset, selection = NULL, variable = NULL) {

  if (methods::is(dataset, "SummarizedExperiment")) {

    # Both arguments are mandatory when using SE input
    stopifnot(is.character(selection), is.character(variable))

    # Extract the binary outcome from the sample-level metadata
    outcome = dataset@colData@listData[[variable]]

    # Subset assay rows to the requested features, transpose, and scale
    mat = scale(
      t(as.matrix(SummarizedExperiment::assay(dataset)[selection, ])),
      center = TRUE, scale = TRUE
    )
    colnames(mat) = selection

    # Bind outcome and feature matrix into a single data.frame
    dataset = as.data.frame(cbind(outcome, mat))
    colnames(dataset)[1] = variable

  } else {
    # Plain data.frame: coerce and pass through without modification
    dataset = as.data.frame(dataset)
  }

  # Sanity check: at least one predictor column must be present
  if (ncol(dataset) < 2)
    stop("Dataset must have at least 2 columns (outcome + one predictor).")

  return(dataset)
}


# ------------------------------------------------------------------------------
# pointsCurve
# ------------------------------------------------------------------------------
# Description:
#   Computes the empirical ROC curve as a sequence of (FPR, TPR) points for a
#   continuous predictor 'y' with respect to a binary class label 'x'.
#
#   The computation follows the same logic as the ROCR package: thresholds are
#   all unique score values ordered descending, flanked by +Inf and -Inf
#   sentinels to guarantee the points (0, 0) and (1, 1) on the curve. TP and
#   FP counts are built by cumulative summation as the threshold is lowered,
#   which handles ties correctly. FPR and TPR are then derived directly from
#   those cumulative counts. Duplicate (FPR, TPR) pairs arising from tied
#   scores are removed, and points are sorted by ascending FPR so that
#   portion_ROC_FPR and portion_ROC_TPR can interpolate precise boundary
#   points for any requested sub-interval, enabling mathematically exact
#   computation of the scaled partial-area indexes.
#
#   Class labels may be coded 0/1 or 1/2.
#
# Parameters:
#   x     - Integer/numeric vector of binary class labels (0/1 or 1/2)
#   y     - Numeric predictor (biomarker) vector, same length as x
#   plot  - Logical; if TRUE prints a ggplot2 ROC curve
#   label - Character string used as the plot title (default: variable name)
#
# Returns:
#   A data.frame with columns FPR and TPR (both in [0, 1]),
#   ordered by ascending FPR, starting at (0, 0) and ending at (1, 1).
# ------------------------------------------------------------------------------

#' @title Compute empirical ROC curve points and optionally plot
#' @description
#'   Derives the (FPR, TPR) pairs for a continuous predictor \code{y} with
#'   respect to a binary outcome \code{x}.
#' @details
#'   The computation follows the same logic as the ROCR package: thresholds are
#'   all unique score values ordered descending, flanked by \code{+Inf} and
#'   \code{-Inf} sentinels to guarantee the points (0, 0) and (1, 1) on the
#'   curve. TP and FP counts are built by cumulative summation as the threshold
#'   is lowered, which handles ties correctly. FPR and TPR are then derived
#'   directly from those cumulative counts. Duplicate (FPR, TPR) pairs arising
#'   from tied scores are removed, and points are sorted by ascending FPR so
#'   that \code{portion_ROC_FPR} and \code{portion_ROC_TPR} can interpolate
#'   precise boundary points for any requested sub-interval, enabling
#'   mathematically exact computation of the scaled partial-area indexes.
#'
#'   Class labels may be coded 0/1 or 1/2.
#' @param x     Integer/numeric vector of binary class labels (0/1 or 1/2).
#' @param y     Numeric predictor vector.
#' @param plot  Logical; if \code{TRUE} prints a ggplot2 ROC plot.
#' @param label Character label for the plot title.
#' @return A \code{data.frame} with columns \code{FPR} and \code{TPR},
#'   ordered by ascending FPR, starting at (0, 0) and ending at (1, 1).
#' @examples
#' library(fission)
#' data(fission)
#' genes  = c("SPNCRNA.1080", "SPAC186.08c")
#' strain = SummarizedExperiment::colData(fission)$strain
#' expr   = t(SummarizedExperiment::assay(fission)[genes, ])
#' Sp_ex  = as.data.frame(cbind(strain = strain, expr))
#'
#' # Input type 1: plain data.frame
#' pointsCurve(Sp_ex[, 1], as.numeric(Sp_ex[, 2]),
#'             plot = TRUE, label = colnames(Sp_ex)[2])
#'
#' # Input type 2: vectors extracted from a SummarizedExperiment
#' outcome = as.numeric(SummarizedExperiment::colData(fission)$strain)
#' pred    = as.numeric(SummarizedExperiment::assay(fission)["SPNCRNA.1080", ])
#' pointsCurve(outcome, pred, plot = TRUE, label = "SPNCRNA.1080 (from SE)")
#' @export
pointsCurve = function(x, y, plot = FALSE, label = deparse(substitute(y))) {

  stopifnot(is.numeric(y) || is.integer(y))

  # Remove observations where x or y is missing
  valid   = which(!is.na(x) & !is.na(y))
  xsample = as.numeric(trimws(as.character(x[valid])))
  ysample = y[valid]

  if (anyNA(xsample))
    stop("Class variable contains non-numeric or missing values.")

  # Harmonise coding: accept 0/1 or 1/2; convert to 0/1
  uvals = sort(unique(xsample))
  if (identical(uvals, c(1, 2))) {
    xsample = ifelse(xsample == 1, 0L, 1L)
  } else if (!identical(uvals, c(0, 1))) {
    stop("Class labels must be coded as 0/1 or 1/2.")
  }
  xsample = as.integer(xsample)

  n_neg = sum(xsample == 0L)
  n_pos = sum(xsample == 1L)
  if (n_neg == 0L) stop("No negative-class observations found.")
  if (n_pos == 0L) stop("No positive-class observations found.")

  # Step 1: unique thresholds (descending) + sentinels
  uniq_desc = sort(unique(ysample), decreasing = TRUE)
  cutoffs   = c(Inf, uniq_desc, -Inf)

  # Steps 2-3: cumulative TP / FP counts → FPR, TPR
  n_pos_at = vapply(uniq_desc,
                    function(th) sum(ysample == th & xsample == 1L),
                    integer(1))
  n_neg_at = vapply(uniq_desc,
                    function(th) sum(ysample == th & xsample == 0L),
                    integer(1))

  tp    = c(0L, cumsum(n_pos_at))   # length == length(cutoffs)
  fp    = c(0L, cumsum(n_neg_at))
  tpr_v = tp / n_pos   # Sensitivity (True Positive Rate)
  fpr_v = fp / n_neg   # 1 - Specificity (False Positive Rate)

  # Step 4: remove duplicate (FPR, TPR) pairs
  keep  = !duplicated(cbind(fpr_v, tpr_v))
  fpr_v = fpr_v[keep]
  tpr_v = tpr_v[keep]

  # Step 5: sort by ascending FPR
  ord   = order(fpr_v)
  fpr_v = fpr_v[ord]
  tpr_v = tpr_v[ord]

  xy = data.frame(FPR = fpr_v, TPR = tpr_v)

  # Optional ggplot2 ROC plot
  if (isTRUE(plot)) {
    p = ggplot2::ggplot(xy, ggplot2::aes(x = FPR, y = TPR)) +
      ggplot2::geom_line(colour = "steelblue", linewidth = 1.1) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", colour = "grey50") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = paste("ROC Curve \u2014", label),
        x     = "False Positive Rate (1 \u2013 Specificity)",
        y     = "True Positive Rate (Sensitivity)"
      ) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
    print(p)
  }

  return(xy)
}


# ------------------------------------------------------------------------------
# portion_ROC_FPR
# ------------------------------------------------------------------------------
# Description:
#   Extracts the segment of the ROC curve within a specified FPR sub-interval
#   [low.limit, up.limit]. If the exact boundary values are not empirical points
#   of the curve, linear interpolation is used to insert precise boundary
#   points, guaranteeing that the partial area is computed over exactly the
#   requested specificity range with no numerical gaps at the endpoints.
#
# Parameters:
#   up.limit  - Upper bound of the FPR interval (e.g., 0.25 for high specificity)
#   low.limit - Lower bound of the FPR interval (typically 0)
#   fpr.roc   - Numeric vector of FPR values of the full ROC curve
#   sen.roc   - Numeric vector of TPR values of the full ROC curve
#
# Returns:
#   A 2-column matrix [FPR | TPR] restricted to [low.limit, up.limit],
#   with interpolated boundary points included.
# ------------------------------------------------------------------------------

portion_ROC_FPR = function(up.limit, low.limit, fpr.roc, sen.roc) {

  # Locate the first empirical point at or above low.limit (lower boundary)
  # and the last empirical point at or below up.limit (upper boundary)
  i.l = min(which(fpr.roc >= low.limit))
  j.l = max(i.l - 1, 1)
  i.u = max(which(fpr.roc <= up.limit))
  j.u = min(1 + i.u, length(fpr.roc))

  # Slice the curve between the two boundary indices
  fpr.p = fpr.roc[i.l:i.u]
  sen.p = sen.roc[i.l:i.u]

  # Interpolate the lower boundary if the first empirical point overshoots low.limit
  # lscale expresses how far low.limit sits between fpr.roc[j.l] and fpr.roc[i.l]
  if (fpr.roc[i.l] > low.limit) {
    fpr.p  = append(fpr.p, low.limit, 0)
    lscale = (fpr.p[1] - fpr.roc[j.l]) / (fpr.roc[i.l] - fpr.roc[j.l])
    sen.p  = append(sen.p,
                    sen.roc[j.l] + (sen.roc[i.l] - sen.roc[j.l]) * lscale,
                    0)
  }

  # Interpolate the upper boundary if the last empirical point falls short of up.limit
  # uscale expresses how far up.limit sits between fpr.roc[i.u] and fpr.roc[j.u]
  if (fpr.roc[i.u] < up.limit) {
    fpr.p  = append(fpr.p, up.limit, length(fpr.p))
    uscale = (fpr.roc[j.u] - fpr.p[length(fpr.p)]) /
               (fpr.roc[j.u] - fpr.roc[i.u])
    sen.p  = append(sen.p,
                    sen.roc[j.u] - (sen.roc[j.u] - sen.roc[i.u]) * uscale,
                    length(sen.p))
  }

  return(cbind(fpr.p, sen.p))
}


# ------------------------------------------------------------------------------
# portion_ROC_TPR
# ------------------------------------------------------------------------------
# Description:
#   Analogous to portion_ROC_FPR but operates on the TPR (sensitivity) axis.
#   Extracts the segment of the ROC curve within [low.limit, up.limit] on the
#   TPR axis, inserting interpolated boundary points where needed.
#   Used by NpAUC and FpAUC (high-sensitivity domain).
#
# Parameters:
#   up.limit  - Upper bound of the TPR interval (typically 1)
#   low.limit - Lower bound of the TPR interval (e.g., 0.9 for high sensitivity)
#   fpr.roc   - Numeric vector of FPR values of the full ROC curve
#   sen.roc   - Numeric vector of TPR values of the full ROC curve
#
# Returns:
#   A 2-column matrix [FPR | TPR] restricted to the requested TPR band,
#   with interpolated boundary points included.
# ------------------------------------------------------------------------------

portion_ROC_TPR = function(up.limit, low.limit, fpr.roc, sen.roc) {

  # Locate the first empirical point at or above low.limit on the TPR axis
  i.low = min(which(sen.roc >= low.limit))
  j.low = max(i.low - 1, 1)

  # Slice from the lower boundary to the end of the curve
  fpr.p = fpr.roc[i.low:length(fpr.roc)]
  sen.p = sen.roc[i.low:length(sen.roc)]

  # Interpolate the lower boundary if the first empirical point overshoots low.limit
  # lscale expresses how far low.limit sits between sen.roc[j.low] and sen.roc[i.low]
  if (sen.roc[i.low] > low.limit && i.low > 1) {
    sen.p  = append(sen.p, low.limit, 0)
    lscale = (sen.p[1] - sen.roc[j.low]) / (sen.roc[i.low] - sen.roc[j.low])
    fpr.p  = append(fpr.p,
                    fpr.roc[j.low] + (fpr.roc[i.low] - fpr.roc[j.low]) * lscale,
                    0)
  }

  return(cbind(fpr.p, sen.p))
}
