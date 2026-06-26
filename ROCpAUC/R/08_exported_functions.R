# ==============================================================================
# 08_exported_functions.R
# ==============================================================================
# All user-facing exported functions. These are the only functions that appear
# in NAMESPACE and are callable by the end user.
#
# Functions in this file:
#   MCpAUC      - McClish standardised pAUC (FPR domain).
#   TpAUC       - Tighter partial area index (FPR domain).
#   NpAUC       - Normalised pAUC (TPR domain).
#   FpAUC       - Fitted partial area index (TPR domain).
#   pAUC        - Raw (unscaled) partial area under the ROC curve.
#   allIndices  - All four indexes in a single wide tibble.
#   MCpAUCboot  - Bootstrap CIs for MCpAUC.
#   TpAUCboot   - Bootstrap CIs for TpAUC.
#   NpAUCboot   - Bootstrap CIs for NpAUC.
#   FpAUCboot   - Bootstrap CIs for FpAUC.
# ==============================================================================


# ==============================================================================
# Point-estimate functions
# ==============================================================================


# ------------------------------------------------------------------------------
# Calculate the partial area under the ROC curve (unscaled) (pAUC)
# ------------------------------------------------------------------------------
# Description:
#   Returns the trapezoidal pAUC for a given TFP or TVP interval without any
#   scaling. Useful as a diagnostic building block or for direct comparison
#   with software that reports raw pAUC values. The four scaled indices
#   (MCpAUC, TpAUC, NpAUC, FpAUC) are derived from this quantity.
#
#   Accepts three types of input:
#     · data.frame           : first column = binary outcome, rest = predictors
#     · SummarizedExperiment : use ‘selection’ and ‘variable’ to specify
#                              features and outcome column in colData
#
#   Two integration domains are supported via the ‘axis’ parameter:
#     · “x” : TFP domain (high specificity) — vertical integration
#     · “y” : TVP domain (high sensitivity) — horizontal integration
#
# Parameters:
#   dataset      - data.frame or SummarizedExperiment
#   low.value    - Lower bound of the interval (TFP or TVP depending on the axis)
#   up.value     - Upper bound of the interval
#   axis         - “x” for TFP domain; “y” for TVP domain
#   selection    - Character vector with feature names (SE input only)
#   variable     - Name of the binary result column in colData (SE input only)
#   output_as_SE - Boolean; if TRUE, returns a SummarizedExperiment
#
# Returns:
#   a tibble with Variable and pAUC columns; or a SummarizedExperiment with
#   pAUC as the only row in the assay.
# ------------------------------------------------------------------------------
#' @title Compute the raw (unscaled) Partial Area Under the ROC Curve (pAUC)
#' @description
#'   Returns the trapezoidal pAUC for a given FPR or TPR interval without any
#'   scaling. Useful as a diagnostic building block or for direct comparison
#'   with software that reports raw pAUC values. All four scaled indexes
#'   (MCpAUC, TpAUC, NpAUC, FpAUC) are derived from this quantity.
#' @param dataset      A \code{data.frame} or \code{SummarizedExperiment}.
#' @param low.value    Lower bound of the interval (FPR or TPR depending on
#'   \code{axis}).
#' @param up.value     Upper bound of the interval.
#' @param axis         \code{"x"} for the FPR domain (high specificity, vertical
#'   integration); \code{"y"} for the TPR domain (high sensitivity, horizontal
#'   integration).
#' @param selection    Character vector of feature names (SE input only).
#' @param variable     Character; outcome column in \code{colData} (SE input only).
#' @param output_as_SE Logical; if \code{TRUE} returns a \code{SummarizedExperiment}.
#' @return \code{tibble} with columns \code{Variable} and \code{pAUC}, or a
#'   \code{SummarizedExperiment} with \code{pAUC} as the single assay row.
#' @examples
#' library(fission)
#' data(fission)
#' genes_of_interest = c("SPNCRNA.1080", "SPAC186.08c",
#'                       "SPNCRNA.1420",  "SPCC70.08c", "SPAC212.04c")
#' strain = SummarizedExperiment::colData(fission)$strain
#' expr   = t(SummarizedExperiment::assay(fission)[genes_of_interest, ])
#' Sp     = as.data.frame(cbind(strain = strain, expr))
#'
#' # Input type 1: SummarizedExperiment, FPR domain [0, 0.25]
#' pAUC(fission, low.value = 0, up.value = 0.25,
#'      axis = "x", selection = genes_of_interest, variable = "strain")
#'
#' # Input type 2: SummarizedExperiment, TPR domain [0.9, 1]
#' pAUC(fission, low.value = 0.9, up.value = 1,
#'      axis = "y", selection = genes_of_interest, variable = "strain")
#'
#' # Input type 3: plain data.frame, FPR domain
#' pAUC(Sp, low.value = 0, up.value = 0.25, axis = "x")
#' @export
pAUC = function(dataset, low.value = 0, up.value = 1,
                axis = c("x", "y"),
                selection = NULL, variable = NULL,
                output_as_SE = FALSE) {
  
  axis    = match.arg(axis)
  dataset = prepare_ROC_dataset(dataset, selection, variable)
  name.var = colnames(dataset)
  
  n_pred    = ncol(dataset) - 1L
  pAUC_list = vector("list", n_pred)
  
  for (i in seq(2L, ncol(dataset))) {
    tmp = cbind(dataset[, 1], dataset[, i])
    roc = pointsCurve(tmp[, 1], tmp[, 2])
    
    val = tryCatch({
      if (axis == "x") {
        portion = portion_ROC_FPR(up.value, low.value, roc[, 1], roc[, 2])
        pA_FPR(portion[, 1], portion[, 2])
      } else {
        portion = portion_ROC_TPR(up.value, low.value, roc[, 1], roc[, 2])
        pA_TPR(portion[, 2], portion[, 1])
      }
    }, error = function(e) NA_real_)
    
    pAUC_list[[i - 1L]] = val
  }
  
  var_names  = name.var[-1]
  clean_val  = function(x) if (length(x) == 1L && is.finite(x)) x else NA_real_
  pAUC_clean = vapply(pAUC_list, clean_val, numeric(1))
  
  if (isTRUE(output_as_SE)) {
    obj = list(pAUC = pAUC_clean)
    return(createSE(obj, var_names))
  }
  
  return(tibble::tibble(Variable = var_names, pAUC = pAUC_clean))
}

# ------------------------------------------------------------------------------
# MCpAUC
# ------------------------------------------------------------------------------
# Description:
#   Computes the standardised partial area under the ROC curve (MCpAUC) in a
#   restricted FPR interval [low.value, up.value]. The raw pAUC is rescaled
#   to [0.5, 1] using the area under the diagonal as the minimum and the
#   rectangle of height 1 as the maximum. Undefined (NA) for improper ROC
#   curves (those that cross the diagonal within the interval of interest).
#
#   Accepts three input types:
#     · data.frame           : first column = binary outcome, remaining = predictors
#     · SummarizedExperiment : use 'selection' and 'variable' to specify
#                              features and outcome column in colData
#
#   Output can be returned as a tidy tibble (default) or as a
#   SummarizedExperiment with metrics as assay rows.
#
# Parameters:
#   dataset        - data.frame or SummarizedExperiment
#   low.value      - Numeric in [0, 1]; lower FPR bound (default 0)
#   up.value       - Numeric in [0, 1]; upper FPR bound (default 1)
#   plot           - Logical; if TRUE generates ROC curve plots with shaded
#                    pAUC region
#   plot_type      - One of "combined", "individual", or "both"
#   plots_per_page - Integer; number of panels per page in individual layout
#   selection      - Character vector of feature names  (SE input only)
#   variable       - Name of the binary outcome column in colData (SE input only)
#   output_as_SE   - Logical; if TRUE returns a SummarizedExperiment
#
# Returns:
#   tibble with columns Variable, MCpAUC, pAUC; or a SummarizedExperiment
#   with those metrics as assay rows. MCpAUC = NA for improper ROC curves.
#
# References:
#   McClish, D.K. (1989). Analyzing a portion of the ROC curve.
#   Medical Decision Making, 9(3), 190-195.
# ------------------------------------------------------------------------------
#' @title Compute the McClish Standardised Partial AUC (MCpAUC)
#' @description
#'   Computes the standardised partial area under the ROC curve (MCpAUC) in a
#'   restricted FPR interval [low.value, up.value]. The raw pAUC is rescaled to
#'   [0.5, 1] using the area under the diagonal as the minimum and the rectangle
#'   of height 1 as the maximum. Undefined (NA) for improper ROC curves.
#' @param dataset      A \code{data.frame} (first column = binary outcome) or
#'   \code{SummarizedExperiment}. For SE inputs use \code{selection} and
#'   \code{variable} to specify predictors and outcome.
#' @param low.value    Numeric in [0, 1]; lower FPR bound (default 0).
#' @param up.value     Numeric in [0, 1]; upper FPR bound (default 1).
#' @param plot         Logical; if \code{TRUE} generates ROC curve plots with
#'   shaded pAUC region.
#' @param plot_type    One of \code{"combined"}, \code{"individual"},
#'   \code{"both"}.
#' @param plots_per_page Integer; panels per page in individual layout.
#' @param selection    Character vector of feature names (SE input only).
#' @param variable     Character; binary outcome column in \code{colData}
#'   (SE input only).
#' @param output_as_SE Logical; if \code{TRUE} returns a
#'   \code{SummarizedExperiment} instead of a \code{tibble}.
#' @return \code{tibble} with columns \code{Variable}, \code{MCpAUC},
#'   \code{pAUC}, or a \code{SummarizedExperiment} with those metrics as
#'   assay rows. \code{MCpAUC = NA} for improper ROC curves.
#' @references
#'   McClish, D.K. (1989). Analyzing a portion of the ROC curve.
#'   \emph{Medical Decision Making}, 9(3), 190--195.
#'
#'   Thompson, M.L. & Zucchini, W. (1989). On the statistical analysis of
#'   ROC curves. \emph{Statistics in Medicine}, 8(10), 1277--1294.
#' @examples
#' library(fission)
#' data(fission)
#' genes_of_interest = c("SPNCRNA.1080", "SPAC186.08c",
#'                       "SPNCRNA.1420",  "SPCC70.08c", "SPAC212.04c")
#' strain = SummarizedExperiment::colData(fission)$strain
#' expr   = t(SummarizedExperiment::assay(fission)[genes_of_interest, ])
#' Sp     = as.data.frame(cbind(strain = strain, expr))
#'
#' # Input type 1: SummarizedExperiment -> tibble (with plots)
#' res_MC = MCpAUC(fission, low.value = 0, up.value = 0.25,
#'                 selection = genes_of_interest, variable = "strain",
#'                 plot = TRUE, plot_type = "both", plots_per_page = 5)
#' res_MC
#'
#' # Input type 2: SummarizedExperiment -> SummarizedExperiment
#' se_MC = MCpAUC(fission, low.value = 0, up.value = 0.25,
#'                selection = genes_of_interest, variable = "strain",
#'                output_as_SE = TRUE)
#' assay(se_MC)
#'
#' # Input type 3: plain data.frame -> tibble
#' MCpAUC(Sp, low.value = 0, up.value = 0.25)
#' @export
MCpAUC = function(dataset, low.value = 0, up.value = 1, plot = FALSE,
                  plot_type = c("combined", "individual", "both"),
                  selection = NULL, variable = NULL,
                  plots_per_page = 2, output_as_SE = FALSE) {

  plot_type = match.arg(plot_type)
  dataset   = prepare_ROC_dataset(dataset, selection, variable)
  name.var  = colnames(dataset)
  n_pred    = ncol(dataset) - 1L

  idx_list  = vector("list", n_pred)
  pAUC_list = vector("list", n_pred)
  roc_list  = vector("list", n_pred)

  for (i in seq(2L, ncol(dataset))) {
    tmp              = cbind(dataset[, 1], dataset[, i])
    roc              = pointsCurve(tmp[, 1], tmp[, 2])
    portion          = portion_ROC_FPR(up.value, low.value, roc[, 1], roc[, 2])
    idx_list[[i-1]]  = MCpA(portion[, 2], portion[, 1])
    pAUC_list[[i-1]] = pA_FPR(portion[, 1], portion[, 2])
    roc_list[[i-1]]  = list(ROC_points = as.matrix(roc))
  }

  if (isTRUE(plot))
    plot_ROC_curves(
      list(name.variable = name.var[-1], results = roc_list),
      low.value = low.value, up.value = up.value,
      plot_type = plot_type, axis = "x", plots_per_page = plots_per_page
    )

  return(format_results(idx_list, pAUC_list, name.var, output_as_SE, "MCpAUC"))
}


# ------------------------------------------------------------------------------
# TpAUC
# ------------------------------------------------------------------------------
# Description:
#   Computes the tighter partial area index (TpAUC) over a restricted FPR
#   interval [low.value, up.value]. Unlike MCpAUC, TpAUC uses bounds derived
#   from the actual TPR values at the interval endpoints rather than the
#   diagonal, making it valid for both proper and improper ROC curves. The
#   index is always in [0.5, 1] and can distinguish crossing curves with
#   equal raw pAUC.
#
#   Accepts three input types:
#     · data.frame           : first column = binary outcome, remaining = predictors
#     · SummarizedExperiment : use 'selection' and 'variable' to specify
#                              features and outcome column in colData
#
# Parameters:
#   dataset        - data.frame or SummarizedExperiment
#   low.value      - Numeric in [0, 1]; lower FPR bound (default 0)
#   up.value       - Numeric in [0, 1]; upper FPR bound (default 1)
#   plot           - Logical; if TRUE generates ROC curve plots with shaded
#                    pAUC region
#   plot_type      - One of "combined", "individual", or "both"
#   plots_per_page - Integer; number of panels per page in individual layout
#   selection      - Character vector of feature names  (SE input only)
#   variable       - Name of the binary outcome column in colData (SE input only)
#   output_as_SE   - Logical; if TRUE returns a SummarizedExperiment
#
# Returns:
#   tibble with columns Variable, TpAUC, pAUC; or a SummarizedExperiment
#   with those metrics as assay rows. Always returns a numeric value in
#   [0.5, 1] — never NA.
#
# References:
#   Vivo, J.-M., Franco, M. & Vicari, D. (2018). Rethinking an ROC partial
#   area index for evaluating the classification performance at a high
#   specificity range. Advances in Data Analysis and Classification,
#   12(3), 683-704.
# ------------------------------------------------------------------------------
#' @title Compute the Tighter Partial Area Index (TpAUC)
#' @description
#'   Computes TpAUC over a restricted FPR interval. Unlike MCpAUC, TpAUC uses
#'   tighter bounds derived from the actual TPR values at the interval endpoints,
#'   making it valid for both proper and improper ROC curves. The index is always
#'   in [0.5, 1] and can distinguish crossing curves with equal raw pAUC.
#' @inheritParams MCpAUC
#' @return \code{tibble} with columns \code{Variable}, \code{TpAUC},
#'   \code{pAUC}, or a \code{SummarizedExperiment} with those metrics as
#'   assay rows. Always returns a numeric value in [0.5, 1] — never \code{NA}.
#' @references
#'   Vivo, J.-M., Franco, M. & Vicari, D. (2018). Rethinking an ROC partial
#'   area index for evaluating the classification performance at a high
#'   specificity range. \emph{Advances in Data Analysis and Classification},
#'   12(3), 683--704.
#' @examples
#' library(fission)
#' data(fission)
#' genes_of_interest = c("SPNCRNA.1080", "SPAC186.08c",
#'                       "SPNCRNA.1420",  "SPCC70.08c", "SPAC212.04c")
#' strain = SummarizedExperiment::colData(fission)$strain
#' expr   = t(SummarizedExperiment::assay(fission)[genes_of_interest, ])
#' Sp     = as.data.frame(cbind(strain = strain, expr))
#'
#' # Input type 1: SummarizedExperiment -> tibble (with plots)
#' res_Tp = TpAUC(fission, low.value = 0, up.value = 0.25,
#'                selection = genes_of_interest, variable = "strain",
#'                plot = TRUE, plot_type = "both", plots_per_page = 5)
#' res_Tp
#'
#' # Input type 2: SummarizedExperiment -> SummarizedExperiment
#' se_Tp = TpAUC(fission, low.value = 0, up.value = 0.25,
#'               selection = genes_of_interest, variable = "strain",
#'               output_as_SE = TRUE)
#' assay(se_Tp)
#'
#' # Input type 3: plain data.frame -> tibble
#' TpAUC(Sp, low.value = 0, up.value = 0.25)
#' @export
TpAUC = function(dataset, low.value = 0, up.value = 1, plot = FALSE,
                 plot_type = c("combined", "individual", "both"),
                 selection = NULL, variable = NULL,
                 plots_per_page = 2, output_as_SE = FALSE) {

  plot_type = match.arg(plot_type)
  dataset   = prepare_ROC_dataset(dataset, selection, variable)
  name.var  = colnames(dataset)
  n_pred    = ncol(dataset) - 1L

  idx_list  = vector("list", n_pred)
  pAUC_list = vector("list", n_pred)
  roc_list  = vector("list", n_pred)

  for (i in seq(2L, ncol(dataset))) {
    tmp              = cbind(dataset[, 1], dataset[, i])
    roc              = pointsCurve(tmp[, 1], tmp[, 2])
    portion          = portion_ROC_FPR(up.value, low.value, roc[, 1], roc[, 2])
    idx_list[[i-1]]  = TpA(portion[, 1], portion[, 2])
    pAUC_list[[i-1]] = pA_FPR(portion[, 1], portion[, 2])
    roc_list[[i-1]]  = list(ROC_points = as.matrix(roc))
  }

  if (isTRUE(plot))
    plot_ROC_curves(
      list(name.variable = name.var[-1], results = roc_list),
      low.value = low.value, up.value = up.value,
      plot_type = plot_type, axis = "x", plots_per_page = plots_per_page
    )

  return(format_results(idx_list, pAUC_list, name.var, output_as_SE, "TpAUC"))
}


# ------------------------------------------------------------------------------
# NpAUC
# ------------------------------------------------------------------------------
# Description:
#   Computes the normalised partial area index (NpAUC) over a restricted TPR
#   interval [low.value, up.value]. The raw horizontal-band pAUC is divided
#   by (1 - TPR0), giving a value in [0.5, 1]. Returns NA when the NLR
#   condition is not satisfied (curve is not NLR-bounded). For those cases,
#   use FpAUC instead.
#
#   Accepts three input types:
#     · data.frame           : first column = binary outcome, remaining = predictors
#     · SummarizedExperiment : use 'selection' and 'variable' to specify
#                              features and outcome column in colData
#
# Parameters:
#   dataset        - data.frame or SummarizedExperiment
#   low.value      - Numeric in [0, 1]; lower TPR bound (e.g., 0.9)
#   up.value       - Numeric in [0, 1]; upper TPR bound (default 1)
#   plot           - Logical; if TRUE generates ROC curve plots
#   plot_type      - One of "combined", "individual", or "both"
#   plots_per_page - Integer; number of panels per page in individual layout
#   selection      - Character vector of feature names  (SE input only)
#   variable       - Name of the binary outcome column in colData (SE input only)
#   output_as_SE   - Logical; if TRUE returns a SummarizedExperiment
#
# Returns:
#   tibble with columns Variable, NpAUC, pAUC; or a SummarizedExperiment
#   with those metrics as assay rows. NpAUC = NA when the NLR condition fails.
#
# References:
#   Jiang, Y., Metz, C.E. & Nishikawa, R.M. (1996). A receiver operating
#   characteristic partial area index for highly sensitive diagnostic tests.
#   Radiology, 201(3), 745-750.
# ------------------------------------------------------------------------------
#' @title Compute the Normalised Partial AUC (NpAUC)
#' @description
#'   Computes NpAUC over a restricted TPR interval [low.value, up.value].
#'   The raw horizontal-band pAUC is divided by (1 - TPR0), giving a value in
#'   [0, 1]. Returns NA when the NLR condition is not satisfied (curve is not
#'   NLR-bounded). For those cases, use FpAUC instead.
#' @param dataset      A \code{data.frame} or \code{SummarizedExperiment}.
#' @param low.value    Numeric in [0, 1]; lower \strong{TPR} bound (e.g., 0.9).
#' @param up.value     Numeric in [0, 1]; upper \strong{TPR} bound (default 1).
#' @param plot         Logical; if \code{TRUE} generates ROC curve plots.
#' @param plot_type    One of \code{"combined"}, \code{"individual"},
#'   \code{"both"}.
#' @param plots_per_page Integer; panels per page.
#' @param selection    Character vector of feature names (SE input only).
#' @param variable     Character; outcome column in \code{colData} (SE input only).
#' @param output_as_SE Logical; if \code{TRUE} returns a \code{SummarizedExperiment}.
#' @return \code{tibble} with columns \code{Variable}, \code{NpAUC},
#'   \code{pAUC}, or a \code{SummarizedExperiment} with those metrics as
#'   assay rows. \code{NpAUC = NA} when the NLR condition fails.
#' @references
#'   Jiang, Y., Metz, C.E. & Nishikawa, R.M. (1996). A receiver operating
#'   characteristic partial area index for highly sensitive diagnostic tests.
#'   \emph{Radiology}, 201(3), 745--750.
#' @examples
#' library(fission)
#' data(fission)
#' genes_of_interest = c("SPNCRNA.1080", "SPAC186.08c",
#'                       "SPNCRNA.1420",  "SPCC70.08c", "SPAC212.04c")
#' strain = SummarizedExperiment::colData(fission)$strain
#' expr   = t(SummarizedExperiment::assay(fission)[genes_of_interest, ])
#' Sp     = as.data.frame(cbind(strain = strain, expr))
#'
#' # Input type 1: SummarizedExperiment -> tibble (with plots)
#' # Note: variables where NpAUC is undefined (NLR not bounded) return NA
#' res_Np = NpAUC(fission, low.value = 0.9, up.value = 1,
#'                selection = genes_of_interest, variable = "strain",
#'                plot = TRUE, plot_type = "both", plots_per_page = 5)
#' res_Np
#'
#' # Input type 2: SummarizedExperiment -> SummarizedExperiment
#' se_Np = NpAUC(fission, low.value = 0.9, up.value = 1,
#'               selection = genes_of_interest, variable = "strain",
#'               output_as_SE = TRUE)
#' assay(se_Np)
#'
#' # Input type 3: plain data.frame -> tibble
#' NpAUC(Sp, low.value = 0.9, up.value = 1)
#' @export
NpAUC = function(dataset, low.value = 0.9, up.value = 1, plot = FALSE,
                 plot_type = c("combined", "individual", "both"),
                 selection = NULL, variable = NULL,
                 plots_per_page = 2, output_as_SE = FALSE) {

  plot_type = match.arg(plot_type)
  dataset   = prepare_ROC_dataset(dataset, selection, variable)
  name.var  = colnames(dataset)
  n_pred    = ncol(dataset) - 1L

  idx_list  = vector("list", n_pred)
  pAUC_list = vector("list", n_pred)
  roc_list  = vector("list", n_pred)

  for (i in seq(2L, ncol(dataset))) {
    tmp              = cbind(dataset[, 1], dataset[, i])
    roc              = pointsCurve(tmp[, 1], tmp[, 2])
    portion          = portion_ROC_TPR(up.value, low.value, roc[, 1], roc[, 2])
    idx_list[[i-1]]  = NpA(portion[, 2], portion[, 1])
    pAUC_list[[i-1]] = pA_TPR(portion[, 2], portion[, 1])
    roc_list[[i-1]]  = list(ROC_points = as.matrix(roc))
  }

  if (isTRUE(plot))
    plot_ROC_curves(
      list(name.variable = name.var[-1], results = roc_list),
      low.value = low.value, up.value = up.value,
      plot_type = plot_type, axis = "y", plots_per_page = plots_per_page
    )

  return(format_results(idx_list, pAUC_list, name.var, output_as_SE, "NpAUC"))
}


# ------------------------------------------------------------------------------
# FpAUC
# ------------------------------------------------------------------------------
# Description:
#   Computes the fitted partial area index (FpAUC) over a restricted TPR
#   interval [low.value, up.value]. Unlike NpAUC, FpAUC uses tighter
#   NLR-based bounds selected according to the shape of the ROC curve in the
#   region of interest (Algorithm 1 of Franco & Vivo, 2021). The index is
#   always in [0.5, 1] and never returns NA, making it valid for any ROC
#   curve shape including improper curves.
#
#   Accepts three input types:
#     · data.frame           : first column = binary outcome, remaining = predictors
#     · SummarizedExperiment : use 'selection' and 'variable' to specify
#                              features and outcome column in colData
#
# Parameters:
#   dataset        - data.frame or SummarizedExperiment
#   low.value      - Numeric in [0, 1]; lower TPR bound (e.g., 0.9)
#   up.value       - Numeric in [0, 1]; upper TPR bound (default 1)
#   plot           - Logical; if TRUE generates ROC curve plots
#   plot_type      - One of "combined", "individual", or "both"
#   plots_per_page - Integer; number of panels per page in individual layout
#   selection      - Character vector of feature names  (SE input only)
#   variable       - Name of the binary outcome column in colData (SE input only)
#   output_as_SE   - Logical; if TRUE returns a SummarizedExperiment
#
# Returns:
#   tibble with columns Variable, FpAUC, pAUC; or a SummarizedExperiment
#   with those metrics as assay rows. Always returns a numeric value in
#   [0.5, 1] — never NA.
#
# References:
#   Franco, M. & Vivo, J.-M. (2021). Evaluating the performances of
#   biomarkers over a restricted domain of high sensitivity. Mathematics,
#   9(21), 2826.
# ------------------------------------------------------------------------------
#' @title Compute the Fitted Partial Area Index (FpAUC)
#' @description
#'   Computes FpAUC over a restricted TPR interval [low.value, up.value].
#'   Unlike NpAUC, FpAUC uses tighter NLR-based bounds selected according to
#'   the shape of the ROC curve in the region (Algorithm 1 of Franco & Vivo,
#'   2021). The index is always in [0.5, 1] and never returns NA.
#' @inheritParams NpAUC
#' @return \code{tibble} with columns \code{Variable}, \code{FpAUC},
#'   \code{pAUC}, or a \code{SummarizedExperiment} with those metrics as
#'   assay rows. Always returns a numeric value in [0.5, 1] — never \code{NA}.
#' @references
#'   Franco, M. & Vivo, J.-M. (2021). Evaluating the performances of
#'   biomarkers over a restricted domain of high sensitivity.
#'   \emph{Mathematics}, 9(21), 2826.
#' @examples
#' library(fission)
#' data(fission)
#' genes_of_interest = c("SPNCRNA.1080", "SPAC186.08c",
#'                       "SPNCRNA.1420",  "SPCC70.08c", "SPAC212.04c")
#' strain = SummarizedExperiment::colData(fission)$strain
#' expr   = t(SummarizedExperiment::assay(fission)[genes_of_interest, ])
#' Sp     = as.data.frame(cbind(strain = strain, expr))
#'
#' # Input type 1: SummarizedExperiment -> tibble (with plots)
#' res_Fp = FpAUC(fission, low.value = 0.9, up.value = 1,
#'                selection = genes_of_interest, variable = "strain",
#'                plot = TRUE, plot_type = "both", plots_per_page = 5)
#' res_Fp
#'
#' # Input type 2: SummarizedExperiment -> SummarizedExperiment
#' se_Fp = FpAUC(fission, low.value = 0.9, up.value = 1,
#'               selection = genes_of_interest, variable = "strain",
#'               output_as_SE = TRUE)
#' assay(se_Fp)
#'
#' # Input type 3: plain data.frame -> tibble
#' FpAUC(Sp, low.value = 0.9, up.value = 1)
#' @export
FpAUC = function(dataset, low.value = 0.9, up.value = 1, plot = FALSE,
                 plot_type = c("combined", "individual", "both"),
                 selection = NULL, variable = NULL,
                 plots_per_page = 2, output_as_SE = FALSE) {

  plot_type = match.arg(plot_type)
  dataset   = prepare_ROC_dataset(dataset, selection, variable)
  name.var  = colnames(dataset)
  n_pred    = ncol(dataset) - 1L

  idx_list  = vector("list", n_pred)
  pAUC_list = vector("list", n_pred)
  roc_list  = vector("list", n_pred)

  for (i in seq(2L, ncol(dataset))) {
    tmp              = cbind(dataset[, 1], dataset[, i])
    roc              = pointsCurve(tmp[, 1], tmp[, 2])
    portion          = portion_ROC_TPR(up.value, low.value, roc[, 1], roc[, 2])
    idx_list[[i-1]]  = FpA(portion[, 2], portion[, 1])
    pAUC_list[[i-1]] = pA_TPR(portion[, 2], portion[, 1])
    roc_list[[i-1]]  = list(ROC_points = as.matrix(roc))
  }

  if (isTRUE(plot))
    plot_ROC_curves(
      list(name.variable = name.var[-1], results = roc_list),
      low.value = low.value, up.value = up.value,
      plot_type = plot_type, axis = "y", plots_per_page = plots_per_page
    )

  return(format_results(idx_list, pAUC_list, name.var, output_as_SE, "FpAUC"))
}

# ------------------------------------------------------------------------------
# allIndices
# ------------------------------------------------------------------------------
# Description:
#   Convenience wrapper that runs MCpAUC, TpAUC, NpAUC, and FpAUC and joins
#   the results into a single wide tibble. FPR-domain indexes (MCpAUC, TpAUC)
#   use fpr.low / fpr.up; TPR-domain indexes (NpAUC, FpAUC) use tpr.low /
#   tpr.up. This is the recommended entry point for exploratory multi-index
#   comparison.
#
#   Accepts two input types:
#     · data.frame           : first column = binary outcome, remaining = predictors
#     · SummarizedExperiment : use 'selection' and 'variable' to specify
#                              features and outcome column in colData
#
#   ROC plots are generated only for MCpAUC (FPR domain) and NpAUC (TPR
#   domain) to avoid duplicate figures; TpAUC and FpAUC share the same
#   regions respectively.
#
# Parameters:
#   dataset        - data.frame or SummarizedExperiment
#   fpr.low        - Lower FPR bound for MCpAUC and TpAUC (default 0)
#   fpr.up         - Upper FPR bound for MCpAUC and TpAUC (default 0.25)
#   tpr.low        - Lower TPR bound for NpAUC and FpAUC (default 0.9)
#   tpr.up         - Upper TPR bound for NpAUC and FpAUC (default 1)
#   plot           - Logical; if TRUE generates ROC curve plots for all indexes
#   plot_type      - One of "combined", "individual", or "both"
#   plots_per_page - Integer; number of panels per page in individual layout
#   selection      - Character vector of feature names  (SE input only)
#   variable       - Name of the binary outcome column in colData (SE input only)
#
# Returns:
#   A wide tibble with columns Variable, pAUC_FPR, MCpAUC, TpAUC,
#   pAUC_TPR, NpAUC, FpAUC.
#
# References:
#   McClish, D.K. (1989). Analyzing a portion of the ROC curve.
#   Medical Decision Making, 9(3), 190-195.
#
#   Jiang, Y., Metz, C.E. & Nishikawa, R.M. (1996). A receiver operating
#   characteristic partial area index for highly sensitive diagnostic tests.
#   Radiology, 201(3), 745-750.
#
#   Vivo, J.-M., Franco, M. & Vicari, D. (2018). Rethinking an ROC partial
#   area index for evaluating the classification performance at a high
#   specificity range. Advances in Data Analysis and Classification,
#   12(3), 683-704.
#
#   Franco, M. & Vivo, J.-M. (2021). Evaluating the performances of
#   biomarkers over a restricted domain of high sensitivity. Mathematics,
#   9(21), 2826.
# ------------------------------------------------------------------------------
#' @title Compute all four pAUC indexes for a shared dataset in one call
#' @description
#'   Convenience wrapper that runs \code{MCpAUC}, \code{TpAUC}, \code{NpAUC},
#'   and \code{FpAUC} and joins the results into a single wide \code{tibble}.
#'   FPR-domain indexes (MCpAUC, TpAUC) use \code{fpr.low} / \code{fpr.up};
#'   TPR-domain indexes (NpAUC, FpAUC) use \code{tpr.low} / \code{tpr.up}.
#'   This is the recommended entry point for exploratory multi-index comparison.
#' @param dataset      A \code{data.frame} or \code{SummarizedExperiment}.
#' @param fpr.low      Lower FPR bound for MCpAUC and TpAUC (default 0).
#' @param fpr.up       Upper FPR bound for MCpAUC and TpAUC (default 0.25).
#' @param tpr.low      Lower TPR bound for NpAUC and FpAUC (default 0.9).
#' @param tpr.up       Upper TPR bound for NpAUC and FpAUC (default 1).
#' @param plot         Logical; if \code{TRUE} generates ROC curve plots for
#'   all indexes.
#' @param plot_type    One of \code{"combined"}, \code{"individual"},
#'   \code{"both"}.
#' @param plots_per_page Integer; panels per page for individual layout.
#' @param selection    Character vector of feature names (SE input only).
#' @param variable     Character; outcome column in \code{colData} (SE input only).
#' @return A wide \code{tibble} with columns \code{Variable}, \code{pAUC_FPR},
#'   \code{MCpAUC}, \code{TpAUC}, \code{pAUC_TPR}, \code{NpAUC}, \code{FpAUC}.
#' @references
#'   McClish, D.K. (1989). Analyzing a portion of the ROC curve.
#'   \emph{Medical Decision Making}, 9(3), 190--195.
#'
#'   Thompson, M.L. & Zucchini, W. (1989). On the statistical analysis of
#'   ROC curves. \emph{Statistics in Medicine}, 8(10), 1277--1294.
#'
#'   Jiang, Y., Metz, C.E. & Nishikawa, R.M. (1996). A receiver operating
#'   characteristic partial area index for highly sensitive diagnostic tests.
#'   \emph{Radiology}, 201(3), 745--750.
#'
#'   Vivo, J.-M., Franco, M. & Vicari, D. (2018). Rethinking an ROC partial
#'   area index for evaluating the classification performance at a high
#'   specificity range. \emph{Advances in Data Analysis and Classification},
#'   12(3), 683--704.
#'
#'   Franco, M. & Vivo, J.-M. (2021). Evaluating the performances of
#'   biomarkers over a restricted domain of high sensitivity.
#'   \emph{Mathematics}, 9(21), 2826.
#' @examples
#' library(fission)
#' data(fission)
#' genes_of_interest = c("SPNCRNA.1080", "SPAC186.08c",
#'                       "SPNCRNA.1420",  "SPCC70.08c", "SPAC212.04c")
#' strain = SummarizedExperiment::colData(fission)$strain
#' expr   = t(SummarizedExperiment::assay(fission)[genes_of_interest, ])
#' Sp     = as.data.frame(cbind(strain = strain, expr))
#'
#' # Input type 1: SummarizedExperiment -> tibble (with plots)
#' allIndices(fission, fpr.low = 0, fpr.up = 0.25,
#'            tpr.low = 0.9, tpr.up = 1,
#'            selection = genes_of_interest, variable = "strain",
#'            plot = TRUE, plot_type = "both", plots_per_page = 5)
#'
#' # Input type 2: plain data.frame -> tibble
#' allIndices(Sp, fpr.low = 0, fpr.up = 0.25, tpr.low = 0.9, tpr.up = 1)
#' @export
allIndices = function(dataset,
                      fpr.low = 0,   fpr.up = 0.25,
                      tpr.low = 0.9, tpr.up = 1,
                      plot = FALSE,
                      plot_type = c("combined", "individual", "both"),
                      plots_per_page = 2,
                      selection = NULL, variable = NULL) {

  plot_type = match.arg(plot_type)

  r_MC = MCpAUC(dataset, low.value = fpr.low, up.value = fpr.up,
                plot = plot, plot_type = plot_type,
                plots_per_page = plots_per_page,
                selection = selection, variable = variable)

  r_Tp = TpAUC(dataset, low.value = fpr.low, up.value = fpr.up,
               plot = FALSE,
               selection = selection, variable = variable)

  r_Np = NpAUC(dataset, low.value = tpr.low, up.value = tpr.up,
               plot = plot, plot_type = plot_type,
               plots_per_page = plots_per_page,
               selection = selection, variable = variable)

  r_Fp = FpAUC(dataset, low.value = tpr.low, up.value = tpr.up,
               plot = FALSE,
               selection = selection, variable = variable)

  return(
    r_MC |>
      dplyr::rename(pAUC_FPR = pAUC) |>
      dplyr::left_join(dplyr::select(r_Tp, Variable, TpAUC),                   by = "Variable") |>
      dplyr::left_join(dplyr::select(r_Np, Variable, NpAUC, pAUC_TPR = pAUC), by = "Variable") |>
      dplyr::left_join(dplyr::select(r_Fp, Variable, FpAUC),                   by = "Variable") |>
      dplyr::select(Variable, pAUC_FPR, MCpAUC, TpAUC, pAUC_TPR, NpAUC, FpAUC)
  )
}


# ==============================================================================
# Bootstrap confidence interval functions
# ==============================================================================


# ------------------------------------------------------------------------------
# MCpAUCboot
# ------------------------------------------------------------------------------
# Description:
#   Estimates confidence intervals for the McClish standardised pAUC (MCpAUC)
#   via nonparametric bootstrap resampling. Variables for which MCpAUC is
#   undefined on the original sample (improper ROC curves) are automatically
#   skipped and return NA in all output columns.
#
#   The full bootstrap lifecycle is handled by run_bootstrap():
#     · Point estimates are computed on the original sample for all predictors.
#     · boot::boot() is called only for variables with a valid estimate.
#     · CIs are extracted via boot::boot.ci() using the requested method.
#     · Optional diagnostic plots show estimate, bootstrap mean, and CI bounds.
#
#   Accepts three input types:
#     · data.frame           : first column = binary outcome, remaining = predictors
#     · SummarizedExperiment : use 'selection' and 'variable' to specify
#                              features and outcome column in colData
#
# Parameters:
#   dataset        - data.frame or SummarizedExperiment
#   low.value      - Numeric in [0, 1]; lower FPR bound (default 0)
#   up.value       - Numeric in [0, 1]; upper FPR bound (default 1)
#   r              - Integer; number of bootstrap replicates (default 2000)
#   level          - Numeric in (0, 1); confidence level (default 0.95)
#   type.interval  - CI method: "norm", "basic", "stud", "perc", or "bca"
#   selection      - Character vector of feature names  (SE input only)
#   variable       - Name of the binary outcome column in colData (SE input only)
#   output_as_SE   - Logical; if TRUE returns a SummarizedExperiment
#   plot           - Logical; if TRUE generates bootstrap diagnostic plots
#   plot_type      - One of "both", "ci", or "distribution"
#   parallel       - Parallelisation method: "no", "multicore", or "snow"
#   ncpus          - Integer; number of CPU cores when parallel != "no"
#
# Returns:
#   tibble with columns Variable, MCpAUC, bias, sd, lwr, upr; or a
#   SummarizedExperiment with those metrics as assay rows.
#   All columns are NA for variables with undefined MCpAUC.
#
# References:
#   McClish, D.K. (1989). Analyzing a portion of the ROC curve.
#   Medical Decision Making, 9(3), 190-195.
#
#   Davison, A.C. & Hinkley, D.V. (1997). Bootstrap Methods and Their
#   Applications. Cambridge: Cambridge University Press.
# ------------------------------------------------------------------------------
#' @title Bootstrap Confidence Intervals for MCpAUC (MCpAUCboot)
#' @description
#'   Estimates confidence intervals for the McClish standardised pAUC (MCpAUC)
#'   via nonparametric bootstrap resampling. Variables for which MCpAUC is
#'   undefined on the original sample (improper ROC curves) are automatically
#'   skipped and return NA in all output columns.
#' @inheritParams MCpAUC
#' @param r             Integer; number of bootstrap replicates (default 2000).
#'   Use at least 1000 for reliable CI estimates; 10000 for publication results.
#' @param level         Numeric in (0, 1); confidence level (default 0.95).
#' @param type.interval Character; CI method passed to \code{boot::boot.ci()}.
#'   One of \code{"norm"} (normal approximation), \code{"basic"},
#'   \code{"stud"} (studentised), \code{"perc"} (percentile, recommended),
#'   or \code{"bca"} (bias-corrected and accelerated).
#' @param parallel      Parallelisation method: \code{"no"} (default, sequential),
#'   \code{"multicore"} (fork-based, Linux/macOS only), or
#'   \code{"snow"} (socket cluster, all platforms).
#' @param ncpus         Integer; number of CPU cores when \code{parallel != "no"}.
#' @param plot          Logical; if \code{TRUE} generates bootstrap diagnostic
#'   plots.
#' @param plot_type     One of \code{"both"} (CI panel + distribution boxplot),
#'   \code{"ci"} (CI panel only), or \code{"distribution"} (boxplot only).
#' @return \code{tibble} with columns \code{Variable}, \code{MCpAUC} (original
#'   estimate), \code{bias}, \code{sd}, \code{lwr}, \code{upr}; or a
#'   \code{SummarizedExperiment} with those metrics as assay rows.
#'   All columns are \code{NA} for variables with undefined MCpAUC.
#' @references
#'   McClish, D.K. (1989). Analyzing a portion of the ROC curve.
#'   \emph{Medical Decision Making}, 9(3), 190--195.
#'
#'   Davison, A.C. & Hinkley, D.V. (1997). \emph{Bootstrap Methods and Their
#'   Applications}. Cambridge: Cambridge University Press.
#' @examples
#' library(fission)
#' data(fission)
#' genes_of_interest = c("SPNCRNA.1080", "SPAC186.08c",
#'                       "SPNCRNA.1420",  "SPCC70.08c", "SPAC212.04c")
#' strain = SummarizedExperiment::colData(fission)$strain
#' expr   = t(SummarizedExperiment::assay(fission)[genes_of_interest, ])
#' Sp     = as.data.frame(cbind(strain = strain, expr))
#'
#' # Input type 1: SummarizedExperiment -> tibble (with diagnostic plots)
#' boot_MC = MCpAUCboot(fission, low.value = 0, up.value = 0.25,
#'                      selection = genes_of_interest, variable = "strain",
#'                      r = 100, plot = TRUE, plot_type = "both")
#' boot_MC
#'
#' # Input type 2: SummarizedExperiment -> SummarizedExperiment
#' boot_MC_se = MCpAUCboot(fission, low.value = 0, up.value = 0.25,
#'                         selection = genes_of_interest, variable = "strain",
#'                         r = 100, output_as_SE = TRUE)
#' assay(boot_MC_se)
#'
#' # Input type 3: plain data.frame -> tibble
#' MCpAUCboot(Sp, low.value = 0, up.value = 0.25, r = 100)
#' @export
MCpAUCboot = function(dataset, low.value = 0, up.value = 1,
                      r = 2000, level = 0.95, type.interval = "perc",
                      selection = NULL, variable = NULL,
                      output_as_SE = FALSE, plot = FALSE,
                      plot_type = c("both", "ci", "distribution"),
                      parallel = c("no", "multicore", "snow"), ncpus = 1L) {

  plot_type = match.arg(plot_type)
  parallel  = match.arg(parallel)
  stopifnot(
    methods::is(dataset, "data.frame") || methods::is(dataset, "SummarizedExperiment"),
    is.numeric(low.value), low.value >= 0, low.value <= 1,
    is.numeric(up.value),  up.value  >= 0, up.value  <= 1,
    low.value < up.value,
    is.numeric(level), level > 0, level < 1,
    is.numeric(r), r > 1,
    type.interval %in% c("norm", "perc", "basic", "stud", "bca")
  )

  dataset  = prepare_ROC_dataset(dataset, selection, variable)
  name.var = colnames(dataset)[-1]

  return(run_bootstrap(
    dataset       = dataset,      boot_fn     = fbootM,
    low.value     = low.value,    up.value    = up.value,
    r             = r,            level       = level,
    type.interval = type.interval,
    name.variable = name.var,     index_label = "MCpAUC",
    output_as_SE  = output_as_SE, plot        = plot, plot_type = plot_type,
    parallel      = parallel,     ncpus       = ncpus
  ))
}


# ------------------------------------------------------------------------------
# TpAUCboot
# ------------------------------------------------------------------------------
# Description:
#   Estimates confidence intervals for the tighter partial area index (TpAUC)
#   via nonparametric bootstrap resampling. Because TpAUC is valid for any
#   ROC curve shape, no variables are skipped and no NA values arise from
#   shape conditions.
#
#   The full bootstrap lifecycle is handled by run_bootstrap():
#     · Point estimates are computed on the original sample for all predictors.
#     · boot::boot() is called for all variables (no skipping).
#     · CIs are extracted via boot::boot.ci() using the requested method.
#     · Optional diagnostic plots show estimate, bootstrap mean, and CI bounds.
#
#   Accepts three input types:
#     · data.frame           : first column = binary outcome, remaining = predictors
#     · SummarizedExperiment : use 'selection' and 'variable' to specify
#                              features and outcome column in colData
#
# Parameters:
#   dataset        - data.frame or SummarizedExperiment
#   low.value      - Numeric in [0, 1]; lower FPR bound (default 0)
#   up.value       - Numeric in [0, 1]; upper FPR bound (default 1)
#   r              - Integer; number of bootstrap replicates (default 2000)
#   level          - Numeric in (0, 1); confidence level (default 0.95)
#   type.interval  - CI method: "norm", "basic", "stud", "perc", or "bca"
#   selection      - Character vector of feature names  (SE input only)
#   variable       - Name of the binary outcome column in colData (SE input only)
#   output_as_SE   - Logical; if TRUE returns a SummarizedExperiment
#   plot           - Logical; if TRUE generates bootstrap diagnostic plots
#   plot_type      - One of "both", "ci", or "distribution"
#   parallel       - Parallelisation method: "no", "multicore", or "snow"
#   ncpus          - Integer; number of CPU cores when parallel != "no"
#
# Returns:
#   tibble with columns Variable, TpAUC, bias, sd, lwr, upr; or a
#   SummarizedExperiment with those metrics as assay rows.
#   Always returns numeric values — never NA from shape conditions.
#
# References:
#   Vivo, J.-M., Franco, M. & Vicari, D. (2018). Rethinking an ROC partial
#   area index for evaluating the classification performance at a high
#   specificity range. Advances in Data Analysis and Classification,
#   12(3), 683-704.
#
#   Davison, A.C. & Hinkley, D.V. (1997). Bootstrap Methods and Their
#   Applications. Cambridge: Cambridge University Press.
# ------------------------------------------------------------------------------
#' @title Bootstrap Confidence Intervals for TpAUC (TpAUCboot)
#' @description
#'   Estimates confidence intervals for the tighter partial area index (TpAUC)
#'   via nonparametric bootstrap resampling. Because TpAUC is valid for any ROC
#'   curve shape, no variables are skipped (no NA from shape conditions).
#' @inheritParams MCpAUCboot
#' @return \code{tibble} with columns \code{Variable}, \code{TpAUC},
#'   \code{bias}, \code{sd}, \code{lwr}, \code{upr}; or a
#'   \code{SummarizedExperiment} with those metrics as assay rows.
#'   Always returns numeric values — never \code{NA} from shape conditions.
#' @references
#'   Vivo, J.-M., Franco, M. & Vicari, D. (2018). Rethinking an ROC partial
#'   area index for evaluating the classification performance at a high
#'   specificity range. \emph{Advances in Data Analysis and Classification},
#'   12(3), 683--704.
#'
#'   Davison, A.C. & Hinkley, D.V. (1997). \emph{Bootstrap Methods and Their
#'   Applications}. Cambridge: Cambridge University Press.
#' @examples
#' library(fission)
#' data(fission)
#' genes_of_interest = c("SPNCRNA.1080", "SPAC186.08c",
#'                       "SPNCRNA.1420",  "SPCC70.08c", "SPAC212.04c")
#' strain = SummarizedExperiment::colData(fission)$strain
#' expr   = t(SummarizedExperiment::assay(fission)[genes_of_interest, ])
#' Sp     = as.data.frame(cbind(strain = strain, expr))
#'
#' # Input type 1: SummarizedExperiment -> tibble (with diagnostic plots)
#' boot_Tp = TpAUCboot(fission, low.value = 0, up.value = 0.25,
#'                     selection = genes_of_interest, variable = "strain",
#'                     r = 100, plot = TRUE, plot_type = "both")
#' boot_Tp
#'
#' # Input type 2: plain data.frame -> tibble
#' TpAUCboot(Sp, low.value = 0, up.value = 0.25, r = 100)
#' @export
TpAUCboot = function(dataset, low.value = 0, up.value = 1,
                     r = 2000, level = 0.95, type.interval = "perc",
                     selection = NULL, variable = NULL,
                     output_as_SE = FALSE, plot = FALSE,
                     plot_type = c("both", "ci", "distribution"),
                     parallel = c("no", "multicore", "snow"), ncpus = 1L) {

  plot_type = match.arg(plot_type)
  parallel  = match.arg(parallel)
  stopifnot(
    methods::is(dataset, "data.frame") || methods::is(dataset, "SummarizedExperiment"),
    is.numeric(low.value), low.value >= 0, low.value <= 1,
    is.numeric(up.value),  up.value  >= 0, up.value  <= 1,
    low.value < up.value,
    is.numeric(level), level > 0, level < 1,
    is.numeric(r), r > 1,
    type.interval %in% c("norm", "perc", "basic", "stud", "bca")
  )

  dataset  = prepare_ROC_dataset(dataset, selection, variable)
  name.var = colnames(dataset)[-1]

  return(run_bootstrap(
    dataset       = dataset,      boot_fn     = fbootT,
    low.value     = low.value,    up.value    = up.value,
    r             = r,            level       = level,
    type.interval = type.interval,
    name.variable = name.var,     index_label = "TpAUC",
    output_as_SE  = output_as_SE, plot        = plot, plot_type = plot_type,
    parallel      = parallel,     ncpus       = ncpus
  ))
}


# ------------------------------------------------------------------------------
# NpAUCboot
# ------------------------------------------------------------------------------
# Description:
#   Estimates confidence intervals for the normalised partial AUC (NpAUC) via
#   nonparametric bootstrap resampling. Variables for which NpAUC is undefined
#   (NLR condition not satisfied) on the original sample are automatically
#   skipped and return NA in all output columns. Consider FpAUCboot for those
#   variables.
#
#   The full bootstrap lifecycle is handled by run_bootstrap():
#     · Point estimates are computed on the original sample for all predictors.
#     · boot::boot() is called only for variables with a valid estimate.
#     · CIs are extracted via boot::boot.ci() using the requested method.
#     · Optional diagnostic plots show estimate, bootstrap mean, and CI bounds.
#
#   Accepts three input types:
#     · data.frame           : first column = binary outcome, remaining = predictors
#     · SummarizedExperiment : use 'selection' and 'variable' to specify
#                              features and outcome column in colData
#
# Parameters:
#   dataset        - data.frame or SummarizedExperiment
#   low.value      - Numeric in [0, 1]; lower TPR bound (e.g., 0.9)
#   up.value       - Numeric in [0, 1]; upper TPR bound (default 1)
#   r              - Integer; number of bootstrap replicates (default 2000)
#   level          - Numeric in (0, 1); confidence level (default 0.95)
#   type.interval  - CI method: "norm", "basic", "stud", "perc", or "bca"
#   selection      - Character vector of feature names  (SE input only)
#   variable       - Name of the binary outcome column in colData (SE input only)
#   output_as_SE   - Logical; if TRUE returns a SummarizedExperiment
#   plot           - Logical; if TRUE generates bootstrap diagnostic plots
#   plot_type      - One of "both", "ci", or "distribution"
#   parallel       - Parallelisation method: "no", "multicore", or "snow"
#   ncpus          - Integer; number of CPU cores when parallel != "no"
#
# Returns:
#   tibble with columns Variable, NpAUC, bias, sd, lwr, upr; or a
#   SummarizedExperiment with those metrics as assay rows.
#   All columns are NA for variables where the NLR condition is not satisfied.
#
# References:
#   Jiang, Y., Metz, C.E. & Nishikawa, R.M. (1996). A receiver operating
#   characteristic partial area index for highly sensitive diagnostic tests.
#   Radiology, 201(3), 745-750.
#
#   Davison, A.C. & Hinkley, D.V. (1997). Bootstrap Methods and Their
#   Applications. Cambridge: Cambridge University Press.
# ------------------------------------------------------------------------------
#' @title Bootstrap Confidence Intervals for NpAUC (NpAUCboot)
#' @description
#'   Estimates confidence intervals for the normalised partial AUC (NpAUC) via
#'   nonparametric bootstrap resampling. Variables for which NpAUC is undefined
#'   (NLR condition not satisfied) on the original sample are automatically
#'   skipped and return NA in all output columns. Consider FpAUCboot for those
#'   variables.
#' @inheritParams MCpAUCboot
#' @param low.value Numeric in [0, 1]; lower bound of the \strong{TPR} interval
#'   (e.g., 0.9).
#' @param up.value  Numeric in [0, 1]; upper bound of the \strong{TPR} interval
#'   (default 1).
#' @return \code{tibble} with columns \code{Variable}, \code{NpAUC},
#'   \code{bias}, \code{sd}, \code{lwr}, \code{upr}; or a
#'   \code{SummarizedExperiment} with those metrics as assay rows.
#'   All columns are \code{NA} for variables where the NLR condition fails.
#' @references
#'   Jiang, Y., Metz, C.E. & Nishikawa, R.M. (1996). A receiver operating
#'   characteristic partial area index for highly sensitive diagnostic tests.
#'   \emph{Radiology}, 201(3), 745--750.
#'
#'   Davison, A.C. & Hinkley, D.V. (1997). \emph{Bootstrap Methods and Their
#'   Applications}. Cambridge: Cambridge University Press.
#' @examples
#' library(fission)
#' data(fission)
#' genes_of_interest = c("SPNCRNA.1080", "SPAC186.08c",
#'                       "SPNCRNA.1420",  "SPCC70.08c", "SPAC212.04c")
#' strain = SummarizedExperiment::colData(fission)$strain
#' expr   = t(SummarizedExperiment::assay(fission)[genes_of_interest, ])
#' Sp     = as.data.frame(cbind(strain = strain, expr))
#'
#' # Input type 1: SummarizedExperiment -> tibble (with diagnostic plots)
#' # Note: variables where NpAUC is undefined (NLR condition not satisfied) return NA
#' boot_Np = NpAUCboot(fission, low.value = 0.9, up.value = 1,
#'                     selection = genes_of_interest, variable = "strain",
#'                     r = 100, plot = TRUE, plot_type = "both")
#' boot_Np
#'
#' # Input type 2: plain data.frame -> tibble
#' NpAUCboot(Sp, low.value = 0.9, up.value = 1, r = 100)
#' @export
NpAUCboot = function(dataset, low.value = 0.9, up.value = 1,
                     r = 2000, level = 0.95, type.interval = "perc",
                     selection = NULL, variable = NULL,
                     output_as_SE = FALSE, plot = FALSE,
                     plot_type = c("both", "ci", "distribution"),
                     parallel = c("no", "multicore", "snow"), ncpus = 1L) {

  plot_type = match.arg(plot_type)
  parallel  = match.arg(parallel)
  stopifnot(
    methods::is(dataset, "data.frame") || methods::is(dataset, "SummarizedExperiment"),
    is.numeric(low.value), low.value >= 0, low.value <= 1,
    is.numeric(up.value),  up.value  >= 0, up.value  <= 1,
    low.value < up.value,
    is.numeric(level), level > 0, level < 1,
    is.numeric(r), r > 1,
    type.interval %in% c("norm", "perc", "basic", "stud", "bca")
  )

  dataset  = prepare_ROC_dataset(dataset, selection, variable)
  name.var = colnames(dataset)[-1]

  return(run_bootstrap(
    dataset       = dataset,      boot_fn     = fbootN,
    low.value     = low.value,    up.value    = up.value,
    r             = r,            level       = level,
    type.interval = type.interval,
    name.variable = name.var,     index_label = "NpAUC",
    output_as_SE  = output_as_SE, plot        = plot, plot_type = plot_type,
    parallel      = parallel,     ncpus       = ncpus
  ))
}


# ------------------------------------------------------------------------------
# FpAUCboot
# ------------------------------------------------------------------------------
# Description:
#   Estimates confidence intervals for the fitted partial area index (FpAUC)
#   via nonparametric bootstrap resampling. Because FpAUC is valid for any ROC
#   curve shape, no variables are skipped. This makes FpAUCboot the most
#   robust bootstrap function in the package, especially useful when NpAUCboot
#   returns NA for some variables.
#
#   The full bootstrap lifecycle is handled by run_bootstrap():
#     · Point estimates are computed on the original sample for all predictors.
#     · boot::boot() is called for all variables (no skipping).
#     · CIs are extracted via boot::boot.ci() using the requested method.
#     · Optional diagnostic plots show estimate, bootstrap mean, and CI bounds.
#
#   Accepts three input types:
#     · data.frame           : first column = binary outcome, remaining = predictors
#     · SummarizedExperiment : use 'selection' and 'variable' to specify
#                              features and outcome column in colData
#
# Parameters:
#   dataset        - data.frame or SummarizedExperiment
#   low.value      - Numeric in [0, 1]; lower TPR bound (e.g., 0.9)
#   up.value       - Numeric in [0, 1]; upper TPR bound (default 1)
#   r              - Integer; number of bootstrap replicates (default 2000)
#   level          - Numeric in (0, 1); confidence level (default 0.95)
#   type.interval  - CI method: "norm", "basic", "stud", "perc", or "bca"
#   selection      - Character vector of feature names  (SE input only)
#   variable       - Name of the binary outcome column in colData (SE input only)
#   output_as_SE   - Logical; if TRUE returns a SummarizedExperiment
#   plot           - Logical; if TRUE generates bootstrap diagnostic plots
#   plot_type      - One of "both", "ci", or "distribution"
#   parallel       - Parallelisation method: "no", "multicore", or "snow"
#   ncpus          - Integer; number of CPU cores when parallel != "no"
#
# Returns:
#   tibble with columns Variable, FpAUC, bias, sd, lwr, upr; or a
#   SummarizedExperiment with those metrics as assay rows.
#   Always returns numeric values — never NA from shape conditions.
#
# References:
#   Franco, M. & Vivo, J.-M. (2021). Evaluating the performances of
#   biomarkers over a restricted domain of high sensitivity. Mathematics,
#   9(21), 2826.
#
#   Davison, A.C. & Hinkley, D.V. (1997). Bootstrap Methods and Their
#   Applications. Cambridge: Cambridge University Press.
# ------------------------------------------------------------------------------
#' @title Bootstrap Confidence Intervals for FpAUC (FpAUCboot)
#' @description
#'   Estimates confidence intervals for the fitted partial area index (FpAUC)
#'   via nonparametric bootstrap resampling. Because FpAUC is valid for any ROC
#'   curve shape, no variables are skipped. This makes FpAUCboot the most
#'   robust bootstrap function in the package, especially useful when NpAUCboot
#'   returns NA for some variables.
#' @inheritParams NpAUCboot
#' @return \code{tibble} with columns \code{Variable}, \code{FpAUC},
#'   \code{bias}, \code{sd}, \code{lwr}, \code{upr}; or a
#'   \code{SummarizedExperiment} with those metrics as assay rows.
#'   Always returns numeric values — never \code{NA} from shape conditions.
#' @references
#'   Franco, M. & Vivo, J.-M. (2021). Evaluating the performances of
#'   biomarkers over a restricted domain of high sensitivity.
#'   \emph{Mathematics}, 9(21), 2826.
#'
#'   Davison, A.C. & Hinkley, D.V. (1997). \emph{Bootstrap Methods and Their
#'   Applications}. Cambridge: Cambridge University Press.
#' @examples
#' library(fission)
#' data(fission)
#' genes_of_interest = c("SPNCRNA.1080", "SPAC186.08c",
#'                       "SPNCRNA.1420",  "SPCC70.08c", "SPAC212.04c")
#' strain = SummarizedExperiment::colData(fission)$strain
#' expr   = t(SummarizedExperiment::assay(fission)[genes_of_interest, ])
#' Sp     = as.data.frame(cbind(strain = strain, expr))
#'
#' # Input type 1: SummarizedExperiment -> tibble (with diagnostic plots)
#' # FpAUCboot is valid for any curve shape; no variables are ever skipped
#' boot_Fp = FpAUCboot(fission, low.value = 0.9, up.value = 1,
#'                     selection = genes_of_interest, variable = "strain",
#'                     r = 100, plot = TRUE, plot_type = "both")
#' boot_Fp
#'
#' # Input type 2: SummarizedExperiment -> SummarizedExperiment
#' boot_Fp_se = FpAUCboot(fission, low.value = 0.9, up.value = 1,
#'                        selection = genes_of_interest, variable = "strain",
#'                        r = 100, output_as_SE = TRUE)
#' assay(boot_Fp_se)
#'
#' # Input type 3: plain data.frame -> tibble
#' FpAUCboot(Sp, low.value = 0.9, up.value = 1, r = 100)
#' @export
FpAUCboot = function(dataset, low.value = 0.9, up.value = 1,
                     r = 2000, level = 0.95, type.interval = "perc",
                     selection = NULL, variable = NULL,
                     output_as_SE = FALSE, plot = FALSE,
                     plot_type = c("both", "ci", "distribution"),
                     parallel = c("no", "multicore", "snow"), ncpus = 1L) {

  plot_type = match.arg(plot_type)
  parallel  = match.arg(parallel)
  stopifnot(
    methods::is(dataset, "data.frame") || methods::is(dataset, "SummarizedExperiment"),
    is.numeric(low.value), low.value >= 0, low.value <= 1,
    is.numeric(up.value),  up.value  >= 0, up.value  <= 1,
    low.value < up.value,
    is.numeric(level), level > 0, level < 1,
    is.numeric(r), r > 1,
    type.interval %in% c("norm", "perc", "basic", "stud", "bca")
  )

  dataset  = prepare_ROC_dataset(dataset, selection, variable)
  name.var = colnames(dataset)[-1]

  return(run_bootstrap(
    dataset       = dataset,      boot_fn     = fbootF,
    low.value     = low.value,    up.value    = up.value,
    r             = r,            level       = level,
    type.interval = type.interval,
    name.variable = name.var,     index_label = "FpAUC",
    output_as_SE  = output_as_SE, plot        = plot, plot_type = plot_type,
    parallel      = parallel,     ncpus       = ncpus
  ))
}
