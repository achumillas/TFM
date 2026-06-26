# ==============================================================================
# 07_format_results.R
# ==============================================================================
# Output formatting utilities.
#
# Functions in this file:
#   createSE           - Build a SummarizedExperiment from a named list of
#                        metric vectors.
#   format_results     - Format point-estimate results as tibble or SE.
#   format_results_boot - Format bootstrap results as tibble or SE.
#
# These three functions are the final step in the computation pipeline: every
# exported function ends by calling either format_results or
# format_results_boot, which in turn may call createSE when output_as_SE = TRUE.
# ==============================================================================


# ------------------------------------------------------------------------------
# createSE
# ------------------------------------------------------------------------------
# Description:
#   Constructs a SummarizedExperiment from a named list of metric vectors.
#   Each list element becomes one row (metric) in the assay matrix; the
#   predictor names become the columns. rowData stores metric names; colData
#   stores predictor names. Used by format_results and format_results_boot
#   when output_as_SE = TRUE.
#
# Parameters:
#   object    - Named list of numeric vectors (one per metric, all same length)
#   var_names - Character vector of predictor/variable names (column labels)
#
# Returns:
#   A SummarizedExperiment with one assay called "metrics".
# ------------------------------------------------------------------------------

createSE = function(object, var_names) {

  # Validate: input must be a named list and all vectors must have the same length
  if (!is.list(object))
    stop("Expected a named list of metric vectors.")
  if (any(lengths(object) != length(var_names)))
    stop("Every metric vector must have length == number of variables.")

  # Stack the list elements as rows of a matrix (metrics x variables)
  mat           = do.call(rbind, object)
  rownames(mat) = names(object)     # metric names as row labels
  colnames(mat) = var_names         # predictor names as column labels

  # Wrap in a SummarizedExperiment with metric and variable metadata
  return(SummarizedExperiment::SummarizedExperiment(
    assays  = list(metrics = mat),
    rowData = S4Vectors::DataFrame(metric   = rownames(mat)),
    colData = S4Vectors::DataFrame(variable = var_names)
  ))
}


# ------------------------------------------------------------------------------
# format_results
# ------------------------------------------------------------------------------
# Description:
#   Formats point-estimate results (index values and raw pAUC) into the
#   output type requested by the user: tibble (default) or SummarizedExperiment.
#
#   The function is robust to whether 'name.variable' includes the outcome
#   column (colnames(dataset), length = n_pred + 1) or only the predictor
#   columns (colnames(dataset)[-1], length = n_pred).
#
# Parameters:
#   index_values  - List/vector of computed index values (one per predictor)
#   pAUC_values   - List/vector of raw pAUC values (one per predictor)
#   name.variable - Column names of the dataset (with or without outcome column)
#   output_as_SE  - Logical; if TRUE return a SummarizedExperiment
#   index_label   - Name of the index column (e.g., "MCpAUC", "TpAUC")
#
# Returns:
#   tibble with columns Variable, <index_label>, pAUC
#   OR SummarizedExperiment with two metric rows.
# ------------------------------------------------------------------------------

format_results = function(index_values, pAUC_values,
                          name.variable, output_as_SE,
                          index_label = "MCpAUC") {

  # Drop the outcome column name if it is included in name.variable
  # (a length mismatch between name.variable and index_values signals its presence)
  var_names = if (length(name.variable) > length(index_values)) {
    name.variable[-1]   # remove the first element (outcome column name)
  } else {
    name.variable       # already contains predictor names only
  }

  # Sanitise values: replace any length-zero or non-finite scalars with NA
  clean_val   = function(x) if (length(x) == 1L && is.finite(x)) x else NA_real_
  index_clean = vapply(index_values, clean_val, numeric(1))
  pAUC_clean  = vapply(pAUC_values,  clean_val, numeric(1))

  # Return a SummarizedExperiment if requested
  if (isTRUE(output_as_SE)) {
    obj        = list(index_clean, pAUC_clean)
    names(obj) = c(index_label, "pAUC")
    return(createSE(obj, var_names))
  }

  # Default: return a tidy tibble with Variable, index and raw pAUC columns
  tbl = tibble::tibble(
    Variable = var_names,
    index_   = index_clean,
    pAUC     = pAUC_clean
  )
  names(tbl)[2] = index_label
  return(tbl)
}


# ------------------------------------------------------------------------------
# format_results_boot
# ------------------------------------------------------------------------------
# Description:
#   Formats bootstrap results (point estimate, bias, standard deviation, and
#   confidence interval bounds) into the output type requested by the user.
#
# Parameters:
#   index_vals    - Numeric vector of original point estimates (t0)
#   bias_vals     - Numeric vector of bootstrap biases (mean(t*) - t0)
#   sd_vals       - Numeric vector of bootstrap standard deviations
#   lwr_vals      - Numeric vector of lower CI bounds
#   upr_vals      - Numeric vector of upper CI bounds
#   name.variable - Character vector of predictor names (predictors only)
#   output_as_SE  - Logical; if TRUE return a SummarizedExperiment
#   index_label   - Name of the index (e.g., "MCpAUC")
#
# Returns:
#   tibble with columns Variable, <index_label>, bias, sd, lwr, upr
#   OR SummarizedExperiment with five metric rows.
# ------------------------------------------------------------------------------

format_results_boot = function(index_vals, bias_vals, sd_vals,
                               lwr_vals, upr_vals,
                               name.variable, output_as_SE = FALSE,
                               index_label = "MCpAUC") {

  # Build a named list with all bootstrap metrics for easy conversion
  obj = list()
  obj[[index_label]] = index_vals   # original point estimate (t0)
  obj[["bias"]]      = bias_vals    # bootstrap bias = mean(t*) - t0
  obj[["sd"]]        = sd_vals      # bootstrap standard deviation
  obj[["lwr"]]       = lwr_vals     # lower bound of the confidence interval
  obj[["upr"]]       = upr_vals     # upper bound of the confidence interval

  # Return a SummarizedExperiment if requested
  if (isTRUE(output_as_SE)) return(createSE(obj, name.variable))

  # Default: return a tidy tibble with all bootstrap summary columns
  tbl = tibble::tibble(
    Variable = name.variable,
    index_   = index_vals,
    bias     = bias_vals,
    sd       = sd_vals,
    lwr      = lwr_vals,
    upr      = upr_vals
  )
  names(tbl)[2] = index_label
  return(tbl)
}
