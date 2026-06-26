# ==============================================================================
# 05_bootstrap_stats.R
# ==============================================================================
# Bootstrap statistic functions and the centralised bootstrap runner.
#
# Functions in this file:
#   fbootM        - Bootstrap statistic for MCpAUC (FPR domain).
#   fbootT        - Bootstrap statistic for TpAUC  (FPR domain).
#   fbootN        - Bootstrap statistic for NpAUC  (TPR domain).
#   fbootF        - Bootstrap statistic for FpAUC  (TPR domain).
#   run_bootstrap - Centralised runner shared by all four *boot exported
#                   functions. Handles NA-skipping, boot::boot(), CI
#                   extraction, plotting, and output formatting.
#
# The four fboot* functions follow the signature required by boot::boot():
#   statistic(data, indices, ...). Each resamples the dataset with the
#   provided indices and returns a numeric vector of index values — one per
#   predictor column.
#
# All four fboot* functions use tryCatch to guard against degenerate bootstrap
# resamples. Common failure modes include: single-class outcome after resampling,
# constant predictor values, FPR/TPR interval that no longer intersects the
# resampled ROC curve, or a constant ROC curve (all TPR values identical).
# Any such replicate is set to NA rather than propagating an invalid value.
#
# run_bootstrap is placed here (rather than in 07_format_results.R) because
# its primary role is orchestrating the bootstrap resampling pipeline; the
# final call to format_results_boot is simply its last output step.
# ==============================================================================


# ------------------------------------------------------------------------------
# fbootM
# ------------------------------------------------------------------------------
# Description:
#   Bootstrap statistic function for MCpAUC. Resamples the dataset using the
#   bootstrap indices provided by boot::boot(), computes the ROC curve for each
#   predictor, extracts the FPR sub-interval, and returns the MCpAUC value for
#   every predictor in the resampled dataset.
#
# Parameters:
#   dataset   - Original data.frame (outcome in col 1, predictors in cols 2+)
#   bssample  - Integer vector of row indices drawn by the bootstrap procedure
#   low.limit - Lower bound of the FPR interval
#   up.limit  - Upper bound of the FPR interval
#
# Returns:
#   Numeric vector of length ncol(dataset) - 1 (one MCpAUC per predictor).
#   NA for any predictor where the index is undefined on the resample.
# ------------------------------------------------------------------------------

fbootM = function(dataset, bssample, low.limit, up.limit) {

  # Resample the dataset using the bootstrap indices provided by boot::boot()
  bsdata = dataset[bssample, ]

  # Loop over each predictor column (skip column 1 = outcome), compute MCpAUC
  # tryCatch guards against degenerate resamples (single-class outcome, constant
  # predictor, or FPR interval not intersecting the resampled ROC curve)
  vapply(seq(2, ncol(bsdata)), function(i) {
    tmp = cbind(bsdata[, 1], bsdata[, i])
    val = tryCatch({
      roc     = pointsCurve(tmp[, 1], tmp[, 2])
      portion = portion_ROC_FPR(up.limit, low.limit, roc[, 1], roc[, 2])
      MCpA(portion[, 2], portion[, 1])
    }, warning = function(w) NA_real_,
       error   = function(e) NA_real_)
    # Additional guard: reject any non-finite value that slipped through
    if (!is.finite(val)) NA_real_ else val
  }, numeric(1))
}


# ------------------------------------------------------------------------------
# fbootT
# ------------------------------------------------------------------------------
# Description:
#   Bootstrap statistic function for TpAUC. Same pipeline as fbootM but
#   computes TpA instead of MCpA on the FPR sub-interval.
#   TpA issues a warning when the bootstrap resample produces a constant ROC
#   curve (all TPR values identical); that warning is captured and the replicate
#   is set to NA to avoid propagating an unstable value.
#
# Parameters:
#   dataset   - Original data.frame (outcome in col 1, predictors in cols 2+)
#   bssample  - Integer vector of row indices
#   low.limit - Lower bound of the FPR interval
#   up.limit  - Upper bound of the FPR interval
#
# Returns:
#   Numeric vector of length ncol(dataset) - 1 (one TpAUC per predictor).
#   NA for any predictor where the index is undefined or unstable on the resample.
# ------------------------------------------------------------------------------

fbootT = function(dataset, bssample, low.limit, up.limit) {

  # Resample the dataset using the bootstrap indices provided by boot::boot()
  bsdata = dataset[bssample, ]

  # Loop over each predictor column, compute TpAUC
  # Warnings from constant ROC curves (all TPR identical) are captured and
  # converted to NA to avoid propagating an unstable value into the distribution
  vapply(seq(2, ncol(bsdata)), function(i) {
    tmp = cbind(bsdata[, 1], bsdata[, i])
    val = tryCatch({
      roc     = pointsCurve(tmp[, 1], tmp[, 2])
      portion = portion_ROC_FPR(up.limit, low.limit, roc[, 1], roc[, 2])
      TpA(portion[, 1], portion[, 2])
    }, warning = function(w) NA_real_,   # constant curve → unstable → NA
       error   = function(e) NA_real_)
    # Additional guard: reject any non-finite value that slipped through
    if (!is.finite(val)) NA_real_ else val
  }, numeric(1))
}


# ------------------------------------------------------------------------------
# fbootN
# ------------------------------------------------------------------------------
# Description:
#   Bootstrap statistic function for NpAUC. Extracts the TPR sub-interval and
#   computes NpA on each resampled dataset. NpA already returns NA when the NLR
#   condition fails; tryCatch additionally guards against errors from degenerate
#   resamples (e.g. TPR interval not intersecting the resampled ROC curve).
#
# Parameters:
#   dataset   - Original data.frame (outcome in col 1, predictors in cols 2+)
#   bssample  - Integer vector of row indices
#   low.limit - Lower bound of the TPR interval
#   up.limit  - Upper bound of the TPR interval
#
# Returns:
#   Numeric vector of length ncol(dataset) - 1 (one NpAUC per predictor).
#   NA for any predictor where the index is undefined or the resample is degenerate.
# ------------------------------------------------------------------------------

fbootN = function(dataset, bssample, low.limit, up.limit) {

  # Resample the dataset using the bootstrap indices provided by boot::boot()
  bsdata = dataset[bssample, ]

  # Loop over each predictor column, compute NpAUC in the TPR domain
  vapply(seq(2, ncol(bsdata)), function(i) {
    tmp = cbind(bsdata[, 1], bsdata[, i])
    val = tryCatch({
      roc     = pointsCurve(tmp[, 1], tmp[, 2])
      portion = portion_ROC_TPR(up.limit, low.limit, roc[, 1], roc[, 2])
      NpA(portion[, 2], portion[, 1])
    }, warning = function(w) NA_real_,
       error   = function(e) NA_real_)
    # Additional guard: reject any non-finite value that slipped through
    if (!is.finite(val)) NA_real_ else val
  }, numeric(1))
}


# ------------------------------------------------------------------------------
# fbootF
# ------------------------------------------------------------------------------
# Description:
#   Bootstrap statistic function for FpAUC. Extracts the TPR sub-interval and
#   computes FpA on each resampled dataset. FpA is valid for any ROC curve shape
#   and never returns NA from shape conditions; tryCatch guards against errors
#   from degenerate resamples (e.g. single-class outcome, TPR interval not
#   intersecting the resampled ROC curve).
#
# Parameters:
#   dataset   - Original data.frame (outcome in col 1, predictors in cols 2+)
#   bssample  - Integer vector of row indices
#   low.limit - Lower bound of the TPR interval
#   up.limit  - Upper bound of the TPR interval
#
# Returns:
#   Numeric vector of length ncol(dataset) - 1 (one FpAUC per predictor).
#   NA only for truly degenerate resamples; never NA from curve shape conditions.
# ------------------------------------------------------------------------------

fbootF = function(dataset, bssample, low.limit, up.limit) {

  # Resample the dataset using the bootstrap indices provided by boot::boot()
  bsdata = dataset[bssample, ]

  # Loop over each predictor column, compute FpAUC in the TPR domain
  vapply(seq(2, ncol(bsdata)), function(i) {
    tmp = cbind(bsdata[, 1], bsdata[, i])
    val = tryCatch({
      roc     = pointsCurve(tmp[, 1], tmp[, 2])
      portion = portion_ROC_TPR(up.limit, low.limit, roc[, 1], roc[, 2])
      FpA(portion[, 2], portion[, 1])
    }, warning = function(w) NA_real_,
       error   = function(e) NA_real_)
    # Additional guard: reject any non-finite value that slipped through
    if (!is.finite(val)) NA_real_ else val
  }, numeric(1))
}


# ------------------------------------------------------------------------------
# run_bootstrap
# ------------------------------------------------------------------------------
# Description:
#   Centralised bootstrap runner shared by all four exported *boot functions.
#   Rather than duplicating the resampling logic in each function, they all
#   delegate to run_bootstrap, which handles the full lifecycle: computing
#   point estimates on the original sample, running boot::boot() only for
#   variables where the index is well-defined, extracting confidence intervals
#   via boot::boot.ci(), generating optional diagnostic plots, and assembling
#   the formatted output table.
#
#   Variables for which the index is undefined on the original sample (t0 = NA,
#   e.g. MCpAUC on improper curves or NpAUC on partially-proper curves) are
#   detected before any resampling begins. All output columns for those
#   variables are set to NA and boot::boot() is never called for them, avoiding
#   wasted computation and boot.ci() failures.
#
# Parameters:
#   dataset        - Preprocessed data.frame (outcome col 1, predictors cols 2+)
#   boot_fn        - One of fbootM, fbootT, fbootN, fbootF
#   low.value      - Lower bound of the FPR or TPR interval
#   up.value       - Upper bound of the FPR or TPR interval
#   r              - Number of bootstrap replicates
#   level          - Confidence level (e.g., 0.95)
#   type.interval  - CI method: "norm", "basic", "stud", "perc", or "bca"
#   name.variable  - Character vector of predictor names (predictors only)
#   index_label    - Name of the index (e.g., "MCpAUC")
#   output_as_SE   - Logical; if TRUE return SummarizedExperiment
#   plot           - Logical; if TRUE generate bootstrap diagnostic plots
#   plot_type      - One of "both", "ci", "distribution"
#   parallel       - Parallelisation: "no", "multicore" (Linux/macOS), "snow"
#   ncpus          - Number of CPU cores to use when parallel != "no"
#
# Returns:
#   tibble or SummarizedExperiment as produced by format_results_boot().
# ------------------------------------------------------------------------------

run_bootstrap = function(dataset, boot_fn, low.value, up.value,
                         r, level, type.interval,
                         name.variable, index_label,
                         output_as_SE, plot, plot_type,
                         parallel = c("no", "multicore", "snow"),
                         ncpus    = 1L) {

  parallel = match.arg(parallel)
  n_vars   = length(name.variable)

  # Evaluate the index on the original (un-resampled) sample for all predictors.
  # Passing the full row index sequence is equivalent to using the original data.
  t0_all = tryCatch(
    boot_fn(dataset, seq_len(nrow(dataset)), low.value, up.value),
    error = function(e) rep(NA_real_, n_vars)
  )

  valid_idx   = which(!is.na(t0_all))   # variables with a valid estimate
  skipped_idx = which( is.na(t0_all))   # variables where the index is undefined

  # Initialise output matrix (all NA); fill in t0 for valid variables
  mat_ci = matrix(NA_real_, nrow = n_vars, ncol = 5,
                  dimnames = list(NULL, c(index_label, "bias", "sd", "lwr", "upr")))
  mat_ci[valid_idx, index_label] = t0_all[valid_idx]

  if (length(skipped_idx) > 0L)
    message(
      sprintf("[%s bootstrap] Skipping %d variable(s) with undefined index (NA): ",
              index_label, length(skipped_idx)),
      paste(name.variable[skipped_idx], collapse = ", "), "."
    )

  # If no variable has a valid estimate, return the all-NA table immediately
  if (length(valid_idx) == 0L) {
    if (isTRUE(plot))
      message("[", index_label, " bootstrap] No valid variables to plot.")
    return(format_results_boot(
      index_vals    = mat_ci[, index_label],
      bias_vals     = mat_ci[, "bias"],
      sd_vals       = mat_ci[, "sd"],
      lwr_vals      = mat_ci[, "lwr"],
      upr_vals      = mat_ci[, "upr"],
      name.variable = name.variable,
      output_as_SE  = output_as_SE,
      index_label   = index_label
    ))
  }

  # Subset to outcome column + valid predictor columns only, then run bootstrap
  dataset_valid = dataset[, c(1L, valid_idx + 1L), drop = FALSE]

  result_boot = boot::boot(
    data      = dataset_valid,
    statistic = boot_fn,
    R         = r,
    low.limit = low.value,
    up.limit  = up.value,
    parallel  = parallel,
    ncpus     = ncpus
  )

  # Map the user-facing type.interval string to the element name in boot.ci output
  ci_slot = switch(type.interval,
    norm  = "normal",
    basic = "basic",
    stud  = "student",
    perc  = "percent",
    bca   = "bca"
  )

  extract_ci = function(ci_obj) {
    if (is.null(ci_obj)) return(c(NA_real_, NA_real_))
    slot = ci_obj[[ci_slot]]
    if (is.null(slot)) return(c(NA_real_, NA_real_))
    tail = as.numeric(slot)
    n    = length(tail)
    return(c(tail[n - 1L], tail[n]))
  }

  # Extract CI bounds, bias and sd for each valid variable
  for (j in seq_along(valid_idx)) {
    i      = valid_idx[j]
    ci     = tryCatch(
      boot::boot.ci(result_boot, type = type.interval, conf = level, index = j),
      error = function(e) NULL
    )
    t0_j   = result_boot$t0[j]
    t_star = result_boot$t[, j]
    bias_j = mean(t_star, na.rm = TRUE) - t0_j
    sd_j   = stats::sd(t_star, na.rm = TRUE)
    bounds = extract_ci(ci)
    mat_ci[i, ] = c(t0_j, bias_j, sd_j, bounds[1], bounds[2])
  }

  # Generate optional diagnostic plots for valid variables only
  if (isTRUE(plot))
    plot_bootstrap_results(
      result_boot   = result_boot,
      name.variable = name.variable[valid_idx],
      index_label   = index_label,
      lwr_vals      = mat_ci[valid_idx, "lwr"],
      upr_vals      = mat_ci[valid_idx, "upr"],
      plot_type     = plot_type
    )

  return(format_results_boot(
    index_vals    = mat_ci[, index_label],
    bias_vals     = mat_ci[, "bias"],
    sd_vals       = mat_ci[, "sd"],
    lwr_vals      = mat_ci[, "lwr"],
    upr_vals      = mat_ci[, "upr"],
    name.variable = name.variable,
    output_as_SE  = output_as_SE,
    index_label   = index_label
  ))
}
