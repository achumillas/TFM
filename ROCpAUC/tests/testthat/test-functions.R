# ==============================================================================
# test-functions.R
# ==============================================================================
# Unit tests for the ROCpAIv2 package using testthat.
# Run with: devtools::test() or testthat::test_package("ROCpAIv2")
# ==============================================================================

library(testthat)
library(ROCpAIv2)
library(SummarizedExperiment)
library(fission)

# ── Shared test data ──────────────────────────────────────────────────────────
data(fission)
genes_test = c("SPNCRNA.1080", "SPAC186.08c", "SPNCRNA.1420")
strain     = SummarizedExperiment::colData(fission)$strain
expr       = t(SummarizedExperiment::assay(fission)[genes_test, ])
Sp_test    = as.data.frame(cbind(strain = strain, expr))


# ==============================================================================
# pointsCurve
# ==============================================================================

test_that("pointsCurve returns a data.frame with FPR and TPR columns", {
  res = pointsCurve(Sp_test[, 1], as.numeric(Sp_test[, 2]))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("FPR", "TPR"))
})

test_that("pointsCurve output starts at (0,0) and ends at (1,1)", {
  res = pointsCurve(Sp_test[, 1], as.numeric(Sp_test[, 2]))
  expect_equal(res$FPR[1], 0)
  expect_equal(res$TPR[1], 0)
  expect_equal(res$FPR[nrow(res)], 1)
  expect_equal(res$TPR[nrow(res)], 1)
})

test_that("pointsCurve FPR is non-decreasing", {
  res = pointsCurve(Sp_test[, 1], as.numeric(Sp_test[, 2]))
  expect_true(all(diff(res$FPR) >= 0))
})

test_that("pointsCurve accepts 0/1 and 1/2 coding", {
  x01 = as.integer(Sp_test[, 1])
  x12 = x01 + 1L
  r01 = pointsCurve(x01, as.numeric(Sp_test[, 2]))
  r12 = pointsCurve(x12, as.numeric(Sp_test[, 2]))
  expect_equal(r01, r12)
})

test_that("pointsCurve errors on invalid class labels", {
  expect_error(pointsCurve(rep(3, nrow(Sp_test)), as.numeric(Sp_test[, 2])))
})


# ==============================================================================
# MCpAUC
# ==============================================================================

test_that("MCpAUC returns a tibble with correct columns (data.frame input)", {
  res = MCpAUC(Sp_test, low.value = 0, up.value = 0.25)
  expect_s3_class(res, "tbl_df")
  expect_true(all(c("Variable", "MCpAUC", "pAUC") %in% names(res)))
  expect_equal(nrow(res), ncol(Sp_test) - 1L)
})

test_that("MCpAUC returns a tibble with correct columns (SE input)", {
  res = MCpAUC(fission, low.value = 0, up.value = 0.25,
               selection = genes_test, variable = "strain")
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), length(genes_test))
})

test_that("MCpAUC values are in [0.5, 1] or NA", {
  res = MCpAUC(Sp_test, low.value = 0, up.value = 0.25)
  vals = res$MCpAUC[!is.na(res$MCpAUC)]
  expect_true(all(vals >= 0.5 & vals <= 1))
})

test_that("MCpAUC output_as_SE returns a SummarizedExperiment", {
  res = MCpAUC(fission, low.value = 0, up.value = 0.25,
               selection = genes_test, variable = "strain",
               output_as_SE = TRUE)
  expect_s4_class(res, "SummarizedExperiment")
  expect_true("MCpAUC" %in% rownames(assay(res)))
})


# ==============================================================================
# TpAUC
# ==============================================================================

test_that("TpAUC values are in [0.5, 1]", {
  res = TpAUC(Sp_test, low.value = 0, up.value = 0.25)
  vals = res$TpAUC[!is.na(res$TpAUC)]
  expect_true(all(vals >= 0.5 & vals <= 1))
})

test_that("TpAUC returns no NA for any curve shape", {
  res = TpAUC(Sp_test, low.value = 0, up.value = 0.25)
  # TpAUC is valid for any curve shape; expect no NA
  expect_true(all(!is.na(res$TpAUC)))
})


# ==============================================================================
# NpAUC
# ==============================================================================

test_that("NpAUC returns a tibble with NpAUC and pAUC columns", {
  res = NpAUC(Sp_test, low.value = 0.9, up.value = 1)
  expect_true(all(c("Variable", "NpAUC", "pAUC") %in% names(res)))
})

test_that("NpAUC values are in [0, 1] or NA", {
  res = NpAUC(Sp_test, low.value = 0.9, up.value = 1)
  vals = res$NpAUC[!is.na(res$NpAUC)]
  expect_true(all(vals >= 0 & vals <= 1))
})


# ==============================================================================
# FpAUC
# ==============================================================================

test_that("FpAUC values are in [0.5, 1] (no NA ever)", {
  res = FpAUC(Sp_test, low.value = 0.9, up.value = 1)
  expect_true(all(!is.na(res$FpAUC)))
  expect_true(all(res$FpAUC >= 0.5 & res$FpAUC <= 1))
})


# ==============================================================================
# pAUC
# ==============================================================================

test_that("pAUC FPR domain returns values in [0, 0.25] range", {
  res = pAUC(Sp_test, low.value = 0, up.value = 0.25, axis = "x")
  vals = res$pAUC[!is.na(res$pAUC)]
  expect_true(all(vals >= 0 & vals <= 0.25))
})

test_that("pAUC accepts axis = 'y'", {
  res = pAUC(Sp_test, low.value = 0.9, up.value = 1, axis = "y")
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("Variable", "pAUC"))
})


# ==============================================================================
# allIndices
# ==============================================================================

test_that("allIndices returns all six index columns", {
  res = allIndices(Sp_test,
                   fpr.low = 0, fpr.up = 0.25,
                   tpr.low = 0.9, tpr.up = 1)
  expected_cols = c("Variable", "pAUC_FPR", "MCpAUC", "TpAUC",
                     "pAUC_TPR", "NpAUC", "FpAUC")
  expect_true(all(expected_cols %in% names(res)))
  expect_equal(nrow(res), ncol(Sp_test) - 1L)
})


# ==============================================================================
# Bootstrap functions (light: r = 10 for speed)
# ==============================================================================

test_that("MCpAUCboot returns correct columns", {
  res = MCpAUCboot(Sp_test, low.value = 0, up.value = 0.25, r = 10)
  expect_true(all(c("Variable", "MCpAUC", "bias", "sd", "lwr", "upr") %in% names(res)))
})

test_that("TpAUCboot returns correct columns", {
  res = TpAUCboot(Sp_test, low.value = 0, up.value = 0.25, r = 10)
  expect_true(all(c("Variable", "TpAUC", "bias", "sd", "lwr", "upr") %in% names(res)))
})

test_that("NpAUCboot returns correct columns", {
  res = NpAUCboot(Sp_test, low.value = 0.9, up.value = 1, r = 10)
  expect_true(all(c("Variable", "NpAUC", "bias", "sd", "lwr", "upr") %in% names(res)))
})

test_that("FpAUCboot returns correct columns with no NA estimates", {
  res = FpAUCboot(Sp_test, low.value = 0.9, up.value = 1, r = 10)
  expect_true(all(c("Variable", "FpAUC", "bias", "sd", "lwr", "upr") %in% names(res)))
  # FpAUC is always defined; no NA in the point estimate column
  expect_true(all(!is.na(res$FpAUC)))
})

test_that("MCpAUCboot output_as_SE works", {
  res = MCpAUCboot(Sp_test, low.value = 0, up.value = 0.25,
                   r = 10, output_as_SE = TRUE)
  expect_s4_class(res, "SummarizedExperiment")
})

test_that("bootstrap CI bounds satisfy lwr <= estimate <= upr when defined", {
  res = FpAUCboot(Sp_test, low.value = 0.9, up.value = 1, r = 50)
  ok  = !is.na(res$lwr) & !is.na(res$upr)
  if (any(ok)) {
    expect_true(all(res$lwr[ok] <= res$FpAUC[ok]))
    expect_true(all(res$FpAUC[ok] <= res$upr[ok]))
  }
})
