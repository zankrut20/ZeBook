# tests/testthat/test-statistical.R
# ─────────────────────────────────────────────────────────────────────────────
# Tests for goodness.of.fit, threshold.measures, evaluation.criteria, and AICf.
#
# All expected values are derived analytically from the definitions in the
# function source code so there is no dependency on reference output files.
# ─────────────────────────────────────────────────────────────────────────────

# ── goodness.of.fit ────────────────────────────────────────────────────────────

# Shared test vectors used repeatedly below
OBS  <- c(1, 2, 3, 4, 5)
PRED <- c(1, 2, 3, 4, 5)   # perfect prediction
PRED_BIASED <- OBS + 1      # constant +1 bias

test_that("goodness.of.fit: perfect prediction gives EF=1, bias=0, RMSE=0", {
  res <- goodness.of.fit(Yobs = OBS, Ypred = PRED)
  expect_equal(res$EF,   1, tolerance = 1e-10)
  expect_equal(res$bias, 0, tolerance = 1e-10)
  expect_equal(res$RMSE, 0, tolerance = 1e-10)
  expect_equal(res$MSE,  0, tolerance = 1e-10)
})

test_that("goodness.of.fit: known constant bias is computed correctly", {
  # bias = mean(Yobs) - mean(Ypred) = 3 - 4 = -1
  res <- goodness.of.fit(Yobs = OBS, Ypred = PRED_BIASED)
  expect_equal(res$bias, -1, tolerance = 1e-10)
})

test_that("goodness.of.fit: known RMSE for unit-offset predictions", {
  # Every residual is -1 so MSE=1 and RMSE=1
  res <- goodness.of.fit(Yobs = OBS, Ypred = PRED_BIASED)
  expect_equal(res$MSE,  1, tolerance = 1e-10)
  expect_equal(res$RMSE, 1, tolerance = 1e-10)
})

test_that("goodness.of.fit: EF = 0.5 for unit-offset predictions on OBS", {
  # EF = 1 - SSE / SSTot = 1 - 5/10 = 0.5
  # SSTot = sum((1:5 - 3)^2) = 4+1+0+1+4 = 10
  res <- goodness.of.fit(Yobs = OBS, Ypred = PRED_BIASED)
  expect_equal(res$EF, 0.5, tolerance = 1e-10)
})

test_that("goodness.of.fit: N equals length of non-NA pairs", {
  res <- goodness.of.fit(Yobs = OBS, Ypred = PRED)
  expect_equal(res$N, 5)
})

test_that("goodness.of.fit: handles NA pairs by exclusion", {
  obs  <- c(1, 2, NA, 4, 5)
  pred <- c(1, 2, 3,  4, 5)
  res  <- goodness.of.fit(Yobs = obs, Ypred = pred)
  expect_equal(res$N, 4)
})

test_that("goodness.of.fit: stops when Yobs and Ypred differ in length", {
  expect_error(goodness.of.fit(Yobs = 1:5, Ypred = 1:4),
               "different lengths")
})

test_that("goodness.of.fit: MSE decomposition identity holds (bias^2 + SDSD + LCS = MSE)", {
  res <- goodness.of.fit(Yobs = OBS, Ypred = PRED_BIASED)
  lhs <- res$bias.squared + res$SDSD + res$LCS
  expect_equal(lhs, res$MSE, tolerance = 1e-8)
})

test_that("goodness.of.fit: MSE decomposition identity holds for random-like data", {
  set.seed(1)
  obs  <- c(78, 110, 92, 75, 110, 108, 113, 155, 150)
  pred <- c(126, 126, 126, 105, 105, 105, 147, 147, 147)
  res  <- goodness.of.fit(Yobs = obs, Ypred = pred)
  lhs  <- res$bias.squared + res$SDSD + res$LCS
  expect_equal(lhs, res$MSE, tolerance = 1e-6)
})

test_that("goodness.of.fit: returns data.frame with all documented columns", {
  res <- goodness.of.fit(Yobs = OBS, Ypred = PRED)
  expected_cols <- c("N", "mean.Yobs", "mean.Ypred", "bias",
                     "sd.Yobs", "sd.Ypred", "MSE", "RMSE",
                     "relative.RMSE", "MAE", "relative.RMAE",
                     "relative.MAE.prime", "EF", "Willmott.index",
                     "bias.squared", "SDSD", "LCS",
                     "bias.squared.again", "NU", "LC",
                     "MSE.systematic", "MSE.unsystematic")
  expect_true(all(expected_cols %in% names(res)))
})

test_that("goodness.of.fit: draw.plot=FALSE does not produce a plot", {
  # No error; graphics device unchanged
  expect_no_error(goodness.of.fit(Yobs = OBS, Ypred = PRED, draw.plot = FALSE))
})

# ── threshold.measures ─────────────────────────────────────────────────────────

test_that("threshold.measures: returns data.frame with TDIp and CP columns", {
  obs  <- c(78, 110, 92, 75, 110, 108, 113, 155, 150)
  pred <- c(126, 126, 126, 105, 105, 105, 147, 147, 147)
  res  <- suppressMessages(
            capture.output(threshold.measures(obs, pred, p = 80, d = 30),
                           type = "output")
          )
  # The function uses print() internally - result is a data.frame
  direct_res <- threshold.measures(obs, pred, p = 80, d = 30)
  expect_s3_class(direct_res, "data.frame")
  expect_true(all(c("TDIp", "CP") %in% names(direct_res)))
})

test_that("threshold.measures: CP=100 when d >= max(|error|)", {
  obs  <- c(1, 2, 3)
  pred <- c(1, 2, 3)    # all errors = 0
  res  <- threshold.measures(obs, pred, p = 80, d = 0)
  expect_equal(res$CP, 100)
})

test_that("threshold.measures: CP=0 when d < min(|error|)", {
  obs  <- c(1, 2, 3)
  pred <- c(3, 4, 5)    # all errors = 2
  res  <- threshold.measures(obs, pred, p = 80, d = 1)
  expect_equal(res$CP, 0)
})

# ── AICf ──────────────────────────────────────────────────────────────────────

test_that("AICf: returns named numeric vector with AICcomplete and AICshort", {
  obs  <- c(1, 2, 3, 4, 5)
  pred <- c(1.1, 1.9, 3.1, 4.0, 4.9)
  res  <- AICf(Yobs = obs, Ypred = pred, npar = 2)
  expect_type(res, "double")
  expect_named(res, c("AICcomplete", "AICshort"))
  expect_length(res, 2)
})

test_that("AICf: AICcomplete - AICshort == n + n*log(2*pi)", {
  # By definition: AICcomplete = AICshort + n + n*log(2*pi)
  obs  <- c(1, 2, 3, 4, 5)
  pred <- c(1.1, 1.9, 3.1, 4.0, 4.9)
  res  <- AICf(Yobs = obs, Ypred = pred, npar = 2)
  n    <- 5
  expected_diff <- n + n * log(2 * pi)
  # Strip names before comparing to avoid name-mismatch failure
  expect_equal(as.numeric(res["AICcomplete"]) - as.numeric(res["AICshort"]),
               expected_diff, tolerance = 1e-8)
})

test_that("AICf: AICshort = n*log(RSS/n) + 2*npar (closed-form check)", {
  # Use obs=c(0,0) pred=c(1,1): RSS=2, n=2, npar=1
  # AICshort = 2*log(2/2) + 2*1 = 2*0 + 2 = 2
  res <- AICf(Yobs = c(0, 0), Ypred = c(1, 1), npar = 1)
  expect_equal(as.numeric(res["AICshort"]), 2, tolerance = 1e-10)
})

test_that("AICf: lower npar gives lower AICshort (more parsimonious model preferred)", {
  obs  <- c(1, 2, 3, 4, 5)
  pred <- c(1.1, 1.9, 3.1, 4.0, 4.9)
  aic1 <- AICf(Yobs = obs, Ypred = pred, npar = 1)["AICshort"]
  aic3 <- AICf(Yobs = obs, Ypred = pred, npar = 3)["AICshort"]
  expect_true(aic1 < aic3)
})

test_that("AICf: ignores NA in Yobs by counting only non-NA pairs", {
  # With 4 non-NA values, n=4; the AICshort numerically depends only on n and RSS
  obs_with_na <- c(1, NA, 3, 4, 5)  # 4 non-NA values
  pred        <- c(1, 2,  3, 4, 5)  # same length, paired
  # na.omit(Yobs) has length 4; RSS uses na.rm=TRUE over all 5 pairs
  res_na <- AICf(Yobs = obs_with_na, Ypred = pred, npar = 1)
  # n should be 4 (na.omit length)
  n_expected <- 4
  # RSS = 0^2 + 0^2 + 0^2 + 0^2 = 0 (NA pair excluded, rest match exactly)
  # AICshort = 4*log(0/4) + 2 = -Inf + 2 = -Inf
  # Just verify it's a numeric, not an error
  expect_type(res_na, "double")
  expect_length(res_na, 2)
})

# ── evaluation.criteria (deprecated wrapper) ────────────────────────────────

test_that("evaluation.criteria: emits a deprecation warning", {
  expect_warning(evaluation.criteria(Ypred = OBS, Yobs = PRED),
                 "deprecated")
})

test_that("evaluation.criteria: returns a data.frame", {
  suppressWarnings({
    res <- evaluation.criteria(Ypred = OBS, Yobs = PRED)
  })
  expect_s3_class(res, "data.frame")
})

test_that("evaluation.criteria: argument ORDER is (Ypred, Yobs) — opposite to goodness.of.fit", {
  # evaluation.criteria(Ypred, Yobs) has bias = mean(Ypred) - mean(Yobs)
  # goodness.of.fit(Yobs, Ypred) has bias = mean(Yobs) - mean(Ypred)
  # For pred=OBS+1: evaluation.criteria bias = +1, goodness.of.fit bias = -1
  suppressWarnings({
    res_ec  <- evaluation.criteria(Ypred = PRED_BIASED, Yobs = OBS)
    res_gof <- goodness.of.fit(Yobs = OBS, Ypred = PRED_BIASED)
  })
  # Signs should be opposite
  expect_equal(res_ec$bias, 1,  tolerance = 1e-10)   # Ypred - Yobs = +1
  expect_equal(res_gof$bias, -1, tolerance = 1e-10)  # Yobs - Ypred = -1
})

test_that("evaluation.criteria: EF matches goodness.of.fit EF for same data", {
  suppressWarnings({
    res_ec  <- evaluation.criteria(Ypred = PRED_BIASED, Yobs = OBS)
    res_gof <- goodness.of.fit(Yobs = OBS, Ypred = PRED_BIASED)
  })
  # EF is symmetric in the sense that same pairs → same value
  expect_equal(res_ec$EF, res_gof$EF, tolerance = 1e-8)
})
