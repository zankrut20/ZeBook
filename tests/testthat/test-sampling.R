# tests/testthat/test-sampling.R
# ─────────────────────────────────────────────────────────────────────────────
# Tests for the parameter-sampling utilities:
#   param.runif, param.rtriangle, q.arg.fast.runif
#
# Because these functions use random number generation, tests use set.seed()
# for reproducibility and focus on structural / distributional properties
# rather than exact values (which would be fragile across RNG versions).
# ─────────────────────────────────────────────────────────────────────────────

# Shared reference parameter matrix used throughout
FACTORS <- weed.define.param()   # 3-row (nominal/binf/bsup), 16-col matrix

# ── param.runif ────────────────────────────────────────────────────────────────

test_that("param.runif: returns a data.frame with N rows and p columns", {
  set.seed(1)
  res <- param.runif(FACTORS, N = 50)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 50)
  expect_equal(ncol(res), ncol(FACTORS))
})

test_that("param.runif: column names match model.factors column names", {
  set.seed(1)
  res <- param.runif(FACTORS, N = 10)
  expect_equal(names(res), colnames(FACTORS))
})

test_that("param.runif: all values lie within [binf, bsup] for each parameter", {
  set.seed(42)
  res <- param.runif(FACTORS, N = 500)
  for (p in colnames(FACTORS)) {
    lo <- FACTORS["binf", p]
    hi <- FACTORS["bsup", p]
    expect_true(all(res[[p]] >= lo), info = paste("param", p, "below binf"))
    expect_true(all(res[[p]] <= hi), info = paste("param", p, "above bsup"))
  }
})

test_that("param.runif: produces different draws on each call (random)", {
  set.seed(1)
  r1 <- param.runif(FACTORS, N = 10)
  r2 <- param.runif(FACTORS, N = 10)
  # Extremely unlikely both draws are identical
  expect_false(identical(r1, r2))
})

test_that("param.runif: same seed gives reproducible results", {
  set.seed(99)
  r1 <- param.runif(FACTORS, N = 20)
  set.seed(99)
  r2 <- param.runif(FACTORS, N = 20)
  expect_identical(r1, r2)
})

test_that("param.runif: works with a single-parameter matrix", {
  single <- FACTORS[, "mu", drop = FALSE]
  set.seed(1)
  res <- param.runif(single, N = 30)
  # Result is a data.frame (or named vector) with 30 values in [binf, bsup]
  vals <- if (is.data.frame(res)) res[[1]] else as.numeric(res)
  expect_length(vals, 30)
  expect_true(all(vals >= FACTORS["binf", "mu"] &
                  vals <= FACTORS["bsup", "mu"]))
})

# ── param.rtriangle ────────────────────────────────────────────────────────────

test_that("param.rtriangle: returns a data.frame with N rows and p columns", {
  set.seed(1)
  res <- param.rtriangle(FACTORS, N = 50)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 50)
  expect_equal(ncol(res), ncol(FACTORS))
})

test_that("param.rtriangle: column names match model.factors column names", {
  set.seed(1)
  res <- param.rtriangle(FACTORS, N = 10)
  expect_equal(names(res), colnames(FACTORS))
})

test_that("param.rtriangle: all values lie within [binf, bsup]", {
  set.seed(42)
  res <- param.rtriangle(FACTORS, N = 500)
  for (p in colnames(FACTORS)) {
    lo <- FACTORS["binf", p]
    hi <- FACTORS["bsup", p]
    expect_true(all(res[[p]] >= lo), info = paste("param", p, "below binf"))
    expect_true(all(res[[p]] <= hi), info = paste("param", p, "above bsup"))
  }
})

test_that("param.rtriangle: same seed gives reproducible results", {
  set.seed(7)
  r1 <- param.rtriangle(FACTORS, N = 20)
  set.seed(7)
  r2 <- param.rtriangle(FACTORS, N = 20)
  expect_identical(r1, r2)
})

test_that("param.rtriangle: mean is closer to nominal than endpoints (triangular mode)", {
  set.seed(123)
  res <- param.rtriangle(FACTORS, N = 5000)
  # For a symmetric triangle the mean equals the mode (nominal)
  # For asymmetric triangles mean ≈ (binf + nominal + bsup) / 3
  for (p in colnames(FACTORS)) {
    lo  <- FACTORS["binf",    p]
    nom <- FACTORS["nominal", p]
    hi  <- FACTORS["bsup",   p]
    theoretical_mean <- (lo + nom + hi) / 3
    sample_mean      <- mean(res[[p]])
    expect_equal(sample_mean, theoretical_mean,
                 tolerance = 0.02 * abs(theoretical_mean) + 0.001,
                 info = paste("param", p))
  }
})

# ── q.arg.fast.runif ──────────────────────────────────────────────────────────

test_that("q.arg.fast.runif: returns a list of length = ncol(model.factors)", {
  res <- q.arg.fast.runif(FACTORS)
  expect_type(res, "list")
  expect_equal(length(res), ncol(FACTORS))
})

test_that("q.arg.fast.runif: each element contains min and max entries", {
  res <- q.arg.fast.runif(FACTORS)
  for (el in res) {
    expect_true(all(c("min", "max") %in% names(el)))
  }
})

test_that("q.arg.fast.runif: min/max values match binf/bsup rows", {
  res <- q.arg.fast.runif(FACTORS)
  for (i in seq_along(res)) {
    p <- colnames(FACTORS)[i]
    expect_equal(res[[i]]$min, FACTORS["binf", p])
    expect_equal(res[[i]]$max, FACTORS["bsup", p])
  }
})
