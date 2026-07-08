# tests/testthat/test-utils-internal.R
# Phase 3 — Tests for internal utility helpers (utils_internal.r)
#
# These functions are unexported, so we access them via the package namespace
# using the ::: operator.
################################################################################

# ── .filter_site_year ─────────────────────────────────────────────────────────

# Build a minimal weather data frame mirroring the column structure after
# .weather() functions have renamed the raw columns.
make_fw <- function() {
  data.frame(
    idsite = c(1L, 1L, 2L, 2L),
    year   = c(2000L, 2001L, 2000L, 2001L),
    day    = c(1L, 1L, 1L, 1L),
    RAIN   = c(0.1, 0.2, 0.3, 0.4),
    stringsAsFactors = FALSE
  )
}

test_that(".filter_site_year: returns all rows when both args are NA", {
  fw  <- make_fw()
  out <- ZeBook:::.filter_site_year(fw, NA, NA)
  expect_equal(nrow(out), nrow(fw))
})

test_that(".filter_site_year: filters by year only", {
  fw  <- make_fw()
  out <- ZeBook:::.filter_site_year(fw, working.year = 2000, working.site = NA)
  expect_equal(nrow(out), 2L)
  expect_true(all(out$year == 2000))
})

test_that(".filter_site_year: filters by site only", {
  fw  <- make_fw()
  out <- ZeBook:::.filter_site_year(fw, working.year = NA, working.site = 1)
  expect_equal(nrow(out), 2L)
  expect_true(all(out$idsite == 1))
})

test_that(".filter_site_year: filters by both year and site", {
  fw  <- make_fw()
  out <- ZeBook:::.filter_site_year(fw, working.year = 2001, working.site = 2)
  expect_equal(nrow(out), 1L)
  expect_equal(out$RAIN, 0.4)
})

test_that(".filter_site_year: returns zero rows when no match", {
  fw  <- make_fw()
  out <- ZeBook:::.filter_site_year(fw, working.year = 9999, working.site = 9)
  expect_equal(nrow(out), 0L)
})

test_that(".filter_site_year: does not modify column structure", {
  fw  <- make_fw()
  out <- ZeBook:::.filter_site_year(fw, working.year = 2000, working.site = NA)
  expect_equal(names(out), names(fw))
})

# ── .make_param_matrix ────────────────────────────────────────────────────────

test_that(".make_param_matrix: returns a matrix", {
  df  <- data.frame(a = c(1, 0.5, 1.5), b = c(2, 1, 3))
  out <- ZeBook:::.make_param_matrix(df)
  expect_true(is.matrix(out))
})

test_that(".make_param_matrix: row names are nominal/binf/bsup", {
  df  <- data.frame(a = c(1, 0.5, 1.5), b = c(2, 1, 3))
  out <- ZeBook:::.make_param_matrix(df)
  expect_equal(rownames(out), c("nominal", "binf", "bsup"))
})

test_that(".make_param_matrix: column names are preserved from data.frame", {
  df  <- data.frame(alpha = c(1, 0.5, 1.5), beta = c(2, 1, 3))
  out <- ZeBook:::.make_param_matrix(df)
  expect_equal(colnames(out), c("alpha", "beta"))
})

test_that(".make_param_matrix: values are preserved exactly", {
  df  <- data.frame(x = c(10, 5, 20))
  out <- ZeBook:::.make_param_matrix(df)
  expect_equal(as.numeric(out[, "x"]), c(10, 5, 20))
})

test_that(".make_param_matrix: roundtrips through real define.param output", {
  # Test via maize.define.param() which now uses .make_param_matrix()
  mp <- maize.define.param()
  expect_equal(rownames(mp), c("nominal", "binf", "bsup"))
  expect_true(is.matrix(mp))
})

# ── .apply_simule ─────────────────────────────────────────────────────────────

# Build a tiny 3-column parameter matrix for testing
make_X <- function() {
  matrix(c(1, 2, 3,
            4, 5, 6),
         nrow = 2, byrow = TRUE,
         dimnames = list(NULL, c("a", "b", "c")))
}

test_that(".apply_simule: returns a matrix with one row per input row", {
  X   <- make_X()
  out <- ZeBook:::.apply_simule(X, fn = function(v) sum(v), output_name = "S", all = FALSE)
  expect_true(is.matrix(out))
  expect_equal(nrow(out), nrow(X))
})

test_that(".apply_simule: scalar function applied correctly (sum of row)", {
  X   <- make_X()
  out <- ZeBook:::.apply_simule(X, fn = function(v) sum(v), output_name = "S", all = FALSE)
  expect_equal(as.numeric(out), c(6, 15))
})

test_that(".apply_simule: all=TRUE appends output column to X", {
  X   <- make_X()
  out <- ZeBook:::.apply_simule(X, fn = function(v) sum(v), output_name = "S", all = TRUE)
  expect_true("S" %in% colnames(out))
  expect_equal(ncol(out), ncol(X) + 1L)
})

test_that(".apply_simule: all=FALSE does not append columns", {
  X   <- make_X()
  out <- ZeBook:::.apply_simule(X, fn = function(v) sum(v), output_name = "S", all = FALSE)
  expect_false("S" %in% colnames(out))
})

test_that(".apply_simule: all=TRUE preserves original column values", {
  X   <- make_X()
  out <- ZeBook:::.apply_simule(X, fn = function(v) sum(v), output_name = "S", all = TRUE)
  # first row of X: a=1, b=2, c=3 → sum = 6
  expect_equal(as.numeric(out[1, "a"]), 1)
  expect_equal(as.numeric(out[1, "S"]), 6)
})

# ── .parse_siteyear ───────────────────────────────────────────────────────────

test_that(".parse_siteyear: returns a named character vector of length 2", {
  out <- ZeBook:::.parse_siteyear("18-2006")
  expect_length(out, 2L)
  expect_named(out, c("site", "year"))
})

test_that(".parse_siteyear: extracts site correctly", {
  expect_equal(ZeBook:::.parse_siteyear("18-2006")[["site"]], "18")
})

test_that(".parse_siteyear: extracts year correctly", {
  expect_equal(ZeBook:::.parse_siteyear("18-2006")[["year"]], "2006")
})

test_that(".parse_siteyear: works with multi-digit site codes", {
  expect_equal(ZeBook:::.parse_siteyear("123-1999")[["site"]], "123")
  expect_equal(ZeBook:::.parse_siteyear("123-1999")[["year"]], "1999")
})

test_that(".parse_siteyear: consistent with legacy strsplit approach", {
  sy   <- "64-2004"
  legacy_site <- strsplit(sy, "-")[[1]][1]
  legacy_year <- strsplit(sy, "-")[[1]][2]
  out  <- ZeBook:::.parse_siteyear(sy)
  expect_equal(out[["site"]], legacy_site)
  expect_equal(out[["year"]], legacy_year)
})
