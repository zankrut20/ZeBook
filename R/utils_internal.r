# R/utils_internal.r
################################################################################
# ZeBook — Internal helper utilities
#
# All functions in this file are INTERNAL (not exported).
# They consolidate code that was previously duplicated across model files.
#
# Do NOT export any function from this file.
# Do NOT call these functions from user-facing vignettes or demos directly.
################################################################################

# ── Weather filtering ──────────────────────────────────────────────────────────

#' Subset a weather data frame by site and/or year.
#'
#' @param weather  A data frame that already contains columns `year` and
#'   `idsite` (columns must have been renamed by the calling `.weather()`
#'   function before this helper is invoked).
#' @param working.year  Integer or character year to filter on, or `NA` for all.
#' @param working.site  Integer or character site id to filter on, or `NA` for
#'   all.
#' @return The filtered (or unfiltered) `weather` data frame.
#' @keywords internal
.filter_site_year <- function(weather, working.year, working.site) {
  if (!is.na(working.year) & !is.na(working.site)) {
    weather <- weather[(weather$year == working.year) &
                         (weather$idsite == working.site), ]
  } else {
    if (!is.na(working.year)) weather <- weather[weather$year == working.year, ]
    if (!is.na(working.site)) weather <- weather[weather$idsite == working.site, ]
  }
  weather
}

# ── Parameter matrix builder ───────────────────────────────────────────────────

#' Convert a parameter data frame into the standard ZeBook parameter matrix.
#'
#' All `*.define.param()` functions produce a matrix with three rows named
#' `"nominal"`, `"binf"`, and `"bsup"`.  This helper enforces that contract
#' and replaces the repetitive two-liner that appeared in every function.
#'
#' @param df  A data frame whose rows represent nominal/binf/bsup values.
#'   The data frame must already have exactly three rows; row names will be
#'   set to `c("nominal", "binf", "bsup")` by this function.
#' @return A named numeric matrix with row names `c("nominal","binf","bsup")`.
#' @keywords internal
.make_param_matrix <- function(df) {
  row.names(df) <- c("nominal", "binf", "bsup")
  as.matrix(df)
}

# ── Simulation wrapper helper ──────────────────────────────────────────────────

#' Apply a scalar-output model function over rows of a parameter matrix.
#'
#' All `*.simule()` wrappers follow the pattern:
#' ```
#'   Y <- apply(X, 1, fn)
#'   if (all) Y <- cbind(X, <name> = Y)
#'   return(as.matrix(Y))
#' ```
#' This helper encapsulates that pattern.
#'
#' @param X           A matrix or data frame; each row is one parameter vector.
#' @param fn          A function of one argument (a row of `X`) that returns a
#'   single numeric scalar.
#' @param output_name Character string used as the column name when `all=TRUE`.
#' @param all         Logical; if `TRUE` return `X` augmented with an output
#'   column named `output_name`.
#' @return A numeric matrix.
#' @keywords internal
.apply_simule <- function(X, fn, output_name = "Y", all = FALSE) {
  Y <- apply(X, 1, fn)
  if (all) Y <- cbind(X, stats::setNames(data.frame(Y), output_name))
  as.matrix(Y)
}

# ── Site-year string parser ────────────────────────────────────────────────────

#' Parse a "site-year" string of the form `"<site>-<year>"`.
#'
#' Used in `maize.multisy()` and `maize.multisy240()` to extract site and year
#' from the conventional `"site-year"` key used in the maize multi-site
#' simulation functions.
#'
#' @param sy  A character string in the format `"<site>-<year>"`, e.g.
#'   `"18-2006"`.
#' @return A named character vector with elements `"site"` and `"year"`.
#' @keywords internal
.parse_siteyear <- function(sy) {
  parts <- strsplit(sy, "-", fixed = TRUE)[[1]]
  c(site = parts[1], year = parts[2])
}
