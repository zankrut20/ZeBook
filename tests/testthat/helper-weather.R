# tests/testthat/helper-weather.R
# ─────────────────────────────────────────────────────────────────────────────
# Synthetic, deterministic weather data factories used across the test suite.
# All values are constant so tests are fully reproducible and need no
# package datasets to be loaded.
# ─────────────────────────────────────────────────────────────────────────────

#' Minimal weather frame for carrot.weevil.model (needs TMIN, TMAX, day)
make_weather_carrot <- function(n = 365) {
  data.frame(
    day  = seq_len(n),
    TMIN = rep(10, n),
    TMAX = rep(20, n)
  )
}

#' Minimal weather frame for maize.model (needs Tmin, Tmax, I, day)
make_weather_maize <- function(n = 300, sdate = 100) {
  data.frame(
    day  = seq_len(n),
    Tmin = rep(12, n),
    Tmax = rep(24, n),
    I    = rep(15, n)   # MJ/m2/day
  )
}

#' Minimal weather frame for watbal.model (needs RAIN, ETr, day)
make_weather_watbal <- function(n = 150) {
  data.frame(
    day  = seq_len(n),
    RAIN = rep(2.0, n),
    ETr  = rep(3.5, n)
  )
}

# Nominal parameter vectors for commonly-tested models (row "nominal" from
# define.param functions, extracted once here to avoid repetition).

MAIZE_PARAM <- c(Tbase = 7, RUE = 1.85, K = 0.7,
                 alpha = 0.00243, LAImax = 7, TTM = 1200, TTL = 700)

WATBAL_PARAM <- c(WHC = 0.15, MUF = 0.1, DC = 0.55, z = 400, CN = 65)
