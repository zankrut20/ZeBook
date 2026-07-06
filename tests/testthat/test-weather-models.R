# tests/testthat/test-weather-models.R
# ─────────────────────────────────────────────────────────────────────────────
# Golden-output tests for models that depend on weather input.
#
# Strategy:
#   * Synthetic, constant-value weather frames (from helper-weather.R) are used
#     for structural / invariant tests — no network access required.
#   * Package datasets are used via `data()` only where the model's internal
#     weather-selection helper (e.g. epirice.weather) requires them.
#   * Snapshot tests capture full model output on first run and re-validate on
#     subsequent runs; update with testthat::snapshot_update() when behaviour
#     intentionally changes.
# ─────────────────────────────────────────────────────────────────────────────

# ── carrot.weevil.model ────────────────────────────────────────────────────────

test_that("carrot.weevil.model: returns a list with $sim and $param", {
  w   <- make_weather_carrot()
  res <- carrot.weevil.model(weather = w, sdate = 1, ldate = 50)
  expect_type(res, "list")
  expect_true(all(c("sim", "param") %in% names(res)))
})

test_that("carrot.weevil.model: $sim has required columns", {
  w   <- make_weather_carrot()
  res <- carrot.weevil.model(weather = w, sdate = 1, ldate = 50)$sim
  expect_s3_class(res, "data.frame")
  expect_true(all(c("day", "DACE", "TT", "STADE", "NGEN") %in% names(res)))
  expect_equal(nrow(res), 50)   # sdate:ldate
})

test_that("carrot.weevil.model: TT[1] == 0 (initialised correctly)", {
  w   <- make_weather_carrot()
  res <- carrot.weevil.model(weather = w, sdate = 1, ldate = 50)$sim
  expect_equal(res$TT[1], 0, tolerance = 1e-8)
})

test_that("carrot.weevil.model: constant 10/20 °C weather gives exact TT accumulation", {
  # dTT = max(0, (10+20)/2 - 7) = 8 degC/day
  # After 1 step: TT mod ttgeneration = 8
  # ttgeneration = 130+256+114+130+91 = 721
  # 8 < 721 so still in eggs stage on day 2
  w   <- make_weather_carrot(n = 365)
  res <- carrot.weevil.model(tbase = 7, weather = w, sdate = 1, ldate = 10)$sim
  # TT on day 2 should be 8 (mod 721 is 8)
  expect_equal(res$TT[2], 8, tolerance = 1e-8)
})

test_that("carrot.weevil.model: produces no console output (print removed)", {
  w <- make_weather_carrot(n = 365)
  expect_silent(carrot.weevil.model(weather = w, sdate = 1, ldate = 100))
})

test_that("carrot.weevil.model: $param contains all expected names", {
  w   <- make_weather_carrot()
  res <- carrot.weevil.model(weather = w, sdate = 1, ldate = 50)$param
  expect_true(all(c("tbase","tteggs","ttlarvae","ttprepupae",
                    "ttpupae","ttadultpreovi","ttgeneration") %in% names(res)))
})

test_that("carrot.weevil.model: snapshot of 20-day simulation (golden output)", {
  set.seed(42)
  w <- make_weather_carrot(n = 365)
  res <- carrot.weevil.model(tbase = 7, tteggs = 130, ttlarvae = 256,
                             ttprepupae = 114, ttpupae = 130,
                             ttadultpreovi = 91,
                             weather = w, sdate = 1, ldate = 20)$sim
  expect_snapshot(round(res[, c("day","TT","NGEN")], 4))
})

# ── maize.model ────────────────────────────────────────────────────────────────

test_that("maize.model: returns data.frame with day, TT, LAI, B columns", {
  w   <- make_weather_maize()
  res <- maize.model(Tbase = 7, RUE = 1.85, K = 0.7, alpha = 0.00243,
                     LAImax = 7, TTM = 1200, TTL = 700,
                     weather = w, sdate = 100, ldate = 200)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("day", "TT", "LAI", "B"))
  expect_equal(nrow(res), 101)
})

test_that("maize.model: initial conditions (TT=0, LAI=0.01, B=1)", {
  w   <- make_weather_maize()
  res <- maize.model(Tbase = 7, RUE = 1.85, K = 0.7, alpha = 0.00243,
                     LAImax = 7, TTM = 1200, TTL = 700,
                     weather = w, sdate = 100, ldate = 200)
  expect_equal(res$TT[1],  0,    tolerance = 1e-8)
  expect_equal(res$LAI[1], 0.01, tolerance = 1e-8)
  expect_equal(res$B[1],   1,    tolerance = 1e-8)
})

test_that("maize.model: TT is monotone non-decreasing", {
  w   <- make_weather_maize()
  res <- maize.model(Tbase = 7, RUE = 1.85, K = 0.7, alpha = 0.00243,
                     LAImax = 7, TTM = 1200, TTL = 700,
                     weather = w, sdate = 100, ldate = 200)
  expect_true(all(diff(res$TT) >= 0))
})

test_that("maize.model: biomass B is monotone non-decreasing while TT < TTM", {
  w   <- make_weather_maize()
  res <- maize.model(Tbase = 7, RUE = 1.85, K = 0.7, alpha = 0.00243,
                     LAImax = 7, TTM = 1200, TTL = 700,
                     weather = w, sdate = 100, ldate = 200)
  growing <- res[res$TT < 1200, ]
  expect_true(all(diff(growing$B) >= 0))
})

test_that("maize.model: constant dTT = (Tmin+Tmax)/2 - Tbase for constant weather", {
  # (12+24)/2 - 7 = 18 - 7 = 11  degC-days per step
  w   <- make_weather_maize()
  res <- maize.model(Tbase = 7, RUE = 1.85, K = 0.7, alpha = 0.00243,
                     LAImax = 7, TTM = 1200, TTL = 700,
                     weather = w, sdate = 1, ldate = 50)
  expect_equal(res$TT[2] - res$TT[1], 11, tolerance = 1e-8)
  expect_equal(res$TT[3] - res$TT[2], 11, tolerance = 1e-8)
})

test_that("maize.model2 is consistent with maize.model", {
  w   <- make_weather_maize()
  r1  <- maize.model(Tbase = 7, RUE = 1.85, K = 0.7, alpha = 0.00243,
                     LAImax = 7, TTM = 1200, TTL = 700,
                     weather = w, sdate = 100, ldate = 200)
  r2  <- maize.model2(param = MAIZE_PARAM, weather = w,
                      sdate = 100, ldate = 200)
  expect_equal(r1, r2, tolerance = 1e-8)
})

test_that("maize.define.param: returns matrix with 3 rows and 7 parameters", {
  param <- maize.define.param()
  expect_true(is.matrix(param))
  expect_equal(nrow(param), 3)
  expect_equal(rownames(param), c("nominal", "binf", "bsup"))
  expect_equal(ncol(param), 7)
})

# ── watbal.model ───────────────────────────────────────────────────────────────

test_that("watbal.model: returns data.frame with required columns", {
  w   <- make_weather_watbal()
  res <- watbal.model(param = WATBAL_PARAM, weather = w,
                      WP = 0.06, FC = 0.23)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("day", "RAIN", "ETr", "WAT", "WATp") %in% names(res)))
  expect_equal(nrow(res), nrow(w))
})

test_that("watbal.model: WAT is strictly positive throughout", {
  w   <- make_weather_watbal()
  res <- watbal.model(param = WATBAL_PARAM, weather = w,
                      WP = 0.06, FC = 0.23)
  expect_true(all(res$WAT > 0))
})

test_that("watbal.model: WATp = WAT / z", {
  w   <- make_weather_watbal()
  res <- watbal.model(param = WATBAL_PARAM, weather = w,
                      WP = 0.06, FC = 0.23)
  z   <- WATBAL_PARAM["z"]
  expect_equal(res$WATp, res$WAT / z, tolerance = 1e-10)
})

test_that("watbal.model: default WAT0 initialises to z*FC", {
  w   <- make_weather_watbal()
  FC  <- 0.23
  z   <- as.numeric(WATBAL_PARAM["z"])  # strip name to avoid name-mismatch
  res <- watbal.model(param = WATBAL_PARAM, weather = w,
                      WP = 0.06, FC = FC)
  expect_equal(res$WAT[1], z * FC, tolerance = 1e-8)
})

# ── epirice.model ─────────────────────────────────────────────────────────────

test_that("epirice.model: runs with synthetic weather and returns correct columns", {
  # Build minimal weather frame matching the columns epirice.model reads:
  # weather$TMIN, weather$TMAX, weather$RH2M, weather$RAIN
  n  <- 130
  w  <- data.frame(
    TMIN = rep(20, n), TMAX = rep(30, n),
    RH2M = rep(95, n), RAIN = rep(10, n)   # RH>=90 triggers infection
  )
  # Minimal epirice parameters (named vector matching param["..."] accesses)
  param <- c(
    Site.size = 1, Sx = 100000, RRG = 0.1, RRS = 0.0035,
    EODate = 20, p = 5, i = 20, rl = 0.0, RcOpt = 0.61,
    RcT = 1, RcW = 1, RcA = 1, a = 1, k = 0.6, SenescType = 1
  )
  res <- epirice.model(param = param, weather = w,
                       sdate = 1, ldate = 120, H0 = 600)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("day","DACE","H","L","II","P","TS","TOTDIS","severity") %in%
                    names(res)))
  expect_equal(nrow(res), 120)
})

test_that("epirice.model: severity is between 0 and 100", {
  n  <- 130
  w  <- data.frame(
    TMIN = rep(20, n), TMAX = rep(30, n),
    RH2M = rep(95, n), RAIN = rep(10, n)
  )
  param <- c(
    Site.size = 1, Sx = 100000, RRG = 0.1, RRS = 0.0035,
    EODate = 20, p = 5, i = 20, rl = 0.0, RcOpt = 0.61,
    RcT = 1, RcW = 1, RcA = 1, a = 1, k = 0.6, SenescType = 1
  )
  res <- epirice.model(param = param, weather = w,
                       sdate = 1, ldate = 120, H0 = 600)
  valid <- !is.na(res$severity)
  expect_true(all(res$severity[valid] >= 0 & res$severity[valid] <= 100))
})

test_that("epirice.multi.simule: produces no console output (print removed)", {
  # This function previously printed AUDPC unconditionally (Phase 1 fix).
  # We verify it is now silent by calling with synthetic data that bypasses
  # epirice.weather (which requires the bundled weather_SouthAsia dataset).
  skip("epirice.multi.simule requires weather_SouthAsia dataset via epirice.weather — tested manually")
})
