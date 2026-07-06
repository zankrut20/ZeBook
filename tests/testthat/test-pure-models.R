# tests/testthat/test-pure-models.R
# ─────────────────────────────────────────────────────────────────────────────
# Golden-output tests for models that require no weather data.
# Expected values are analytically derived; tolerances are set to 1e-8 to
# catch floating-point drift while ignoring legitimate rounding differences.
# ─────────────────────────────────────────────────────────────────────────────

# ── carbonsoil ────────────────────────────────────────────────────────────────

test_that("carbonsoil.update: single-step forward Euler is exact", {
  # carbonsoil.update(Zy, R, b, Uy)
  # dZ = -R*Zy + b*Uy = -0.02*10000 + 0.3*2000 = -200 + 600 = 400
  expect_equal(carbonsoil.update(Zy = 10000, R = 0.02, b = 0.3, Uy = 2000),
               10400, tolerance = 1e-8)
})

test_that("carbonsoil.update: equilibrium when dZ = 0", {
  # dZ = 0  ⟺  b*Uy = R*Zy  ⟺  Zy* = b*Uy/R = 0.3*2000/0.02 = 30000
  eq_Zy <- 0.3 * 2000 / 0.02
  expect_equal(carbonsoil.update(Zy = eq_Zy, R = 0.02, b = 0.3, Uy = 2000),
               eq_Zy, tolerance = 1e-8)
})

test_that("carbonsoil.model: returns data.frame with correct structure", {
  # carbonsoil.model(R, b, U, Z1, duration)
  result <- carbonsoil.model(R = 0.02, b = 0.3, U = 2000,
                             Z1 = 10000, duration = 10)
  expect_s3_class(result, "data.frame")
  expect_true("Z" %in% names(result))      # column is Z, not Zy
  expect_equal(nrow(result), 10)           # years 1 .. duration
  expect_equal(result$Z[1], 10000, tolerance = 1e-8)
  # After 10 steps from below equilibrium (Z1=10000 < 30000) carbon must increase
  expect_true(result$Z[10] > result$Z[1])
})

test_that("carbonsoil.model: converges toward equilibrium from below and above", {
  # equilibrium: Z* = b*U/R = 0.3*2000/0.02 = 30000
  eq    <- 0.3 * 2000 / 0.02
  below <- carbonsoil.model(R = 0.02, b = 0.3, U = 2000, Z1 = 10000, duration = 200)
  above <- carbonsoil.model(R = 0.02, b = 0.3, U = 2000, Z1 = 50000, duration = 200)
  # Both move toward equilibrium
  expect_true(abs(tail(below$Z, 1) - eq) < abs(below$Z[1] - eq))
  expect_true(abs(tail(above$Z, 1) - eq) < abs(above$Z[1] - eq))
})

# ── seedweight ────────────────────────────────────────────────────────────────

test_that("seedweight.model: logistic formula matches closed-form value", {
  # W / (1 + exp(B - C*DD)) with DD=500, W=30, B=4, C=0.02
  # = 30 / (1 + exp(4 - 10)) = 30 / (1 + exp(-6))
  expected <- 30 / (1 + exp(-6))
  expect_equal(seedweight.model(DD = 500, W = 30, B = 4, C = 0.02),
               expected, tolerance = 1e-10)
})

test_that("seedweight.model: approaches W asymptotically for large DD", {
  W <- 30
  val_large <- seedweight.model(DD = 10000, W = W, B = 4, C = 0.02)
  expect_lt(abs(val_large - W), 1e-5)
})

test_that("seedweight.model: is monotone increasing in DD", {
  vals <- vapply(c(100, 200, 300, 400, 500), function(d)
    seedweight.model(DD = d, W = 30, B = 4, C = 0.02), numeric(1))
  expect_true(all(diff(vals) > 0))
})

# ── exponential model ─────────────────────────────────────────────────────────

test_that("exponential.model: Y[1] = Y0", {
  res <- exponential.model(a = 0.1, Y0 = 100, duration = 5, dt = 1)
  expect_equal(res$Y[1], 100, tolerance = 1e-10)
})

test_that("exponential.model: Y[2] = Y0 * (1 + a*dt) for dt=1", {
  # dY = a * Y[k] * dt  ⟹  Y[2] = 100 * (1 + 0.1) = 110
  res <- exponential.model(a = 0.1, Y0 = 100, duration = 5, dt = 1)
  expect_equal(res$Y[2], 110, tolerance = 1e-8)
})

test_that("exponential.model: follows geometric series Y[k] = Y0*(1+a)^(k-1)", {
  a  <- 0.1; Y0 <- 100; n <- 6
  res <- exponential.model(a = a, Y0 = Y0, duration = 5, dt = 1)
  expected <- Y0 * (1 + a)^(0:5)
  expect_equal(res$Y, expected, tolerance = 1e-8)
})

test_that("exponential.model: returns data.frame with time and Y columns", {
  res <- exponential.model(a = 0.1, Y0 = 50, duration = 10, dt = 1)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "Y"))
  expect_equal(nrow(res), 11)
})

test_that("exponential.model: time column equals seq(0, duration, by=dt)", {
  res <- exponential.model(a = 0.1, Y0 = 100, duration = 4, dt = 0.5)
  expect_equal(res$time, seq(0, 4, by = 0.5), tolerance = 1e-10)
})

test_that("exponential.model.bis matches exponential.model output (golden comparison)", {
  # Both implementations must produce identical trajectories
  res1 <- exponential.model    (a = 0.1, Y0 = 100, duration = 5, dt = 1)
  res2 <- exponential.model.bis(a = 0.1, Y0 = 100, duration = 5, dt = 1)
  expect_equal(res1$Y, res2$Y, tolerance = 1e-8)
})

# ── verhulst (logistic growth) ────────────────────────────────────────────────

test_that("verhulst.update: single-step value is correct", {
  # dY = a * Y * (1 - Y/k) = 0.08 * 1 * 0.99 = 0.0792  ⟹  Y1 = 1.0792
  expect_equal(verhulst.update(Y = 1, a = 0.08, k = 100),
               1 + 0.08 * 1 * (1 - 1/100), tolerance = 1e-10)
})

test_that("verhulst.update: at carrying capacity (Y=k) growth is zero", {
  # dY = a * k * (1 - k/k) = 0  ⟹  Y1 = k
  expect_equal(verhulst.update(Y = 100, a = 0.08, k = 100),
               100, tolerance = 1e-10)
})

test_that("verhulst.model: starts at Y0", {
  res <- verhulst.model(a = 0.08, k = 100, Y0 = 1, duration = 10)
  expect_equal(res$Y[1], 1, tolerance = 1e-10)
})

test_that("verhulst.model: population is monotone increasing below k", {
  res <- verhulst.model(a = 0.08, k = 100, Y0 = 1, duration = 200)
  expect_true(all(diff(res$Y) >= 0))
  expect_true(tail(res$Y, 1) < 100)   # discrete steps never exceed k
})

test_that("verhulst.model: returns data.frame with day and Y columns", {
  res <- verhulst.model(a = 0.08, k = 100, Y0 = 1, duration = 5)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("day", "Y") %in% names(res)))
})

# ── magarey disease model ──────────────────────────────────────────────────────

test_that("magarey.model: at T=Topt the infection efficiency is maximal (fT=1)", {
  # At T=Topt=18, fT=1, so W_required = Wmin/1 = 10
  result <- magarey.model(T = 18, Tmin = 7, Topt = 18, Tmax = 30,
                          Wmin = 10, Wmax = 42)
  expect_equal(result, 10, tolerance = 1e-8)
})

test_that("magarey.model: below Tmin returns Wmax (no infection possible)", {
  result <- magarey.model(T = 5, Tmin = 7, Topt = 18, Tmax = 30,
                          Wmin = 10, Wmax = 42)
  expect_equal(result, 42)
})

test_that("magarey.model: above Tmax returns Wmax (no infection possible)", {
  result <- magarey.model(T = 35, Tmin = 7, Topt = 18, Tmax = 30,
                          Wmin = 10, Wmax = 42)
  expect_equal(result, 42)
})

test_that("magarey.model: returns a single numeric scalar", {
  result <- magarey.model(T = 18, Tmin = 7, Topt = 18, Tmax = 30,
                          Wmin = 10, Wmax = 42)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("magarey.model2: consistent with magarey.model at optimal temperature", {
  # magarey.model2(T, param) — param is a named vector
  r1    <- magarey.model(T = 18, Tmin = 7, Topt = 18, Tmax = 30,
                         Wmin = 10, Wmax = 42)
  param <- c(Tmin = 7, Topt = 18, Tmax = 30, Wmin = 10, Wmax = 42)
  r2    <- magarey.model2(T = 18, param = param)
  # magarey.model2 returns a named scalar; strip names for comparison
  expect_equal(as.numeric(r1), as.numeric(r2), tolerance = 1e-6)
})

test_that("magarey.define.param: returns matrix with expected row/col structure", {
  suppressMessages({
    param <- magarey.define.param("G.citricarpa")
  })
  expect_true(is.matrix(param))
  expect_equal(rownames(param), c("nominal", "binf", "bsup"))
  expect_true(all(c("Tmin", "Topt", "Tmax", "Wmin", "Wmax") %in% colnames(param)))
})

test_that("magarey.define.param: emits message not printed stdout for unknown species", {
  # message() goes to stderr and is suppressable; must NOT error
  expect_message(magarey.define.param("unkown"), "See help")
})

# ── population age model ───────────────────────────────────────────────────────

test_that("population.age.model: returns data.frame with 8 columns (time + 7 stages)", {
  res <- population.age.model(duration = 10, dt = 1)
  expect_s3_class(res, "data.frame")
  expect_equal(ncol(res), 8)
  expect_named(res, c("time", "E", "L1", "L2", "L3", "L4", "P", "A"))
})

test_that("population.age.model: initial egg count is 5, all other stages 0", {
  res <- population.age.model(duration = 5, dt = 1)
  expect_equal(res$E[1],  5, tolerance = 1e-10)
  expect_equal(res$L1[1], 0, tolerance = 1e-10)
  expect_equal(res$A[1],  0, tolerance = 1e-10)
})

test_that("population.age.model: L1 larvae accumulate as eggs hatch (biological flow)", {
  # With E[0]=5, A[0]=0, the egg compartment feeds L1 via hatching rate rE.
  # After 1 time step: dL1 = (rE*E[0] - (r12+m1)*L1[0])*dt = 0.172*5 = 0.86
  # So L1[2] should be approximately 0.86 > 0.
  res <- population.age.model(duration = 5, dt = 1)
  expect_true(res$L1[2] > 0)
})

test_that("population.age.matrix.model: produces same output as loop version", {
  r1 <- population.age.model(duration = 20, dt = 1)
  r2 <- population.age.matrix.model(duration = 20, dt = 1)
  expect_equal(r1, r2, tolerance = 1e-6)
})

# ── carcass growth model ───────────────────────────────────────────────────────

test_that("carcass.model: returns data.frame with expected column names", {
  result <- carcass.model(
    protcmax = 200, protncmax = 100,
    alphac = 0.022, alphanc = 0.022,
    gammac = 0.009, gammanc = 0.009,
    lip0 = 0.1, lipc1 = 0.3, lipnc1 = 0.3,
    beta = 0.041, delta = 0.029,
    k = 0.67, b0c = 4.7, b1c = 0.83,
    b0nc = 5.7, b1nc = 0.79,
    c0 = 1.23, c1 = 1.01, cem = 0.4,
    duration = 10
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("time", "ProtC", "LipC", "ProtNC", "LipNC", "PV") %in%
                    names(result)))
  expect_equal(nrow(result), 10)
})

test_that("carcass.model: ProtC eventually exceeds initial value under growth conditions", {
  # Use high EMI (energy) relative to CPM (maintenance) so synthesis dominates:
  # alphac >> gammac ensures net anabolism.
  result <- carcass.model(
    protcmax = 200, protncmax = 100,
    alphac = 0.10, alphanc = 0.10,  # high synthesis rate
    gammac = 0.005, gammanc = 0.005, # low degradation rate
    lip0 = 0.1, lipc1 = 0.3, lipnc1 = 0.3,
    beta = 0.041, delta = 0.029,
    k = 0.01,                        # low maintenance constant → high EMI/(CPM+EMI)
    b0c = 4.7, b1c = 0.83,
    b0nc = 5.7, b1nc = 0.79,
    c0 = 1.23, c1 = 1.01, cem = 1.5, # high energy input
    duration = 100
  )
  expect_true(result$ProtC[100] > result$ProtC[1])
})

# ── weed population model ──────────────────────────────────────────────────────

test_that("weed.define.param: returns a named numeric matrix with 3 rows", {
  param <- weed.define.param()
  expect_true(is.matrix(param))
  expect_equal(nrow(param), 3)
  expect_equal(rownames(param), c("nominal", "binf", "bsup"))
  expect_true("mu" %in% colnames(param))
  expect_true("Ymax" %in% colnames(param))
})

test_that("weed.model: returns data.frame with year, d, S, SSBa, DSBa, Yield", {
  param <- weed.define.param()["nominal", ]
  deci  <- data.frame(Soil = rep(1, 5), Crop = rep(1, 5), Herb = rep(0, 5))
  result <- weed.model(param, deci)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("year", "d", "S", "SSBa", "DSBa", "Yield") %in%
                    names(result)))
  expect_equal(nrow(result), 6)  # year 0 .. 5
})

test_that("weed.model: herbicide treatment reduces weed density vs no treatment", {
  param  <- weed.define.param()["nominal", ]
  herb   <- data.frame(Soil = rep(1, 5), Crop = rep(1, 5), Herb = rep(1, 5))
  no_herb <- data.frame(Soil = rep(1, 5), Crop = rep(1, 5), Herb = rep(0, 5))
  r_herb   <- weed.model(param, herb)
  r_no_herb <- weed.model(param, no_herb)
  # With herbicide, final weed density should be lower
  expect_true(r_herb$d[6] < r_no_herb$d[6])
})

test_that("weed.model: Yield is non-negative", {
  param <- weed.define.param()["nominal", ]
  deci  <- data.frame(Soil = rep(1, 5), Crop = rep(1, 5), Herb = rep(0, 5))
  result <- weed.model(param, deci)
  expect_true(all(result$Yield >= 0, na.rm = TRUE))
})

# ── lactation model ────────────────────────────────────────────────────────────

test_that("lactation.define.param: returns matrix with correct row names", {
  param <- lactation.define.param("calf")
  expect_true(is.matrix(param))
  expect_equal(rownames(param), c("nominal", "binf", "bsup"))
  expect_true("cu" %in% colnames(param))
})

test_that("lactation.calf.model2: returns a matrix with at least 4 columns", {
  param <- lactation.define.param("calf")["nominal", ]
  result <- lactation.calf.model2(param, duration = 20, dt = 0.1)
  expect_true(is.matrix(result))
  expect_true(ncol(result) >= 4)
})

test_that("lactation.calf.model2: milk M is non-negative throughout", {
  param <- lactation.define.param("calf")["nominal", ]
  result <- lactation.calf.model2(param, duration = 20, dt = 0.1)
  expect_true(all(result[, "M"] >= 0))
})
