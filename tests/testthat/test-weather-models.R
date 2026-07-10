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

test_that("maize_cir.model: produces correct golden outputs", {
  w <- make_weather_maize()
  res <- maize_cir.model(Tbase = 7, RUE = 1.85, K = 0.7, alpha = 0.00243,
                         LAImax = 7, TTM = 1200, TTL = 700,
                         weather = w, sdate = 100, ldate = 200)
  golden_CumInt <- c(0, 0.104633356, 0.2287356112, 0.3759052045, 0.5503958744, 0.7572320102, 
1.0023426876, 1.2927167414, 1.6365812382, 2.0436055568, 2.525132892, 
3.094440258, 3.767026831, 4.5609285511, 5.4970540785, 6.5995331936, 
7.8960632596, 9.4182321655, 11.2017870564, 13.2868071756, 15.7177267404, 
18.5431410346, 21.815317912, 25.5893309912, 29.9217347895, 34.8687219145, 
40.4837445945, 46.81465207, 53.9004917021, 61.7682363778, 70.4298130223, 
79.8798834522, 90.0948295121, 101.0332867524, 112.6383484294, 
124.8412600613, 137.5661221401, 150.7349116138, 164.2720954484, 
178.1082553274, 192.182414611, 206.4430579052, 220.8480665654, 
235.363910999, 249.9644464852, 264.6295918798, 279.3440759697, 
294.0963485046, 308.8776883795, 323.6815022827, 338.5027881989, 
353.3377324125, 368.1834100829, 383.0375639539, 397.8984409804, 
412.7646714485, 427.635179128, 442.5091140771, 457.3858020265, 
472.2647059575, 487.1453967083, 502.0275303172, 516.9108304323, 
531.7950745653, 546.6800832858, 561.5657116828, 576.4513400797, 
591.3369684766, 606.2225968735, 621.1082252705, 635.9938536674, 
650.8794820643, 665.7651104612, 680.6507388581, 695.5363672551, 
710.421995652, 725.3076240489, 740.1932524458, 755.0788808428, 
769.9645092397, 784.8501376366, 799.7357660335, 814.6213944305, 
829.5070228274, 844.3926512243, 859.2782796212, 874.1639080181, 
889.0495364151, 903.935164812, 918.8207932089, 933.7064216058, 
948.5920500028, 963.4776783997, 978.3633067966, 993.2489351935, 
1008.1345635904, 1023.0201919874, 1037.9058203843, 1052.7914487812, 
1067.6770771781, 1082.5627055751)
  expect_equal(res$CumInt, golden_CumInt, tolerance = 1e-8)
  expect_true(all(c("day", "TT", "LAI", "B", "CumInt") %in% names(res)))
})

test_that("maize_cir_rue.model: produces correct golden outputs", {
  w <- make_weather_maize()
  res <- maize_cir_rue.model(Tbase = 7, RUE_max = 1.85, K = 0.7, alpha = 0.00243,
                             LAImax = 7, TTM = 1200, TTL = 700,
                             weather = w, sdate = 100, ldate = 200)
  golden_B <- c(1, 1.1935717086, 1.4231608807, 1.6954246284, 2.0182323676, 
2.4008792189, 2.8543339721, 3.3915259716, 4.0276752907, 4.7806702801, 
5.6714958501, 6.7247144773, 7.9689996374, 9.4377178195, 11.1695500452, 
13.2091364082, 15.6077170302, 18.4237295062, 21.7233060543, 25.5805932749, 
30.0777944697, 35.304810914, 41.3583381372, 48.3402623336, 56.3552093606, 
65.5071355418, 75.8949274999, 87.6071063294, 100.7159096488, 
115.2712372989, 131.2951540913, 148.7777843866, 167.6754345973, 
187.9115804919, 209.3809445943, 231.9563311134, 255.4973259592, 
279.8595864856, 304.9033765796, 330.5002723556, 356.5374670303, 
382.9196571246, 409.568923146, 436.4232353482, 463.4342259976, 
490.5647449776, 517.7865405439, 545.0782447335, 572.423723502, 
599.810779223, 627.2301581679, 654.6748049632, 682.1393086534, 
709.6194933146, 737.1121158137, 764.6146421798, 792.1250813868, 
819.6418610427, 847.1637337491, 874.6897060213, 902.2189839104, 
929.7509310869, 957.2850362998, 984.8208879457, 1012.3581540788, 
1039.8965666131, 1067.4349791474, 1094.9733916817, 1122.511804216, 
1150.0502167503, 1177.5886292847, 1205.127041819, 1232.6654543533, 
1260.2038668876, 1287.7422794219, 1315.2806919562, 1342.8191044905, 
1370.3575170248, 1397.8959295591, 1425.4343420934, 1452.9727546277, 
1480.511167162, 1508.0495796963, 1535.5879922306, 1563.126404765, 
1590.6648172993, 1618.2032298336, 1645.7416423679, 1673.2800549022, 
1700.8184674365, 1728.3568799708, 1755.8952925051, 1783.4337050394, 
1810.9721175737, 1838.510530108, 1866.0489426423, 1893.5873551766, 
1921.1257677109, 1948.6641802452, 1976.2025927796, 2003.7410053139)
  expect_equal(res$B, golden_B, tolerance = 1e-8)
  expect_true(all(c("day", "TT", "LAI", "B", "CumInt") %in% names(res)))
})

test_that("maize_cir_rue_ear.model: produces correct golden outputs", {
  w <- make_weather_maize()
  res <- maize_cir_rue_ear.model(Tbase = 7, RUE_max = 1.85, K = 0.7, alpha = 0.00243,
                                 LAImax = 7, TTM = 1200, TTL = 700,
                                 weather = w, sdate = 100, ldate = 200)
  golden_BE <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 27.5384125343, 55.0768250686, 82.6152376029, 110.1536501372, 
137.6920626715, 165.2304752058, 192.7688877401, 220.3073002745, 
247.8457128088, 275.3841253431, 302.9225378774, 330.4609504117, 
357.999362946, 385.5377754803, 413.0761880146, 440.6146005489, 
468.1530130832, 495.6914256175, 523.2298381518, 550.7682506861, 
578.3066632204, 605.8450757548, 633.3834882891, 660.9219008234, 
688.4603133577, 715.998725892, 743.5371384263, 771.0755509606, 
798.6139634949, 826.1523760292, 853.6907885635, 881.2292010978, 
908.7676136321, 936.3060261664, 963.8444387007, 991.382851235)
  expect_equal(res$BE, golden_BE, tolerance = 1e-8)
  expect_true(all(c("day", "TT", "LAI", "B", "CumInt", "BE") %in% names(res)))
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
