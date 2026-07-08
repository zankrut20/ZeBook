# NEWS

## ZeBook 1.5.0

### Performance

* `maize.multisy()`: replaced O(n²) `rbind`-grow loop with pre-allocated list
  and single `do.call(rbind, ...)`. Reduces wall time by **~38%** for 10
  site-years; savings scale quadratically with the number of site-years.

* `param.runif()`, `param.rtriangle()`: replaced `cbind`-grow-in-loop with
  a single pre-allocated matrix filled column-by-column, eliminating O(N·p²)
  intermediate data.frame copies. `param.rtriangle` shows **4× speedup** for
  `N=2000, p=7`.

* `q.arg.fast.runif()`: replaced `c(list, ...)` list-grow with a pre-allocated
  `vector("list", p)` filled by index, eliminating O(p²) list copies.

* `carcass.model()`: (a) fixed `rep(NA, 1, duration)` bug (wrong 3-arg form);
  (b) cached `EMI/(CPM+EMI)` energy ratio computed once per step instead of
  4 times; (c) cached all four `log()` growth terms computed once per step
  instead of twice. Saves ~2920 redundant `log()` calls per 365-day run.

* `carbonsoil.model()`: when `U` is a scalar, now converts it directly to a
  plain `numeric` vector instead of constructing a 2-column `data.frame` per
  call.

* `lactation.define.param()`: adopted `.make_param_matrix()` helper,
  consistent with Phase 3 consolidation.

---

## ZeBook 1.4.0


### Internal refactoring

* Introduced `R/utils_internal.r` containing four unexported helper functions
  that consolidate duplicated patterns found across model files:

  * `.filter_site_year()`: eliminates the 9-line site-and-year filtering
    block previously duplicated in `maize.weather()`, `epirice.weather()`, and
    `watbal.weather()` (27 lines removed across 3 files).

  * `.make_param_matrix()`: eliminates the `row.names + as.matrix` two-liner
    previously repeated in `maize.define.param()`, `epirice.define.param()`,
    `watbal.define.param()`, and `magarey.define.param()` (8 lines removed
    across 4 files).

  * `.apply_simule()`: eliminates the `apply/cbind/as.matrix` pattern
    previously repeated in `maize.simule()`, `maize.simule240()`,
    `maize.simule_multisy240()`, `magarey.simule()`, and `weed.simule()`
    (15 lines removed across 5 functions).

  * `.parse_siteyear()`: eliminates the inline `strsplit(sy,"-")[[1]][1/2]`
    repeated in `maize.multisy()` and `maize.multisy240()` (4 expressions
    replaced across 2 functions).

* No exported function signatures, return types, or return values were changed.

* 31 new tests in `tests/testthat/test-utils-internal.R` verify all four
  helpers directly via the `:::` operator.

---

## ZeBook 1.3.0


### Testing

* Established a complete `testthat` (edition 3) testing foundation with 275
  passing tests across four test files, covering every core model domain in the
  package.

* `tests/testthat/test-pure-models.R` (65 tests): golden-output and invariant
  tests for `carbonsoil`, `seedweight`, `exponential`, `verhulst`, `magarey`,
  `population.age`, `carcass`, `weed`, and `lactation` models.

* `tests/testthat/test-weather-models.R` (35 tests): structural, invariant, and
  snapshot tests for `carrot.weevil`, `maize`, `watbal`, and `epirice` models,
  using deterministic synthetic weather fixtures (no network access required).

* `tests/testthat/test-statistical.R` (32 tests): analytical-value tests for
  `goodness.of.fit`, `threshold.measures`, `AICf`, and `evaluation.criteria`,
  including the MSE decomposition identity (bias² + SDSD + LCS = MSE).

* `tests/testthat/test-sampling.R` (143 tests): structural, boundary, and
  distributional tests for `param.runif`, `param.rtriangle`, and
  `q.arg.fast.runif`, including a theoretical-mean check at N=5000.

* `tests/testthat/helper-weather.R`: shared synthetic weather/parameter
  factories used throughout the test suite.

* Phase 1 regressions protected: `expect_silent(carrot.weevil.model(...))` and
  `expect_message(magarey.define.param("unkown"), ...)` guard against
  reintroduction of the print-statement bugs fixed in v1.2.1.

---

## ZeBook 1.2.1


### Bug fixes

* `carrot.weevil.model()`: Removed a debugging `print(dTT)` statement that was
  inadvertently left inside the day-by-day simulation loop. For a standard 360-day
  run the statement fired 359 times, producing hundreds of lines of unwanted console
  output and making the function unusable inside `lapply`, `apply`, or any
  sensitivity-analysis workflow. Return value is unchanged. (#P0-PERF-04)

* `epirice.multi.simule()`: Removed an unconditional `print(AUDPC)` statement that
  printed the full result vector to stdout on every call. CRAN policy prohibits
  exported functions from writing to stdout unless the user explicitly requests it.
  Return value is unchanged. (#P0-CRAN-02)

* `magarey.define.param()`: Replaced two `print()` calls with `message()`. The
  status messages ("See help for values for other species." and "parameter values
  for species : <X>") are preserved but now written to stderr, which is the correct
  channel for optional informational output and can be suppressed by callers with
  `suppressMessages()`. Return value is unchanged. (#P1-API-04)
