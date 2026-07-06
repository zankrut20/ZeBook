# NEWS

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
