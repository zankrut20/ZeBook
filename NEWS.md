# NEWS

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
