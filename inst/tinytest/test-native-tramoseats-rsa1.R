source("helpers/fixtures.R", local = TRUE)
dir <- system.file("tinytest/fixtures/native/tramo-seats-bridge", package = "rkalos")
cases <- fixture_named_rsa0_cases(file.path(dir, "rsa1-x86_64.txt"))
additional_cases <- function(file, prefix) {
  x <- fixture_named_rsa0_cases(file.path(dir, file))
  stats::setNames(x, paste0(prefix, names(x)))
}
cases <- c(cases, additional_cases("rsa1-level-effects-x86_64.txt", "level_"),
           additional_cases("rsa1-level-tc-x86_64.txt", "level_tc_"),
           additional_cases("rsa1-log-tc-x86_64.txt", "log_tc_"))
tokens <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]
keys <- c(Trend = "finalTrend", SeasonallyAdjusted = "finalSA", Seasonal = "finalSeasonal", Irregular = "finalIrregular")
nrows <- points <- 0L
for (id in names(cases)) {
  lines <- cases[[id]]
  input <- tokens(lines[startsWith(lines, "input ")])
  series <- stats::ts(fixture_hex_doubles(input[-seq_len(4L)]), start = as.integer(input[2:3]), frequency = 12)
  fit <- rkalos::native_tramoseats_fit(series, "RSA1")
  pre <- rkalos::native_tramo_preprocess(series, "TR1")
  expect_identical(fit$model, pre$model, info = id)
  expect_identical(fit$likelihood, pre$likelihood, info = id)
  expect_identical(fit$residuals, pre$residuals, info = id)
  expect_identical(fit$template, "RSA1", info = id)
  expect_identical(fit$status, "ok", info = id)
  expect_identical(fit$mode, if (pre$model$log) "multiplicative" else "additive", info = id)
  product <- tokens(lines[startsWith(lines, "product ")])
  expect_identical(fit$model_changed, product[4] == "true", info = id)
  expect_identical(fit$series_forecast_stdev_status, "not_qualified", info = id)
  expect_identical(fit$series_forecast_stdev, numeric(), info = id)
  oracle <- lines[startsWith(lines, "merged ")]
  expect_identical(nrow(fit$merged), length(oracle), info = id)
  for (i in seq_along(oracle)) {
    row <- tokens(oracle[i])
    values <- fit$merged$values[[i]]
    expect_identical(c(fit$merged$component[i], fit$merged$information[i]), row[2:3], info = id)
    expect_identical(as.integer(stats::start(values)), as.integer(row[4:5]), info = id)
    expect_identical(stats::frequency(values), 12, info = id)
    expect_identical(length(values), as.integer(row[6L]), info = id)
    expect_identical(writeBin(as.numeric(values), raw()), writeBin(fixture_hex_doubles(row[-seq_len(6L)]), raw()), info = id)
    if (row[3L] == "Value" && row[2L] %in% names(keys)) {
      final <- fit[[keys[[row[2L]]]]]
      expect_identical(stats::tsp(final), stats::tsp(series), info = id)
      expect_identical(as.numeric(final), as.numeric(values), info = id)
    }
    points <- points + length(values)
  }
  nrows <- nrows + nrow(fit$merged)
  if (id == "fit_multi" || startsWith(id, "level_") || startsWith(id, "log_tc_")) {
    # The raw route and independently constructed fitted bridge receive the same
    # estimated outlier model and custom forecast and backcast month counts.
    model <- pre$model
    fit <- rkalos::native_tramoseats_fit(series, "tr1", forecast = 6L, backcast = 3L)
    bridge <- rkalos::seats_from_tramo(series, pre$likelihood$ssq, pre$likelihood$effective,
      pre$likelihood$estimated, d = model$orders[["d"]], bd = model$orders[["bd"]],
      theta = model$theta, btheta = model$btheta, log = model$log, mean_correction = model$mean_correction,
      outliers = model$outliers[c("type", "position", "coefficient")], forecast = 6L, backcast = 3L)
    for (key in c("status", "initial", "final", "merged")) expect_identical(fit[[key]], bridge[[key]], info = key)
  }
}
expect_identical(c(length(cases), nrows, points), c(12L, 120L, 7560L))
zero <- rkalos::native_tramoseats_fit(AirPassengers, "RSA1", forecast = 0L)
expect_identical(unique(zero$merged$information), "Value")
expect_identical(zero$series_forecast_stdev_status, "no_requested_predictions")
scope <- rkalos::capabilities()$tramo_seats_fit$by_template
expect_identical(scope$RSA1$residuals, "tsres")
expect_identical(scope$RSA1$outlier_types, c("AO", "LS", "TC"))
expect_identical(scope$RSA1$tc_rate, 0.7)
expect_false(scope$RSA1$series_forecast_stdev)
expect_true(scope$RSA0$series_forecast_stdev)
expect_error(rkalos::native_tramoseats_fit(stats::ts(1:23, frequency = 12), "RSA1"), "two years")
