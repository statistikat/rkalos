source("helpers/fixtures.R", local = TRUE)
dir <- system.file("tinytest/fixtures/native/tramo-seats-bridge", package = "rkalos")
inputs <- fixture_named_rsa0_cases(file.path(dir, "rsa0-x86_64.txt"))
errors <- fixture_named_rsa0_cases(file.path(dir, "rsa0-forecast-errors-x86_64.txt"))
n <- 0L
for (id in names(errors)) {
  fields <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]
  input <- fields(inputs[[id]][2L])
  x <- stats::ts(fixture_hex_doubles(input[-seq_len(4L)]), start = as.integer(input[2:3]), frequency = 12)
  fit <- rkalos::native_tramoseats_fit(x)
  pre <- fields(errors[[id]][2L])
  want <- fixture_hex_doubles(pre[-seq_len(4L)])
  actual <- fit$series_forecast_stdev
  expect_identical(fit$series_forecast_stdev_status, "available", info = id)
  expect_identical(as.integer(stats::start(actual)), as.integer(pre[2:3]), info = id)
  expect_identical(stats::frequency(actual), 12, info = id)
  expect_identical(length(actual), as.integer(pre[4L]), info = id)
  expect_true(all(abs(as.numeric(actual) - want) <= 2e-8 * (1 + abs(want))), info = id)
  wk <- fields(errors[[id]][3L])
  row <- fit$final[fit$final$component == "Series" & fit$final$information == "StdevForecast", ]
  expect_identical(nrow(row), 1L, info = id)
  se <- row$values[[1L]]
  expect_identical(as.integer(stats::start(se)), as.integer(wk[2:3]), info = id)
  expect_identical(writeBin(as.numeric(se), raw()), writeBin(fixture_hex_doubles(wk[-seq_len(4L)]), raw()), info = id)
  expect_false(identical(as.numeric(actual), as.numeric(se)), info = id)
  n <- n + length(want)
}
expect_identical(n, 144L)
short <- rkalos::native_tramoseats_fit(x, forecast = 6L, backcast = 3L)
expect_identical(short$series_forecast_stdev_status, "available")
expect_identical(length(short$series_forecast_stdev), 6L)
expect_identical(stats::start(short$series_forecast_stdev), stats::start(actual))
expect_true(all(abs(as.numeric(short$series_forecast_stdev) - want[1:6]) <= 2e-8 * (1 + abs(want[1:6]))))
zero <- rkalos::native_tramoseats_fit(x, forecast = 0L)
expect_identical(zero$series_forecast_stdev_status, "no_requested_predictions")
expect_identical(zero$series_forecast_stdev, numeric())
expect_identical(rkalos::capabilities()$tramo_seats_fit$series_forecast_stdev, TRUE)
