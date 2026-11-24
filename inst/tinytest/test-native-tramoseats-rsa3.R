source("helpers/fixtures.R", local = TRUE)
dir <- system.file("tinytest/fixtures/native/tramo-seats-bridge", package = "rkalos")
files <- c("rsa3-x86_64.txt", "rsa3-additional-x86_64.txt", "rsa3-mean-retry-x86_64.txt",
           "rsa3-ref-retry-x86_64.txt")
cases <- do.call(c, lapply(files, function(file) fixture_named_rsa0_cases(file.path(dir, file))))
diagnostics <- c(fixture_named_rsa0_cases(file.path(dir, "rsa3-residual-diagnostics-x86_64.txt")),
                 fixture_named_rsa0_cases(file.path(dir, "rsa3-mean-retry-diagnostics-x86_64.txt")),
                 fixture_named_rsa0_cases(file.path(dir, "rsa3-ref-retry-diagnostics-x86_64.txt")))
expect_identical(sort(names(cases)), sort(names(diagnostics)))
forecast_errors <- c(fixture_named_rsa0_cases(file.path(dir, "rsa3-forecast-errors-x86_64.txt")),
                     fixture_named_rsa0_cases(file.path(dir, "rsa3-ref-retry-forecast-errors-x86_64.txt")))
expect_identical(sort(names(cases)), sort(names(forecast_errors)))
wk_cases <- fixture_named_rsa0_cases(file.path(dir, "rsa3-ref-retry-wk-x86_64.txt"))
names(wk_cases) <- sub(" .*", "", names(wk_cases))
tokens <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]
exact <- function(x, expected, id) {
  expect_identical(writeBin(as.numeric(x), raw()), writeBin(fixture_hex_doubles(expected), raw()), info = id)
}
residual_count <- points <- rows <- 0L
for (id in names(cases)) {
  lines <- cases[[id]]
  field <- function(key) tokens(lines[startsWith(lines, paste0(key, " "))])
  input <- field("input")
  series <- stats::ts(fixture_hex_doubles(input[-seq_len(4L)]), start = as.integer(input[2:3]), frequency = 12)
  pre <- rkalos::native_tramo_preprocess(series, "rsa3")
  fit <- rkalos::native_tramoseats_fit(series, "Tr3")
  expect_identical(c(pre$template, fit$template), c("TR3", "RSA3"), info = id)
  expect_identical(c(pre$status, fit$status), c("available", "ok"), info = id)
  expect_identical(fit$model, pre$model, info = id)
  expect_identical(fit$likelihood, pre$likelihood, info = id)
  expect_identical(fit$residuals, pre$residuals, info = id)
  expect_identical(names(fit$residuals), c("tsres", "ser", "lb", "doornikhansen"), info = id)
  oracle <- diagnostics[[id]]
  test <- function(key) tokens(oracle[startsWith(oracle, paste0("residuals.", key, " "))])
  exact(fit$residuals$ser, test("ser")[2L], id)
  for (key in c("lb", "doornikhansen")) {
    result <- fit$residuals[[key]]
    exact(c(result$value, result$pvalue), test(key)[2:3], paste(id, key))
  }
  expect_identical(fit$residuals$lb$df, 24 - length(fit$model$parameters), info = id)
  expect_identical(fit$residuals$lb$lags, 24L, info = id)
  expect_identical(fit$residuals$doornikhansen$df, 2, info = id)
  forecast_reference <- forecast_errors[[id]]
  errors <- tokens(forecast_reference[startsWith(forecast_reference, "y_ef ")])
  values <- fixture_hex_doubles(errors[-seq_len(4L)])
  actual <- fit$series_forecast_stdev
  expect_identical(fit$series_forecast_stdev_status, "available", info = id)
  expect_identical(as.integer(stats::start(actual)), as.integer(errors[2:3]), info = id)
  expect_identical(length(actual), as.integer(errors[4L]), info = id)
  expect_true(all(abs(as.numeric(actual) - values) <= 2e-8 * (1 + abs(values))), info = id)
  wk <- tokens(forecast_reference[startsWith(forecast_reference, "wk_series_stdev_forecast ")])
  index <- which(fit$final$component == "Series" & fit$final$information == "StdevForecast")
  exact(fit$final$values[[index]], wk[-seq_len(4L)], id)
  if (id %in% names(wk_cases)) {
    wk <- wk_cases[[id]]
    expected_rows <- wk[startsWith(wk, "wk ")]
    expect_identical(nrow(fit$final), length(expected_rows), info = id)
    for (line in expected_rows) {
      row <- tokens(line)
      index <- which(fit$final$component == row[2L] & fit$final$information == row[3L])
      expect_identical(length(index), 1L, info = paste(id, row[2L], row[3L]))
      output <- fit$final$values[[index]]
      expect_identical(as.integer(stats::start(output)), as.integer(row[4:5]), info = id)
      exact(output, row[-seq_len(6L)], paste(id, row[2L], row[3L]))
    }
  }
  custom <- rkalos::native_tramoseats_fit(series, "RSA3", forecast = 7L)
  expect_identical(stats::start(custom$series_forecast_stdev), stats::start(actual), info = id)
  expect_identical(length(custom$series_forecast_stdev), 7L, info = id)
  expect_true(all(abs(as.numeric(custom$series_forecast_stdev) - head(values, 7L)) <=
                    2e-8 * (1 + abs(head(values, 7L)))), info = id)
  expect_identical(fit$model_changed, field("product")[4L] == "true", info = id)
  model <- field("model")
  expect_identical(unname(fit$model$orders), as.integer(model[7:12]), info = id)
  expect_identical(c(fit$model$log, fit$model$mean_correction), model[c(3L, 5L)] == "true", info = id)
  params <- field("parameters")[-1L]
  exact(fit$model$parameters, params, id)
  exact(c(fit$model$phi, fit$model$bphi, fit$model$theta, fit$model$btheta), params, id)
  mean <- lines[startsWith(lines, "variable const ")]
  if (length(mean)) exact(fit$model$mean_coefficient, tokens(mean)[3L], id)
  else expect_identical(fit$model$mean_coefficient, NA_real_, info = id)
  outliers <- lines[startsWith(lines, "variable ") & !startsWith(lines, "variable const ")]
  estimated <- fit$model$outliers
  expect_identical(nrow(estimated), length(outliers), info = id)
  for (i in seq_along(outliers)) {
    row <- tokens(outliers[i])
    date <- as.integer(strsplit(gsub("[()]", "", row[3L]), "-", fixed = TRUE)[[1L]])
    start <- stats::start(series)
    expect_identical(estimated$name[i], paste(row[2:3], collapse = " "), info = id)
    expect_identical(estimated$type[i], row[2L], info = id)
    expect_equal(estimated$position[i], 1 + (date[1L] - start[1L]) * 12 + date[2L] - start[2L], info = id)
    exact(estimated$coefficient[i], row[4L], id)
  }
  ll <- field("likelihood")
  exact(unlist(fit$likelihood[c("log_likelihood", "transformation_adjustment", "adjusted_log_likelihood", "ssq")]),
        ll[2:5], id)
  expect_identical(unname(unlist(fit$likelihood[c("effective", "estimated")])), as.integer(ll[6:7]), info = id)
  residual <- field("residuals")
  tsres <- fit$residuals$tsres
  expect_identical(as.integer(stats::start(tsres)), as.integer(residual[2:3]), info = id)
  expect_identical(stats::frequency(tsres), 12, info = id)
  exact(tsres, residual[-seq_len(4L)], id)
  residual_count <- residual_count + length(tsres)
  expected <- lines[startsWith(lines, "merged ")]
  expect_identical(nrow(fit$merged), length(expected), info = id)
  for (i in seq_along(expected)) {
    row <- tokens(expected[i])
    value <- fit$merged$values[[i]]
    expect_identical(c(fit$merged$component[i], fit$merged$information[i]), row[2:3], info = id)
    expect_identical(as.integer(stats::start(value)), as.integer(row[4:5]), info = id)
    expect_identical(stats::frequency(value), 12, info = id)
    expect_identical(length(value), as.integer(row[6L]), info = id)
    exact(value, row[-seq_len(6L)], id)
    points <- points + length(value)
  }
  rows <- rows + length(expected)
}
expect_identical(c(length(cases), residual_count, rows, points), c(20L, 2105L, 200L, 13980L))
scope <- rkalos::capabilities()$tramo_seats_fit$by_template$RSA3
expect_identical(scope$residuals, c("tsres", "ser", "lb", "doornikhansen"))
expect_true(scope$series_forecast_stdev)
expect_true(scope$seasonal_restart)
expect_true(scope$fast_mean_refit)
expect_identical(scope$identification, "reference_retries")
zero <- rkalos::native_tramoseats_fit(AirPassengers, "RSA3", forecast = 0L)
expect_identical(zero$series_forecast_stdev_status, "no_requested_predictions")
expect_identical(zero$series_forecast_stdev, numeric())
custom <- rkalos::native_tramoseats_fit(AirPassengers, "RSA3", forecast = 6L, backcast = 3L)
expect_identical(length(custom$merged$values[[which(custom$merged$component == "Trend" & custom$merged$information == "Forecast")]]), 6L)
expect_identical(length(custom$merged$values[[which(custom$merged$component == "Trend" & custom$merged$information == "Backcast")]]), 3L)
i <- seq_len(36L)
retry <- stats::ts(round(2000 + i^2 / 100 + 70 * sin(2 * pi * i / 12) + 10 * cos(i)), frequency = 12)
expect_identical(rkalos::native_tramo_preprocess(retry, "TR3")$status, "available")
expect_error(rkalos::native_tramoseats_fit(retry, "RSA3"), "approximation_required")
short <- rkalos::native_tramoseats_fit(stats::ts(as.numeric(AirPassengers)[1:37], frequency = 12), "RSA3")
expect_identical(length(short$residuals$tsres), 24L)
expect_true(is.nan(short$residuals$lb$value))
expect_true(is.nan(short$residuals$lb$pvalue))
expect_true(is.finite(short$residuals$doornikhansen$pvalue))
