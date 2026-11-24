source("helpers/fixtures.R", local = TRUE)
directory <- system.file("tinytest/fixtures/native/tramo-seats-bridge", package = "rkalos", mustWork = TRUE)
panel <- function(suffix) fixture_named_rsa0_cases(file.path(directory, paste0("rsafull", suffix, "-x86_64.txt")))
cases <- c(panel(""), panel("-late-calendar"))
calendars <- c(panel("-calendar-effect"), panel("-late-calendar-effect"))
diagnostics <- c(panel("-residual-diagnostics"), panel("-late-residual-diagnostics"))
forecasts <- c(panel("-forecast-errors"), panel("-late-forecast-errors"))
for (oracle in list(calendars, diagnostics, forecasts)) expect_identical(sort(names(oracle)), sort(names(cases)))
tokens <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]
row <- function(lines, key) tokens(lines[startsWith(lines, paste0(key, " "))])
exact <- function(x, hex, label) {
  expect_identical(writeBin(as.numeric(x), raw()), writeBin(fixture_hex_doubles(hex), raw()), info = label)
}
check_ts <- function(x, fields, label) {
  expect_true(stats::is.ts(x), info = label)
  expect_identical(as.integer(stats::start(x)), as.integer(fields[2:3]), info = label)
  expect_identical(stats::frequency(x), 12, info = label)
  expect_identical(length(x), as.integer(fields[4L]), info = label)
  exact(x, fields[-seq_len(4L)], label)
}
innovations <- points <- 0L
for (id in names(cases)) {
  lines <- cases[[id]]
  input <- row(lines, "input")
  series <- stats::ts(fixture_hex_doubles(input[-seq_len(4L)]), frequency = 12, start = as.integer(input[2:3]))
  pre <- rkalos::native_tramo_preprocess(series, "rsafull")
  fit <- rkalos::native_tramoseats_fit(series, "TRFULL")
  expect_identical(c(pre$template, fit$template), c("TRfull", "RSAfull"), info = id)
  expect_identical(c(pre$status, fit$status), c("available", "ok"), info = id)
  for (key in c("model", "residuals", "likelihood")) expect_identical(pre[[key]], fit[[key]], info = paste(id, key))
  model <- row(lines, "model")
  expect_identical(c(fit$model$log, fit$model$mean_correction), model[c(3, 5)] == "true", info = id)
  expect_identical(unname(fit$model$orders), as.integer(model[7:12]), info = id)
  exact(fit$model$parameters, row(lines, "parameters")[-1L], id)
  variables <- lapply(lines[startsWith(lines, "variable ")], tokens)
  labels <- vapply(variables, `[[`, character(1L), 2L)
  estimates <- vapply(variables, function(fields) tail(fields, 1L), character(1L), USE.NAMES = FALSE)
  if ("const" %in% labels) exact(fit$model$mean_coefficient, estimates[labels == "const"], id)
  else expect_identical(fit$model$mean_coefficient, NA_real_, info = id)
  calendar <- fit$model$calendar
  selected <- startsWith(labels, "td") | labels %in% c("lp", "easter")
  expect_identical(names(calendar$coefficients), unname(labels[selected]), info = id)
  exact(calendar$coefficients, estimates[selected], id)
  days <- if ("td" %in% labels) "working_days" else if ("td[0]" %in% labels) "trading_days" else "none"
  expect_identical(calendar$trading_days, days, info = id)
  expect_identical(calendar$leap_year, "lp" %in% labels, info = id)
  expect_identical(calendar$easter, "easter" %in% labels, info = id)
  outliers <- !selected & labels != "const"
  expect_identical(fit$model$outliers$name, unname(vapply(variables[outliers], function(v) paste(v[2:3], collapse = " "), "")), info = id)
  exact(fit$model$outliers$coefficient, estimates[outliers], id)
  likelihood <- row(lines, "likelihood")
  exact(unlist(fit$likelihood[1:4]), likelihood[2:5], id)
  expect_identical(unname(unlist(fit$likelihood[5:6])), as.integer(likelihood[6:7]), info = id)
  check_ts(fit$residuals$tsres, row(lines, "residuals"), id)
  innovations <- innovations + length(fit$residuals$tsres)
  diag <- diagnostics[[id]]
  exact(fit$residuals$ser, row(diag, "residuals.ser")[2L], id)
  for (key in c("lb", "doornikhansen")) {
    test <- fit$residuals[[key]]
    exact(c(test$value, test$pvalue), row(diag, paste0("residuals.", key))[2:3], id)
  }
  expect_identical(fit$model_changed, row(lines, "product")[4L] == "true", info = id)
  merged <- lapply(lines[startsWith(lines, "merged ")], tokens)
  expect_identical(nrow(fit$merged), length(merged), info = id)
  for (i in seq_along(merged)) {
    ref <- merged[[i]]
    expect_identical(c(fit$merged$component[i], fit$merged$information[i]), ref[2:3], info = id)
    check_ts(fit$merged$values[[i]], ref[-(1:2)], paste(id, ref[2L], ref[3L]))
    points <- points + length(fit$merged$values[[i]])
  }
  check_ts(fit$calendar_effect, row(calendars[[id]], "cal"), id)
  check_ts(fit$calendar_effect_forecast, row(calendars[[id]], "cal_f"), id)
  expect_identical(stats::tsp(fit$calendar_effect), stats::tsp(series), info = id)
  errors <- row(forecasts[[id]], "y_ef")
  expected <- fixture_hex_doubles(errors[-seq_len(4L)])
  expect_identical(fit$series_forecast_stdev_status, "available", info = id)
  expect_identical(as.integer(stats::start(fit$series_forecast_stdev)), as.integer(errors[2:3]), info = id)
  expect_true(all(abs(as.numeric(fit$series_forecast_stdev) - expected) <= 2e-8 * (1 + abs(expected))), info = id)
  index <- which(fit$final$component == "Series" & fit$final$information == "StdevForecast")
  check_ts(fit$final$values[[index]], row(forecasts[[id]], "wk_series_stdev_forecast"), id)
}
expect_identical(c(length(cases), innovations, points), c(29L, 4165L, 26160L))
zero <- rkalos::native_tramoseats_fit(AirPassengers, "RSAfull", forecast = 0L)
expect_identical(zero$calendar_effect_forecast, numeric())
expect_identical(zero$series_forecast_stdev_status, "no_requested_predictions")
short <- rkalos::native_tramoseats_fit(AirPassengers, "RSAfull", forecast = 7L, backcast = 2L)
expect_identical(length(short$calendar_effect_forecast), 7L)
expect_identical(as.integer(stats::start(short$calendar_effect_forecast)), c(1961L, 1L))
expect_identical(rkalos::capabilities()$tramo_seats_fit$by_template$RSAfull$calendar_effect, c("observed", "forecast"))
