source("helpers/fixtures.R", local = TRUE)
directory <- system.file("tinytest/fixtures/native", package = "rkalos", mustWork = TRUE)
cases <- fixture_named_rsa0_cases(file.path(directory, "x13-estimated-templates/templates-x86_64.txt"))
diagnostics <- fixture_named_rsa0_cases(file.path(directory, "x13-estimated-templates/diagnostics-x86_64.txt"))
tokens <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]
rows <- innovations <- points <- 0L
for (header in names(cases)) {
  parts <- tokens(header)
  input <- readLines(file.path(directory, "x13-refusal-jvm-profiles", paste0(parts[2L], ".csv")))
  start <- as.integer(strsplit(input[1L], ",", fixed = TRUE)[[1L]])
  series <- stats::ts(fixture_doubles(input[-1L]), start = start, frequency = 12)
  fit <- rkalos::native_x13_estimated_fit(series, parts[1L])
  expect_identical(tolower(fit$template), parts[1L], info = header)
  expect_identical(fit$status, "available", info = header)
  expect_identical(class(fit), c("native_x13_estimated_fit", "list"), info = header)
  expect_identical(fit$series, series, info = header)
  expect_identical(c(fit$effective_backcast, fit$effective_forecast), c(0L, 12L), info = header)
  diagnostic <- diagnostics[[header]]
  item <- function(name) tokens(diagnostic[startsWith(diagnostic, paste0(name, " "))])
  expect_identical(fit$residuals$ser, fixture_hex_doubles(item("residuals.ser")[2L]), info = header)
  for (name in c("lb", "doornikhansen")) {
    ref <- item(paste0("residuals.", name))
    actual <- fit$residuals[[name]]
    expect_identical(c(actual$value, actual$pvalue), fixture_hex_doubles(ref[2:3]), info = paste(header, name))
    expect_identical(as.numeric(actual$df), as.numeric(ref[6L]), info = paste(header, name))
  }
  ref <- item("y_ef_12")
  errors <- fit$series_forecast_stdev
  expected <- fixture_hex_doubles(ref[-seq_len(4L)])
  expect_identical(fit$series_forecast_stdev_status, "available", info = header)
  expect_identical(as.integer(stats::start(errors)), as.integer(ref[2:3]), info = header)
  expect_identical(stats::tsp(errors), stats::tsp(fit$preadjustment$tables$a1a), info = header)
  expect_identical(length(errors), length(expected), info = header)
  expect_true(all(abs(as.numeric(errors) - expected) <= 2e-8 * (1 + abs(expected))), info = header)
  lines <- cases[[header]]
  field <- function(key) tokens(lines[startsWith(lines, paste0(key, " "))])
  model <- field("model")
  expect_identical(c(fit$model$log, fit$model$mean), model[2:3] == "true", info = header)
  expect_identical(unname(unlist(fit$model[c("p", "d", "q", "bp", "bd", "bq")])), as.integer(model[4:9]), info = header)
  expect_identical(unname(unlist(fit$model$coefficients)), fixture_hex_doubles(field("arma")[-1L]), info = header)
  vars <- lapply(lines[startsWith(lines, "variable ")], tokens)
  labels <- unname(vapply(vars, `[[`, character(1), 2L))
  is_calendar <- labels %in% c("td", "lp", "easter")
  if (parts[1L] == "rsa2c") {
    calendar <- fit$model$calendar
    expect_identical(names(calendar$coefficients), labels[is_calendar], info = header)
    expect_identical(calendar$working_days, "td" %in% labels, info = header)
    expect_identical(calendar$easter_duration > 0, "easter" %in% labels, info = header)
  }
  vars <- vars[!is_calendar]
  fitted <- fit$model$variables
  expect_identical(fitted$name, unname(vapply(vars, function(v) paste(v[2:3], collapse = " "), "")), info = header)
  expected <- fixture_hex_doubles(tokens(lines[grepl("^coefficients( |$)", lines)])[-1L])
  values <- c(unname(fit$model$calendar$coefficients), fitted$coefficient)
  expect_true(all(abs(values - expected) <= 1e-8 * pmax(1, abs(expected))), info = header)
  ll <- field("likelihood")
  expect_identical(unname(unlist(fit$likelihood[c("effective", "estimated")])), as.numeric(ll[2:3]), info = header)
  expect_identical(unname(unlist(fit$likelihood[c("log_likelihood", "ssq")])), fixture_hex_doubles(ll[4:5]), info = header)
  residual <- field("residuals")
  res <- fit$residuals$tsres
  expect_identical(as.integer(stats::start(res)), as.integer(residual[2:3]), info = header)
  expect_identical(writeBin(as.numeric(res), raw()), writeBin(fixture_hex_doubles(residual[-seq_len(4L)]), raw()), info = header)
  innovations <- innovations + length(res)
  selection <- field("selected")
  expect_identical(fit$mode_effective, tolower(selection[2L]), info = header)
  expect_identical(c(fit$selection$stages$d7$length, fit$selection$stages$final$length), as.integer(selection[3:4]), info = header)
  expect_identical(as.integer(fit$selection$d$seas_filter), as.integer(sub(".*X([0-9]+).*", "\\1", selection[6:17])), info = header)
  tables <- list(b1 = fit$preadjustment$tables$out, y_f = fit$preadjustment$tables$a1a,
                 b7 = fit$selection$b$tables$B7, c7 = fit$selection$c$tables$C7, d7 = fit$selection$d$tables$D7,
                 d10 = fit$results$d10, d11 = fit$results$d11, d12 = fit$results$d12, d13 = fit$results$d13,
                 final_sa = fit$finalSA, final_trend = fit$finalTrend, final_irregular = fit$finalIrregular,
                 forecast_sa = fit$finals$tables$d11a, forecast_trend = fit$finals$tables$d12a)
  for (key in names(tables)) {
    row <- field(key)
    actual <- tables[[key]]
    expected <- fixture_hex_doubles(row[-seq_len(4L)])
    expect_identical(as.integer(stats::start(actual)), as.integer(row[2:3]), info = paste(header, key))
    expect_identical(stats::frequency(actual), 12, info = paste(header, key))
    expect_identical(length(actual), length(expected), info = paste(header, key))
    expect_true(all(abs(as.numeric(actual) - expected) <= 1e-8 * pmax(1, abs(expected))), info = paste(header, key))
    points <- points + length(actual)
  }
  rows <- rows + length(tables)
}
expect_identical(c(rows, innovations), c(126L, 891L))
expect_identical(rkalos::capabilities()$estimated_x13_fit$templates, c("RSA0", "RSA1", "RSA2c"))
expect_error(rkalos::native_x13_estimated_fit(UKgas), "monthly")
expect_error(rkalos::native_x13_estimated_fit(replace(AirPassengers, 2L, NA_real_)), "finite")
expect_error(rkalos::native_x13_estimated_fit(AirPassengers, "RSA3"), "RSA0")
