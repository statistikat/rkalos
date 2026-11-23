source("helpers/fixtures.R", local = TRUE)
dir <- system.file("tinytest/fixtures/native/tramo-seats-bridge", package = "rkalos")
cases <- fixture_named_rsa0_cases(file.path(dir, "rsa1-x86_64.txt"))
tokens <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]
exact <- function(actual, expected, id) {
  expect_identical(writeBin(as.numeric(actual), raw()), writeBin(fixture_hex_doubles(expected), raw()), info = id)
}
residual_count <- variables <- 0L
for (id in names(cases)) {
  lines <- cases[[id]]
  field <- function(key) tokens(lines[startsWith(lines, paste0(key, " "))])
  input <- field("input")
  start <- as.integer(input[2:3])
  series <- stats::ts(fixture_hex_doubles(input[-seq_len(4L)]), start = start, frequency = 12)
  fit <- rkalos::native_tramo_preprocess(series)
  expect_identical(class(fit), c("native_tramo_preprocess", "list"), info = id)
  expect_identical(fit$status, "available", info = id)
  expect_identical(fit$template, "TR1", info = id)
  model <- field("model")
  expect_identical(unname(fit$model$orders), as.integer(model[7:12]), info = id)
  expect_identical(fit$model$log, model[3L] == "true", info = id)
  expect_identical(fit$model$mean_correction, model[5L] == "true", info = id)
  params <- field("parameters")[-1L]
  exact(fit$model$parameters, params, id)
  exact(c(fit$model$phi, fit$model$bphi, fit$model$theta, fit$model$btheta), params, id)
  regression <- lines[startsWith(lines, "variable ")]
  mean <- tokens(regression[1L])
  expect_identical(mean[2L], "const", info = id)
  exact(fit$model$mean_coefficient, mean[3L], id)
  outliers <- fit$model$outliers
  expect_identical(nrow(outliers), length(regression) - 1L, info = id)
  for (j in seq_len(nrow(outliers))) {
    var <- tokens(regression[j + 1L])
    date <- as.integer(strsplit(gsub("[()]", "", var[3L]), "-", fixed = TRUE)[[1L]])
    expect_identical(outliers$name[j], paste(var[2:3], collapse = " "), info = id)
    expect_identical(outliers$type[j], var[2L], info = id)
    expect_identical(outliers$position[j], 1L + (date[1L] - start[1L]) * 12L + date[2L] - start[2L], info = id)
    exact(outliers$coefficient[j], var[4L], id)
  }
  variables <- variables + nrow(outliers)
  likelihood <- field("likelihood")
  exact(unlist(fit$likelihood[c("log_likelihood", "transformation_adjustment", "adjusted_log_likelihood", "ssq")]),
        likelihood[2:5], id)
  expect_identical(unname(unlist(fit$likelihood[c("effective", "estimated")])), as.integer(likelihood[6:7]), info = id)
  residuals <- field("residuals")
  tsres <- fit$residuals$tsres
  expect_identical(stats::frequency(tsres), 12, info = id)
  expect_identical(as.integer(stats::start(tsres)), as.integer(residuals[2:3]), info = id)
  expect_identical(length(tsres), as.integer(residuals[4L]), info = id)
  exact(tsres, residuals[-seq_len(4L)], id)
  residual_count <- residual_count + length(tsres)
}
expect_identical(c(length(cases), residual_count, variables), c(6L, 570L, 5L))
alias <- rkalos::native_tramo_preprocess(AirPassengers, "rsa1")
expect_identical(alias$template, "TR1")
expect_identical(as.integer(stats::start(alias$residuals$tsres)), c(1950L, 2L))
expect_identical(rkalos::capabilities()$tramo_preprocessing$templates, c("TR1", "TR3"))
expect_error(rkalos::native_tramo_preprocess(AirPassengers, "TR2"), "TR1")
expect_error(rkalos::native_tramo_preprocess(stats::ts(1:30, frequency = 4)), "monthly")
expect_error(rkalos::native_tramo_preprocess(stats::ts(1:23, frequency = 12)), "24")
expect_error(rkalos::native_tramo_preprocess(replace(AirPassengers, 2L, NA_real_)), "finite")
