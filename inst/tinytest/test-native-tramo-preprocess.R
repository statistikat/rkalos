source("helpers/fixtures.R", local = TRUE)
dir <- system.file("tinytest/fixtures/native/tramo-seats-bridge", package = "rkalos")
files <- c("rsa3-x86_64.txt", "rsa3-additional-x86_64.txt", "rsa3-mean-retry-x86_64.txt", "rsa3-ref-retry-x86_64.txt")
cases <- list(TR1 = fixture_named_rsa0_cases(file.path(dir, "rsa1-x86_64.txt")),
              TR3 = do.call(c, lapply(files, function(file) fixture_named_rsa0_cases(file.path(dir, file)))))
diagnostic_files <- c("rsa3-residual-diagnostics-x86_64.txt", "rsa3-mean-retry-diagnostics-x86_64.txt",
                      "rsa3-ref-retry-diagnostics-x86_64.txt")
diagnostics <- do.call(c, lapply(diagnostic_files, function(file) fixture_named_rsa0_cases(file.path(dir, file))))
expect_identical(sort(names(cases$TR3)), sort(names(diagnostics)))
tokens <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]
exact <- function(actual, expected, id) {
  expect_identical(writeBin(as.numeric(actual), raw()), writeBin(fixture_hex_doubles(expected), raw()), info = id)
}
residual_count <- variables <- 0L
for (template in names(cases)) for (id in names(cases[[template]])) {
  lines <- cases[[template]][[id]]
  field <- function(key) tokens(lines[startsWith(lines, paste0(key, " "))])
  input <- field("input")
  start <- as.integer(input[2:3])
  series <- stats::ts(fixture_hex_doubles(input[-seq_len(4L)]), start = start, frequency = 12)
  fit <- rkalos::native_tramo_preprocess(series, template)
  expect_identical(class(fit), c("native_tramo_preprocess", "list"), info = id)
  expect_identical(fit$status, "available", info = id)
  expect_identical(fit$template, template, info = id)
  model <- field("model")
  expect_identical(unname(fit$model$orders), as.integer(model[7:12]), info = id)
  expect_identical(fit$model$log, model[3L] == "true", info = id)
  expect_identical(fit$model$mean_correction, model[5L] == "true", info = id)
  params <- field("parameters")[-1L]
  exact(fit$model$parameters, params, id)
  exact(c(fit$model$phi, fit$model$bphi, fit$model$theta, fit$model$btheta), params, id)
  constant <- lines[startsWith(lines, "variable const ")]
  if (length(constant)) exact(fit$model$mean_coefficient, tokens(constant)[3L], id)
  else expect_identical(fit$model$mean_coefficient, NA_real_, info = id)
  regression <- lines[startsWith(lines, "variable ") & !startsWith(lines, "variable const ")]
  outliers <- fit$model$outliers
  expect_identical(nrow(outliers), length(regression), info = id)
  for (j in seq_len(nrow(outliers))) {
    var <- tokens(regression[j])
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
  if (template == "TR3") {
    oracle <- diagnostics[[id]]
    test <- function(key) tokens(oracle[startsWith(oracle, paste0("residuals.", key, " "))])
    exact(fit$residuals$ser, test("ser")[2L], id)
    for (key in c("lb", "doornikhansen")) {
      result <- fit$residuals[[key]]
      exact(c(result$value, result$pvalue), test(key)[2:3], paste(id, key))
    }
    expect_identical(fit$residuals$lb$df, 24 - length(fit$model$parameters), info = id)
  }
  residual_count <- residual_count + length(tsres)
}
expect_identical(c(sum(lengths(cases)), residual_count, variables), c(26L, 2675L, 24L))
alias <- rkalos::native_tramo_preprocess(AirPassengers, "rsa1")
expect_identical(alias$template, "TR1")
expect_identical(as.integer(stats::start(alias$residuals$tsres)), c(1950L, 2L))
expect_identical(rkalos::capabilities()$tramo_preprocessing$templates, c("TR1", "TR3", "TRfull"))
expect_error(rkalos::native_tramo_preprocess(AirPassengers, "TR2"), "TR1")
expect_error(rkalos::native_tramo_preprocess(stats::ts(1:30, frequency = 4)), "monthly")
expect_error(rkalos::native_tramo_preprocess(stats::ts(1:23, frequency = 12)), "24")
expect_error(rkalos::native_tramo_preprocess(replace(AirPassengers, 2L, NA_real_)), "finite")
