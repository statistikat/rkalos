source("helpers/fixtures.R", local = TRUE)
oracle <- system.file("tinytest/fixtures/native/tramo-seats-bridge/rsa0-x86_64.txt", package = "rkalos")
lines <- readLines(oracle)
rows <- split(lines, cumsum(startsWith(lines, "case ")))
points <- cases <- 0L
keys <- c(SeasonallyAdjusted = "finalSA", Trend = "finalTrend",
          Irregular = "finalIrregular", Seasonal = "finalSeasonal")
for (records in rows) {
  id <- sub("case ", "", records[1L], fixed = TRUE)
  tokens <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]
  input <- tokens(records[2L])
  x <- stats::ts(fixture_hex_doubles(input[-seq_len(4)]), start = as.integer(input[2:3]), frequency = 12)
  fit <- rkalos::native_tramoseats_fit(x)
  expect_identical(class(fit), c("native_tramoseats_fit", "list"), info = id)
  expect_identical(fit$mode, "additive", info = id)
  expect_identical(fit$model$orders, c(p = 0L, d = 1L, q = 1L, bp = 0L, bd = 1L, bq = 1L), info = id)
  expect_identical(fit$model$phi, numeric(), info = id)
  expect_identical(fit$model$bphi, numeric(), info = id)
  expect_identical(fit$status, "ok", info = id)
  expect_identical(fit$template, "RSA0", info = id)
  expect_identical(fit$method, "tramoseats", info = id)
  expect_identical(fit$series, x, info = id)
  fitted <- tokens(records[3L])
  actual <- c(fit$model$theta, fit$model$btheta, fit$model$mean_coefficient,
              fit$likelihood$log_likelihood, fit$likelihood$ssq)
  expect_identical(writeBin(actual, raw()), writeBin(fixture_hex_doubles(fitted[2:6]), raw()), info = id)
  expect_identical(c(fit$likelihood$effective, fit$likelihood$estimated), as.integer(fitted[7:8]), info = id)
  expect_identical(fit$model$innovation_variance,
    fit$likelihood$ssq / (fit$likelihood$effective - fit$likelihood$estimated), info = id)
  expect_identical(c(fit$model$d, fit$model$bd), c(1L, 1L), info = id)
  expect_false(fit$model$log, info = id)
  expect_true(fit$model$mean_correction, info = id)
  expected <- records[startsWith(records, "merged ")]
  expect_identical(nrow(fit$merged), length(expected), info = id)
  for (i in seq_along(expected)) {
    row <- tokens(expected[i])
    table <- fit$merged
    expect_identical(c(table$component[i], table$information[i]), row[2:3], info = id)
    value <- table$values[[i]]
    expect_identical(as.integer(stats::start(value)), as.integer(row[4:5]), info = id)
    expect_identical(stats::frequency(value), 12, info = id)
    expect_identical(length(value), as.integer(row[6L]), info = id)
    expect_identical(writeBin(as.numeric(value), raw()),
                     writeBin(fixture_hex_doubles(row[-seq_len(6)]), raw()), info = id)
    if (row[3L] == "Value" && row[2L] %in% names(keys)) {
      observed <- fit[[keys[[row[2L]]]]]
      expect_identical(as.numeric(observed), as.numeric(value), info = id)
      expect_identical(stats::tsp(observed), stats::tsp(x), info = id)
    }
    points <- points + length(value)
  }
  if (id == "fit_tr0_air") {
    wk <- readLines(file.path(dirname(oracle), "merge-x86_64.txt"))
    wk <- wk[seq_len(which(startsWith(wk, "case "))[2L] - 1L)]
    uncertainty <- wk[grepl("^wk [^ ]+ Stdev", wk)]
    table <- fit$final[startsWith(fit$final$information, "Stdev"), , drop = FALSE]
    expect_identical(nrow(table), length(uncertainty))
    nvalues <- 0L
    for (j in seq_along(uncertainty)) {
      row <- tokens(uncertainty[j])
      expected <- fixture_hex_doubles(row[-seq_len(6L)])
      value <- table$values[[j]]
      expect_identical(c(table$component[j], table$information[j]), row[2:3])
      expect_identical(as.integer(stats::start(value)), as.integer(row[4:5]))
      expect_identical(stats::frequency(value), 12)
      expect_identical(writeBin(as.numeric(value), raw()), writeBin(expected, raw()))
      nvalues <- nvalues + length(value)
    }
    expect_identical(c(nrow(table), nvalues), c(9L, 696L))
  }
  cases <- cases + 1L
}
expect_identical(c(cases, points), c(6L, 3960L))
x <- stats::ts(as.numeric(AirPassengers), start = c(2010, 3), frequency = 12)
# Some ts producers preserve an equivalent month endpoint in different low bits.
attr(x, "tsp")[2L] <- attr(x, "tsp")[2L] + 2^-41
fit <- rkalos::native_tramoseats_fit(x, template = "tr0", forecast = 6L, backcast = 3L)
expect_identical(fit$status, "ok")
for (key in unname(keys)) expect_identical(stats::tsp(fit[[key]]), stats::tsp(x))
expect_identical(as.integer(stats::start(fit$merged$values[[which(fit$merged$component == "Series" &
  fit$merged$information == "Backcast")]])), c(2009L, 12L))
expect_identical(as.integer(stats::start(fit$merged$values[[which(fit$merged$component == "Series" &
  fit$merged$information == "Forecast")]])), c(2022L, 3L))
zero <- rkalos::native_tramoseats_fit(x, forecast = 0L)
expect_identical(unique(zero$merged$information), "Value")
expect_error(rkalos::native_tramoseats_fit(x, "RSA3"), "RSA0")
expect_error(rkalos::native_tramoseats_fit(stats::ts(1:30, frequency = 4)), "monthly")
expect_error(rkalos::native_tramoseats_fit(replace(x, 6L, NA_real_)), "finite")
expect_error(rkalos::native_tramoseats_fit(stats::ts(1:8, frequency = 12)), "degrees of freedom")
expect_error(rkalos::native_tramoseats_fit(x, forecast = 1.5), "integer")
expect_error(rkalos::native_tramoseats_fit(x, backcast = -.Machine$integer.max), "range")
expect_identical(rkalos::capabilities()$tramo_seats_fit$templates, c("RSA0", "RSA1"))
