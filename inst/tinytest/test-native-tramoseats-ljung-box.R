source("helpers/fixtures.R", local = TRUE)
dir <- system.file("tinytest/fixtures/native/tramo-seats-bridge", package = "rkalos")
inputs <- fixture_named_rsa0_cases(file.path(dir, "rsa0-x86_64.txt"))
oracle <- fixture_named_rsa0_cases(file.path(dir, "rsa0-ljung-box-x86_64.txt"))
cases <- edges <- 0L
check <- function(actual, tokens, id) {
  expected <- fixture_hex_doubles(tokens)
  expect_identical(actual$df, 22, info = id)
  expect_identical(actual$lags, 24L, info = id)
  for (i in seq_along(expected)) {
    got <- actual[[c("value", "pvalue")[i]]]
    if (is.nan(expected[i])) {
      expect_true(is.nan(got), info = id)
    } else {
      expect_identical(writeBin(got, raw()), writeBin(expected[i], raw()), info = id)
    }
  }
}
for (id in names(oracle)) {
  fields <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]
  input <- fields(inputs[[id]][2L])
  x <- stats::ts(fixture_hex_doubles(input[-seq_len(4L)]), frequency = 12, start = as.integer(input[2:3]))
  fit <- rkalos::native_tramoseats_fit(x)
  want <- fields(oracle[[id]][2L])
  check(fit$residuals$lb, want[2:3], id)
  cases <- cases + 1L
  for (line in oracle[[id]][startsWith(oracle[[id]], "short ")]) {
    edge <- fields(line)
    n <- as.integer(edge[2L])
    short <- stats::ts(as.numeric(x)[seq_len(n)], frequency = 12, start = stats::start(x))
    result <- rkalos::native_tramoseats_fit(short)$residuals
    expect_identical(length(result$tsres), as.integer(edge[3L]), info = id)
    check(result$lb, edge[4:5], paste(id, n))
    edges <- edges + 1L
  }
}
expect_identical(c(cases, edges), c(6L, 3L))
expect_true("lb" %in% rkalos::capabilities()$tramo_seats_fit$residuals)
