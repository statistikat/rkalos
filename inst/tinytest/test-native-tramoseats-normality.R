source("helpers/fixtures.R", local = TRUE)
dir <- system.file("tinytest/fixtures/native/tramo-seats-bridge", package = "rkalos")
inputs <- fixture_named_rsa0_cases(file.path(dir, "rsa0-x86_64.txt"))
oracle <- fixture_named_rsa0_cases(file.path(dir, "rsa0-doornik-hansen-x86_64.txt"))
check <- function(actual, tokens, id) {
  expect_identical(names(actual), c("value", "pvalue", "df"), info = id)
  expect_identical(actual$df, 2, info = id)
  expect_identical(writeBin(c(actual$value, actual$pvalue), raw()),
                   writeBin(fixture_hex_doubles(tokens), raw()), info = id)
}
fields <- function(line) strsplit(line, " ", fixed = TRUE)[[1L]]
cases <- edges <- 0L
for (id in names(oracle)) {
  input <- fields(inputs[[id]][2L])
  series <- stats::ts(fixture_hex_doubles(input[-seq_len(4L)]), frequency = 12, start = as.integer(input[2:3]))
  fit <- rkalos::native_tramoseats_fit(series)
  expected <- fields(oracle[[id]][2L])
  check(fit$residuals$doornikhansen, expected[2:3], id)
  cases <- cases + 1L
  for (line in oracle[[id]][startsWith(oracle[[id]], "short ")]) {
    edge <- fields(line)
    n <- as.integer(edge[2L])
    short <- stats::ts(as.numeric(series)[seq_len(n)], frequency = 12, start = stats::start(series))
    residuals <- rkalos::native_tramoseats_fit(short)$residuals
    expect_identical(length(residuals$tsres), as.integer(edge[3L]), info = id)
    check(residuals$doornikhansen, edge[4:5], paste(id, n))
    edges <- edges + 1L
  }
}
expect_identical(c(cases, edges), c(6L, 3L))
expect_true("doornikhansen" %in% rkalos::capabilities()$tramo_seats_fit$residuals)
