source("helpers/fixtures.R", local = TRUE)
dir <- system.file("tinytest/fixtures/native", package = "rkalos", mustWork = TRUE)
oracle <- fixture_named_rsa0_cases(file.path(dir, "adjusted-seasonality", "public-x86_64.txt"))
tokens <- function(x) strsplit(x, " ", fixed = TRUE)[[1L]]
keys <- list(qs = "diagnostics.seas-sa-qs", seasonal_f = "diagnostics.seas-sa-f")
count <- numeric_count <- 0L
for (header in names(oracle)) {
  record <- tokens(header)
  id <- record[1L]
  lines <- oracle[[header]]
  series <- fixture_hex_doubles(tokens(lines[startsWith(lines, "model_scale_sa ")])[-1L])
  result <- monthly_adjusted_seasonality(series, as.integer(record[4L]), record[5L] == "true")
  expect_identical(names(result), c("qs", "seasonal_f"), info = id)
  for (name in names(keys)) {
    line <- tokens(lines[startsWith(lines, keys[[name]])][1L])
    test <- result[[name]]
    if (line[2L] == "unavailable") {
      expect_null(test, info = paste(id, name))
      next
    }
    expect_identical(names(test), c("value", "p_value"), info = paste(id, name))
    expected <- fixture_hex_doubles(line[-1L])
    actual <- unlist(test, use.names = FALSE)
    agree <- vapply(seq_along(expected), function(i)
      (is.nan(expected[i]) && is.nan(actual[i])) || actual[i] == expected[i], logical(1))
    expect_true(all(agree), info = paste(id, name))
    numeric_count <- numeric_count + length(actual)
  }
  count <- count + 1L
}
expect_identical(c(count, numeric_count), c(32L, 124L))
absent <- monthly_adjusted_seasonality(1, 1L)
expect_null(absent$qs)
expect_null(absent$seasonal_f)
expect_error(monthly_adjusted_seasonality(c(1, NaN), 1L))
expect_error(monthly_adjusted_seasonality(c(1, 2), -1))
expect_error(monthly_adjusted_seasonality(c(1, 2), 4L))
expect_error(monthly_adjusted_seasonality(c(1, 2), 1, mean = NA))
cap <- rkalos::capabilities()$adjusted_seasonality
expect_identical(cap$tests, c("seas-sa-qs", "seas-sa-f"))
expect_identical(cap$frequency, 12L)
expect_identical(cap$entry_point, "monthly_adjusted_seasonality")
