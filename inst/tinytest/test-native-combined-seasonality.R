source("helpers/fixtures.R", local = TRUE)
dir <- system.file("tinytest/fixtures/native", package = "rkalos", mustWork = TRUE)
read_panel <- function(directory, id) fixture_named_rsa0_cases(file.path(dir, directory, paste0(id, "-x86_64.txt")))
oracle <- c(read_panel("combined-seasonality", "public"), read_panel("combined-seasonality", "fixed-x13"))
input <- c(read_panel("tramo-seats-bridge", "rsa3"), read_panel("tramo-seats-bridge", "rsa3-ref-retry"))
tokens <- function(x) strsplit(x, " ", fixed = TRUE)[[1L]]
count <- numeric_count <- 0L
for (header in names(oracle)) {
  record <- tokens(header)
  id <- record[1L]
  kind <- record[2L]
  if (!startsWith(kind, "wk_") && kind != "fixed_d8") next
  lines <- oracle[[header]]
  input_id <- substring(id, regexpr("fit_", id, fixed = TRUE))
  fields <- tokens(input[[input_id]][startsWith(input[[input_id]], "input ")])
  series <- stats::ts(fixture_hex_doubles(fields[-seq_len(4L)]), start = as.integer(fields[2:3]), frequency = 12)
  if (kind == "fixed_d8") {
    mode <- c(Additive = "additive", Multiplicative = "multiplicative", LogAdditive = "log_additive",
              PseudoAdditive = "pseudo_additive")[[strsplit(id, "_", fixed = TRUE)[[1L]][2L]]]
    fit <- native_x13_fit(series, log = record[3L] == "true", seasonal = rep("s3x5", 12L),
                          trend_length = 13L, mode = mode)
    expect_identical(fit$mode_effective, mode, info = id)
  } else {
    fit <- native_tramoseats_fit(series, sub("_.*", "", id))
  }
  test <- fit$diagnostics$combined_seasonality
  summary <- tokens(lines[startsWith(lines, "summary ")])[2L]
  expect_identical(test$summary, summary, info = id)
  expect_identical(names(test), c("summary", "stable", "moving", "kruskal_wallis"), info = id)
  expected <- fixture_hex_doubles(tokens(lines[startsWith(lines, "statistics ")])[-1L])
  actual <- unlist(test[-1L], use.names = FALSE)
  expect_true(all(abs(actual - expected) <= 1e-12 * (1 + abs(expected))), info = id)
  count <- count + 1L
  numeric_count <- numeric_count + length(actual)
}
expect_identical(c(count, numeric_count), c(24L, 144L))
short <- stats::ts(as.numeric(AirPassengers)[1:35], start = c(1980, 1), frequency = 12)
expect_null(native_tramoseats_fit(short, "RSA0")$diagnostics$combined_seasonality)
cap <- rkalos::capabilities()$combined_seasonality
expect_identical(cap$getter, "diagnostics.seas-si-combined")
expect_identical(cap$tramo_seats_fit, c("RSA0", "RSA1", "RSA3", "RSAfull"))
expect_true(cap$restricted_x13_fit)
expect_true(cap$automatic_x13_fit)
expect_identical(cap$estimated_x13_fit, c("RSA0", "RSA1", "RSA2c"))

stock <- c(read_panel("combined-seasonality", "stock"), read_panel("combined-seasonality", "stock-rsa2c"))
full_inputs <- c(read_panel("tramo-seats-bridge", "rsafull"), read_panel("tramo-seats-bridge", "rsafull-late-calendar"))
expect_identical(length(stock), 38L)
for (header in names(stock)) {
  id <- tokens(header)[1L]
  if (startsWith(id, "x13_")) {
    key <- strsplit(id, "_", fixed = TRUE)[[1L]]
    rows <- readLines(file.path(dir, "x13-refusal-jvm-profiles", paste0(key[3L], ".csv")))
    series <- stats::ts(fixture_doubles(rows[-1L]), start = as.integer(strsplit(rows[1L], ",")[[1L]]), frequency = 12)
    fit <- native_x13_estimated_fit(series, key[2L])
  } else {
    lines <- full_inputs[[sub("^rsafull_", "", id)]]
    field <- tokens(lines[startsWith(lines, "input ")])
    series <- stats::ts(fixture_hex_doubles(field[-seq_len(4L)]), start = as.integer(field[2:3]), frequency = 12)
    fit <- native_tramoseats_fit(series, "RSAfull")
  }
  lines <- stock[[header]]
  test <- fit$diagnostics$combined_seasonality
  reference <- fixture_hex_doubles(tokens(lines[startsWith(lines, "statistics ")])[-1L])
  actual <- unlist(test[-1L], use.names = FALSE)
  expect_identical(test$summary, tokens(lines[startsWith(lines, "summary ")])[2L], info = id)
  expect_identical(length(actual), length(reference), info = id)
  expect_true(all(abs(actual - reference) <= 1e-12 * (1 + abs(reference))), info = id)
}
