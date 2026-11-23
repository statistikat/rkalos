source("helpers/fixtures.R", local = TRUE)
dir <- system.file("tinytest/fixtures/native/tramo-differencing", package = "rkalos")
read <- function(file) utils::read.csv(file.path(dir, file), colClasses = "character")
cases <- read("cases.csv")
inputs <- read("inputs.csv")
facts <- read("diff-facts.csv")
events <- read("diff-events.csv")
orders <- c("p", "d", "q", "bp", "bd", "bq")
for (id in cases$id[cases$mode %in% c("DIFF", "BAT")]) {
  lines <- inputs$line[inputs$id == id]
  header <- strsplit(lines[1L], " ", fixed = TRUE)[[1L]]
  battery <- header[1L] == "BAT"
  offset <- as.integer(battery)
  values <- fixture_doubles(strsplit(sub("^values ", "", lines[2L]), " ", fixed = TRUE)[[1L]])
  series <- stats::ts(values, frequency = as.integer(header[3L]))
  entry <- stats::setNames(as.integer(header[6:11 + offset]), c("d", "bd", "p", "q", "bp", "bq"))
  fit <- rkalos::tramo_differencing_decision(
    series, entry, mean = header[12L + offset] == "y", log = header[5L + offset] == "y",
    seasonal = header[4L] == "y", initial_round = header[17L + 2L * offset] == "y")
  fact <- function(key) strsplit(facts$value[facts$id == id & facts$key == key], " ", fixed = TRUE)[[1L]]
  if (!battery) {
    after <- fact("after")
    expected <- as.integer(strsplit(sub("@.*$", "", after[3L]), ",", fixed = TRUE)[[1L]])
    expect_identical(fit$orders, stats::setNames(expected, orders), info = id)
    expect_identical(fit$mean, after[5L] == "1", info = id)
  }
  expect_identical(fit$decision, fact("diff_result")[2L], info = id)
  info <- fact("diff_info")
  expect_identical(fit$mean_correction, info[7L] == "1", info = id)
  expect_identical(writeBin(fit$tmean, raw()), writeBin(fixture_doubles(info[9L]), raw()), info = id)
  stage <- events[events$id == id, ]
  stage <- stage[seq_len(match("pop", stage$kind)), ]
  expected_events <- stage$kind
  expected_events[stage$kind == "push"] <- paste("push", stage$msg[stage$kind == "push"])
  expected_events[stage$kind == "info"] <- paste("info", sub(" .*", "", stage$detail[stage$kind == "info"]))
  expect_identical(fit$events, expected_events, info = id)
  expect_false(any(fit$default_model, fit$selection_failed, fit$airline_applied), info = id)
}
short <- stats::ts(1:4, frequency = 12)
fallback <- rkalos::tramo_differencing_decision(short)
expect_true(fallback$default_model)
expect_true(fallback$airline_applied)
expect_false(fallback$selection_failed)
expect_identical(fallback$decision, "Changed")
expect_identical(fallback$orders, stats::setNames(c(0L, 1L, 1L, 0L, 1L, 1L), orders))
expect_identical(rkalos::tramo_differencing_decision(short, fallback$orders)$decision, "Unprocessed")
expect_error(rkalos::tramo_differencing_decision(short, unname(fallback$orders)), "orders")
expect_error(rkalos::tramo_differencing_decision(short, replace(fallback$orders, "bd", 2L)), "orders")
expect_error(rkalos::tramo_differencing_decision(-short, log = TRUE), "positive")
expect_error(rkalos::tramo_differencing_decision(replace(short, 1L, NA_real_)), "complete finite")
expect_error(rkalos::tramo_differencing_decision(short, seasonal = NA), "logical")
expect_true("tramo_differencing_decision" %in% rkalos::capabilities()$primitives)
