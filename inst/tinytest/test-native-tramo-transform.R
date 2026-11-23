source("helpers/fixtures.R", local = TRUE)
fixture <- system.file("tinytest/fixtures/native/tramo-transform", package = "rkalos")
read_capture <- function(file) utils::read.csv(file.path(fixture, file), colClasses = "character")
cases <- read_capture("cases.csv")
inputs <- read_capture("inputs.csv")
facts <- read_capture("transform-facts.csv")

# The evidence uses round-trip binary64 strings, including the degenerate test's NaN.
numbers <- function(text) {
  tokens <- strsplit(text, " ", fixed = TRUE)[[1L]]
  values <- fixture_doubles(tokens)
  values[tokens == "NaN"] <- NaN
  values
}
for (id in cases$id) {
  fact <- with(facts[facts$id == id, ], stats::setNames(value, key))
  kind <- cases$kind[cases$id == id]
  if (kind == "calccv") {
    fields <- strsplit(fact[["calccv"]], " ", fixed = TRUE)[[1L]]
    expect_identical(rkalos::tramo_critical_value(as.integer(fields[2L])),
                     numbers(fields[4L]), info = id)
    next
  }
  records <- inputs$line[inputs$id == id]
  x <- stats::ts(numbers(sub("^values ", "", records[startsWith(records, "values ")])),
                 start = c(2001, 12), frequency = 12)
  if (kind == "season") {
    result <- rkalos::tramo_seasonality(x)
    expect_identical(paste(result$verdict, result$to_int), fact[["season_verdict"]], info = id)
    for (key in c(qs = "qs", friedman = "np")) {
      test <- if (key == "qs") result$qs else result$friedman
      expect_identical(test$value, numbers(fact[[paste0(key, "_val")]]), info = id)
      expect_identical(test$p_value, numbers(fact[[paste0(key, "_p")]]), info = id)
      expect_identical(test$valid, fact[[paste0(key, "_valid")]] == "true", info = id)
    }
    expect_identical(as.numeric(result$differenced),
                     numbers(sub("^n [0-9]+ ", "", fact[["diff"]])), info = id)
    expect_identical(stats::tsp(result$differenced), stats::tsp(diff(x)), info = id)
    next
  }
  fields <- strsplit(fact[["in"]], " ", fixed = TRUE)[[1L]]
  request <- stats::setNames(fields[seq(2, length(fields), 2)], fields[seq(1, length(fields), 2)])
  decide <- function() rkalos::tramo_transform_decision(
    x, numbers(request[["logpref"]]), numbers(request[["prec"]]),
    mean = request[["mean"]] == "1", detect_seasonality = request[["detect"]] == "1",
    spec_seasonality = request[["specseas"]] == "1", robust = request[["robust"]] == "1")
  if (request[["robust"]] == "1") {
    expect_error(decide(), "robust AO/LS", info = id)
    next
  }
  result <- decide()
  for (key in c("spec_pre", "loglevel_status")) {
    expect_identical(result[[key]], fact[[key]], info = id)
  }
  for (key in c("seasonal", "log_transform", "choosing_log")) {
    expect_identical(result[[key]], fact[[key]] == "1", info = id)
  }
  for (key in c("log_ll", "level_ll", "log_correction")) {
    expect_identical(result[[key]], numbers(fact[[key]]), info = id)
  }
  expect_identical(result$log_preference, numbers(fact[["log_pref"]]), info = id)
  expect_identical(result$nonpositive, any(x <= 0), info = id)
  if (request[["detect"]] == "1") {
    expect_identical(result$seasonality, rkalos::tramo_seasonality(x), info = id)
    expect_identical(paste(result$seasonality$verdict, result$seasonality$to_int),
                     fact[["season_verdict"]], info = id)
  } else {
    expect_null(result$seasonality, info = id)
  }
  if (id == "tr_neg") {
    expect_identical(result$remarks, c("No identifiable seasonality", "negative values, levels are chosen"))
  }
}

expect_identical(nrow(cases), 29L)
expect_true(all(c("tramo_critical_value", "tramo_seasonality", "tramo_transform_decision") %in%
                rkalos::capabilities()$primitives))
expect_false("tramoseats" %in% rkalos::capabilities()$methods)
expect_false(rkalos::capabilities()$fit)
for (n in list(-1L, 1.3, NA, 2^31, matrix(1))) expect_error(rkalos::tramo_critical_value(n))
for (x in list(1:24, stats::ts(1:23, frequency = 12), stats::ts(1:24, frequency = 4),
               stats::ts(matrix(1:48, ncol = 2), frequency = 12),
               stats::ts(c(1:23, Inf), frequency = 12))) expect_error(rkalos::tramo_seasonality(x))
expect_error(rkalos::tramo_transform_decision(AirPassengers, -0.05, 0), "precision")
expect_error(rkalos::tramo_transform_decision(AirPassengers, Inf, 1e-5), "finite numeric")
expect_error(rkalos::tramo_transform_decision(AirPassengers, -0.05, 1e-5, mean = 1), "logical")
expect_error(rkalos::tramo_transform_decision(AirPassengers, -0.05, 1e-5, robust = NA), "logical")
