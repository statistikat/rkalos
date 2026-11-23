source("helpers/fixtures.R", local = TRUE)
dir <- system.file("tinytest/fixtures/native/tramo-model-builder", package = "rkalos")
read <- function(file) utils::read.csv(file.path(dir, file), colClasses = "character")
cases <- read("cases.csv")
inputs <- read("inputs.csv")
facts <- read("builder-facts.csv")
ids <- c(paste0("tr", 0:5), "trfull", "missing_lead", "missing_interior", "missing_trail", "feb_start")
for (id in ids) {
  ca <- cases[cases$id == id, ]
  request <- inputs$line[inputs$id == id]
  template <- sub("req template ", "", request[startsWith(request, "req template ")])
  tokens <- sub("values ", "", request[startsWith(request, "values ")])
  tokens <- strsplit(tokens, ";", fixed = TRUE)[[1]]
  tokens[tokens == "nan"] <- NA_character_
  values <- fixture_doubles(tokens)
  start <- as.integer(strsplit(ca$start, " ", fixed = TRUE)[[1]])
  x <- stats::ts(values, start = start, frequency = as.integer(ca$freq))
  result <- rkalos::tramo_prepare_model(x, template)
  expect_identical(result, rkalos::tramo_prepare_model(x, tolower(template)), info = id)
  expect_identical(result$template, template, info = id)
  expect_identical(result$automatic_model, ca$auto == "y", info = id)
  expect_identical(result$log, ca$log == "y", info = id)
  expect_identical(result$mean, ca$mean == "y", info = id)
  orders <- c("p", "d", "q", "bp", "bd", "bq")
  expect_identical(result$orders, stats::setNames(as.integer(unlist(ca[orders])), orders), info = id)
  expect_identical(result$period, as.integer(ca$period), info = id)
  expect_identical(length(result$variables), as.integer(ca$nvar), info = id)
  expect_identical(result$warnings, character(), info = id)
  get <- function(key) facts$value[facts$id == id & facts$key == key]
  missing <- if (get("missing") == "none") integer() else
    1L + as.integer(strsplit(get("missing"), ";", fixed=TRUE)[[1]])
  expect_identical(result$missing, missing, info = id)
  if (length(missing)) {
    replacements <- strsplit(get("interp"), ";", fixed=TRUE)[[1]]
    values[missing] <- fixture_doubles(sub("^[0-9]+=", "", replacements))
  }
  expect_identical(stats::tsp(result$interpolated), stats::tsp(x), info = id)
  expect_identical(writeBin(as.numeric(result$interpolated), raw()), writeBin(values, raw()), info = id)
}
expect_error(rkalos::tramo_prepare_model(AirPassengers, "TRFULL"), "Unknown TRAMO template")
expect_error(rkalos::tramo_prepare_model(AirPassengers, "rsa3"), "Unknown TRAMO template")
expect_error(rkalos::tramo_prepare_model(stats::ts(rep(NA_real_,24), frequency=12), "tr0"), "all-missing")
expect_error(rkalos::tramo_prepare_model(stats::ts(c(1:23,Inf),frequency=12), "tr0"), "infinite")
expect_error(rkalos::tramo_prepare_model(stats::ts(1:24,frequency=4), "tr0"), "D11")
expect_true("tramo_prepare_model" %in% rkalos::capabilities()$primitives)
expect_false("tramoseats" %in% rkalos::capabilities()$methods)
