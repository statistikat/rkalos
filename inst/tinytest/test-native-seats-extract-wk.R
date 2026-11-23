source("helpers/fixtures.R", local = TRUE)
dir <- system.file("tinytest/fixtures/native/seats-wk-extraction", package = "rkalos")
read <- function(file) utils::read.csv(file.path(dir, file), colClasses = "character")
cases <- read("cases.csv")
inputs <- read("inputs.csv")
facts <- read("model-facts.csv")
refusals <- c(x_airline = "approximation_required", x_crash_fixma = "approximation_required",
              x_quarterly = "D11", x_var0 = "innovation variance 0",
              x_kalman = "KalmanSmoother", x_kalman_wn = "KalmanSmoother")
count <- points <- 0L
for (i in seq_len(nrow(cases))) {
  ca <- cases[i, ]
  id <- ca$id
  records <- inputs$line[inputs$id == id]
  fields <- function(key) {
    prefix <- paste0(key, " ")
    line <- records[startsWith(records, prefix)]
    stopifnot(length(line) == 1L)
    strsplit(substring(line, nchar(prefix)+1L), " ", fixed = TRUE)[[1]]
  }
  order <- as.integer(fields("req spec"))
  data <- fields("req log")
  knobs <- fields("req knobs")
  series <- stats::ts(fixture_doubles(fields("values")), frequency = as.integer(ca$freq),
                      start = as.integer(strsplit(ca$start, " ", fixed=TRUE)[[1]]))
  invoke <- function() rkalos::seats_extract_wk(series,
    innovation_variance = fixture_doubles(data[5]), d = order[2], bd = order[5],
    phi = fixture_doubles(fields("req phi")), theta = fixture_doubles(fields("req theta")),
    bphi = fixture_doubles(fields("req bphi")), btheta = fixture_doubles(fields("req btheta")),
    log = data[1] == "1", mean_correction = data[3] == "1",
    xl = fixture_doubles(knobs[1]), epsphi = fixture_doubles(knobs[2]), rmod = fixture_doubles(knobs[3]),
    smod = fixture_doubles(knobs[4]), smod_pi = fixture_doubles(knobs[5]),
    approximation = tolower(knobs[6]), forecast = as.integer(knobs[7]), backcast = as.integer(knobs[8]),
    method = if (knobs[9] == "Burman") "burman" else "kalman_smoother", bias = tolower(knobs[10]))
  if (id %in% names(refusals)) {
    expect_error(invoke(), refusals[[id]], info = id)
    next
  }
  result <- invoke()
  expected <- facts[facts$id == id, ]
  expect_identical(result$status, expected$value[expected$key == "status"], info = id)
  for (phase in c("initial", "final")) {
    table <- result[[phase]]
    key <- if (phase == "initial") "init" else "final"
    rows <- if (result$status == "ok") expected$value[expected$key == key] else character()
    expect_identical(nrow(table), length(rows), info = paste(id, phase))
    for (j in seq_along(rows)) {
      tokens <- strsplit(rows[j], " ", fixed=TRUE)[[1]]
      n <- as.integer(tokens[4])
      start <- as.integer(tokens[6:7]) + c(0L, 1L)
      values <- table$values[[j]]
      expect_identical(table$component[j], tokens[1], info = id)
      expect_identical(table$information[j], tokens[2], info = id)
      expect_identical(as.integer(c(table$start_year[j], table$start_month[j])), start, info = id)
      expect_identical(length(values), n, info = id)
      if (n > 0L) {
        expect_identical(as.integer(stats::start(values)), start, info = id)
        expect_identical(stats::frequency(values), 12, info = id)
      }
      want <- fixture_doubles(tokens[-seq_len(8)])
      # Compare the actual finite bit patterns, including signed zero.
      expect_identical(writeBin(as.numeric(values), raw()), writeBin(want, raw()), info = paste(id, key, j))
      count <- count + 1L
      points <- points + n
    }
  }
}
expect_identical(nrow(cases), 36L)
expect_identical(count, 1055L)
expect_identical(points, 55812L)
expect_true("seats_extract_wk" %in% rkalos::capabilities()$primitives)
expect_false(rkalos::capabilities()$fit)
expect_error(rkalos::seats_extract_wk(AirPassengers, 1, forecast = -.Machine$integer.max), "index range")
expect_error(rkalos::seats_extract_wk(AirPassengers, 1, phi = NA_real_), "finite numeric vector")
expect_error(rkalos::seats_extract_wk(AirPassengers, 1, d = 0.5), "integer")
expect_error(rkalos::seats_extract_wk(c(AirPassengers), 1), "univariate numeric ts")
expect_error(rkalos::seats_extract_wk(-AirPassengers, 1, log = TRUE), "positive observations")
