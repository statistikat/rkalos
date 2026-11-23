source("helpers/fixtures.R", local = TRUE)
dir <- system.file("tinytest/fixtures/native/tramo-seats-bridge", package = "rkalos")
cases <- function(file) {
  lines <- readLines(file.path(dir, file))
  split(lines, cumsum(startsWith(lines, "case ")))
}
fitted <- list()
for (lines in cases("bridge.txt")) {
  fields <- function(key) strsplit(lines[startsWith(lines, paste0(key, " "))], " ", fixed = TRUE)[[1L]]
  domain <- fields("interpolated")
  start <- as.integer(domain[2:3])
  model <- fields("model")
  likelihood <- fields("likelihood")
  request <- list(series = stats::ts(fixture_hex_doubles(domain[-seq_len(4)]), start = start, frequency = 12),
    ssq = fixture_hex_doubles(likelihood[3L]), effective = as.integer(likelihood[5L]),
    estimated = as.integer(likelihood[7L]), d = as.integer(model[9L]), bd = as.integer(model[11L]),
    mean_correction = model[5L] == "true", log = model[3L] == "true")
  for (name in c("phi", "bphi", "theta", "btheta")) {
    params <- fields(name)[-c(1L, 2L)]
    request[[name]] <- if (length(params)) fixture_hex_doubles(sub("/Fixed$", "", params)) else numeric()
  }
  for (line in lines[startsWith(lines, "variable AO ") | startsWith(lines, "variable LS ")]) {
    tokens <- strsplit(line, " ", fixed = TRUE)[[1L]]
    date <- as.integer(strsplit(gsub("[()]", "", tokens[3L]), "-", fixed = TRUE)[[1L]])
    row <- data.frame(type = tokens[2L], position = 1L + (date[1L] - start[1L]) * 12L + date[2L] - start[2L],
                      coefficient = fixture_hex_doubles(tail(tokens, 1L)))
    request$outliers <- rbind(request$outliers, row)
  }
  fitted[[fields("case")[2L]]] <- request
}
accepted <- rows <- points <- 0L
for (lines in cases("merge-x86_64.txt")) {
  header <- strsplit(lines[1L], " ", fixed = TRUE)[[1L]]
  id <- header[2L]
  request <- fitted[[id]]
  request$forecast <- as.integer(header[7L])
  request$backcast <- as.integer(header[5L])
  result <- do.call(rkalos::seats_from_tramo, request)
  expect_identical(result$status, "ok", info = id)
  for (phase in c("wk", "merged")) {
    table <- if (phase == "wk") result$final else result$merged
    oracle <- lines[startsWith(lines, paste0(phase, " "))]
    expect_identical(nrow(table), length(oracle), info = id)
    rows <- rows + nrow(table)
    for (i in seq_along(oracle)) {
      row <- strsplit(oracle[i], " ", fixed = TRUE)[[1L]]
      expect_identical(c(table$component[i], table$information[i]), row[2:3], info = id)
      expect_identical(as.integer(c(table$start_year[i], table$start_month[i])), as.integer(row[4:5]), info = id)
      values <- table$values[[i]]
      expect_identical(length(values), as.integer(row[6L]), info = id)
      if (length(values)) {
        expect_identical(as.integer(stats::start(values)), as.integer(row[4:5]), info = id)
        expect_identical(stats::frequency(values), 12, info = id)
        expected <- fixture_hex_doubles(row[-seq_len(6)])
        expect_identical(writeBin(as.numeric(values), raw()), writeBin(expected, raw()), info = id)
        points <- points + length(values)
      }
    }
  }
  accepted <- accepted + 1L
}
expect_identical(c(accepted, rows, points), c(6L, 189L, 11343L))
# The seasonal AR-to-MA refit is qualified only without a mean term.
expect_error(do.call(rkalos::seats_from_tramo, within(fitted$fit_multi, mean_correction <- TRUE)),
             "outside the qualified monthly retry envelope")
x <- fitted[[1L]]
expect_error(do.call(rkalos::seats_from_tramo, within(x, effective <- estimated)), "likelihood")
expect_error(do.call(rkalos::seats_from_tramo, within(x, ssq <- 0)), "likelihood")
bad <- data.frame(type = "AO", position = 0L, coefficient = 1)
x$outliers <- bad
expect_error(do.call(rkalos::seats_from_tramo, x), "position")
x$outliers$position <- 1L
expect_identical(do.call(rkalos::seats_from_tramo, x)$status, "ok")
x$outliers$type <- "TC"
expect_error(do.call(rkalos::seats_from_tramo, x), "AO or LS")
expect_true("seats_from_tramo" %in% rkalos::capabilities()$primitives)
