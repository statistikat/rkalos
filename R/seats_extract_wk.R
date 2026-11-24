seats_extract_wk <- function(series, innovation_variance, d = 0L, bd = 0L,
    phi = numeric(), theta = numeric(), bphi = numeric(), btheta = numeric(),
    log = FALSE, mean_correction = TRUE, forecast = -2L, backcast = 0L,
    approximation = "legacy", method = "burman", bias = "legacy", xl = 0.95,
    epsphi = 2, rmod = 0.5, smod = 0.8, smod_pi = 0.8) {
  spec <- as.list(environment())[setdiff(names(formals(seats_extract_wk)), "series")]
  domain <- seats_model_domain(series, spec)
  for (key in c("innovation_variance", "xl", "epsphi", "rmod", "smod", "smod_pi")) {
    x <- spec[[key]]
    if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || length(x) != 1L || !is.finite(x))
      stop(key, " must be a finite number.", call. = FALSE)
  }
  if (innovation_variance < 0) stop("Innovation variance must be nonnegative.", call. = FALSE)
  choices <- list(approximation = c("none", "legacy", "noisy"), method = c("burman", "kalman_smoother"),
                  bias = c("none", "legacy"))
  for (key in names(choices)) {
    x <- spec[[key]]
    if (!is.character(x) || !is.null(dim(x)) || length(x) != 1L || !x %in% choices[[key]])
      stop(key, " must be one of ", paste(choices[[key]], collapse = ", "), ".", call. = FALSE)
  }
  result <- .native_seats_extract_wk(spec, domain, as.numeric(series))
  seats_date_components(result, domain$frequency, c("initial", "final"))
}

seats_date_components <- function(result, frequency, phases) {
  for (phase in phases) {
    table <- result[[phase]]
    for (i in seq_len(nrow(table))) {
      if (length(table$values[[i]])) {
        table$values[[i]] <- stats::ts(table$values[[i]],
          start = c(table$start_year[i], table$start_month[i]), frequency = frequency)
      }
    }
    result[[phase]] <- table
  }
  result
}

seats_model_domain <- function(series, spec) {
  domain <- native_regular_domain(series)
  if (any(!is.finite(series))) stop("series must be complete and finite.", call. = FALSE)
  if (abs(domain$start[1]) > .Machine$integer.max)
    stop("SEATS start year is outside the native integer range.", call. = FALSE)
  for (key in c("phi", "theta", "bphi", "btheta")) {
    x <- spec[[key]]
    if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || any(!is.finite(x)))
      stop(key, " must be a finite numeric vector.", call. = FALSE)
  }
  for (key in c("d", "bd")) {
    x <- spec[[key]]
    if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || length(x) != 1L || !is.finite(x) ||
        x != trunc(x) || abs(x) > .Machine$integer.max)
      stop(key, " must be an integer in the native index range.", call. = FALSE)
  }
  if (spec$d < 0 || spec$bd < 0 || spec$d + as.double(spec$bd) * domain$frequency >= length(series))
    stop("Differencing must leave observed data.", call. = FALSE)
  seats_validate_horizons(spec$backcast, spec$forecast, domain$frequency, length(series))
  for (key in c("log", "mean_correction")) native_regression_flag(spec[[key]], key)
  if (spec$log && any(series <= 0))
    stop("SEATS log transformation requires positive observations.", call. = FALSE)
  domain
}

seats_validate_horizons <- function(backcast, forecast, frequency, observations) {
  for (x in list(backcast, forecast)) {
    if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || length(x) != 1L || !is.finite(x) ||
        x != trunc(x) || abs(x) > .Machine$integer.max)
      stop("SEATS horizon must be an integer in the native index range.", call. = FALSE)
  }
  counts <- c(as.double(backcast), as.double(forecast))
  counts[counts < 0] <- -counts[counts < 0] * frequency
  if (observations + sum(counts) > .Machine$integer.max)
    stop("The extended SEATS series exceeds the native index range.", call. = FALSE)
  as.integer(counts)
}
