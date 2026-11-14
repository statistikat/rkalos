# Serialization compatibility from Persephone 8450d6e R/forecasts.R.
# Public forecast values use a backend-independent, serializable contract.
native_forecast_component_names <- function() c("y", "t", "sa", "s", "i", "cal")

native_forecast_series <- function(value, component, observed = NULL, error = FALSE) {
  if (is.null(value) || !length(value)) return(NULL)
  if (!stats::is.ts(value) || !is.numeric(value) || !is.null(dim(value))) {
    stop("Forecast '", component, "' must be a univariate numeric ts.", call. = FALSE)
  }
  if (any(is.infinite(value)) || (error && any(value < 0, na.rm = TRUE))) {
    stop("Forecast '", component, "' has invalid ",
         if (error) "standard errors." else "nonfinite values.", call. = FALSE)
  }
  if (!is.null(observed)) {
    if (stats::frequency(value) != stats::frequency(observed)) {
      stop("Forecast frequency must match the observed series.", call. = FALSE)
    }
    expected <- stats::tsp(observed)[2] + 1 / stats::frequency(observed)
    if (abs((stats::tsp(value)[1] - expected) * stats::frequency(value)) > 1e-6) {
      stop("Forecasts must begin immediately after the observed series.", call. = FALSE)
    }
  }
  value
}

native_forecast_state <- function(value, reason = NULL) {
  if (is.null(value)) return(if (is.null(reason)) "not_provided" else "unavailable")
  if (all(is.na(value))) return("missing")
  if (anyNA(value)) return("partially_missing")
  "available"
}

native_new_forecast_contract <- function(points, standard_errors = list(), observed = NULL,
                                  mode = "unknown", aggregation = "direct",
                                  point_reasons = list(), error_reasons = list()) {
  components <- native_forecast_component_names()
  validate_fields <- function(values, name) {
    if (!is.list(values) || (length(values) && (is.null(names(values)) ||
        anyNA(names(values)) || anyDuplicated(names(values)) ||
        any(!names(values) %in% components)))) {
      stop(name, " must be a list named by forecast component.", call. = FALSE)
    }
  }
  validate_fields(points, "points")
  validate_fields(standard_errors, "standard_errors")
  points <- stats::setNames(lapply(components, function(name) {
    native_forecast_series(points[[name]], name, observed)
  }), components)
  if (all(vapply(points, is.null, logical(1)))) return(NULL)
  present <- points[!vapply(points, is.null, logical(1))]
  frequencies <- vapply(present, stats::frequency, numeric(1))
  if (any(frequencies != frequencies[1])) {
    stop("Forecast component frequencies must match.", call. = FALSE)
  }
  standard_errors <- stats::setNames(lapply(components, function(name) {
    value <- native_forecast_series(standard_errors[[name]], name, observed, error = TRUE)
    if (!is.null(value) && (is.null(points[[name]]) ||
        stats::frequency(value) != stats::frequency(points[[name]]) ||
        abs(stats::tsp(value)[1] - stats::tsp(points[[name]])[1]) > 1e-9 ||
        length(value) > length(points[[name]]))) {
      stop("Forecast standard errors must cover an aligned prefix of their point forecasts.",
           call. = FALSE)
    }
    value
  }), components)
  if (!is.character(mode) || length(mode) != 1L || is.na(mode)) mode <- "unknown"
  scales <- stats::setNames(c(rep("level", 3), rep(switch(mode,
    Additive = "level", Multiplicative = "factor", "unknown"), 3)), components)
  for (name in components) {
    if (!is.null(standard_errors[[name]]) &&
        length(standard_errors[[name]]) < length(points[[name]])) {
      error_reasons[[name]] <- paste0("Standard errors cover only the first ",
        length(standard_errors[[name]]), " of ", length(points[[name]]), " forecast periods.")
    }
  }
  availability <- data.frame(
    component = components,
    points = vapply(components, function(name) native_forecast_state(points[[name]], point_reasons[[name]]),
                     character(1)),
    standard_errors = vapply(components, function(name) {
      state <- native_forecast_state(standard_errors[[name]], error_reasons[[name]])
      if (state %in% c("available", "partially_missing") &&
          length(standard_errors[[name]]) < length(points[[name]])) "partial_horizon" else state
    }, character(1)),
    point_reason = vapply(components, function(name) if (is.null(point_reasons[[name]])) "" else point_reasons[[name]], character(1)),
    uncertainty_reason = vapply(components, function(name) if (is.null(error_reasons[[name]])) "" else error_reasons[[name]], character(1)),
    stringsAsFactors = FALSE, row.names = NULL)
  structure(list(schema_version = 1L, points = points,
                 standard_errors = standard_errors,
                 domains = lapply(points, native_result_domain),
                 standard_error_domains = lapply(standard_errors, native_result_domain), mode = mode,
                 scales = scales, aggregation = aggregation,
                 availability = availability), class = "persephone_forecasts")
}

