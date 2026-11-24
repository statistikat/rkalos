native_tramoseats_fit <- function(series, template = "RSA0", forecast = -2L, backcast = 0L) {
  templates <- c(RSA0 = "RSA0", TR0 = "RSA0", RSA1 = "RSA1", TR1 = "RSA1", RSA3 = "RSA3", TR3 = "RSA3",
                 RSAFULL = "RSAfull", TRFULL = "RSAfull")
  if (!is.character(template) || !is.null(dim(template)) || length(template) != 1L ||
      is.na(template) || !toupper(template) %in% names(templates))
    stop("The native TRAMO/SEATS fit supports RSA0 (TR0), RSA1 (TR1), RSA3 (TR3) and RSAfull (TRfull).", call. = FALSE)
  template <- unname(templates[[toupper(template)]])
  domain <- native_regular_domain(series)
  if (domain$frequency != 12L || any(!is.finite(series)))
    stop("Native ", template, " requires complete finite monthly observations.", call. = FALSE)
  if (abs(domain$start[1L]) > .Machine$integer.max)
    stop("The start year exceeds the supported native range.", call. = FALSE)
  months <- seats_validate_horizons(backcast, forecast, domain$frequency, length(series))
  if (template != "RSA0") {
    fit <- switch(template, RSA1 = .native_tramoseats_rsa1, RSA3 = .native_tramoseats_rsa3, RSAfull = .native_tramoseats_rsafull)
    result <- fit(as.numeric(series), domain$start[1L], domain$start[2L], months[1L], months[2L])
    result <- format_tramo_preprocessing(result, series, domain)
  } else {
    result <- .native_tramoseats_rsa0(as.numeric(series), domain$start[1L], domain$start[2L], months[1L], months[2L])
    result$residuals$tsres <- stats::ts(result$residuals$tsres, start = result$residuals$start,
                                       frequency = domain$frequency)
    result$residuals$start <- NULL
  }
  result <- seats_date_components(result, domain$frequency, c("initial", "final", "merged"))
  first <- domain$start_ordinal + length(series)
  future_start <- c(first %/% domain$frequency, first %% domain$frequency + 1)
  for (name in c("series_forecast_stdev", "calendar_effect_forecast")) {
    if (length(result[[name]]))
      result[[name]] <- stats::ts(result[[name]], start = future_start, frequency = domain$frequency)
  }
  if (length(result$calendar_effect)) attributes(result$calendar_effect) <- list(tsp = stats::tsp(series), class = "ts")
  result$template <- template
  result$method <- "tramoseats"
  result$mode <- if (isTRUE(result$model$log)) "multiplicative" else "additive"
  result$series <- series
  if (identical(result$status, "ok")) {
    components <- c(finalSA = "SeasonallyAdjusted", finalTrend = "Trend",
               finalIrregular = "Irregular", finalSeasonal = "Seasonal")
    observed <- result$merged[result$merged$information == "Value", , drop = FALSE]
    for (key in names(components)) {
      value <- observed$values[[match(components[[key]], observed$component)]]
      attributes(value) <- list(tsp = stats::tsp(series), class = "ts")
      result[[key]] <- value
    }
  }
  class(result) <- c("native_tramoseats_fit", "list")
  result
}
