native_tramoseats_fit <- function(series, template = "RSA0", forecast = -2L, backcast = 0L) {
  if (!is.character(template) || !is.null(dim(template)) || length(template) != 1L ||
      is.na(template) || !toupper(template) %in% c("RSA0", "TR0"))
    stop("The native TRAMO/SEATS fit supports the RSA0 (TR0) template.", call. = FALSE)
  domain <- native_regular_domain(series)
  if (domain$frequency != 12L || any(!is.finite(series)))
    stop("Native RSA0 requires complete finite monthly observations.", call. = FALSE)
  if (abs(domain$start[1L]) > .Machine$integer.max)
    stop("The RSA0 start year exceeds the supported native range.", call. = FALSE)
  seats_validate_horizons(backcast, forecast, domain$frequency, length(series))
  result <- .native_tramoseats_rsa0(as.numeric(series), domain$start[1L], domain$start[2L], backcast, forecast)
  result <- seats_date_components(result, domain$frequency, c("initial", "final", "merged"))
  result$template <- "RSA0"
  result$method <- "tramoseats"
  result$mode <- "additive"
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
