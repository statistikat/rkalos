#' Seasonally adjust a monthly series with TRAMO–SEATS
#'
#' Fit a TRAMO model to the observations, use SEATS to separate trend,
#' seasonal, and irregular variation, and restore estimated deterministic
#' effects to the appropriate components.
#'
#' @param series Complete finite univariate monthly \code{ts}. \code{RSA0}
#'   requires enough observations to fit the airline model. \code{RSA1},
#'   \code{RSA3}, and \code{RSAfull} require at least 24.
#' @param template \code{"RSA0"}, \code{"RSA1"}, \code{"RSA3"}, or
#'   \code{"RSAfull"}, also accepting \code{"TR0"}, \code{"TR1"},
#'   \code{"TR3"}, and \code{"TRfull"}, respectively. Matching is
#'   case-insensitive.
#' @param forecast,backcast Numbers of months to extend after and before
#'   the series. Negative counts denote years, so \code{forecast = -2L}
#'   requests 24 forecast months. Use zero for no extension.
#' @details
#' \code{RSA0} fits the airline model to levels. \code{RSA1} also selects
#' the transformation and outliers. \code{RSA3} additionally selects the
#' model orders, and \code{RSAfull} also considers trading-day, leap-year,
#' and Easter effects. Decomposition is additive for level fits and
#' multiplicative for log fits.
#'
#' The combined test for seasonality in the seasonal-irregular series is
#' in \code{diagnostics$combined_seasonality}. Its summary is
#' \code{"Present"}, \code{"ProbablyNone"}, or \code{"None"} when at least
#' three years are available. Residual trading-day results, when available,
#' are in \code{diagnostics[["td.ftest.sa"]]} and
#' \code{diagnostics[["td.ftest.i"]]}.
#' @return An object of class \code{native_tramoseats_fit}. A successful
#'   fit has \code{status = "ok"} and returns \code{finalSA},
#'   \code{finalTrend}, \code{finalSeasonal}, and \code{finalIrregular}
#'   as monthly series on the input calendar. \code{initial} and \code{final}
#'   contain dated SEATS component values, extensions, and standard
#'   deviations. \code{merged} contains the final component values and
#'   extensions on the original scale.
#'
#'   \code{model} and \code{likelihood} describe the fitted model.
#'   \code{model_changed}, when returned, indicates whether SEATS adjusted
#'   the model before decomposition. \code{residuals} contains dated
#'   innovations (\code{tsres}), their standard error (\code{ser}), the
#'   24-lag Ljung–Box test (\code{lb}), and the Doornik–Hansen normality
#'   test (\code{doornikhansen}). The tests report \code{value},
#'   \code{pvalue}, and \code{df} when available. If
#'   \code{series_forecast_stdev_status} is \code{"available"},
#'   \code{series_forecast_stdev} gives dated standard errors for forecasts
#'   of the original series. A zero forecast horizon reports
#'   \code{"no_requested_predictions"}. \code{RSAfull} also returns selected calendar
#'   terms in \code{model$calendar}, and their observed and future effects
#'   in \code{calendar_effect} and \code{calendar_effect_forecast}.
#' @seealso \code{\link{native_tramo_preprocess}}, \code{\link{seats_from_tramo}}
#' @examples
#' fit <- native_tramoseats_fit(AirPassengers, forecast = 12)
#' fit$status
#' head(fit$finalSA)
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
    # Public residual trading-day F-tests (M7 cutover): reference-parity input
    # selection (Persephone's M7 root cause) — the JVM's GenericSaTests feeds
    # the SEATS INITIAL (stochastic) decomposition on the model's working
    # scale, and the bridge publishes exactly that block: the initial
    # SeasonallyAdjusted/Irregular Value rows with their stamped domain start
    # replace the final merged rows for every template. The two blocks coincide
    # for additive fits (transform.fun None); for multiplicative fits the
    # initial rows carry the estimation-scale series the reference machine
    # tests. Missing rows or disagreeing starts are a decline: the attach then
    # omits the fields, never fabricates a component.
    initial_sa <- .seats_initial_td_components(result$initial, "SeasonallyAdjusted")
    initial_irregular <- .seats_initial_td_components(result$initial, "Irregular")
    if (!is.null(initial_sa) && !is.null(initial_irregular) &&
        identical(initial_sa$start, initial_irregular$start))
      result$diagnostics <- .attach_td_tests(result$diagnostics, initial_sa$values,
                                             initial_irregular$values, initial_sa$start)
  }
  class(result) <- c("native_tramoseats_fit", "list")
  result
}
