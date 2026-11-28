#' Seasonally adjust a monthly series with a standard X-13 specification
#'
#' Estimate the regular and seasonal moving-average coefficients of the
#' airline model, then use X-11 to separate trend, seasonal, and irregular
#' variation. X-11 chooses the trend and seasonal filters automatically.
#'
#' @param series Complete finite univariate monthly \code{ts} with start
#'   year from 1 to 9999.
#' @param template \code{"RSA0"}, \code{"RSA1"}, or \code{"RSA2c"},
#'   matched without regard to case.
#' @details
#' \code{RSA0} models the series in levels. \code{RSA1} also selects the
#' transformation and detects outliers. \code{RSA2c} also tests for
#' working-day, leap-year, and Easter effects. All three templates request
#' 12 months of forecasts and zero months of backcasts.
#' @return An object of class \code{native_x13_estimated_fit}. On
#'   \code{status = "available"}, \code{series}, \code{finalSA},
#'   \code{finalTrend}, and \code{finalIrregular} are dated monthly series
#'   on the input calendar. \code{model} and \code{likelihood} describe the
#'   estimated model. \code{model$calendar} reports the selected calendar
#'   terms and coefficients for \code{RSA2c}. \code{residuals} contains
#'   the monthly innovations, their standard error, and the Ljung–Box and
#'   Doornik–Hansen tests when available.
#'
#'   \code{preadjustment}, \code{selection}, \code{results}, and
#'   \code{finals} contain the X-11 tables, including the original-series
#'   forecasts in \code{preadjustment$tables$a1a}. When
#'   \code{series_forecast_stdev_status} is \code{"available"},
#'   \code{series_forecast_stdev} gives dated forecast standard errors.
#'   \code{x11_mstats} reports the X-11 diagnostic battery for an
#'   available fit: the published statistics \code{m1} to \code{m11},
#'   \code{q} and \code{qm2}, the quality grades \code{q_grade} and
#'   \code{qm2_grade}, and the battery's own \code{status}.
#'   \code{diagnostics} contains the combined seasonality and residual
#'   trading-day tests when available.
#' @seealso \code{\link{native_x13_fit}}, \code{\link{native_tramoseats_fit}}
#' @examples
#' fit <- native_x13_estimated_fit(AirPassengers)
#' head(fit$finalSA)
native_x13_estimated_fit <- function(series, template = "RSA0") {
  templates <- c(RSA0 = "RSA0", RSA1 = "RSA1", RSA2C = "RSA2c")
  if (!is.character(template) || length(template) != 1L || !is.null(dim(template)) ||
      is.na(template) || !toupper(template) %in% names(templates))
    stop("The supported stock estimated X-13 templates are RSA0, RSA1 and RSA2c.", call. = FALSE)
  domain <- native_regular_domain(series)
  if (domain$frequency != 12L || any(!is.finite(series)))
    stop("Stock estimated X-13 fits require complete finite monthly observations.", call. = FALSE)
  if (domain$start[1L] < 1L || domain$start[1L] > 9999L)
    stop("Stock estimated X-13 fits require a start year between 1 and 9999.", call. = FALSE)
  fit <- .native_x13_estimated_fit(as.numeric(series), domain$start[1L], domain$start[2L], templates[[toupper(template)]])
  fit$sigma_provenance <- "stock_template"
  fit <- date_estimated_x13_result(fit, series, domain$start, domain$frequency)
  if (length(fit$residuals$tsres)) {
    first <- domain$start_ordinal + fit$model$d + domain$frequency * fit$model$bd
    fit$residuals$tsres <- stats::ts(fit$residuals$tsres,
      start = c(first %/% domain$frequency, first %% domain$frequency + 1), frequency = domain$frequency)
  }
  structure(fit, class = c("native_x13_estimated_fit", "list"))
}
