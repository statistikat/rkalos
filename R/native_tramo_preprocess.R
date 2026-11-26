#' Fit a monthly TRAMO preprocessing model
#'
#' Fit a seasonal ARIMA model with automatic transformation and outlier
#' selection. Use the fitted model and residual diagnostics to inspect the
#' series before seasonal adjustment.
#'
#' @param series Complete finite univariate monthly \code{ts} with at least
#'   24 observations.
#' @param template \code{"TR1"}, \code{"TR3"}, or \code{"TRfull"}.
#'   \code{"RSA1"}, \code{"RSA3"}, and \code{"RSAfull"} are equivalent
#'   names. Matching is case-insensitive.
#' @details
#' \code{TR1} fits the airline model, ARIMA(0,1,1)(0,1,1) with seasonal
#' period 12. \code{TR3} also selects the model orders. \code{TRfull} also
#' considers trading-day, leap-year, and Easter effects.
#' @return An object of class \code{native_tramo_preprocess}. \code{model}
#'   contains the transformation, orders, coefficients, and detected outliers
#'   with one-based positions. For \code{TRfull}, \code{model$calendar}
#'   reports selected terms and their coefficients. \code{likelihood} gives
#'   the log likelihood, transformation adjustment, residual sum of squares,
#'   and observation and parameter counts. \code{residuals} contains the
#'   dated innovations (\code{tsres}), standard error (\code{ser}), 24-lag
#'   Ljung–Box test (\code{lb}), and Doornik–Hansen normality test
#'   (\code{doornikhansen}).
#' @seealso \code{\link{native_tramoseats_fit}}
#' @examples
#' fit <- native_tramo_preprocess(AirPassengers, "TR1")
#' fit$model$orders
#' fit$residuals$lb
native_tramo_preprocess <- function(series, template = "TR1") {
  templates <- c(TR1 = "TR1", RSA1 = "TR1", TR3 = "TR3", RSA3 = "TR3", TRFULL = "TRfull", RSAFULL = "TRfull")
  if (!is.character(template) || length(template) != 1L || !is.null(dim(template)) ||
      is.na(template) || !toupper(template) %in% names(templates))
    stop("The native preprocessing templates are TR1 (RSA1), TR3 (RSA3) and TRfull (RSAfull).", call. = FALSE)
  template <- unname(templates[[toupper(template)]])
  domain <- native_regular_domain(series)
  if (domain$frequency != 12L || any(!is.finite(series)) || length(series) < 24L)
    stop(template, " requires at least 24 complete finite monthly observations.", call. = FALSE)
  if (abs(domain$start[1L]) > .Machine$integer.max)
    stop("The ", template, " start year exceeds the supported native range.", call. = FALSE)
  fit <- switch(template, TR1 = .native_tramo_tr1, TR3 = .native_tramo_tr3, TRfull = .native_tramo_trfull)
  result <- fit(as.numeric(series), domain$start[1L], domain$start[2L])
  result$template <- template
  result <- format_tramo_preprocessing(result, series, domain)
  class(result) <- c("native_tramo_preprocess", "list")
  result
}

format_tramo_preprocessing <- function(result, series, domain) {
  orders <- result$model$orders[c("p", "bp", "q", "bq")]
  blocks <- rep(c("phi", "bphi", "theta", "btheta"), orders)
  for (name in c("phi", "bphi", "theta", "btheta"))
    result$model[[name]] <- result$model$parameters[blocks == name]
  first <- domain$start_ordinal + length(series) - result$likelihood$effective
  result$residuals$tsres <- stats::ts(result$residuals$tsres, frequency = 12,
    start = c(first %/% 12, first %% 12 + 1))
  result
}
