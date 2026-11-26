#' Inspect TRAMO outlier and transformation decisions
#'
#' Get the outlier threshold for a series length, assess monthly seasonality,
#' or compare a log transformation with a model fitted in levels.
#'
#' @param observations Number of observations, as a nonnegative integer.
#' @param series Complete finite univariate monthly \code{ts} with at least
#'   24 observations.
#' @param log_preference Additive adjustment to the log-model criterion.
#'   Smaller values favor the log model.
#' @param precision Positive numerical tolerance for the log/level comparison.
#' @param mean Logical mean flag of the input specification. The comparison
#'   itself uses mean correction.
#' @param detect_seasonality Whether to run the seasonality tests before the
#'   transformation comparison.
#' @param spec_seasonality Logical seasonal flag in the input specification,
#'   used when seasonality detection is disabled.
#' @param robust Logical flag. Use \code{FALSE}. Robust correction for
#'   outliers before the comparison is unavailable.
#' @return \code{tramo_critical_value()} returns the numeric outlier
#'   threshold. \code{tramo_seasonality()} returns a verdict, an integer score
#'   (\code{to_int}), Qs and Friedman tests, and the dated centered first
#'   differences. \code{tramo_transform_decision()} returns the criteria
#'   for the two models when available, the adjustment, the resulting
#'   \code{log_transform} flag, and any seasonality report. A log comparison
#'   is unavailable if any observation is nonpositive.
#' @examples
#' tramo_critical_value(length(AirPassengers))
#' tramo_seasonality(AirPassengers)$verdict
#' tramo_transform_decision(AirPassengers, log(0.95), 1e-5)$log_transform
#' @rdname tramo_transform_decision
tramo_critical_value <- function(observations) {
  if (!is.null(dim(observations))) stop("observations must be scalar.", call. = FALSE)
  .native_tramo_critical_value(native_horizon(observations, "observations"))
}

tramo_monthly_values <- function(series) {
  domain <- native_regular_domain(series)
  if (domain$frequency != 12L || length(series) < 24L ||
      length(series) > .Machine$integer.max || any(!is.finite(series))) {
    stop("series must contain at least two full years of complete finite monthly data within the supported index range.",
         call. = FALSE)
  }
  as.numeric(series)
}

tramo_date_differences <- function(report, series) {
  report$differenced <- stats::ts(report$differenced, end = stats::tsp(series)[2L], frequency = 12)
  report
}

#' @rdname tramo_transform_decision
tramo_seasonality <- function(series) {
  report <- .native_tramo_seasonality(tramo_monthly_values(series))
  tramo_date_differences(report, series)
}

#' @rdname tramo_transform_decision
tramo_transform_decision <- function(series, log_preference, precision, mean = TRUE,
                                     detect_seasonality = TRUE, spec_seasonality = TRUE,
                                     robust = FALSE) {
  if (native_regression_flag(robust, "robust")) {
    stop("The TRAMO transform decision requires a series without robust AO/LS pre-correction.",
         call. = FALSE)
  }
  values <- tramo_monthly_values(series)
  for (value in list(log_preference, precision)) {
    if (!is.numeric(value) || is.complex(value) || !is.null(dim(value)) ||
        length(value) != 1L || !is.finite(value)) {
      stop("log_preference and precision must be finite numeric scalars.", call. = FALSE)
    }
  }
  result <- .native_tramo_transform_decision(
    values, native_regression_flag(mean, "mean"),
    native_regression_flag(detect_seasonality, "detect_seasonality"),
    native_regression_flag(spec_seasonality, "spec_seasonality"), log_preference, precision)
  if (detect_seasonality) {
    result$seasonality <- tramo_date_differences(result$seasonality, series)
  }
  result
}
