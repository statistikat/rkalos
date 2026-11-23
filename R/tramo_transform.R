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

tramo_seasonality <- function(series) {
  report <- .native_tramo_seasonality(tramo_monthly_values(series))
  tramo_date_differences(report, series)
}

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
