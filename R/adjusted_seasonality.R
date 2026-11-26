# Residual seasonality of the model-scale adjusted series. The fitted
# preprocessing model supplies the ordinary difference count and its mean flag.
#' Test for residual seasonality in an adjusted monthly series
#'
#' Apply the Qs test and a seasonal F-test to seasonally adjusted
#' observations on the fitted model's scale, using the specified differencing
#' order and mean setting.
#'
#' @param series Nonempty finite numeric vector of monthly adjusted values in
#'   chronological order. Supply the observed adjusted series on the model's
#'   scale, for example log values for a log-scale model.
#' @param differencing_order Number of ordinary differences to take,
#'   from zero to three. Use the fitted model's \code{d + bd}.
#' @param mean Whether to remove the mean after differencing. Use the
#'   mean flag of the fitted preprocessing model.
#' @details Qs measures autocorrelation at seasonal lags. The F-test
#'   tests for month-of-year effects after accounting for the previous
#'   observation.
#' @return A list with \code{qs} and \code{seasonal_f}, each containing
#'   the test statistic \code{value} and upper-tail \code{p_value}, or
#'   \code{NULL} if the test is unavailable. The F-test is absent below
#'   14 observations after differencing. With exactly 14 and a full-rank
#'   design, it returns undefined values.
#' @examples
#' fit <- native_x13_estimated_fit(AirPassengers)
#' monthly_adjusted_seasonality(fit$finalSA, fit$model$d + fit$model$bd)
monthly_adjusted_seasonality <- function(series, differencing_order, mean = TRUE) {
  if (!is.numeric(series) || is.complex(series) || !is.null(dim(series)))
    stop("series must be a finite numeric vector.", call. = FALSE)
  values <- as.numeric(series)
  if (!length(values) || any(!is.finite(values)))
    stop("series must contain at least one finite observation.", call. = FALSE)
  .native_adjusted_seasonality(values, native_horizon(differencing_order, "differencing_order"),
                               native_regression_flag(mean, "mean"))
}
