# Residual seasonality of the model-scale adjusted series. The fitted
# preprocessing model supplies the ordinary difference count and its mean flag.
monthly_adjusted_seasonality <- function(series, differencing_order, mean = TRUE) {
  if (!is.numeric(series) || is.complex(series) || !is.null(dim(series)))
    stop("series must be a finite numeric vector.", call. = FALSE)
  values <- as.numeric(series)
  if (!length(values) || any(!is.finite(values)))
    stop("series must contain at least one finite observation.", call. = FALSE)
  .native_adjusted_seasonality(values, native_horizon(differencing_order, "differencing_order"),
                               native_regression_flag(mean, "mean"))
}
