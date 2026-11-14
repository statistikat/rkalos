# Capabilities refer to the bundled, qualified surface, not planned features.
capabilities <- function() {
    list(backend = "native", schema_version = 1L,
                status = "unimplemented", fit = FALSE,
                primitives = c("time_domain", "difference_polynomial", "convolve", "finite_filter",
                               "cholesky_solve", "qr_least_squares",
                               "transform", "inverse_transform", "prepare_inputs",
                               "stationary_likelihood", "arma_likelihood",
                               "diffuse_likelihood", "sarima_likelihood",
                               "fixed_likelihood", "regression_likelihood",
                               "arma_parameter_map", "arma_objective", "arma_optimization", "arma_identification",
                               "arma_information", "arma_residual_information", "arma_legacy_covariance",
                               "regression_information", "joint_regression_information", "arma_residuals", "arma_legacy_residuals",
                               "interpolation_information", "interpolation_arma_information", "forecast_information", "forecast_arma_information", "prediction_moments", "x11_symmetric_filters", "x11_musgrave_endpoints", "x11_mean_endpoints", "x11_seasonal_endpoints", "x11_seasonal_processing", "x11_trend_selection", "x11_seasonal_stage_policy", "x11_msr_tables", "x11_msr_pilot", "x11_msr_selection", "x11_additive_normalization", "x11_additive_initial_tables"),
                methods = character(), specification_blocks = character(),
                packages_available = NA, runtime_checked = FALSE,
                reason = "Native seasonal-adjustment fitting is not implemented.")
}

core_version <- function() {
  result <- .native_core_version()
  result$parallel <- .native_parallel_info()
  result$R <- R.version.string
  result$BLAS <- unname(extSoftVersion()[["BLAS"]])
  result$LAPACK <- list(version = La_version(), library = La_library())
  result
}
