# Capabilities refer to the bundled, qualified surface, not planned features.
# schema_version 2 adds the restricted_fit block and the x13 primitive names;
# every pre-existing field keeps its schema-1 meaning. In particular fit
# stays FALSE and keeps its schema-1 meaning: UNRESTRICTED fitting (the ARIMA
# route: estimation, missing values, nonzero horizons, the length-of-period
# transform, benchmarking) remains unimplemented. What is qualified is the
# restricted fixed-coefficient X-13 fit described by restricted_fit below;
# gate on native_x13_fit_restricted (or restricted_fit$...) for it, never on
# fit. Field semantics are documented in IMPLEMENTATION_STATUS.md.
capabilities <- function() {
    list(backend = "native", schema_version = 2L,
                status = "unimplemented", fit = FALSE,
                native_x13_fit_restricted = TRUE,
                restricted_fit = list(
                  schema_version = 1L,
                  native_x13_fit_restricted = TRUE,
                  fixed_coefficients_only = TRUE,    # variables carry pinned coefficients; nothing is estimated
                  fixed_seasonal_requests_only = TRUE, # one pinned filter per period; x11default/msr are refused
                  fixed_trend_length_only = TRUE,    # odd Henderson >= 5; the automatic 0 is refused
                  sigma_policy_required = TRUE,      # "legacy"|"accurate"; NULL answers the core's fail-closed status
                  zero_horizon_only = TRUE,          # backcast = forecast = 0 is the only delivered domain
                  monthly_only = TRUE,               # qualified on the frequency-12 corpus only
                  frequency = 12L,
                  modes = c("additive", "multiplicative", "log_additive", "pseudo_additive"),
                  refusals = "verbatim unsupported: statuses from the core; see IMPLEMENTATION_STATUS.md",
                  entry_points = "native_x13_fit(), native_x13_preadjust(), native_x13_restore()"),
                primitives = c("time_domain", "difference_polynomial", "convolve", "finite_filter",
                               "cholesky_solve", "qr_least_squares",
                               "transform", "inverse_transform", "prepare_inputs",
                               "stationary_likelihood", "arma_likelihood",
                               "diffuse_likelihood", "sarima_likelihood",
                               "fixed_likelihood", "regression_likelihood",
                               "arma_parameter_map", "arma_objective", "arma_optimization", "arma_identification",
                               "arma_information", "arma_residual_information", "arma_legacy_covariance",
                               "regression_information", "joint_regression_information", "arma_residuals", "arma_legacy_residuals",
                               "interpolation_information", "interpolation_arma_information", "forecast_information", "forecast_arma_information", "prediction_moments", "x11_symmetric_filters", "x11_musgrave_endpoints", "x11_mean_endpoints", "x11_seasonal_endpoints", "x11_seasonal_processing", "x11_trend_selection", "x11_seasonal_stage_policy", "x11_msr_tables", "x11_msr_pilot", "x11_msr_selection", "x11_additive_normalization", "x11_additive_initial_tables", "x11_pass_level_selection", "x11_decomposition_modes", "x13_deterministic_allocation", "x13_restricted_fixed_fit"),
                methods = character(), specification_blocks = character(),
                packages_available = NA, runtime_checked = FALSE,
                reason = "Native seasonal-adjustment fitting is not implemented; only the restricted fixed-coefficient X-13 fit of restricted_fit is qualified.")
}

core_version <- function() {
  result <- .native_core_version()
  result$parallel <- .native_parallel_info()
  result$R <- R.version.string
  result$BLAS <- unname(extSoftVersion()[["BLAS"]])
  result$LAPACK <- list(version = La_version(), library = La_library())
  result
}
