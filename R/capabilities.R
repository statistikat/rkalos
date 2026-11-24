# Capabilities refer to the bundled, qualified surface, not planned features.
# The legacy fit field describes unrestricted fitting. Use restricted_fit
# and automatic_fit for qualified X-13 routes and tramo_seats_fit for RSA0.
capabilities <- function() {
    list(backend = "native", schema_version = 2L,
                status = "unimplemented", fit = FALSE,
                native_x13_fit_restricted = TRUE,
                native_x13_fit_automatic = TRUE,
                tramo_preprocessing = list(available = TRUE, templates = c("TR1", "TR3", "TRfull"), aliases = c("RSA1", "RSA3", "RSAfull"), frequency = 12L,
                  minimum_observations = 24L, complete_finite_input_only = TRUE,
                  by_template = list(TR1 = list(residuals = c("tsres", "ser", "lb", "doornikhansen")),
                    TR3 = list(residuals = c("tsres", "ser", "lb", "doornikhansen"), identification = "reference_retries",
                               seasonal_restart = TRUE, fast_mean_refit = TRUE),
                    TRfull = list(residuals = c("tsres", "ser", "lb", "doornikhansen"),
                                  calendar = c("working_days", "trading_days", "leap_year", "easter"))),
                  residuals = c("tsres", "ser", "lb", "doornikhansen"),
                  entry_point = "native_tramo_preprocess"),
                tramo_seats_fit = list(available = TRUE, templates = c("RSA0", "RSA1", "RSA3", "RSAfull"), aliases = c("TR0", "TR1", "TR3", "TRfull"),
                  by_template = list(
                    RSA0 = list(modes = "additive", residuals = c("tsres", "ser", "lb", "doornikhansen"),
                                series_forecast_stdev = TRUE),
                    RSA1 = list(modes = c("additive", "multiplicative"), residuals = c("tsres", "ser", "lb", "doornikhansen"),
                                outlier_types = c("AO", "LS", "TC"), tc_rate = 0.7,
                                series_forecast_stdev = TRUE),
                    RSA3 = list(modes = c("additive", "multiplicative"), residuals = c("tsres", "ser", "lb", "doornikhansen"),
                                identification = "reference_retries", seasonal_restart = TRUE,
                                fast_mean_refit = TRUE, series_forecast_stdev = TRUE),
                    RSAfull = list(modes = c("additive", "multiplicative"), residuals = c("tsres", "ser", "lb", "doornikhansen"),
                                   calendar = c("working_days", "trading_days", "leap_year", "easter"),
                                   calendar_effect = c("observed", "forecast"), series_forecast_stdev = TRUE)),
                  series_forecast_stdev = TRUE, residuals = c("tsres", "ser", "lb", "doornikhansen"),
                  frequency = 12L, complete_finite_input_only = TRUE, mode = c("additive", "multiplicative"),
                  entry_point = "native_tramoseats_fit"),
                combined_seasonality = list(getter = "diagnostics.seas-si-combined",
                  tramo_seats_fit = c("RSA0", "RSA1", "RSA3", "RSAfull"), restricted_x13_fit = TRUE,
                  automatic_x13_fit = FALSE, minimum_years = 3L),
                seats_approximation = list(mode = "legacy", frequency = 12L,
                  families = c("airline seasonal MA replacement", "mean-free (2,1,0)(1,1,0) seasonal AR-to-MA refit")),
                restricted_fit = list(
                  schema_version = 2L,     # 2L: sigma_policy_required -> sigma_policy_default (maintainer 2026-11-17)
                  native_x13_fit_restricted = TRUE,
                  fixed_coefficients_only = TRUE,    # variables carry pinned coefficients; nothing is estimated
                  fixed_seasonal_requests_only = TRUE, # one pinned filter per period; x11default/msr are refused
                  fixed_trend_length_only = TRUE,    # odd Henderson >= 5; the automatic 0 is refused
                  sigma_policy_default = "legacy",   # maintainer 2026-11-17: NULL/omitted sigma is "legacy" at this interface; the core keeps its explicitness contract (raw entries refuse a missing policy)
                  zero_horizon_only = TRUE,          # backcast = forecast = 0 is the only delivered domain
                  monthly_only = TRUE,               # qualified on the frequency-12 corpus only
                  frequency = 12L,
                  modes = c("additive", "multiplicative", "log_additive", "pseudo_additive"),
                  refusals = "verbatim unsupported: statuses from the core; see IMPLEMENTATION_STATUS.md",
                  entry_points = "native_x13_fit(), native_x13_preadjust(), native_x13_restore()"),
                # The automatic-spec fit (item 21-H stage Z), reported in the
                # rsa-template module shape the request arrives in. Each
                # module lists what is delivered; every module scalar is
                # required request data — the descriptor promotes no default
                # except the maintainer sigma-rule default below.
                automatic_fit = list(
                  schema_version = 3L,     # 3L: x11 zero_horizon_only -> horizon_envelope (item 21-H.4 R horizon qualification)
                  native_x13_fit_automatic = TRUE,
                  transform = list(
                    auto_log_test = TRUE,            # LogLevelModule transform.function = Auto is the only wired branch (decoder)
                    aicc_correction_required = TRUE  # aicDiff arrives as request data (JVM corpus: -2.0)
                  ),
                  automodel = list(
                    differencing = "required module state: max_d, max_bd, ub1, ub2, cancel, precision",
                    arma_grid = "required module state: models_count, max_p/max_q/max_bp/max_bq, balanced, mixed, accept_white_noise, eps"
                  ),
                  outliers = list(
                    types = c("AO", "LS", "TC", "SO"), # the pinned factory order subset
                    span = "full estimation domain only (a narrower span is refused)",
                    max_round_max_outliers_tcrate = "required request data"
                  ),
                  ami = list(kernel_options = "required request data: precision, mean_precision, va, reduce_va, ljung_box_limit, check_mu, accept_airline, ur_limit"),
                  regression = list(
                    mean_tests = "initial/final cutoffs required request data (decoder 1.96/1.0)",
                    final_tsig = "required request data (cascade ratioLimit, decoder 1.0)",
                    user_detection = FALSE           # deferred to 21-I; the flag answers the explicit unsupported status
                  ),
                  x11 = list(
                    msr_requests_pinned = TRUE,      # one "msr" per calendar period: the decoder G4 surface, not an option
                    automatic_henderson = TRUE,      # trend_length 0 (the restricted fit refuses it)
                    sigma_policy_default = "legacy", # same maintainer 2026-11-17 decision as restricted_fit; $sigma_provenance reports "legacy (default)" vs "(explicit)"
                    horizon_tables = c("a1a", "a1b", "d11a", "d12a", "d16a", "d18a",
                                       "d11b", "d12b", "d16b", "d18b"),
                    horizon_envelope = "backcast/forecast month counts, a negative count meaning |h| annual frequencies (the core normalizes); the 21-H.4 initialStep route extends the downstream composition only \u2014 decisions, landing and the preadjustment stay horizon-invariant",
                                         # replaces zero_horizon_only (schema 3L): R-qualified by the 4-case JVM horizon capture in test-native-x13-horizon.R; the RESTRICTED fit stays zero-horizon-only
                    modes = c("additive", "multiplicative", "log_additive", "pseudo_additive")
                  ),
                  monthly_only = TRUE,               # qualified on the frequency-12 corpus only
                  frequency = 12L,
                  complete_finite_input_only = TRUE, # missing-value interpolation is not part of the automatic description build
                  deferred = "calendar, easter-detection, automatic-trading-days and user-variable detection answer 'unsupported:x13 automatic fit defers calendar/easter/TD/user detection (21-I)'",
                  qualification = "13-case pinned JVM oracle, decision-bit exact end to end at the core tier (21-H.3); re-pinned at the R boundary in test-native-x13-automatic.R; horizons qualified on the 4-case 21-H.4 horizon capture in test-native-x13-horizon.R (bit decisions, transcribed table bands)",
                  refusals = "verbatim unsupported:/failed: statuses from the core; NEVER a fallback to the restricted fit; see IMPLEMENTATION_STATUS.md",
                  entry_points = "native_x13_fit_automatic()"),
                primitives = c("time_domain", "difference_polynomial", "convolve", "finite_filter",
                               "cholesky_solve", "qr_least_squares",
                               "transform", "inverse_transform", "prepare_inputs",
                               "stationary_likelihood", "arma_likelihood",
                               "diffuse_likelihood", "sarima_likelihood",
                               "fixed_likelihood", "regression_likelihood",
                               "arma_parameter_map", "arma_objective", "arma_optimization", "arma_identification",
                               "arma_information", "arma_residual_information", "arma_legacy_covariance",
                               "regression_information", "joint_regression_information", "arma_residuals", "arma_legacy_residuals",
                               "interpolation_information", "interpolation_arma_information", "forecast_information", "forecast_arma_information", "prediction_moments", "x11_symmetric_filters", "x11_musgrave_endpoints", "x11_mean_endpoints", "x11_seasonal_endpoints", "x11_seasonal_processing", "x11_trend_selection", "x11_seasonal_stage_policy", "x11_msr_tables", "x11_msr_pilot", "x11_msr_selection", "x11_additive_normalization", "x11_additive_initial_tables", "x11_pass_level_selection", "x11_decomposition_modes", "x13_deterministic_allocation", "x13_restricted_fixed_fit", "x13_automatic_spec_fit",
                               "tramo_critical_value", "tramo_seasonality", "tramo_transform_decision", "seats_extract_wk", "tramo_prepare_model",
                               "tramo_differencing_decision", "seats_from_tramo"),
                methods = character(), specification_blocks = character(),
                packages_available = NA, runtime_checked = FALSE,
                reason = "Unrestricted fitting is unavailable. Qualified methods are described by restricted_fit, automatic_fit and tramo_seats_fit.")
}

core_version <- function() {
  result <- .native_core_version()
  result$parallel <- .native_parallel_info()
  result$R <- R.version.string
  result$BLAS <- unname(extSoftVersion()[["BLAS"]])
  result$LAPACK <- list(version = La_version(), library = La_library())
  result
}
