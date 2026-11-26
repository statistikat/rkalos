# Public engine entry points; native_* names remain internal for migration tests.

#' Seasonal adjustment and time-series modelling
#'
#' Fit seasonal-adjustment models, inspect seasonal and calendar effects, or
#' work directly with seasonal ARIMA models and X-11 filters.
#'
#' @details
#' Use \code{\link{native_x13_estimated_fit}} for X-13 with a standard monthly
#' specification, \code{\link{native_x13_fit}} for a custom X-13 specification,
#' or \code{\link{native_tramoseats_fit}} for monthly TRAMO–SEATS.
#'
#' For modelling outside a seasonal-adjustment fit, start with
#' \code{\link{fixed_likelihood}} or \code{\link{arma_objective_setup}}.
#' \code{\link{x11_filter}} applies individual X-11 filters to a monthly or
#' quarterly series. \code{\link{capabilities}} lists the available operations.
#' @name rkalos
NULL

#' Prepare a series and aligned regressors
#'
#' Set up the observed time span and any periods requested before and after
#' it. Align regressors by date and, optionally, take logs of the response.
#'
#' @param ... Supply \code{series}, a univariate monthly or quarterly
#'   \code{ts}. Optional arguments are \code{regressors} (a named list of
#'   univariate series or a named multivariate \code{ts}),
#'   \code{transformation} (\code{"none"} or \code{"log"}), and nonnegative
#'   integer \code{forecast_horizon} and \code{backcast_horizon} (both zero
#'   by default, in observations).
#' @details
#' Missing response values retain their dates. A log transformation requires
#' positive observed values. Each regressor must have the same frequency as
#' the response and supply finite values throughout the observed and requested
#' extension periods. Values outside this span are ignored.
#' @return A list with the original and transformed response, the missing
#'   observation mask, the log-Jacobian terms, the observed and required time
#'   spans, and the aligned regressor matrix. \code{observed_rows} indexes the
#'   response dates within the required span.
#' @examples
#' inputs <- prepare_inputs(AirPassengers, transformation = "log", forecast_horizon = 12)
#' inputs$required_domain
prepare_inputs <- function(...) native_prepare_inputs(...)

#' Evaluate a seasonal ARIMA likelihood
#'
#' \code{fixed_likelihood()} evaluates the likelihood for supplied seasonal
#' ARIMA coefficients. \code{regression_likelihood()} also fits any free
#' regression coefficients with the seasonal ARIMA coefficients held fixed.
#'
#' @param ... Supply the arguments described below.
#' @details
#' Both functions take these arguments.
#' \describe{
#'   \item{\code{series}}{A univariate monthly or quarterly \code{ts}.
#'     Missing observations retain their dates within the selected span.}
#'   \item{\code{ar}, \code{ma}, \code{seasonal_ar}, \code{seasonal_ma}}{
#'     Coefficients in increasing lag order within each factor. The AR
#'     polynomial has signs \eqn{1 - \phi_1 B - \cdots} and the MA polynomial
#'     has signs \eqn{1 + \theta_1 B + \cdots}. Seasonal lags are multiples of
#'     the series frequency. Omit a vector to omit that factor.}
#'   \item{\code{d}, \code{seasonal_d}}{Regular and seasonal differencing
#'     orders from zero to two, both zero by default.}
#'   \item{\code{mean}, \code{variance}}{Mean term on the differenced model
#'     scale and positive innovation variance. Defaults are zero and one.}
#'   \item{\code{transformation}}{\code{"none"} (default) or \code{"log"}.
#'     The latter requires positive observed values in the selected span.}
#'   \item{\code{estimation_start}, \code{estimation_end}}{Inclusive endpoints
#'     as \code{c(year, period)}. \code{NULL} uses the corresponding end of
#'     the series. Only observations within the span condition the fit.}
#' }
#' \code{regression_likelihood()} also takes these arguments.
#' \describe{
#'   \item{\code{regressors}, \code{coefficients}}{Regressors as in
#'     \code{\link{prepare_inputs}} and a named vector of their coefficients.
#'     A finite value fixes a coefficient. \code{NA} or an unspecified
#'     coefficient leaves it free to estimate.}
#'   \item{\code{include_mean}, \code{include_drift}}{Logical flags, both
#'     defaulting to \code{FALSE}. A stationary model can add a \code{mean}
#'     regressor with \code{include_mean = TRUE} and \code{mean = 0}.
#'     A differenced model can add a \code{drift} regressor whose differences
#'     are constant. Use those names in \code{coefficients} to fix the
#'     corresponding terms.}
#'   \item{\code{forecast_horizon}, \code{backcast_horizon}}{Nonnegative
#'     observation counts determining the dates the regressors must cover.
#'     Both default to zero.}
#'   \item{\code{rank_tolerance}}{Numerical cutoff for deciding whether the
#'     free regressors are independently estimable. Defaults to \code{1e-10}.}
#' }
#' @return A list containing the estimation span, likelihood values using the
#'   supplied variance and a profiled variance when available, transformation
#'   adjustment, and dated predictions and innovations when available.
#'   \code{regression_likelihood()} also reports regression coefficients,
#'   their estimability and covariance, and fitted effects on the model scale.
#' @seealso \code{\link{arma_objective_setup}} to estimate free ARMA coefficients.
#' @examples
#' set.seed(42)
#' series <- ts(as.numeric(arima.sim(list(ma = 0.4), n = 96)), frequency = 12)
#' result <- fixed_likelihood(series, ma = 0.4)
#' result$likelihood
fixed_likelihood <- function(...) native_fixed_likelihood(...)

#' @rdname fixed_likelihood
regression_likelihood <- function(...) native_regression_likelihood(...)

#' Specify fixed and free ARMA coefficients
#'
#' Build a map between model coefficients and the coordinates used to vary
#' the free coefficients while enforcing stationarity and invertibility.
#'
#' @param ... For \code{arma_parameter_map()}, supply the coefficient vectors
#'   \code{ar}, \code{ma}, \code{seasonal_ar}, and \code{seasonal_ma} as in
#'   \code{\link{fixed_likelihood}}, with \code{NA} marking free coefficients.
#'   \code{initial} optionally gives a named list of starting vectors of
#'   matching lengths. \code{reflection_margin} (default \code{1e-8}) is
#'   the margin from the reflection-coefficient bounds. For
#'   \code{arma_decode()}, supply the resulting \code{parameters} map and
#'   optionally a vector of \code{coordinates} in map order.
#' @return \code{arma_parameter_map()} returns the fixed and free
#'   specification and the initial coordinates. \code{arma_decode()} returns
#'   the corresponding ordinary ARMA coefficients and feasibility information.
#' @examples
#' map <- arma_parameter_map(ar = NA_real_, initial = list(ar = 0.3))
#' arma_decode(map)
arma_parameter_map <- function(...) native_arma_parameter_map(...)

#' @rdname arma_parameter_map
arma_decode <- function(...) native_arma_decode(...)

#' Set up, evaluate, and minimize an ARMA objective
#'
#' Fit a seasonal ARIMA model with optional regression effects. Regular and
#' seasonal differencing are supplied. ARMA and regression coefficients may
#' be fixed or estimated.
#'
#' @param ... Supply the arguments described below.
#' @details
#' \code{arma_objective_setup()} takes the arguments of
#' \code{\link{regression_likelihood}} except \code{variance}. Mark free
#' ARMA coefficients with \code{NA} and optionally supply \code{initial} and
#' \code{reflection_margin} as in \code{\link{arma_parameter_map}}.
#'
#' \code{arma_objective()} and \code{arma_optimize()} take the prepared
#' \code{objective} and optional \code{coordinates}, defaulting to the starting
#' coordinates. The value minimized is the negative log likelihood after
#' estimating free regression effects and profiling out the innovation
#' variance, with the adjustment for any log transformation.
#'
#' \code{arma_optimize()} also accepts a named \code{control} list. It can set
#' \code{max_iterations}, \code{max_evaluations}, \code{max_line_search},
#' \code{gradient_tolerance}, \code{step_tolerance},
#' \code{function_tolerance}, \code{difference_step},
#' \code{difference_refinements}, \code{armijo} and \code{max_step}.
#' @return \code{arma_objective_setup()} returns the prepared data and
#'   parameter specification. \code{arma_objective()} gives the objective
#'   value, status, decoded coefficients and likelihood at the supplied
#'   coordinates. \code{arma_optimize()} returns the starting and final
#'   evaluations, final coordinates, convergence status and parameter
#'   identification diagnostics.
#' @examples
#' set.seed(42)
#' series <- ts(as.numeric(arima.sim(list(ma = 0.4), n = 96)), frequency = 12)
#' objective <- arma_objective_setup(series, ma = NA_real_)
#' arma_objective(objective)$value
#' fit <- arma_optimize(objective)
#' fit$status
arma_objective_setup <- function(...) native_arma_objective_setup(...)

#' @rdname arma_objective_setup
arma_objective <- function(...) native_arma_objective(...)

#' @rdname arma_objective_setup
arma_optimize <- function(...) native_arma_optimize(...)

#' Check local identification of ARMA coefficients
#'
#' Assess how well changes in the free coefficients can be distinguished by
#' their effect on the ARMA model. Close cancellation of AR and MA factors can
#' make coefficients weakly identified even when the model is legal.
#'
#' @param ... Supply \code{parameters} from \code{\link{arma_parameter_map}}
#'   and the seasonal \code{period} (for example 12 for monthly data).
#'   Optional \code{coordinates} default to
#'   the starting coordinates. \code{rank_tolerance} (default \code{1e-10})
#'   and \code{near_tolerance} (default \code{1e-4}) set the cutoffs for
#'   numerical rank and weak identification.
#' @return A list with the identification status, numerical rank, singular
#'   values and diagnostics for the free coefficients.
arma_identification <- function(...) native_arma_identification(...)

#' Estimate uncertainty for ARMA and regression coefficients
#'
#' Calculate local covariance and standard errors at a specified point in a
#' fitted seasonal ARIMA model with optional regression effects.
#'
#' @param ... Supply an \code{objective} from
#'   \code{\link{arma_objective_setup}} and optionally \code{coordinates}.
#'   Use the final coordinates from \code{\link{arma_optimize}} for inference
#'   at the estimated model. Omission uses the starting coordinates. The
#'   methods and other arguments are described below.
#' @details
#' For \code{arma_information()}, \code{method = "profile_ml"} uses the
#' curvature of the profiled log likelihood. \code{"conditional_residual_df"}
#' treats the fitted regression coefficients as fixed and uses residual
#' degrees of freedom. \code{"tramo_residual"} uses the derivatives of the
#' profiled residuals with residual degrees of freedom.
#'
#' For \code{regression_information()}, \code{method = "conditional_residual_df"}
#' estimates regression covariance conditional on the ARMA coefficients.
#' \code{"joint_profile_ml"} also accounts for their uncertainty and estimates
#' cross-covariance. \code{arma_method} chooses the ARMA information method
#' and must be \code{"profile_ml"} for joint inference.
#'
#' \code{control} is a named list accepting \code{difference_step},
#' \code{difference_refinements}, \code{max_evaluations},
#' \code{relative_tolerance}, \code{absolute_tolerance},
#' \code{gradient_tolerance}, \code{rank_tolerance}, and \code{weak_tolerance}.
#' \code{covariance = "accurate"} is the default. The comparison option
#' \code{"legacy"} uses a supplied \code{legacy_state} with
#' \code{final_coefficients}, \code{working_coefficients},
#' \code{difference_steps}, \code{damping}, \code{normalized_sum_squares},
#' and \code{effective_count}. Joint inference requires the default
#' covariance and a \code{NULL} legacy state.
#' @return A list containing estimated uncertainty and its availability
#'   status. ARMA information reports coefficient-space derivatives and
#'   covariance where estimable. Regression information reports covariance
#'   for the free regressors and, for joint inference, results for the
#'   combined regression and ARMA coefficients.
arma_information <- function(...) native_arma_information(...)

#' @rdname arma_information
regression_information <- function(...) native_regression_information(...)

#' Calculate standardized ARMA residuals
#'
#' Compute dated innovations and standardized residuals for evaluating model
#' fit. When the series has gaps, the result can also report interpolated
#' values and residuals for the completed series.
#'
#' @param ... Supply an \code{objective} from
#'   \code{\link{arma_objective_setup}}. Optional \code{coordinates} default
#'   to its starting coordinates. Set \code{scale} to \code{"ml"} (default)
#'   for the maximum-likelihood variance, or \code{"residual_df"} to scale
#'   by residual degrees of freedom. \code{reporting} selects the
#'   \code{"accurate"} (default) or \code{"legacy"} residual convention.
#' @return A list with residuals before and, when available, after filling
#'   gaps, including the dates of missing values, interpolated points,
#'   innovations, standardized residuals and availability statuses.
arma_residuals <- function(...) native_arma_residuals(...)

#' Interpolate, forecast, and backcast a seasonal ARIMA series
#'
#' \code{interpolation_information()} predicts missing values within the
#' estimation span. \code{arma_forecasts()} predicts values before and after
#' it. Both report prediction uncertainty when it can be estimated.
#'
#' @param ... Supply the arguments described below.
#' @details
#' Both functions accept an \code{objective} from
#' \code{\link{arma_objective_setup}} and optional \code{coordinates}. Omission
#' uses the starting coordinates. Use \code{scale = "ml"} (the default) for
#' the maximum-likelihood variance, or \code{"residual_df"} to scale by the
#' residual degrees of freedom. \code{method = "conditional"} holds the ARMA
#' coefficients fixed for calculating uncertainty. \code{"delta_profile_ml"}
#' also propagates local ARMA coefficient uncertainty and requires
#' \code{scale = "ml"}. With that method, \code{control} accepts the
#' information settings described in \code{\link{arma_information}}.
#'
#' \code{arma_forecasts()} also accepts nonnegative observation counts
#' \code{forecast_horizon} and \code{backcast_horizon}; by default they are
#' the counts set on the objective, which are also their upper bounds.
#' \code{forecast_error_horizon} and \code{backcast_error_horizon}
#' default to the corresponding point horizon
#' and can request uncertainty at fewer dates. Regressors must cover the
#' requested periods.
#' @return A list with prediction dates, points, standard errors, covariance
#'   and availability statuses. For a log model, exponentiated points are
#'   medians. Original-scale means and covariance, when available, use
#'   lognormal moments. See \code{\link{prediction_moments}}.
interpolation_information <- function(...) native_interpolation_information(...)

#' @rdname interpolation_information
arma_forecasts <- function(...) native_arma_forecasts(...)

#' Summarize a Gaussian prediction on its original scale
#'
#' Compute means, medians, covariance and standard errors of a joint Gaussian
#' prediction or of its elementwise exponential for a log-transformed series.
#'
#' @param ... Supply \code{location}, a numeric vector of predictive locations,
#'   and \code{covariance}, a symmetric covariance matrix of matching size.
#'   Use \code{transformation = "none"} (default) or \code{"log"}.
#' @details For a log model, the median is \eqn{\exp(\mu)} and the mean is
#'   \eqn{\exp(\mu + \sigma^2/2)}. To get a marginal prediction interval,
#'   exponentiate the corresponding Gaussian quantiles on the log scale.
#' @return A list with \code{medians}, \code{means}, \code{covariance},
#'   \code{standard_errors} and indicators of which quantities are available.
prediction_moments <- function(...) native_prediction_moments(...)

#' Get or apply X-11 moving-average filters
#'
#' Get filter weights and their offsets, or apply a filter to a univariate
#' monthly or quarterly time series.
#'
#' @param ... Supply the arguments described below.
#' @details
#' \code{x11_filter_spec()} takes \code{type} and \code{frequency} (4 or 12).
#' The types are \code{"trend"} for the centered annual moving average,
#' \code{"henderson"}, and the seasonal filters \code{"s3x1"},
#' \code{"s3x3"}, \code{"s3x5"}, \code{"s3x9"}, and \code{"s3x15"}.
#' Supply an odd \code{length} for a Henderson filter.
#'
#' \code{x11_musgrave_spec()} takes \code{frequency}, Henderson \code{length},
#' the number of \code{available} observations on the truncated side, and
#' \code{ic_ratio}, the irregular-to-trend ratio used for the endpoint weights.
#' \code{x11_seasonal_spec()} takes a seasonal \code{type}, \code{frequency},
#' and \code{available}, counted within the same month or quarter. Both take
#' \code{side = "right"} (default) or \code{"left"} for the end being filtered.
#'
#' \code{x11_filter()} takes a \code{series} as a \code{ts} and the same
#' \code{type} and optional \code{length} as \code{x11_filter_spec()}.
#' \code{endpoints = "none"} leaves values unavailable where there is
#' insufficient support. Use \code{"musgrave"} with a Henderson filter and
#' supply \code{ic_ratio}, \code{"mean"} with the centered trend filter, or
#' \code{"seasonal"} with a supported seasonal filter.
#' @return A filter description with \code{weights} and \code{offsets}.
#'   \code{x11_filter()} also returns the dated filtered \code{values}, a
#'   missing-value mask and indices for the available values and endpoints.
#' @examples
#' x11_filter_spec("henderson", frequency = 12, length = 13)$weights
#' smoothed <- x11_filter(AirPassengers, "trend", endpoints = "mean")
#' head(smoothed$values)
x11_filter_spec <- function(...) native_x11_filter_spec(...)

#' @rdname x11_filter_spec
x11_musgrave_spec <- function(...) native_x11_musgrave_spec(...)

#' @rdname x11_filter_spec
x11_seasonal_spec <- function(...) native_x11_seasonal_spec(...)

#' @rdname x11_filter_spec
x11_filter <- function(...) native_x11_filter(...)

#' Smooth, center, and extend X-11 seasonal estimates
#'
#' Smooth seasonal-irregular values within each calendar period, center
#' additive seasonal estimates, or calculate the opening X-11 tables.
#'
#' @param ... Supply the arguments described below.
#' @details
#' All three functions take a monthly or quarterly \code{series} as a
#' \code{ts} and fixed \code{filters}. Use one filter for all calendar
#' periods, or one per period in calendar order. If named, use the period
#' numbers as names. The choices are \code{"stable"} and \code{"s3x1"},
#' \code{"s3x3"}, \code{"s3x5"}, \code{"s3x9"}, and \code{"s3x15"}.
#'
#' \code{x11_seasonal_process()} smooths the subseries for each month or
#' quarter. \code{short_series = "x11"} allows a shorter or stable filter
#' where the requested filter has insufficient data. \code{"reject"}
#' requires enough observations for the requested filter.
#'
#' \code{x11_seasonal_normalize()} subtracts a local annual average to
#' center seasonal estimates around zero. \code{extend} is the number of
#' observations to add at each end using the nearest annual cycle. It defaults
#' to zero. The supported \code{mode} is \code{"additive"}.
#'
#' \code{x11_initial_additive()} takes the series to decompose and computes
#' the first additive B tables, through the initial irregular. Its supported
#' \code{mode} is \code{"additive"}.
#' @return A list with dated seasonal estimates and the applied filters,
#'   endpoint or short-span decisions, and availability information for
#'   seasonal processing and normalization. The initial decomposition
#'   returns dated \code{tables} (\code{B1}, \code{B2}, \code{B3}, \code{B4a},
#'   \code{B4anorm}, and \code{B4d}) and the indices into the original series.
#' @seealso \code{\link{x11_seasonal_plan}} for choosing filters by stage.
x11_seasonal_process <- function(...) native_x11_seasonal_process(...)

#' Select a Henderson trend filter
#'
#' Choose a Henderson filter length from the irregular-to-trend change ratio
#' (I/C), or compute that ratio from a seasonally adjusted series.
#'
#' @param ... For \code{x11_trend_policy()}, supply \code{frequency} (4 or 12)
#'   and exactly one of \code{ic_ratio} on the comparison scale or a fixed
#'   odd Henderson \code{length}. The quarterly comparison ratio is three
#'   times the raw I/C ratio. For \code{x11_trend_select()}, supply the
#'   seasonally adjusted monthly or quarterly \code{series} on its working
#'   scale. Optional \code{backcasts} and \code{forecasts} count observations
#'   to exclude at either end. Set \code{mode} to \code{"additive"} (default)
#'   or \code{"log_additive"} for absolute changes, or \code{"multiplicative"}
#'   or \code{"pseudo_additive"} for relative changes.
#' @return \code{x11_trend_policy()} returns the filter length and the I/C
#'   ratio used for the endpoint weights. \code{x11_trend_select()} reports
#'   the preliminary trend and irregular change magnitudes, the I/C ratio,
#'   selected length when available, status and spans used in the calculation.
x11_trend_policy <- function(...) native_x11_trend_policy(...)

#' @rdname x11_trend_policy
x11_trend_select <- function(...) native_x11_trend_select(...)

#' Plan seasonal filters for an X-11 stage
#'
#' Translate seasonal filter requests into the filters applied at a given
#' stage. A final D-stage plan can require an MSR result before its filters
#' can be applied.
#'
#' @param ... For \code{x11_seasonal_plan()}, supply \code{frequency} (4
#'   or 12), \code{filters}, and \code{stage}. For filters, use one request
#'   or one per calendar period, with period-number names if named. Choices
#'   are \code{"stable"}, \code{"s3x1"}, \code{"s3x3"}, \code{"s3x5"},
#'   \code{"s3x9"}, \code{"s3x15"}, \code{"x11default"}, and \code{"msr"}.
#'   Stages are \code{"b_initial"}, \code{"b_final"}, \code{"c_initial"},
#'   \code{"c_final"}, \code{"d_initial"}, and \code{"d_final"}.
#'   For \code{x11_msr_resolve()}, pass a pending D-final \code{plan} and a
#'   \code{selection} from \code{\link{x11_msr_select}} at the same frequency.
#' @return A list with the requested and resolved filters and reasons for
#'   the choices. \code{requires_msr} indicates a pending decision. A
#'   selection that cannot produce a filter leaves the plan pending.
x11_seasonal_plan <- function(...) native_x11_seasonal_plan(...)

#' Measure seasonal and irregular change for filter selection
#'
#' The moving seasonality ratio (MSR) compares irregular change with seasonal
#' change between observations one year apart. Use it to assess how much
#' smoothing to apply to the seasonal estimates.
#'
#' @param ... Supply the arguments described below.
#' @details
#' \code{x11_msr_table()} takes \code{seasonal} and \code{irregular} monthly
#' or quarterly \code{ts} series on the same calendar. \code{x11_msr_pilot()}
#' and \code{x11_msr_select()} take a seasonal-irregular \code{series} as
#' a \code{ts}. Supply data on the working scale in all three cases. The
#' \code{mode} is \code{"multiplicative"} for relative changes,
#' or \code{"additive"} (default), \code{"log_additive"}, or
#' \code{"pseudo_additive"} for absolute changes.
#'
#' \code{x11_msr_policy()} takes the numeric \code{ratio} and identifies the
#' filter band or a gap between bands.
#' \code{x11_msr_select()} can revisit the decision on a shorter span when
#' the ratio falls in a gap. Its \code{backcasts} and \code{forecasts}
#' arguments identify extension observations to exclude. \code{max_iterations}
#' bounds the number of attempts. \code{iteration_policy = "jdplus"} is
#' the default. \code{"x11_five_year_cap"} limits the span reduction to
#' five years.
#' @return The table reports change magnitudes and ratios by calendar period
#'   and \code{global_msr} for the full span. The policy reports a band and
#'   filter when defined. The pilot returns preliminary seasonal and irregular
#'   series when available. Selection reports \code{selected_filter},
#'   \code{status}, the span used, attempted ratios and whether it used the
#'   default filter.
x11_msr_table <- function(...) native_x11_msr_table(...)

#' @rdname x11_msr_table
x11_msr_policy <- function(...) native_x11_msr_policy(...)

#' @rdname x11_msr_table
x11_msr_pilot <- function(...) native_x11_msr_pilot(...)

#' @rdname x11_msr_table
x11_msr_select <- function(...) native_x11_msr_select(...)

#' @rdname x11_seasonal_plan
x11_msr_resolve <- function(...) native_x11_msr_resolve(...)

#' @rdname x11_seasonal_process
x11_seasonal_normalize <- function(...) native_x11_seasonal_normalize(...)

#' @rdname x11_seasonal_process
x11_initial_additive <- function(...) native_x11_initial_additive(...)
