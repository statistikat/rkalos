# Item 22c M3-finale: the restricted fixed-model X-13 fit — one composition
# of the qualified pieces in the pinned X13Kernel.process order: the
# fixed-coefficient preadjust (native_x13_preadjust's seam), the pass-level
# X-11 selection run over `out` with the updateSpec-forced mode, the
# log-additive results-builder stage where the mode is log_additive, and the
# finals restore (native_x13_restore's seam). The delivered surface is
# restricted by the pinned contract, and the boundary mirrors the core's
# fail-closed answers verbatim instead of inventing R-side refusals: only
# fixed seasonal requests (a pinned SeasonalFilterOption per calendar period;
# "x11default"/"msr" answer "unsupported:x13 fit fixed seasonal requests
# only"), a fixed Henderson trend_length (an odd length >= 5; the automatic 0
# answers "unsupported:x13 fit fixed Henderson trend length"), the pinned
# fixed-path horizon domain backcast=forecast=0 (anything else answers
# "unsupported:x13 fit requires zero horizons"), and the sigma policy: NULL
# (or omitted) receives the maintainer default "legacy" — the Kalos decision
# of 2026-11-17 promotes legacy (bit-exact against the pinned
# DefaultExtremeValuesCorrector / JDemetra behavior) to the product default,
# applied here at the interface, while the core keeps its explicitness
# contract and the raw entry still answers the pinned "unsupported:x13 fit
# requires an explicit sigma policy" for a missing policy. An explicitly
# stated policy always wins. Every refusal still echoes the full provenance
# of the request.
# Estimation, missing values, the length-of-period transform, benchmarking,
# and the ARIMA route stay outside. The completed fit reports combined SI
# seasonality; capabilities() reports the
# restricted_fit scope for callers to gate on.
#' Seasonally adjust a monthly series with a custom X-13 specification
#'
#' \code{native_x13_fit()} adjusts for supplied regression effects and
#' applies X-11 decomposition with specified filters.
#' \code{native_x13_fit_automatic()} selects the transformation, outliers,
#' and seasonal ARIMA model, then estimates the effects and applies X-11
#' with automatic seasonal and Henderson filter selection.
#'
#' @param series Complete finite univariate monthly \code{ts} or numeric
#'   vector in observation order.
#' @param variables List of fixed regression effects. Each element has
#'   \code{class}, \code{effect}, numeric \code{coefficient}, and
#'   \code{columns}, a list with one observed-length numeric vector per
#'   coefficient. Classes are \code{"outlier"}, \code{"days_related"},
#'   \code{"moving_holiday"}, and \code{"user"}. The effect assigns the
#'   variable to \code{"trend"}, \code{"seasonal"}, \code{"irregular"},
#'   \code{"seasonally_adjusted"}, \code{"series"}, or \code{"calendar_effect"}.
#'   Use \code{"undefined"} to leave the assignment unspecified.
#' @param log Whether to use logs for preadjustment. Requires positive data.
#' @param seasonal Fixed seasonal filter names, one for each of the 12
#'   calendar months, chosen from \code{"stable"}, \code{"s3x1"},
#'   \code{"s3x3"}, \code{"s3x5"}, \code{"s3x9"}, and \code{"s3x15"}.
#'   If named, use \code{"1"} through \code{"12"}.
#' @param trend_length Odd Henderson filter length of at least five.
#' @param sigma Extreme-value weighting rule, \code{"legacy"} or
#'   \code{"accurate"}. \code{NULL} uses \code{"legacy"}.
#' @param lower_sigma,upper_sigma Lower and upper sigma cutoffs for
#'   downweighting unusually large irregular values.
#' @param msr_policy Policy for reconsidering the moving seasonality ratio
#'   on a shorter span, \code{"jdplus"} or \code{"x11_five_year_cap"}. Applies
#'   when seasonal filter selection uses the moving seasonality ratio.
#' @param mode Requested X-11 decomposition mode, one of \code{"additive"},
#'   \code{"multiplicative"}, \code{"log_additive"}, and
#'   \code{"pseudo_additive"}.
#' @param backcast,forecast Extension counts before and after the data.
#'   \code{native_x13_fit()} requires zero. For the automatic fit, positive
#'   counts denote months and negative counts denote years.
#' @param frequency,start_period,first_year Calendar frequency (12) and
#'   first month (1–12) and year. For a numeric vector, supply
#'   \code{frequency} and \code{start_period}, and supply \code{first_year}
#'   for the automatic fit. The fixed fit dates a vector from year 1.
#'   For a \code{ts}, these can be read from the series. Explicit calendar
#'   arguments for the automatic fit must agree with its calendar.
#' @param transform Named list with \code{aicc}, the adjustment to the
#'   log-model AICc for the log-versus-level comparison.
#' @param automodel Named list controlling the ARIMA order search.
#'   See Details for required entries.
#' @param outliers Named list controlling automatic outlier detection.
#'   See Details for required entries.
#' @param ami Named list controlling estimation, model checks, and
#'   outlier thresholds. See Details for required entries.
#' @param regression Named list with \code{mean_test_initial} and
#'   \code{mean_test_final} (initial and final mean-test cutoffs, with
#'   zero disabling the corresponding test), and \code{tsig} (the
#'   significance cutoff used to check final ARMA terms).
#' @param calendar,easter,trading_days,user_variables Logical switches
#'   for automatic detection of those effects. Use \code{FALSE}.
#' @details
#' For \code{native_x13_fit_automatic()}, supply the following named
#' controls in addition to \code{transform} and \code{regression}.
#'
#' \describe{
#'   \item{\code{automodel}}{\code{max_d} and \code{max_bd} cap regular
#'     and seasonal differencing. \code{ub1} and \code{ub2} are unit-root
#'     cutoffs and \code{cancel} is the cancellation tolerance.
#'     \code{precision} controls estimation in the differencing search.
#'     \code{models_count} sets the number of ARMA candidates to retain.
#'     \code{max_p}, \code{max_q}, \code{max_bp}, and \code{max_bq}
#'     cap regular AR, regular MA, seasonal AR, and seasonal MA orders.
#'     \code{balanced}, \code{mixed}, and \code{accept_white_noise} are
#'     logical switches for the search. \code{eps} is the estimation
#'     tolerance for ARMA selection.}
#'   \item{\code{outliers}}{\code{used} switches detection on.
#'     \code{types} is the character vector of candidate outliers
#'     (\code{"AO"}, \code{"LS"}, \code{"TC"}, \code{"SO"}), and may be
#'     omitted when detection is off. Use \code{span = NULL} or
#'     \code{c(1L, length(series))} for the full span.
#'     \code{max_round} and \code{max_outliers} limit the search.
#'     \code{tcrate} is the decay rate for temporary changes.}
#'   \item{\code{ami}}{\code{precision} and \code{mean_precision}
#'     control estimation and the mean check, respectively.
#'     \code{va} sets the outlier critical value (zero derives it from the
#'     sample size) and \code{reduce_va} controls threshold reduction.
#'     \code{ljung_box_limit} is the model-check cutoff.
#'     \code{check_mu} and \code{accept_airline} are logical flags for
#'     checking the mean and accepting an airline model.
#'     \code{ur_limit} is a near-unit-root cutoff.}
#' }
#' @return A named list. Check \code{status} or \code{fitted} before
#'   using component values. Successful fits (\code{status = "available"})
#'   give \code{finalSA}, \code{finalTrend}, and \code{finalIrregular}
#'   on the observed calendar. \code{preadjustment}, \code{selection},
#'   \code{results}, and \code{finals} contain the preadjusted, filter,
#'   and decomposition results. \code{diagnostics$combined_seasonality}
#'   contains the combined seasonality verdict and its tests when available.
#'   Residual trading-day tests appear as \code{td.ftest.sa} and
#'   \code{td.ftest.i} in \code{diagnostics} when they can be computed.
#'
#'   Automatic fits also report the model, decisions, and residuals.
#'   For a requested extension, \code{preadjustment$tables$a1a} and
#'   \code{a1b} are forecasts and backcasts of the original series.
#'   \code{finals$tables} holds dated extensions of the adjusted series
#'   (\code{d11}), trend (\code{d12}), combined adjustment factor
#'   (\code{d16}), and calendar effects (\code{d18}), with suffix \code{a}
#'   for forecasts and \code{b} for backcasts. Zero-length extensions
#'   omit their tables. When available, \code{series_forecast_stdev}
#'   reports standard errors of the original-series forecasts.
#'   \code{x11_mstats} reports the X-11 diagnostic battery for an
#'   available fit: the published statistics \code{m1} to \code{m11},
#'   \code{q} and \code{qm2}, the quality grades \code{q_grade} and
#'   \code{qm2_grade}, and the battery's own \code{status}.
#'   \code{irregular_forecast} reports the published irregular-component
#'   forecast for an available fit, as the reference publishes it: the
#'   constant placeholder over the finals trend-forecast domain, as the
#'   pair \code{status} and \code{values} — 1.0 on multiplicative-domain
#'   fits (including log fits) and 0.0 on additive and pseudo-additive
#'   fits, with the zero-horizon status \code{"no_requested_predictions"}
#'   and an empty vector.
#' @seealso \code{\link{native_x13_estimated_fit}}, \code{\link{capabilities}}
#' @examples
#' fit <- native_x13_fit(AirPassengers, seasonal = rep("s3x3", 12),
#'   trend_length = 13L)
#' fit$status
#' head(fit$finalSA)
native_x13_fit <- function(series,variables=list(),log=FALSE,seasonal,trend_length,
    sigma=NULL,msr_policy="jdplus",frequency=NULL,start_period=NULL,lower_sigma=1.5,
    upper_sigma=2.5,mode="additive",backcast=0,forecast=0) {
  if(missing(seasonal)||missing(trend_length))
    stop("The per-period seasonal requests and the fixed trend filter length must be stated explicitly.",call.=FALSE)
  if(is.null(sigma))sigma<-"legacy" # the maintainer 2026-11-17 product default
  if(is.null(frequency)||is.null(start_period)) {
    span<-stats::tsp(series)
    if(is.null(span))stop("A plain numeric series needs explicit frequency and start_period.",call.=FALSE)
    if(is.null(frequency))frequency<-span[3L]
    if(is.null(start_period))start_period<-round((span[1L]-floor(span[1L]))*frequency)+1L
  }
  f<-as.integer(frequency);periods<-as.character(seq_len(f))
  if(!is.null(names(seasonal))) {
    if(length(seasonal)!=f || anyNA(names(seasonal)) || anyDuplicated(names(seasonal)) || !setequal(names(seasonal),periods))
      stop("Named seasonal requests must identify every calendar period exactly once.",call.=FALSE)
    seasonal<-seasonal[match(periods,names(seasonal))]
  }
  if(!is.list(variables))
    stop("The x13 model variables must be a list of variable specifications.",call.=FALSE)
  tsp_orig<-stats::tsp(series)
  y<-if(!is.null(dim(series)))as.vector(series) else series
  if(!is.numeric(y))y<-as.numeric(y)
  result<-.native_x13_fit(y,variables,log,f,start_period,as.character(seasonal),
    trend_length,sigma,lower_sigma,upper_sigma,msr_policy,mode,backcast,forecast)
  ordinal<-if(!is.null(tsp_orig))round(tsp_orig[1L]*f) else f+as.integer(start_period)-1L
  start<-if(!is.null(tsp_orig))
    c(floor(tsp_orig[1L]),as.integer(round((tsp_orig[1L]-floor(tsp_orig[1L]))*f))+1L)
    else c(1L,as.integer(start_period))
  if(!is.null(result[["preadjustment"]])) { # past the fail-closed guards
    result$preadjustment$tables<-.x13_decorate_tables(result$preadjustment$tables,start,f)
    result$preadjustment$row_basis<-rep(list(seq_along(y)),length(result$preadjustment$tables))
    names(result$preadjustment$row_basis)<-names(result$preadjustment$tables)
    result$preadjustment$available_rows<-lapply(names(result$preadjustment$tables),
      function(name)result$preadjustment$row_basis[[name]][!is.na(result$preadjustment$tables[[name]])])
    names(result$preadjustment$available_rows)<-names(result$preadjustment$tables)
  }
  if(!is.null(result[["selection"]])) { # the same record the selection wrapper delivers
    result$selection<-.x11_decorate_selection_record(result$selection,ordinal,f,length(y))
  }
  if(identical(result$status,"available")) {
    for(name in c("d10","d11","d12","d13","c17"))
      result$results[[name]]<-stats::ts(result$results[[name]],start=start,frequency=f)
    result$finals$tables<-.x13_decorate_tables(result$finals$tables,start,f)
    for(name in c("finalSA","finalTrend","finalIrregular"))
      result[[name]]<-stats::ts(result[[name]],start=start,frequency=f)
    # Public residual trading-day F-tests (M7 cutover): computed from this
    # fit's own finalSA/finalIrregular on their published observed domain —
    # the dating the wrapper just stamped (ts calendar, else year 1) — on the
    # reference working scale: the fit's preadjustment log flag (the delivered
    # envelope's preadjustment$log) makes the attach log() the components for
    # log fits and leaves additive fits on the delivered levels (no-op).
    result$diagnostics<-.attach_td_tests(result$diagnostics,result$finalSA,result$finalIrregular,start,isTRUE(result$preadjustment$log))
  } else {
    result$records_only<-"the provenance echo is complete; the pinned stage refuses any run before that"
  }
  result$fitted<-identical(result$status,"available")
  structure(result,class="persephone_native_x13_fit")
}
