\name{native_x13_fit}
\alias{native_x13_fit}
\alias{native_x13_fit_automatic}
\title{Fit a qualified monthly X-13 specification}
\description{Run the restricted fixed-coefficient or automatic X-13 route.
Inspect \code{capabilities()$restricted_fit} and \code{capabilities()$automatic_fit}
for the respective support envelopes.}
\usage{
native_x13_fit(series, variables = list(), log = FALSE, seasonal, trend_length,
  sigma = NULL, msr_policy = "jdplus", frequency = NULL, start_period = NULL,
  lower_sigma = 1.5, upper_sigma = 2.5, mode = "additive",
  backcast = 0, forecast = 0)

native_x13_fit_automatic(series, transform, automodel, outliers, ami, regression,
  sigma = NULL, lower_sigma = 1.5, upper_sigma = 2.5, msr_policy = "jdplus",
  mode = "additive", backcast = 0, forecast = 0, frequency = NULL,
  start_period = NULL, first_year = NULL, calendar = FALSE, easter = FALSE,
  trading_days = FALSE, user_variables = FALSE)
}
\arguments{
\item{series}{Univariate numeric time series or numeric vector.}
\item{frequency,start_period,first_year}{Calendar frequency, one-based initial period,
and initial year. Supply the calendar for a plain vector. A time series supplies
its own calendar. Explicit values on the automatic route must agree with it.}
\item{variables,log}{List of restricted fixed-coefficient variables, and logical
log-transformation flag. Each variable supplies its class, effect, coefficient
vector, and a list of corresponding observed-length columns.}
\item{seasonal,trend_length}{One fixed seasonal filter per calendar period and
a fixed odd Henderson length of at least five for the restricted route.}
\item{transform,automodel,outliers,ami,regression}{Named lists of required automatic
module parameters. The automatic route takes decoder-equivalent module state.
See the automatic pipeline contract for the numerical meaning of each module.}
\item{sigma,lower_sigma,upper_sigma}{Weight policy (\code{"legacy"} by default,
or \code{"accurate"}) and sigma limits.}
\item{msr_policy,mode}{MSR iteration policy and requested decomposition mode.}
\item{backcast,forecast}{On the automatic route, months to extend before and
after the input. Negative values count years. The restricted route requires zero.}
\item{calendar,easter,trading_days,user_variables}{Logical detection flags.
Enabling a deferred automatic detection module returns an unsupported status.}
}
\value{A named result list carrying the core status, request provenance and,
when available, the preadjustment, selection, results, finals and in-sample
seasonally adjusted series, trend and irregular. Automatic fits additionally
return their estimated model and decision chain. For positive effective
horizons, \code{preadjustment$tables} includes the original-scale raw
forecasts (\code{a1a}) or backcasts (\code{a1b}). The corresponding
\code{finals$tables} entries are \code{d11a}, \code{d12a}, \code{d16a},
\code{d18a} for forecasts, and the same roots with suffix \code{b} for
backcasts. Each is a monthly \code{ts} on its own domain. Zero-length windows
omit these optional names. Inspect \code{status} before
using fitted quantities. Invalid input shapes raise errors.}
\seealso{\code{\link{capabilities}}}
\examples{
x <- native_x13_fit(AirPassengers, seasonal = rep("s3x3", 12), trend_length = 13L)
x$status
}
