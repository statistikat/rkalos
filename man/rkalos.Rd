\name{rkalos}
\alias{rkalos}
\alias{capabilities}
\alias{core_version}
\alias{prepare_inputs}
\alias{fixed_likelihood}
\alias{regression_likelihood}
\alias{arma_parameter_map}
\alias{arma_decode}
\alias{arma_objective_setup}
\alias{arma_objective}
\alias{arma_optimize}
\alias{arma_identification}
\alias{arma_information}
\alias{regression_information}
\alias{arma_residuals}
\alias{interpolation_information}
\alias{arma_forecasts}
\alias{prediction_moments}
\alias{x11_filter_spec}
\alias{x11_musgrave_spec}
\alias{x11_seasonal_spec}
\alias{x11_filter}
\alias{x11_seasonal_process}
\alias{x11_trend_policy}
\alias{x11_trend_select}
\alias{x11_seasonal_plan}
\alias{x11_msr_table}
\alias{x11_msr_policy}
\alias{x11_msr_pilot}
\alias{x11_msr_select}
\alias{x11_msr_resolve}
\alias{x11_seasonal_normalize}
\alias{x11_initial_additive}
\title{Experimental Kalos Engine Interfaces}
\description{Independent, Java-free time-series numerical interfaces. The exported \code{native_x13_fit} accepts a fixed monthly specification and \code{native_x13_fit_automatic} accepts decoder-equivalent automatic module blocks. Consult \code{capabilities()} for route requirements and refusals.}
\usage{
capabilities()
core_version()
prepare_inputs(...)
fixed_likelihood(...)
regression_likelihood(...)
arma_parameter_map(...)
arma_decode(...)
arma_objective_setup(...)
arma_objective(...)
arma_optimize(...)
arma_identification(...)
arma_information(...)
regression_information(...)
arma_residuals(...)
interpolation_information(...)
arma_forecasts(...)
prediction_moments(...)
x11_filter_spec(...)
x11_musgrave_spec(...)
x11_seasonal_spec(...)
x11_filter(...)
x11_seasonal_process(...)
x11_trend_policy(...)
x11_trend_select(...)
x11_seasonal_plan(...)
x11_msr_table(...)
x11_msr_policy(...)
x11_msr_pilot(...)
x11_msr_select(...)
x11_msr_resolve(...)
x11_seasonal_normalize(...)
x11_initial_additive(...)
}
\arguments{\item{...}{Arguments forwarded to the corresponding engine operation. See the getting-started guide in the package README and the numerical contracts in the Kalos repository.}}
\details{Inputs use monthly or quarterly univariate \code{ts} series. Regressors are named series or a named multivariate \code{ts}. Prepared objects and result lists preserve explicit unavailable statuses, domain metadata, missing masks and fixed/free parameter states. The accurate final-fit information policy remains the default; legacy comparison modes are explicitly requested through their operation's controls. Result classes retain their historical \code{persephone_native_*} names for serialization compatibility.}
\value{\code{capabilities()} returns the qualified operation inventory and seasonal-adjustment availability. Its \code{tramo_seats_fit} block lists the monthly templates accepted by \code{\link{native_tramoseats_fit}}. Its \code{seats_approximation} block lists the supported monthly legacy retry families. \code{core_version()} reports the embedded revision, provider and thread configuration. Numerical operations return named result lists, including their status and domain metadata; specification helpers return filter or model descriptions.}
\examples{
core_version()
y <- ts(100 + sin(seq_len(48)/3), start=c(2000,1), frequency=12)
inputs <- prepare_inputs(y)
fit <- fixed_likelihood(y, ar=0.3)
fit$likelihood$status
}
