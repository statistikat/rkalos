\name{seats_extract_wk}
\alias{seats_extract_wk}
\title{Extract SEATS components for a fixed monthly SARIMA model}
\description{
Runs the Burman/Wiener--Kolmogorov component extractor on a supplied fixed
model and innovation variance. The caller provides the fitted model parameters.
}
\usage{
seats_extract_wk(series, innovation_variance, d = 0L, bd = 0L,
  phi = numeric(), theta = numeric(), bphi = numeric(), btheta = numeric(),
  log = FALSE, mean_correction = TRUE, forecast = -2L, backcast = 0L,
  approximation = "legacy", method = "burman", bias = "legacy", xl = 0.95,
  epsphi = 2, rmod = 0.5, smod = 0.8, smod_pi = 0.8)
}
\arguments{
\item{series}{A complete, finite, univariate monthly \code{ts}, in original units.}
\item{innovation_variance}{The supplied positive innovation variance.}
\item{d, bd}{Regular and seasonal differencing orders.}
\item{phi, theta, bphi, btheta}{Fixed regular AR and MA coefficients and seasonal AR and MA coefficients, respectively, in Kalos order and sign convention.}
\item{log, mean_correction}{Logical flags for the log-transformed model and mean correction.}
\item{forecast, backcast}{Month counts. A negative count denotes that many years in absolute value.}
\item{approximation}{\code{"none"} or \code{"legacy"}. Legacy monthly retries cover seasonal MA replacement for
(0,1,q)(0,1,BQ) models with q at most one and seasonal AR-to-MA
re-estimation for mean-free (2,1,0)(1,1,0) models. Other required
approximations raise an error.}
\item{method}{The qualified component estimator is \code{"burman"}.}
\item{bias}{\code{"legacy"} or \code{"none"}.}
\item{xl}{Parameter cutoff used by the SEATS model validator.}
\item{epsphi, rmod, smod, smod_pi}{The SEATS decomposition specification's
seasonal tolerance, trend boundary, seasonal boundary and boundary at pi.}
}
\value{A list with the core status and two data frames, \code{initial} and
\code{final}. Each row identifies its component and information (value,
forecast, backcast or corresponding standard deviation), start year and
one-based month, and holds the numeric values in the \code{values} list column.
Nonempty values are monthly \code{ts} objects. A present empty component
retains its dated row and an empty vector. An unavailable component has no row.
Documented core decomposition errors carry their status and empty tables.
Unsupported requests raise an error.}
