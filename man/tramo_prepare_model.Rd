\name{tramo_prepare_model}
\alias{tramo_prepare_model}
\title{Prepare a monthly TRAMO template model}
\description{Initializes a native TRAMO model using a standard template and the
supplied time series, before selection, estimation and decomposition.}
\usage{tramo_prepare_model(series, template)}
\arguments{
\item{series}{A univariate monthly numeric \code{ts}. Missing observations
are interpolated by the core's average interpolator. At least one finite
observation is required.}
\item{template}{One of \code{TR0}, \code{TR1}, \code{TR2}, \code{TR3},
\code{TR4}, \code{TR5} and \code{TRfull}, or its lowercase spelling.}
}
\value{A list containing the canonical template name, automatic-model flag,
initial log and mean flags, SARIMA orders and period, one-based missing
observation positions, interpolated monthly series, registered variable names
and model-builder warnings. The initial flags describe the model at
construction. Automatic transform and model decisions occur later.}
\examples{
x <- AirPassengers
x[20] <- NA_real_
model <- tramo_prepare_model(x, "TR1")
model$missing
model$orders
}
