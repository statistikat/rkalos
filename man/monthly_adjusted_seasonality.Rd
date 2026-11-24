\name{monthly_adjusted_seasonality}
\alias{monthly_adjusted_seasonality}
\title{Residual seasonality tests on a model-scale adjusted series}
\description{
Computes the public all-years monthly Qs and seasonal F tests on the
model-scale seasonally adjusted observations of a fitted monthly product.
The fitted preprocessing model supplies the ordinary difference count
\code{d + bd} and its mean flag. X-13 products use the observed portion of
public D11, logged in multiplicative and log-additive modes; TRAMO-SEATS
products use the stochastic adjusted series from the initial
Wiener--Kolmogorov extraction.
}
\usage{
monthly_adjusted_seasonality(series, differencing_order, mean = TRUE)
}
\arguments{
\item{series}{A nonempty finite numeric vector of the model-scale adjusted
observations in chronological order.}
\item{differencing_order}{Nonnegative integer ordinary difference count
\code{d + bd} of the fitted preprocessing model.}
\item{mean}{Logical mean flag of the fitted preprocessing model.}
}
\value{
A list with the fields \code{qs} and \code{seasonal_f}. Each field reports
the corresponding test independently as a list with the statistic
\code{value} and the upper-tail \code{p_value}, or is \code{NULL} when the
test is absent. Fewer than fourteen stationary observations give an absent
F test. At exactly fourteen the saturated model returns undefined
statistics.
}
\examples{
monthly_adjusted_seasonality(as.numeric(AirPassengers), 2L)
}
