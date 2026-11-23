\name{tramo_transform_decision}
\alias{tramo_critical_value}
\alias{tramo_seasonality}
\alias{tramo_transform_decision}
\title{TRAMO critical value, seasonality tests and log/level decision}
\description{
Inspect the early decisions of the native TRAMO processor. These primitives
work independently of fitting or decomposing a time series.
}
\usage{
tramo_critical_value(observations)

tramo_seasonality(series)

tramo_transform_decision(series, log_preference, precision, mean = TRUE,
  detect_seasonality = TRUE, spec_seasonality = TRUE, robust = FALSE)
}
\arguments{
\item{observations}{A nonnegative integer observation count in R's integer range.}
\item{series}{A complete, finite, univariate monthly \code{ts} with at least
24 observations. The detector uses the order of the observations.}
\item{log_preference}{Finite offset added to the log-model criterion for the
comparison with the level-model criterion. Supply the value from the specification.}
\item{precision}{Positive finite estimation precision for the log/level module.}
\item{mean}{Logical mean flag from the specification. The log/level module
uses mean correction internally.}
\item{detect_seasonality}{Logical. Run the detector before the log/level comparison.}
\item{spec_seasonality}{Logical seasonal flag supplied by the specification.}
\item{robust}{Logical. This primitive accepts \code{FALSE}. Requests for
robust outlier pre-correction raise an error.}
}
\value{
\code{tramo_critical_value} returns the TRAMO outlier critical value as a double.
\code{tramo_seasonality} returns a list with the verdict, integer score
(\code{to_int}), Qs and Friedman test readouts (\code{value}, \code{p_value}
and \code{valid}), and the mean-removed first differences as a monthly series.
\code{tramo_transform_decision} returns the detector report when requested,
the pre-stage model string, seasonal and nonpositive flags, log/level status,
log and level criteria, correction and preference, the module's log choice,
the resulting transform flag, and chronological remarks.
}
\examples{
tramo_critical_value(length(AirPassengers))
tramo_seasonality(AirPassengers)$verdict
tramo_transform_decision(AirPassengers, log(0.95), 1e-5)$log_transform
}
