\name{native_x13_estimated_fit}
\alias{native_x13_estimated_fit}
\title{Estimate a stock monthly X-13 model}
\usage{
native_x13_estimated_fit(series, template = "RSA0")
}
\arguments{
\item{series}{A complete finite univariate monthly time series.}
\item{template}{The stock RSA0, RSA1 or RSA2c template, case-insensitively.}
}
\description{Estimates the regular and seasonal moving-average coefficients
of the stock airline model from the original observations. RSA0 uses level
data; RSA1 selects the transformation and outliers. RSA2c also tests working days, leap year and Easter. X-11 automatically
selects the Henderson and seasonal filters. The stock extension is twelve
forecast months and zero backcast months.}
\value{A \code{native_x13_estimated_fit} list with the template, core
status, fitted ARIMA model and estimated regression variables. The
\code{likelihood} list reports the log likelihood, residual sum of squares
and observation and parameter counts; \code{residuals$tsres} contains the
dated monthly innovations. RSA2c reports its selected calendar and named coefficients in \code{model$calendar}. Successful
fits expose \code{series}, \code{finalSA}, \code{finalTrend} and
\code{finalIrregular} on the input calendar. The \code{preadjustment},
\code{selection}, \code{results} and \code{finals} records contain the
X-11 decisions and tables, including the original-series forecast
\code{preadjustment$tables$a1a} and the dated final forecasts.}
\details{The combined SI diagnostic in \code{diagnostics$combined_seasonality} reports its verdict and three tests on the observed D8 domain.}
\seealso{\code{\link{native_x13_fit}}, \code{\link{native_x13_fit_automatic}}}
\examples{
fit <- native_x13_estimated_fit(AirPassengers, "RSA0")
fit$finalSA
}
