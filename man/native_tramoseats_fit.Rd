\name{native_tramoseats_fit}
\alias{native_tramoseats_fit}
\title{Fit a monthly RSA0 or RSA1 TRAMO--SEATS model from observations}
\usage{
native_tramoseats_fit(series, template = "RSA0", forecast = -2L, backcast = 0L)
}
\arguments{
\item{series}{Complete finite monthly univariate time series. RSA0 requires sufficient residual degrees of freedom. RSA1 requires at least 24 observations.}
\item{template}{Stock RSA0 or RSA1 template, also accepting TR0 and TR1 respectively, case-insensitively.}
\item{forecast, backcast}{Extension counts in months. Negative counts denote years.}
}
\description{Fits the airline model and mean from the observed series, extracts
the Wiener--Kolmogorov components and restores deterministic effects. RSA0
uses the level series. RSA1 selects the transformation and outliers using TR1.}
\value{An object of class \code{native_tramoseats_fit} with native \code{status},
\code{template}, \code{method}, \code{series}, fitted \code{model},
\code{likelihood} and \code{mode}. The \code{initial} and \code{final}
component tables include conditional extraction uncertainty. The \code{merged}
table contains the final product values, forecasts and backcasts, with one
dated monthly time series per row. On success \code{finalSA}, \code{finalTrend},
\code{finalIrregular} and \code{finalSeasonal} retain the input calendar.
Both templates expose the dated preprocessing residuals in \code{residuals$tsres}.

RSA0 also exposes \code{residuals$ser}, the residual standard error, the
24-lag Ljung--Box test in \code{residuals$lb}, and the Doornik--Hansen
normality test in \code{residuals$doornikhansen}. Both tests report \code{value},
\code{pvalue} and \code{df}. The Ljung--Box result is undefined for at most
24 residuals, and normality is absent for at most seven residuals. The dated
\code{series_forecast_stdev} gives RSA0 preprocessing forecast errors when
\code{series_forecast_stdev_status} is \code{"available"}.

RSA1 returns the model, likelihood and residuals of \code{native_tramo_preprocess},
plus \code{model_changed} indicating a SEATS model adjustment. Its fitted
mode is multiplicative for log and additive for level transformations. The
separate RSA1 original-series forecast error projection reports \code{"not_qualified"}
and an empty vector. Zero forecast months report \code{"no_requested_predictions"}.}
\seealso{\code{\link{native_tramo_preprocess}}, \code{\link{seats_from_tramo}}}
