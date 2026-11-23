\name{native_tramoseats_fit}
\alias{native_tramoseats_fit}
\title{Estimate and decompose a monthly series using native TRAMO/SEATS RSA0}
\usage{
native_tramoseats_fit(series, template = "RSA0", forecast = -2L, backcast = 0L)
}
\arguments{
\item{series}{Complete finite monthly univariate time series with sufficient residual degrees of freedom for airline estimation.}
\item{template}{The RSA0 template, also called TR0. Case is ignored.}
\item{forecast, backcast}{Extension lengths in months. A negative count denotes years. The default requests two years of forecasts and no backcasts.}
}
\description{Fits the fixed-structure airline model with estimated regular and seasonal MA coefficients and mean correction, then extracts and merges the SEATS components. The fit begins with only the observed series and estimation dates. Approximation-required models and unsupported templates raise errors.}
\value{An object of class \code{native_tramoseats_fit} carrying additive \code{mode}, native \code{status}, canonical \code{template}, \code{method}, input \code{series}, fitted \code{model} and \code{likelihood}. The \code{initial} and \code{final} component tables include conditional WK uncertainty, and \code{merged} contains product values by component and information type. The preprocessing \code{residuals} list contains the dated monthly residual series \code{tsres} and scalar residual standard error \code{ser}. The dated \code{series_forecast_stdev} gives preprocessing forecast standard deviations for the original series. Its separate \code{series_forecast_stdev_status} reports availability. An empty vector accompanies zero horizons or an unavailable projection. On success, \code{finalSA}, \code{finalTrend}, \code{finalIrregular}, and \code{finalSeasonal} are observed-domain time series.}
\seealso{\code{\link{seats_from_tramo}}}
