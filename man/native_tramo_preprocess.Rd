\name{native_tramo_preprocess}
\alias{native_tramo_preprocess}
\title{Fit a monthly TR1 or TR3 preprocessing model}
\usage{
native_tramo_preprocess(series, template = "TR1")
}
\arguments{
\item{series}{A complete, finite, univariate monthly time series with at least 24 observations.}
\item{template}{TR1 or TR3, also accepting RSA1 or RSA3 case-insensitively.}
}
\description{Selects log or level and estimates coefficients from the original observations. TR1 fits the airline orders. TR3 identifies the orders automatically through the first reference comparison and subsequent controller sweep. The seasonal controller can restart identification and the fast mean test can request a fourth estimation pass. An additional round after the first reference comparison raises an explicit error.}
\value{An object of class \code{native_tramo_preprocess} containing \code{status},
\code{template}, \code{model}, \code{likelihood} and a dated monthly
residual series in \code{residuals$tsres}. For both templates the residuals list also contains
the standard error (\code{ser}), the 24-lag Ljung--Box test (\code{lb})
and the Doornik--Hansen test (\code{doornikhansen}). The model contains its six orders,
fitted AR and MA parameters in regular AR, seasonal AR, regular MA and
seasonal MA order, mean coefficient and an outlier data frame. Outlier
positions are one-based indices on the input calendar. Likelihood fields
include the log likelihood, transformation adjustment, adjusted log
likelihood, residual sum of squares and effective observation and parameter
counts.}
\seealso{\code{\link{native_tramoseats_fit}}}
