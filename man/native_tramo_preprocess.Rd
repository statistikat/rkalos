\name{native_tramo_preprocess}
\alias{native_tramo_preprocess}
\title{Fit the monthly TR1 preprocessing model}
\usage{
native_tramo_preprocess(series, template = "TR1")
}
\arguments{
\item{series}{A complete, finite, univariate monthly time series with at least 24 observations.}
\item{template}{The stock TR1 template. The case-insensitive alias RSA1 is accepted.}
}
\description{Fits the TR1 airline model, selects the log or level transformation
and estimates its mean and detected outlier coefficients from the original observations.}
\value{An object of class \code{native_tramo_preprocess} containing \code{status},
\code{template}, \code{model}, \code{likelihood} and a dated monthly
residual series in \code{residuals$tsres}. The residuals list also contains
the standard error (\code{ser}), the 24-lag Ljung--Box test (\code{lb})
and the Doornik--Hansen test (\code{doornikhansen}). The model contains its six orders,
fitted AR and MA parameters in regular AR, seasonal AR, regular MA and
seasonal MA order, mean coefficient and an outlier data frame. Outlier
positions are one-based indices on the input calendar. Likelihood fields
include the log likelihood, transformation adjustment, adjusted log
likelihood, residual sum of squares and effective observation and parameter
counts.}
\seealso{\code{\link{native_tramoseats_fit}}}
