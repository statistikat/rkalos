\name{seats_from_tramo}
\alias{seats_from_tramo}
\title{Extract and merge SEATS components from supplied fitted TRAMO coefficients}
\usage{
seats_from_tramo(series, ssq, effective, estimated, d = 0L, bd = 0L,
    phi = numeric(), theta = numeric(), bphi = numeric(), btheta = numeric(),
    log = FALSE, mean_correction = TRUE, outliers = NULL,
    forecast = -2L, backcast = 0L)
}
\arguments{
\item{series}{Complete original monthly time series covering the full estimation span.}
\item{ssq}{Positive fitted likelihood residual sum of squares.}
\item{effective, estimated}{Effective fitted observation count and estimated parameter count as recorded in the likelihood statistics.}
\item{d, bd}{Fitted regular and seasonal differencing orders.}
\item{phi, theta, bphi, btheta}{Fitted regular and seasonal AR and MA coefficients in Kalos convention.}
\item{log, mean_correction}{Fitted model's log transformation and mean flag.}
\item{outliers}{Optional data frame with character \code{type} (AO or LS), one-based observation \code{position}, and numeric fitted \code{coefficient}, in fitted model order. Effects are restored on the fitted level or log scale.}
\item{forecast, backcast}{Extension counts. Negative counts denote years.}
}
\description{Runs the fixed fitted-model bridge, Wiener--Kolmogorov extraction and deterministic-effects merge. The envelope requires full-span estimation with the supplied outliers and no additional calendar or user regressors, missing observations or length-of-period adjustment. Monthly legacy retries use the model families described in \code{\link{seats_extract_wk}}.}
\value{A list with native \code{status} and three data frames of dated components. \code{initial} and \code{final} include conditional standard deviations from SEATS. \code{merged} contains value and extension series on the original scale after restoring fitted outlier effects.}
\seealso{\code{\link{seats_extract_wk}}}
