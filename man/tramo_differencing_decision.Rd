\name{tramo_differencing_decision}
\alias{tramo_differencing_decision}
\title{Run the first TRAMO differencing decision}
\usage{
tramo_differencing_decision(series,
  orders = c(p = 0L, d = 0L, q = 0L, bp = 0L, bd = 0L, bq = 0L),
  mean = TRUE, log = FALSE, seasonal = TRUE, initial_round = TRUE)
}
\arguments{
\item{series}{Complete finite monthly or quarterly time series after interpolation and removal of outlier effects. Supply the original scale when requesting a log transformation.}
\item{orders}{Named initial SARIMA orders. Regular orders range from zero through two and seasonal orders from zero through one. Accepts the named orders returned by \code{tramo_prepare_model()}.}
\item{mean}{Whether the initial model includes a mean.}
\item{log}{Whether to transform the input with the natural logarithm. Requires positive data.}
\item{seasonal}{Whether to consider seasonal differencing.}
\item{initial_round}{Whether the ladder starts at zero differencing orders. A subsequent round starts at the supplied orders.}
}
\description{Runs the pinned TRAMO differencing module with its standard threshold settings. It returns the decision before the separate controller battery.}
\value{A list containing the decision, resulting SARIMA orders and mean flag, mean-correction test and statistic, default-model, selection-failure and airline flags, and processing events. Unsupported native estimation states raise an error.}
\seealso{\code{\link{tramo_prepare_model}}}
