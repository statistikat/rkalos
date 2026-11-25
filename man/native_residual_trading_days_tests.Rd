\name{native_residual_trading_days_tests}
\alias{native_residual_trading_days_tests}
\title{Residual trading-day F-tests for the SA and irregular components}
\description{
Serves the four public residual trading-day diagnostics of the JVM dashboard,
\code{diagnostics.td-sa-all}, \code{td-sa-last}, \code{td-i-all} and
\code{td-i-last}, natively. Each component is tested on its own: the
\code{sa} and \code{irregular} series run an ARIMA \code{(0,1,1)(0,0,0)[12]}
(one regular difference, one free MA(1) coefficient) fitted over the six
GenericTradingDays TD7 contrast columns as the only design (no mean), and the
pinned JointTest evaluates the concentrated likelihood as
\code{F = (ssq(b)/6)/(ssq/(dim-7))} with the \code{F(6, dim-7)} upper
probability. The pinned capture option is \code{flast = 0}: the whole
component is tested, so the \code{-all} and \code{-last} names of a component
share one fit and agree bit for bit. The product fit routes surface this test
as \code{diagnostics[["td.ftest.sa"]]} / \code{diagnostics[["td.ftest.i"]]}
from the components on the model's working (estimation) scale, transform-aware
— the input convention of the JVM's GenericSaTests, whose live
\code{diagnostics.td-*} fields these mirror: the X-13 routes apply \code{log}
to the decorated final SA/irregular when the fit's preadjustment log flag is
set (additive fits stay on the delivered levels), and the SEATS routes feed
the initial decomposition's working-scale Value rows. Bit-exactness layering:
the primitive is bit-exact against the 24-case double-JVM oracle for the
oracle's own bit-frozen component series
(fixtures/native/residual-trading-days/public-x86_64.txt, pinned by
\code{test-native-residual-trading-days.R}); component series composed by a
product fit are tolerance-qualified, so a product fit's scalars are
wiring-exact (identical to applying this primitive to the components the
route selects) rather than fixture-bit-exact. The oracle was captured under
the bridge flag profile.
}
\usage{
native_residual_trading_days_tests(sa, irregular, start_year, start_month)
}
\arguments{
\item{sa, irregular}{Nonempty finite numeric vectors: the seasonally-adjusted
and irregular component series, each on its own extraction domain. Both
series need more than eight observations for positive test degrees of
freedom.}
\item{start_year}{Whole-number calendar year (1-9999) of the first
observation; both series start at the same monthly domain start.}
\item{start_month}{Whole-number one-based month (1-12) of the first
observation.}
}
\value{A named list with the four entries \code{td-sa-all},
\code{td-sa-last}, \code{td-i-all} and \code{td-i-last}, each a
\code{list(value, pvalue)} carrying the F statistic and the upper-tail
p-value on \code{F(6, dim-7)} degrees of freedom, \code{dim} being the
differenced length \code{n-1}. With the pinned option \code{flast = 0} the
\code{-last} pair of a component is identical to its \code{-all} pair.
Invalid domains, non-finite or too-short series raise an error (the core
fails closed; there is no fallback).}
\examples{
air <- as.numeric(datasets::AirPassengers)
# Degenerate stand-in (both names get the same test here): real callers
# pass fitted SA and irregular components on their shared observed domain.
res <- native_residual_trading_days_tests(air, air, 1949, 1)
names(res)
res[["td-sa-all"]]
identical(res[["td-sa-all"]], res[["td-sa-last"]])
}
