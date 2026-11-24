#pragma once

#include <Rcpp.h>
#include <kalos/tramo_fit.hpp>

Rcpp::List wrap_tramo_preprocessing(const kalos::TramoFitResult& fit);
void attach_tramo_calendar(Rcpp::List& result, const kalos::TramoFittedCalendar& calendar);
Rcpp::List wrap_monthly_tramo_residuals(const std::vector<double>& values, double standard_error, int arma_parameters);
