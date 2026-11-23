#pragma once

#include <Rcpp.h>
#include <kalos/tramo_fit.hpp>

Rcpp::List wrap_tramo_preprocessing(const kalos::TramoFitResult& fit);
Rcpp::List wrap_monthly_airline_residuals(const std::vector<double>& values, double standard_error);
