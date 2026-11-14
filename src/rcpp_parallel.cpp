#include <Rcpp.h>
#include "kalos/parallel.hpp"

// [[Rcpp::export(name = ".native_parallel_info")]]
Rcpp::List native_parallel_info() {
    const auto info = kalos::parallel_info();
    const auto actual = kalos::parallel_probe_threads(info);
    return Rcpp::List::create(
        Rcpp::Named("compiled") = info.compiled,
        Rcpp::Named("openmp_version") = info.version,
        Rcpp::Named("requested_threads") = info.requested_threads,
        Rcpp::Named("runtime_max_threads") = info.runtime_max_threads,
        Rcpp::Named("selected_threads") = info.selected_threads,
        Rcpp::Named("probe_threads") = actual,
        Rcpp::Named("parallel_enabled") = actual > 1,
        Rcpp::Named("work_threshold") = static_cast<double>(kalos::parallel_work_threshold),
        Rcpp::Named("chunk_work") = static_cast<double>(kalos::parallel_chunk_work));
}
