#include <Rcpp.h>
#include "kalos/residual_trading_days_tests.hpp"
#include <stdexcept>
#include <vector>
// The qualified oracle surface runs BLAS kernels single-threaded: OpenBLAS
// reduces wide kernels in schedule-dependent order, and inside an R process
// a second OpenBLAS (R's own BLAS) can oversubscribe the pool; the pinned
// JVM oracle and C++ qualification are both measured per-kernel serial.
extern "C" int openblas_get_num_threads(void);
extern "C" void openblas_set_num_threads(int);
namespace {
struct ScopedSingleBlas {
    const int previous = openblas_get_num_threads();
    ScopedSingleBlas() { openblas_set_num_threads(1); }
    ~ScopedSingleBlas() { openblas_set_num_threads(previous); }
};
Rcpp::List statistic(const kalos::TradingDaysFTest& test) {
    return Rcpp::List::create(Rcpp::_["value"] = test.statistic, Rcpp::_["pvalue"] = test.p_value);
}
}
// The four public residual trading-day F-tests (M7 item 27 cutover gap):
// diagnostics.td-sa-all / td-sa-last / td-i-all / td-i-last of the JVM
// dashboard, served natively by kalos::residual_trading_days_tests from each
// component series and its monthly domain start. The pinned capture option
// is flast = 0, where the "-all" and "-last" names of a component share one
// fit and agree bit for bit (the core returns the same pair under both
// names). The R wrapper validates the marshalling; every domain and
// degrees-of-freedom guard fails closed from the core, never a fallback.
// [[Rcpp::export(name = ".native_residual_trading_days_tests", rng = false)]]
Rcpp::List native_residual_trading_days_tests(const std::vector<double>& sa,
        const std::vector<double>& irregular, int start_year, int start_month) {
    if (sa.empty() || irregular.empty())
        throw std::invalid_argument("the residual trading-day test needs both component series");
    const ScopedSingleBlas worker;
    const auto tests = kalos::residual_trading_days_tests(sa, start_year, start_month,
                                                           irregular, start_year, start_month);
    return Rcpp::List::create(Rcpp::_["td-sa-all"] = statistic(tests.sa_all),
        Rcpp::_["td-sa-last"] = statistic(tests.sa_last),
        Rcpp::_["td-i-all"] = statistic(tests.irregular_all),
        Rcpp::_["td-i-last"] = statistic(tests.irregular_last));
}
