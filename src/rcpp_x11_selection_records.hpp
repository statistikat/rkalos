#ifndef KALOS_RCPP_X11_SELECTION_RECORDS_HPP
#define KALOS_RCPP_X11_SELECTION_RECORDS_HPP
// Selection-result serialization shared by the pass-level selection boundary
// and the restricted X-13 fit composition: nested pass records keep the exact
// shapes of the fixed-pass boundaries in both, so a result's $b/$c chain into
// .native_x11_c_pass/.native_x11_d_pass unchanged wherever they arrive from.
#include "rcpp_x11_table_support.hpp"
#include "kalos/x11_selection_pass.hpp"
#include <string>
#include <vector>
namespace kalos_x11_rcpp {
// One pinned SeasonalFilterOption request per calendar period, named exactly
// as in the qualified seasonal-plan boundary.
std::vector<kalos::X11SeasonalRequest> requests(SEXP value);
const char* policy_name(kalos::X11SigmaPolicy policy);
kalos::X11SigmaPolicy parse_policy(const std::string& name,const char* label);
kalos::X11MsrIterationPolicy parse_msr_policy(const std::string& name);
double real(double value);
Rcpp::NumericVector reals(const std::vector<double>& values);
Rcpp::IntegerVector integers(const std::vector<int>& values);
// The B and C correction-round records share their fields, so one converter
// serves both (the 1-based row convention matches the fixed-pass boundaries).
template<class Round>
Rcpp::List correction_round(const Round& round) {
    auto plan=reals(round.replacements.replacements);
    int count=0;
    for(R_xlen_t i=0;i<plan.size();++i) if(!ISNAN(plan[i])) ++count;
    return Rcpp::List::create(Rcpp::_["detection"]=values(round.detection),
        Rcpp::_["analysis_position"]=static_cast<int>(round.analysis_position)+1, // 1-based rows
        Rcpp::_["policy"]=policy_name(round.policy),
        Rcpp::_["sigma_weights"]=reals(round.sigma.weights),Rcpp::_["replacements"]=plan,
        Rcpp::_["corrections"]=reals(round.replacements.correction_factors),
        Rcpp::_["replacement_count"]=count);
}
Rcpp::List msr_table(const kalos::X11MsrDiagnostics& msr);
Rcpp::List stage(const kalos::X11TrendStage& stage);
Rcpp::List plan(const kalos::X11SeasonalPlan& plan);
Rcpp::List b_result(const kalos::X11BPass& b,const std::string& policy,double lower,
    double upper,int frequency,int first);
Rcpp::List c_result(const kalos::X11CPass& c,const std::string& policy,double lower,
    double upper,int frequency);
Rcpp::List d_result(const kalos::X11DPass& d,const std::string& policy,int frequency,int first);
// The full nested record of one selection run; the embedding boundary supplies
// the request echoes (policy/msr names, mode, horizons, request names) so the
// record reads as the request that produced it, and its own scope label.
Rcpp::List selection_result(const kalos::X11SelectionPass& result,const std::string& policy,
    const std::string& msr_policy,const std::string& mode,const std::string& scope,
    double lower,double upper,int period,int first,int requested_trend_length,
    int backcast,int forecast,Rcpp::CharacterVector requests);
}
#endif
