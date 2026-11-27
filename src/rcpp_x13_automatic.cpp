#include "rcpp_x13_shared.hpp"
#include "rcpp_combined_seasonality.hpp"
#include "rcpp_monthly_residuals.hpp"
#include "rcpp_x11_selection_records.hpp"
#include "kalos/x13_automatic.hpp"
#include <cmath>
#include <limits>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
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

// The REQUIRED sigma policy arrives as NULL (unset — the fit then answers the
// pinned fail-closed status, never a default) or the policy name string.
std::optional<kalos::X11SigmaPolicy> sigma_policy(SEXP value) {
    if(Rf_isNull(value))return std::nullopt;
    if(TYPEOF(value)!=STRSXP || Rf_xlength(value)!=1 || STRING_ELT(value,0)==NA_STRING ||
        !Rf_isNull(Rf_getAttrib(value,R_DimSymbol)))
        throw std::invalid_argument("x13 automatic fit sigma weight policy must be legacy or accurate");
    return kalos_x11_rcpp::parse_policy(Rcpp::as<std::string>(value),
        "x13 automatic fit sigma weight policy");
}
int integer_at(SEXP value, R_xlen_t index, const char* name, int minimum) {
    if ((TYPEOF(value) != INTSXP && TYPEOF(value) != REALSXP) || Rf_isFactor(value) ||
        !Rf_isNull(Rf_getAttrib(value, R_DimSymbol)))
        Rcpp::stop("%s must contain whole numbers", name);
    const double number = TYPEOF(value) == REALSXP ? REAL(value)[index] : INTEGER(value)[index];
    if (!std::isfinite(number) || std::floor(number) != number ||
        number < minimum || number > std::numeric_limits<int>::max())
        Rcpp::stop("%s must contain whole numbers in the supported range", name);
    return static_cast<int>(number);
}
int integer_scalar(SEXP value, const char* name, int minimum = -std::numeric_limits<int>::max()) {
    if (Rf_xlength(value) != 1) Rcpp::stop("%s must be a scalar", name);
    return integer_at(value, 0, name, minimum);
}
bool flag(SEXP value,const char* label) {
    if(TYPEOF(value)!=LGLSXP || Rf_xlength(value)!=1 || LOGICAL(value)[0]==NA_LOGICAL ||
        !Rf_isNull(Rf_getAttrib(value,R_DimSymbol)))
        throw std::invalid_argument(label);
    return LOGICAL(value)[0]==TRUE;
}
// The pinned factory codes of the outlier type set.
kalos::OutlierType outlier_type(SEXP value) {
    const auto name=Rcpp::as<std::string>(value);
    if(name=="AO")return kalos::OutlierType::Additive;
    if(name=="LS")return kalos::OutlierType::LevelShift;
    if(name=="TC")return kalos::OutlierType::TransitoryChange;
    if(name=="SO")return kalos::OutlierType::SeasonalOutlier;
    throw std::invalid_argument("x13 automatic fit outlier type must be AO, LS, TC or SO");
}
std::string outlier_type_name(kalos::OutlierType type) {
    switch(type) {
        case kalos::OutlierType::Additive:return "AO";
        case kalos::OutlierType::LevelShift:return "LS";
        case kalos::OutlierType::TransitoryChange:return "TC";
        case kalos::OutlierType::SeasonalOutlier:return "SO";
    }
    return "?";
}
Rcpp::LogicalVector flags(const std::vector<bool>& values) {
    Rcpp::LogicalVector out(static_cast<R_xlen_t>(values.size()));
    for(std::size_t i=0;i<values.size();++i)
        out[static_cast<R_xlen_t>(i)]=values[i] ? TRUE : FALSE;
    return out;
}
// The pinned scalar echo pairs: a named numeric vector, insertion order.
Rcpp::NumericVector scalars_vector(const std::vector<std::pair<std::string,double>>& scalars) {
    Rcpp::NumericVector values(static_cast<R_xlen_t>(scalars.size()));
    Rcpp::CharacterVector names(static_cast<R_xlen_t>(scalars.size()));
    for(std::size_t i=0;i<scalars.size();++i) {
        values[static_cast<R_xlen_t>(i)]=scalars[i].second;
        names[static_cast<R_xlen_t>(i)]=scalars[i].first;
    }
    values.attr("names")=names;
    return values;
}
Rcpp::List branch_report(const kalos::TransformationBranchReport& branch) {
    return Rcpp::List::create(Rcpp::_["available"]=branch.available,
        Rcpp::_["log_likelihood"]=branch.log_likelihood,Rcpp::_["jacobian"]=branch.jacobian,
        Rcpp::_["adjusted_log_likelihood"]=branch.adjusted_log_likelihood,
        Rcpp::_["observations"]=static_cast<double>(branch.observations),
        Rcpp::_["parameters"]=static_cast<double>(branch.parameters),
        Rcpp::_["aicc"]=branch.aicc);
}
Rcpp::List transformation_payload(const kalos::TransformationDecisionResult& r) {
    return Rcpp::List::create(Rcpp::_["outcome"]=r.outcome,Rcpp::_["log_chosen"]=r.log_chosen,
        Rcpp::_["aicc_level"]=r.aicc_level,Rcpp::_["aicc_log"]=r.aicc_log,
        Rcpp::_["level_branch"]=branch_report(r.level),
        Rcpp::_["log_branch"]=branch_report(r.log));
}
Rcpp::List mean_payload(const kalos::MeanControllerResult& r) {
    return Rcpp::List::create(Rcpp::_["changed"]=r.changed,Rcpp::_["mean"]=r.mean,
        Rcpp::_["tstat"]=r.tstat);
}
Rcpp::List outliers_payload(const kalos::OutlierScoringResult& r) {
    const R_xlen_t n=static_cast<R_xlen_t>(r.outliers.size());
    Rcpp::IntegerVector position(n),type(n);
    Rcpp::CharacterVector name(n);
    Rcpp::NumericVector t_stat(n),coef(n);
    for(std::size_t i=0;i<r.outliers.size();++i) {
        const auto& d=r.outliers[i];
        position[static_cast<R_xlen_t>(i)]=static_cast<int>(d.position)+1; // 1-based
        type[static_cast<R_xlen_t>(i)]=static_cast<int>(d.type);
        name[static_cast<R_xlen_t>(i)]=outlier_type_name(d.type);
        t_stat[static_cast<R_xlen_t>(i)]=d.t_statistic;
        coef[static_cast<R_xlen_t>(i)]=d.coefficient;
    }
    const R_xlen_t m=static_cast<R_xlen_t>(r.rejections.size());
    Rcpp::IntegerVector reject_position(m),reject_type(m),reject_round(m);
    Rcpp::CharacterVector reject_type_name(m),reject_origin(m);
    Rcpp::NumericVector reject_abs_t(m);
    for(std::size_t i=0;i<r.rejections.size();++i) {
        const auto& j=r.rejections[i];
        const R_xlen_t k=static_cast<R_xlen_t>(i);
        reject_position[k]=static_cast<int>(j.position)+1; // 1-based
        reject_type[k]=static_cast<int>(j.type);
        reject_type_name[k]=outlier_type_name(j.type);
        reject_abs_t[k]=j.abs_t_statistic;
        reject_round[k]=static_cast<int>(j.round);
        reject_origin[k]=j.origin==kalos::OutlierRejection::Origin::AddRollback ?
            "add_rollback" : "verify_model";
    }
    const std::size_t sentinel=kalos::OutlierScanResult::max_type_sentinel;
    const R_xlen_t s=static_cast<R_xlen_t>(r.scans.size());
    Rcpp::IntegerVector scan_position(s),scan_type(s);
    Rcpp::NumericVector scan_t(s);
    Rcpp::LogicalVector scan_beyond(s);
    for(std::size_t i=0;i<r.scans.size();++i) {
        const auto& x=r.scans[i];
        const R_xlen_t k=static_cast<R_xlen_t>(i);
        scan_position[k]=static_cast<int>(x.position)+1; // 1-based
        scan_type[k]=x.type==sentinel ? NA_INTEGER : static_cast<int>(x.type);
        scan_t[k]=x.t_stat;
        scan_beyond[k]=x.beyond_cv ? TRUE : FALSE;
    }
    return Rcpp::List::create(Rcpp::_["changed"]=r.changed,Rcpp::_["stop_reason"]=r.stop_reason,
        Rcpp::_["rounds"]=static_cast<double>(r.rounds),
        Rcpp::_["decisions"]=Rcpp::List::create(Rcpp::_["position"]=position,
            Rcpp::_["type_index"]=type,Rcpp::_["type"]=name,
            Rcpp::_["t_statistic"]=t_stat,Rcpp::_["coefficient"]=coef),
        Rcpp::_["rejections"]=Rcpp::List::create(Rcpp::_["position"]=reject_position,
            Rcpp::_["type_index"]=reject_type,Rcpp::_["type"]=reject_type_name,
            Rcpp::_["abs_t_statistic"]=reject_abs_t,Rcpp::_["round"]=reject_round,
            Rcpp::_["origin"]=reject_origin),
        Rcpp::_["scans"]=Rcpp::List::create(Rcpp::_["position"]=scan_position,
            Rcpp::_["type_index"]=scan_type,Rcpp::_["t_stat"]=scan_t,
            Rcpp::_["beyond_cv"]=scan_beyond),
        Rcpp::_["last_scan_max_t"]=r.last_scan_max_t);
}
Rcpp::List regtest_payload(const kalos::RegressionVariableTestResult& r) {
    return Rcpp::List::create(Rcpp::_["changed"]=r.changed,Rcpp::_["mean"]=r.mean,
        Rcpp::_["preadjustment"]=r.preadjustment,Rcpp::_["removed"]=flags(r.removed),
        Rcpp::_["tstats"]=kalos_x13_rcpp::reals_out(r.tstats),
        Rcpp::_["tested"]=flags(r.tested),Rcpp::_["accepted"]=flags(r.accepted),
        Rcpp::_["derived_tested"]=flags(r.derived_tested),
        Rcpp::_["derived_tstats"]=kalos_x13_rcpp::reals_out(r.derived_tstats));
}
Rcpp::List accept_payload(const kalos::ModelControllerResult& r) {
    return Rcpp::List::create(Rcpp::_["accepted"]=r.accepted,Rcpp::_["lb_built"]=r.lb_built,
        Rcpp::_["lb_length"]=r.lb_length,Rcpp::_["lb_value"]=r.lb_value,Rcpp::_["lb_p"]=r.lb_p,
        Rcpp::_["rvr"]=r.rvr,Rcpp::_["rtval"]=r.rtval);
}
Rcpp::List refine_payload(const kalos::ArmaRefinement& r) {
    const R_xlen_t n=static_cast<R_xlen_t>(r.preferred_models.size());
    Rcpp::IntegerVector p(n),q(n),bp(n),bq(n);
    Rcpp::NumericVector bic(n),loglik(n);
    for(std::size_t i=0;i<r.preferred_models.size();++i) {
        const auto& m=r.preferred_models[i];
        p[static_cast<R_xlen_t>(i)]=m.p;q[static_cast<R_xlen_t>(i)]=m.q;
        bp[static_cast<R_xlen_t>(i)]=m.bp;bq[static_cast<R_xlen_t>(i)]=m.bq;
        bic[static_cast<R_xlen_t>(i)]=m.bic;loglik[static_cast<R_xlen_t>(i)]=m.log_likelihood;
    }
    Rcpp::Nullable<Rcpp::List> differencing=Rcpp::Nullable<Rcpp::List>(R_NilValue);
    if(r.differencing) {
        const auto& d=*r.differencing;
        differencing=Rcpp::Nullable<Rcpp::List>(Rcpp::List::create(Rcpp::_["d"]=d.d,
            Rcpp::_["seasonal_d"]=d.seasonal_d,Rcpp::_["mean_correction"]=d.mean_correction,
            Rcpp::_["t_mean"]=d.t_mean,Rcpp::_["ml_used"]=d.ml_used,
            Rcpp::_["nextstep_calls"]=d.nextstep_calls));
    }
    return Rcpp::List::create(Rcpp::_["p"]=r.p,Rcpp::_["q"]=r.q,Rcpp::_["bp"]=r.bp,
        Rcpp::_["bq"]=r.bq,Rcpp::_["d"]=r.d,Rcpp::_["seasonal_d"]=r.seasonal_d,
        Rcpp::_["mean"]=r.mean,Rcpp::_["specification_changed"]=r.specification_changed,
        Rcpp::_["mean_changed"]=r.mean_changed,Rcpp::_["changed"]=r.changed,
        Rcpp::_["fallback"]=r.fallback,Rcpp::_["fallback_stage"]=r.fallback_stage,
        Rcpp::_["differencing"]=differencing,
        Rcpp::_["preferred_models"]=Rcpp::List::create(Rcpp::_["p"]=p,Rcpp::_["q"]=q,
            Rcpp::_["bp"]=bp,Rcpp::_["bq"]=bq,Rcpp::_["bic"]=bic,
            Rcpp::_["log_likelihood"]=loglik));
}
Rcpp::List pass2_payload(const kalos::AutoPass2Record& r) {
    return Rcpp::List::create(Rcpp::_["filter"]=r.filter,Rcpp::_["restored"]=r.restored,
        Rcpp::_["continue_round"]=r.continue_round,Rcpp::_["plbox"]=r.plbox,
        Rcpp::_["plbox0"]=r.plbox0,Rcpp::_["rvr"]=r.rvr,Rcpp::_["rtval"]=r.rtval,
        Rcpp::_["pcr_after"]=r.pcr_after,Rcpp::_["loop_after"]=r.loop_after,
        Rcpp::_["curva_after"]=r.curva_after,Rcpp::_["curva_reduced"]=r.curva_reduced,
        Rcpp::_["need_outliers_after"]=r.need_outliers_after,
        Rcpp::_["need_auto_modelling_after"]=r.need_auto_modelling_after,
        Rcpp::_["last_solution"]=r.last_solution);
}
Rcpp::List unitroots_payload(const kalos::AutoUnitRootsRecord& r) {
    return Rcpp::List::create(Rcpp::_["regular_count"]=r.regular_count,
        Rcpp::_["seasonal_count"]=r.seasonal_count,Rcpp::_["regular_fired"]=r.regular_fired,
        Rcpp::_["seasonal_fired"]=r.seasonal_fired,Rcpp::_["p_after"]=r.p_after,
        Rcpp::_["d_after"]=r.d_after,Rcpp::_["bp_after"]=r.bp_after,
        Rcpp::_["bd_after"]=r.bd_after,Rcpp::_["regular_margin"]=r.regular_margin,
        Rcpp::_["seasonal_margin"]=r.seasonal_margin,
        Rcpp::_["redo_estimation"]=r.redo_estimation);
}
Rcpp::List matest_payload(const kalos::AutoMaTestRecord& r) {
    return Rcpp::List::create(Rcpp::_["scanned"]=r.scanned,
        Rcpp::_["abs_theta_at_1"]=r.abs_theta_at_1,Rcpp::_["fired"]=r.fired,
        Rcpp::_["q_after"]=r.q_after,Rcpp::_["d_after"]=r.d_after,
        Rcpp::_["mean_set"]=r.mean_set,Rcpp::_["redo_estimation"]=r.redo_estimation);
}
Rcpp::List meanrule_payload(const kalos::AutoMeanRuleRecord& r) {
    return Rcpp::List::create(Rcpp::_["rtval0"]=r.rtval0,Rcpp::_["mean_added"]=r.mean_added);
}
Rcpp::List cascade_round(const kalos::AutoCascadeRoundRecord& r) {
    return Rcpp::List::create(Rcpp::_["round"]=r.round,
        Rcpp::_["estimation_status"]=r.estimation_status,
        Rcpp::_["estimation_ok"]=r.estimation_ok,Rcpp::_["fi_status"]=r.fi_status,
        Rcpp::_["fx"]=r.fx,Rcpp::_["nobs"]=r.nobs,Rcpp::_["p_before"]=r.p_before,
        Rcpp::_["bp_before"]=r.bp_before,Rcpp::_["q_before"]=r.q_before,
        Rcpp::_["bq_before"]=r.bq_before,Rcpp::_["p_after"]=r.p_after,
        Rcpp::_["bp_after"]=r.bp_after,Rcpp::_["q_after"]=r.q_after,
        Rcpp::_["bq_after"]=r.bq_after,Rcpp::_["status"]=r.status,
        Rcpp::_["asy_diagonal"]=kalos_x13_rcpp::reals_out(r.asy_diagonal),
        Rcpp::_["cascade_ok"]=r.cascade_ok);
}
// One decision-chain entry: the stage/round tags, the pinned scalar echoes and
// the consumed unit's qualified record (R NULL when the stage carries no
// payload). The vocabulary mirrors the chain struct of
// docs/contracts/x13-automatic-pipeline-contract.md verbatim.
Rcpp::List decision_record(const kalos::AutoDecisionRecord& entry) {
    Rcpp::Nullable<Rcpp::List> transformation=Rcpp::Nullable<Rcpp::List>(R_NilValue);
    if(entry.transformation)transformation=Rcpp::Nullable<Rcpp::List>(
        transformation_payload(*entry.transformation));
    Rcpp::Nullable<Rcpp::List> mean=Rcpp::Nullable<Rcpp::List>(R_NilValue);
    if(entry.mean)mean=Rcpp::Nullable<Rcpp::List>(mean_payload(*entry.mean));
    Rcpp::Nullable<Rcpp::List> outliers=Rcpp::Nullable<Rcpp::List>(R_NilValue);
    if(entry.outliers)outliers=Rcpp::Nullable<Rcpp::List>(outliers_payload(*entry.outliers));
    Rcpp::Nullable<Rcpp::List> regtest=Rcpp::Nullable<Rcpp::List>(R_NilValue);
    if(entry.regtest)regtest=Rcpp::Nullable<Rcpp::List>(regtest_payload(*entry.regtest));
    Rcpp::Nullable<Rcpp::List> accept=Rcpp::Nullable<Rcpp::List>(R_NilValue);
    if(entry.accept)accept=Rcpp::Nullable<Rcpp::List>(accept_payload(*entry.accept));
    Rcpp::Nullable<Rcpp::List> refine=Rcpp::Nullable<Rcpp::List>(R_NilValue);
    if(entry.refine)refine=Rcpp::Nullable<Rcpp::List>(refine_payload(*entry.refine));
    Rcpp::Nullable<Rcpp::List> pass2=Rcpp::Nullable<Rcpp::List>(R_NilValue);
    if(entry.pass2)pass2=Rcpp::Nullable<Rcpp::List>(pass2_payload(*entry.pass2));
    Rcpp::Nullable<Rcpp::List> unitroots=Rcpp::Nullable<Rcpp::List>(R_NilValue);
    if(entry.unitroots)unitroots=Rcpp::Nullable<Rcpp::List>(unitroots_payload(*entry.unitroots));
    Rcpp::Nullable<Rcpp::List> matest=Rcpp::Nullable<Rcpp::List>(R_NilValue);
    if(entry.matest)matest=Rcpp::Nullable<Rcpp::List>(matest_payload(*entry.matest));
    Rcpp::Nullable<Rcpp::List> meanrule=Rcpp::Nullable<Rcpp::List>(R_NilValue);
    if(entry.meanrule)meanrule=Rcpp::Nullable<Rcpp::List>(meanrule_payload(*entry.meanrule));
    Rcpp::List cascade;
    for(const auto& round : entry.cascade)cascade.push_back(cascade_round(round));
    return Rcpp::List::create(Rcpp::_["stage"]=entry.stage,Rcpp::_["round"]=entry.round,
        Rcpp::_["scalars"]=scalars_vector(entry.scalars),
        Rcpp::_["notes"]=Rcpp::CharacterVector(entry.notes.begin(),entry.notes.end()),
        Rcpp::_["transformation"]=transformation,Rcpp::_["mean"]=mean,
        Rcpp::_["outliers"]=outliers,Rcpp::_["regtest"]=regtest,Rcpp::_["accept"]=accept,
        Rcpp::_["refine"]=refine,Rcpp::_["pass2"]=pass2,Rcpp::_["unitroots"]=unitroots,
        Rcpp::_["matest"]=matest,Rcpp::_["meanrule"]=meanrule,Rcpp::_["cascade"]=cascade);
}

Rcpp::List wrap_estimated_x13(Rcpp::List out, const kalos::X13AutomaticFitResult& result,
        const Rcpp::CharacterVector& seasonal, const std::string& msr_name, const std::string& record_scope) {
    out["status"]=result.status;
    out["effective_backcast"] = result.fit.backcast;
    out["effective_forecast"] = result.fit.forecast;
    out["mode_requested"]=kalos_x11_rcpp::mode_name(result.fit.mode_requested);
    if(!result.chain.empty()) { // the pinned fail-closed matrix refuses before any state
        Rcpp::List chain;
        for(const auto& entry : result.chain)chain.push_back(decision_record(entry));
        out["chain"]=chain;
    }
    if(!result.fit.preadjustment.status.empty()) { // past the fail-closed guards
        out["mode_effective"]=kalos_x11_rcpp::mode_name(result.fit.mode_effective);
        out["sigma_set"]=result.fit.sigma_set;
        out["sigma"]=Rcpp::wrap(std::string(kalos_x11_rcpp::policy_name(result.fit.sigma)));
        out["lower_sigma"]=result.fit.lower_sigma;
        out["upper_sigma"]=result.fit.upper_sigma;
        out["preadjustment"]=kalos_x13_rcpp::preadjustment_list(result.fit.preadjustment,result.fit.log);
    }
    if(!result.fit.selection.status.empty()) { // the run over out, same record shape
        const int f = result.fit.frequency;
        const int extended_period = ((result.fit.first_period - 1 - result.fit.backcast) % f + f) % f + 1;
        out["selection"]=kalos_x11_rcpp::selection_result(result.fit.selection,
            kalos_x11_rcpp::policy_name(result.fit.sigma),msr_name,
            kalos_x11_rcpp::mode_name(result.fit.mode_effective),
            record_scope + "_selection_run",result.fit.lower_sigma,result.fit.upper_sigma,
            f,extended_period,result.fit.trend_length,
            result.fit.backcast,result.fit.forecast,
            seasonal);
    }
    if(result.status=="available") { // the final description and the tables
        out["diagnostics"] = wrap_si_diagnostics(kalos::x13_si_seasonality(result.fit));
        Rcpp::List block;
        for(const auto& coefficients : result.coefficients)
            block.push_back(kalos_x13_rcpp::reals_out(coefficients));
        block.attr("names")=Rcpp::CharacterVector({"ar","ma","sar","sma"});
        std::size_t n_coefficients=0;
        for(const auto& coefficients : result.coefficients)
            n_coefficients+=coefficients.size();
        const R_xlen_t nvars=static_cast<R_xlen_t>(result.variables.size());
        Rcpp::CharacterVector var_names(nvars),var_types(nvars),var_effects(nvars);
        Rcpp::IntegerVector var_types_index(nvars),var_positions(nvars);
        Rcpp::NumericVector var_coefficients(nvars);
        Rcpp::List columns(nvars);
        for(std::size_t i=0;i<result.variables.size();++i) {
            const auto& variable=result.variables[i];
            const R_xlen_t k=static_cast<R_xlen_t>(i);
            var_names[k]=variable.name;
            var_types_index[k]=static_cast<int>(variable.type);
            var_types[k]=outlier_type_name(variable.type);
            var_effects[k]=variable.regeffect;
            var_positions[k]=static_cast<int>(variable.position)+1; // 1-based
            var_coefficients[k]=variable.coefficient;
            columns[k]=kalos_x13_rcpp::reals_out(variable.column);
        }
        out["model"]=Rcpp::List::create(Rcpp::_["log"]=result.log,Rcpp::_["mean"]=result.mean,
            Rcpp::_["p"]=result.p,Rcpp::_["d"]=result.d,Rcpp::_["q"]=result.q,
            Rcpp::_["bp"]=result.bp,Rcpp::_["bd"]=result.bd,Rcpp::_["bq"]=result.bq,
            Rcpp::_["coefficients"]=block,
            Rcpp::_["n_coefficients"]=static_cast<int>(n_coefficients),
            Rcpp::_["variables"]=Rcpp::List::create(Rcpp::_["name"]=var_names,
                Rcpp::_["type_index"]=var_types_index,Rcpp::_["type"]=var_types,
                Rcpp::_["regeffect"]=var_effects,Rcpp::_["position"]=var_positions,
                Rcpp::_["coefficient"]=var_coefficients,Rcpp::_["columns"]=columns),
            Rcpp::_["estimation_status"]=result.estimation_status,
            Rcpp::_["log_likelihood"]=result.log_likelihood,Rcpp::_["aicc"]=result.aicc,
            Rcpp::_["effective_observations"]=static_cast<double>(result.effective_observations));
        out["results"]=Rcpp::List::create(
            Rcpp::_["d10"]=kalos_x13_rcpp::reals_out(result.fit.d10),
            Rcpp::_["d11"]=kalos_x13_rcpp::reals_out(result.fit.d11),
            Rcpp::_["d12"]=kalos_x13_rcpp::reals_out(result.fit.d12),
            Rcpp::_["d13"]=kalos_x13_rcpp::reals_out(result.fit.d13),
            Rcpp::_["c17"]=kalos_x13_rcpp::reals_out(result.fit.c17),
            Rcpp::_["note"]=result.fit.mode_effective==kalos::X11DecompositionMode::log_additive
                ? "the X11Results tables as finals consumed them (log-additive "
                    "results-builder stage, Ratio bias default)"
                : "the X11Results tables as finals consumed them (raw step tables)");
        out["finals"]=Rcpp::List::create(Rcpp::_["schema_version"]=1,
            Rcpp::_["type"]="x13_finals",Rcpp::_["scope"]=record_scope + "_finals",
            Rcpp::_["tables"]=kalos_x13_rcpp::finals_list(result.fit.finals));
        out["finalSA"]=kalos_x13_rcpp::reals_out(result.fit.finalSA);
        out["finalTrend"]=kalos_x13_rcpp::reals_out(result.fit.finalTrend);
        out["finalIrregular"]=kalos_x13_rcpp::reals_out(result.fit.finalIrregular);
    }
    return out;
}
}
// [[Rcpp::export(name = ".native_x13_fit_automatic")]]
Rcpp::List native_x13_fit_automatic(SEXP values,SEXP frequency,SEXP start_period,
    SEXP first_year,SEXP aicc_log_correction,SEXP max_d,SEXP max_bd,SEXP ub1,SEXP ub2,SEXP cancel,
    SEXP differencing_precision,SEXP models_count,SEXP max_p,SEXP max_q,SEXP max_bp,
    SEXP max_bq,SEXP balanced,SEXP mixed,SEXP accept_white_noise,SEXP arma_eps,
    SEXP outliers_used,SEXP otypes,SEXP span,SEXP outlier_max_round,
    SEXP outlier_max_outliers,SEXP tcrate,SEXP precision,SEXP mean_precision,SEXP va,
    SEXP reduce_va,SEXP ljung_box_limit,SEXP check_mu,SEXP accept_airline,SEXP ur_limit,
    SEXP mean_test_initial,SEXP mean_test_final,SEXP tsig,SEXP calendar,SEXP easter,
    SEXP trading_days,SEXP user_variables,SEXP seasonal,SEXP sigma,SEXP lower_sigma,
    SEXP upper_sigma,SEXP msr_policy,SEXP mode,SEXP backcast,SEXP forecast) {
    kalos::X13AutomaticRequest request;
    request.input=kalos_x13_rcpp::reals_in(values,"x13 automatic fit input");
    request.frequency=integer_scalar(frequency,"frequency");
    request.first_period=integer_scalar(start_period,"start_period");
    request.first_year=integer_scalar(first_year,"first_year");
    request.aicc_log_correction=Rcpp::as<double>(aicc_log_correction);
    request.max_d=integer_scalar(max_d,"max_d");
    request.max_bd=integer_scalar(max_bd,"max_bd");
    request.ub1=Rcpp::as<double>(ub1);
    request.ub2=Rcpp::as<double>(ub2);
    request.cancel=Rcpp::as<double>(cancel);
    request.differencing_precision=Rcpp::as<double>(differencing_precision);
    request.models_count=integer_scalar(models_count,"models_count");
    request.max_p=integer_scalar(max_p,"max_p");
    request.max_q=integer_scalar(max_q,"max_q");
    request.max_bp=integer_scalar(max_bp,"max_bp");
    request.max_bq=integer_scalar(max_bq,"max_bq");
    request.balanced=flag(balanced,"x13 automatic fit balanced must be TRUE or FALSE");
    request.mixed=flag(mixed,"x13 automatic fit mixed must be TRUE or FALSE");
    request.arma_accept_white_noise=flag(accept_white_noise,
        "x13 automatic fit accept_white_noise must be TRUE or FALSE");
    request.arma_eps=Rcpp::as<double>(arma_eps);
    request.outliers_used=flag(outliers_used,"x13 automatic fit outliers used must be TRUE or FALSE");
    if(TYPEOF(otypes)!=STRSXP || !Rf_isNull(Rf_getAttrib(otypes,R_DimSymbol)))
        throw std::invalid_argument("x13 automatic fit outlier types must be a character vector of AO, LS, TC, SO");
    for(R_xlen_t i=0;i<Rf_xlength(otypes);++i)
        request.outlier_types.push_back(outlier_type(STRING_ELT(otypes,i)));
    if(Rf_isNull(span))request.outlier_span=std::nullopt;
    else {
        if (Rf_xlength(span) != 2)
            throw std::invalid_argument("x13 automatic fit outlier span must be NULL or a length-2 integer vector");
        const int first = integer_at(span, 0, "span", 1);
        const int last = integer_at(span, 1, "span", first);
        if (static_cast<std::size_t>(last) > request.input.size())
            throw std::invalid_argument("outlier span exceeds the observed series");
        request.outlier_span = std::make_pair(first - 1, last - 1);
    }
    request.outlier_max_round = integer_scalar(outlier_max_round, "outlier_max_round", 0);
    request.outlier_max_outliers = integer_scalar(outlier_max_outliers, "outlier_max_outliers", 0);
    request.tcrate=Rcpp::as<double>(tcrate);
    request.precision=Rcpp::as<double>(precision);
    request.mean_precision=Rcpp::as<double>(mean_precision);
    request.va=Rcpp::as<double>(va);
    request.reduce_va=Rcpp::as<double>(reduce_va);
    request.ljung_box_limit=Rcpp::as<double>(ljung_box_limit);
    request.check_mu=flag(check_mu,"x13 automatic fit check_mu must be TRUE or FALSE");
    request.accept_airline=flag(accept_airline,"x13 automatic fit accept_airline must be TRUE or FALSE");
    request.ur_limit=Rcpp::as<double>(ur_limit);
    request.mean_test_initial=Rcpp::as<double>(mean_test_initial);
    request.mean_test_final=Rcpp::as<double>(mean_test_final);
    request.tsig=Rcpp::as<double>(tsig);
    request.use_calendar_detection=flag(calendar,"x13 automatic fit calendar flag must be TRUE or FALSE");
    request.use_easter_detection=flag(easter,"x13 automatic fit easter flag must be TRUE or FALSE");
    request.automatic_trading_days=flag(trading_days,"x13 automatic fit trading_days flag must be TRUE or FALSE");
    request.use_user_variables=flag(user_variables,"x13 automatic fit user_variables flag must be TRUE or FALSE");
    request.seasonal=kalos_x11_rcpp::requests(seasonal);
    request.trend_length=0; // the decoder automatic surface: automatic Henderson
    request.sigma=sigma_policy(sigma);
    request.lower_sigma=Rcpp::as<double>(lower_sigma);
    request.upper_sigma=Rcpp::as<double>(upper_sigma);
    if(TYPEOF(msr_policy)!=STRSXP || Rf_xlength(msr_policy)!=1 ||
        STRING_ELT(msr_policy,0)==NA_STRING || !Rf_isNull(Rf_getAttrib(msr_policy,R_DimSymbol)))
        throw std::invalid_argument("MSR iteration policy must be jdplus or x11_five_year_cap");
    const auto msr_name=Rcpp::as<std::string>(msr_policy);
    request.msr_policy=kalos_x11_rcpp::parse_msr_policy(msr_name);
    request.mode=kalos_x11_rcpp::parse_mode(mode,"x13 automatic fit decomposition mode");
    request.backcast=integer_scalar(backcast,"backcast");
    request.forecast=integer_scalar(forecast,"forecast");
    // Also guard the installed snapshot, which predates checked core normalization.
    if (request.frequency == 12) {
        if (!request.input.empty() && request.first_period >= 1 && request.first_period <= 12) {
            const auto last = Rf_xlength(values) - 1 + static_cast<long long>(request.first_period - 1);
            if (last > std::numeric_limits<int>::max() ||
                request.first_year + last / 12 > std::numeric_limits<int>::max())
                throw std::invalid_argument("observed calendar exceeds the supported year range");
        }
        const auto months = [](int n) { return n < 0 ? -12LL * n : static_cast<long long>(n); };
        if (months(request.backcast) + months(request.forecast) + Rf_xlength(values) >
            std::numeric_limits<int>::max())
            throw std::invalid_argument("x13 automatic fit extended series exceeds the supported index range");
    }
    // The provenance echo is complete whatever the outcome: the module state
    // exactly as the boundary received it (1-based span kept 1-based).
    Rcpp::List out=Rcpp::List::create(Rcpp::_["schema_version"]=2,
        Rcpp::_["type"]="x13_automatic_fit",Rcpp::_["status"]=std::string(),
        Rcpp::_["request"]=Rcpp::List::create(
            Rcpp::_["frequency"]=request.frequency,Rcpp::_["start_period"]=request.first_period,
            Rcpp::_["first_year"]=request.first_year,
            Rcpp::_["transform"]=Rcpp::List::create(Rcpp::_["aicc"]=request.aicc_log_correction),
            Rcpp::_["automodel"]=Rcpp::List::create(Rcpp::_["max_d"]=request.max_d,
                Rcpp::_["max_bd"]=request.max_bd,Rcpp::_["ub1"]=request.ub1,
                Rcpp::_["ub2"]=request.ub2,Rcpp::_["cancel"]=request.cancel,
                Rcpp::_["precision"]=request.differencing_precision,
                Rcpp::_["models_count"]=request.models_count,Rcpp::_["max_p"]=request.max_p,
                Rcpp::_["max_q"]=request.max_q,Rcpp::_["max_bp"]=request.max_bp,
                Rcpp::_["max_bq"]=request.max_bq,Rcpp::_["balanced"]=request.balanced,
                Rcpp::_["mixed"]=request.mixed,
                Rcpp::_["accept_white_noise"]=request.arma_accept_white_noise,
                Rcpp::_["eps"]=request.arma_eps),
            Rcpp::_["outliers"]=Rcpp::List::create(Rcpp::_["used"]=request.outliers_used,
                Rcpp::_["types"]=Rcpp::CharacterVector(otypes),
                Rcpp::_["span"]=span,Rcpp::_["max_round"]=outlier_max_round,
                Rcpp::_["max_outliers"]=outlier_max_outliers,Rcpp::_["tcrate"]=request.tcrate),
            Rcpp::_["ami"]=Rcpp::List::create(Rcpp::_["precision"]=request.precision,
                Rcpp::_["mean_precision"]=request.mean_precision,Rcpp::_["va"]=request.va,
                Rcpp::_["reduce_va"]=request.reduce_va,
                Rcpp::_["ljung_box_limit"]=request.ljung_box_limit,
                Rcpp::_["check_mu"]=request.check_mu,Rcpp::_["accept_airline"]=request.accept_airline,
                Rcpp::_["ur_limit"]=request.ur_limit),
            Rcpp::_["regression"]=Rcpp::List::create(
                Rcpp::_["mean_test_initial"]=request.mean_test_initial,
                Rcpp::_["mean_test_final"]=request.mean_test_final,
                Rcpp::_["tsig"]=request.tsig),
            Rcpp::_["deferred"]=Rcpp::List::create(Rcpp::_["calendar"]=request.use_calendar_detection,
                Rcpp::_["easter"]=request.use_easter_detection,
                Rcpp::_["trading_days"]=request.automatic_trading_days,
                Rcpp::_["user_variables"]=request.use_user_variables),
            Rcpp::_["x11"]=Rcpp::List::create(
                Rcpp::_["seasonal"]=Rcpp::CharacterVector(seasonal),
                Rcpp::_["trend_length"]=0,Rcpp::_["sigma"]=Rf_isNull(sigma) ? R_NilValue : sigma,
                Rcpp::_["lower_sigma"]=request.lower_sigma,Rcpp::_["upper_sigma"]=request.upper_sigma,
                Rcpp::_["msr_policy"]=msr_name,Rcpp::_["mode"]=mode,
                Rcpp::_["backcast"]=request.backcast,Rcpp::_["forecast"]=request.forecast)),
        Rcpp::_["scope"]="automatic_spec_x13_fit_only",Rcpp::_["backend"]="native");
    const ScopedSingleBlas worker;
    const auto result=kalos::x13_automatic_fit(request);
    auto product=wrap_estimated_x13(out, result, Rcpp::CharacterVector(seasonal), msr_name, "automatic_fit");
    if(result.stock_likelihood) {
        // The public automatic-route residual diagnostics (M7 item 27 cutover): the
        // residuals.ser / residuals.lb / residuals.doornikhansen surface over the
        // final-model likelihood this route populates on its success path, in the
        // same wrap_monthly_residuals shape the estimated entry delivers (the
        // Doornik-Hansen arm is the bridge capture domain pinned by the core's
        // x13_automatic_residual_diagnostics). Refusals never carry the surface;
        // the accessor fails closed on a partial one, never a fallback.
        const auto diagnostic=kalos::x13_automatic_residual_diagnostics(result);
        product["residuals"]=wrap_monthly_residuals(result.stock_likelihood->full_residuals,
            diagnostic.standard_error,diagnostic.ljung_box,diagnostic.normality);
    }
    if(result.status=="available" && result.stock_likelihood) {
        // Public original-scale preprocessing forecast uncertainty (y_ef) on the
        // automatic route (M13-P5): the estimated entry's series_forecast_stdev
        // pair mirrored name-for-name, fed by kalos::x13_automatic_forecast_errors
        // over this fit's own landing horizon. The core entry's domain is an
        // available fit with the complete final-model surface; outside it the
        // fields stay ABSENT (never an estimated-on-the-side value), and inside it
        // the entry's own statuses ("available", "no_requested_predictions",
        // "covariance_range") pass through verbatim, its shape guards thrown.
        const auto errors=kalos::x13_automatic_forecast_errors(result, result.fit.forecast);
        product["series_forecast_stdev_status"]=errors.status;
        product["series_forecast_stdev"]=errors.standard_errors;
    }
    return product;
}

// [[Rcpp::export(name = ".native_x13_estimated_fit", rng = false)]]
Rcpp::List native_x13_estimated_fit(const std::vector<double>& observations, int year, int month,
                                    const std::string& preset) {
    const ScopedSingleBlas worker;
    kalos::X13AutomaticFitResult result;
    if (preset == "RSA0") result = kalos::x13_estimated_rsa0(observations, year, month);
    else if (preset == "RSA1") result = kalos::x13_estimated_rsa1(observations, year, month);
    else if (preset == "RSA2c") result = kalos::x13_estimated_rsa2c(observations, year, month);
    else Rcpp::stop("Unsupported stock X-13 template.");
    Rcpp::List out = Rcpp::List::create(Rcpp::_["schema_version"] = 1,
        Rcpp::_["type"] = "x13_estimated_fit", Rcpp::_["scope"] = "stock_estimated_x13_fit",
        Rcpp::_["template"] = preset, Rcpp::_["backend"] = "native");
    if (result.stock_likelihood) {
        const auto& likelihood = *result.stock_likelihood;
        out["likelihood"] = Rcpp::List::create(Rcpp::_["effective"] = static_cast<double>(result.effective_observations),
            Rcpp::_["estimated"] = static_cast<double>(likelihood.estimated_parameters),
            Rcpp::_["log_likelihood"] = result.log_likelihood, Rcpp::_["ssq"] = likelihood.residual_ssq);
        const auto diagnostic = kalos::x13_stock_residual_diagnostics(result);
        out["residuals"] = wrap_monthly_residuals(likelihood.full_residuals, diagnostic.standard_error,
                                                  diagnostic.ljung_box, diagnostic.normality);
        const auto errors = kalos::x13_stock_forecast_errors(result, result.fit.forecast);
        out["series_forecast_stdev_status"] = errors.status;
        out["series_forecast_stdev"] = errors.standard_errors;
    }
    auto product = wrap_estimated_x13(out, result, Rcpp::CharacterVector(12, "msr"), "jdplus", "estimated_fit");
    if (result.status == "available" && result.stock_calendar) {
        const auto& calendar = *result.stock_calendar;
        Rcpp::NumericVector coefficients = Rcpp::wrap(calendar.coefficients);
        coefficients.attr("names") = calendar.names;
        Rcpp::List model = product["model"];
        model["calendar"] = Rcpp::List::create(Rcpp::_["working_days"] = calendar.selection.working_days,
            Rcpp::_["easter_duration"] = calendar.selection.easter_duration, Rcpp::_["coefficients"] = coefficients);
        product["model"] = model;
    }
    return product;
}
