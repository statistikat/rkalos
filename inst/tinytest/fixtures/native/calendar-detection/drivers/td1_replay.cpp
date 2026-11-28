// m17-td1 Tier-1 replay probe (21-I stage 1): run the native ATD machine on
// the exact corpora the Td1Probe JVM battery consumed and print the replay
// surface in the JVM capture's line shapes: one "det ATD" line per firing
// (names/best/arms at %.17g) plus the final facts (log, ll/aicc, landed
// variable list). Templates: td2noe | td7noe | td2noadj.
#include "kalos/x13_automatic.hpp"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

const double kAir[144] = {
    112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118, 115, 126, 141,
    135, 125, 149, 170, 170, 158, 133, 114, 140, 145, 150, 178, 163, 172, 178,
    199, 199, 184, 162, 146, 166, 171, 180, 193, 181, 183, 218, 230, 242, 209,
    191, 172, 194, 196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
    204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229, 242, 233, 267,
    269, 270, 315, 364, 347, 312, 274, 237, 278, 284, 277, 317, 313, 318, 374,
    413, 405, 355, 306, 271, 306, 315, 301, 356, 348, 355, 422, 465, 467, 404,
    347, 305, 336, 340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
    360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405, 417, 391, 419,
    461, 472, 535, 622, 606, 508, 461, 390, 432};

std::vector<double> read_series(const std::string& path, const std::string& kind) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::string line;
  std::vector<double> out;
  bool active = false;
  while (std::getline(in, line)) {
    std::istringstream head(line);
    std::string tag;
    if (head >> tag && tag == "series") {
      std::string name;
      head >> name;
      active = name == kind;
      continue;
    }
    if (active) out.push_back(std::stod(line));
  }
  return out;
}

kalos::X13AutomaticRequest shape(std::vector<double> input, int year) {
  kalos::X13AutomaticRequest r{};
  r.input = std::move(input);
  r.frequency = 12;
  r.first_period = 1;
  r.first_year = year;
  r.aicc_log_correction = -2.;
  r.max_d = 2;
  r.max_bd = 1;
  r.ub1 = 1. / 1.0309278350515465;
  r.ub2 = 1. / 1.1363636363636365;
  r.cancel = 0.1;
  r.differencing_precision = 1e-7;
  r.models_count = 5;
  r.max_p = 2;
  r.max_q = 2;
  r.max_bp = 1;
  r.max_bq = 1;
  r.balanced = false;
  r.mixed = true;
  r.arma_accept_white_noise = false;
  r.arma_eps = 1e-7;
  r.outliers_used = true;
  r.outlier_types = {kalos::OutlierType::Additive, kalos::OutlierType::LevelShift,
                     kalos::OutlierType::TransitoryChange};
  r.outlier_max_round = 30;
  r.outlier_max_outliers = 30;
  r.tcrate = 0.7;
  r.precision = 1e-7;
  r.mean_precision = 1e-5;
  r.va = 0.;
  r.reduce_va = 0.14286;
  r.ljung_box_limit = 0.95;
  r.check_mu = true;
  r.accept_airline = false;
  r.ur_limit = 0.95;
  r.mean_test_initial = 1.96;
  r.mean_test_final = 1.;
  r.tsig = 1.;
  r.seasonal.assign(12, kalos::X11SeasonalRequest::msr);
  r.trend_length = 0;
  r.sigma = kalos::X11SigmaPolicy::legacy;
  r.lower_sigma = 1.5;
  r.upper_sigma = 2.5;
  r.msr_policy = kalos::X11MsrIterationPolicy::jdplus;
  r.mode = kalos::X11DecompositionMode::additive;
  r.backcast = 0;
  r.forecast = -1;
  r.automatic_trading_days = true;
  r.td_leap_year = true;
  r.td_aicc_diff = 0.;
  r.td_estimation_eps = 1e-5;
  return r;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: td1_replay <td2noe|td7noe|td2noadj> "
                         "<clean|tdeff|air> [corpora-file]\n");
    return 2;
  }
  const std::string tpl = argv[1], kind = argv[2];
  std::vector<double> input;
  int year = 1949;
  if (kind == "air") {
    input.assign(kAir, kAir + 144);
  } else {
    if (argc < 4) {
      std::fprintf(stderr, "corpora file required for %s\n", kind.c_str());
      return 2;
    }
    input = read_series(argv[3], kind);
    year = 1990;
  }
  auto req = shape(std::move(input), year);
  if (tpl == "rsa3") {
    req.automatic_trading_days = false;  // baseline: the pinned no-td shape
  } else if (tpl == "td2noe") {
    req.td_clustering = 0;  // td2
  } else if (tpl == "td7noe") {
    req.td_clustering = 8;  // td7
  } else if (tpl == "td2noadj") {
    req.td_clustering = 0;
    req.td_auto_adjust = false;
  } else {
    std::fprintf(stderr, "unknown template %s\n", tpl.c_str());
    return 2;
  }

  kalos::X13AutomaticFitResult out = kalos::x13_automatic_fit(req);
  std::printf("specinfo %s %s n=%zu\n", tpl.c_str(), kind.c_str(),
              req.input.size());
  std::printf("status %s\n", out.status.c_str());
  for (const auto& e : out.chain) {
    std::printf("chain %s r%d", e.stage.c_str(), e.round);
    for (const auto& [key, value] : e.scalars)
      std::printf(" %s=%.17g", key.c_str(), value);
    if (e.accept)
      std::printf(" acc[rtval=%.17g rvr=%.17g lb_p=%.17g accepted=%d]",
                  e.accept->rtval, e.accept->rvr, e.accept->lb_p,
                  e.accept->accepted ? 1 : 0);
    if (e.pass2)
      std::printf(" p2[filter=%d plbox=%.17g plbox0=%.17g pcr=%.17g cont=%d]",
                  e.pass2->filter, e.pass2->plbox, e.pass2->plbox0,
                  e.pass2->pcr_after, e.pass2->continue_round ? 1 : 0);
    if (e.regtest) {
      std::printf(" removed=[");
      for (std::size_t j = 0; j < e.regtest->removed.size(); ++j)
        std::printf("%s%d", j ? "," : "", e.regtest->removed[j] ? 1 : 0);
      std::printf("] tstats=[");
      for (std::size_t j = 0; j < e.regtest->tstats.size(); ++j)
        std::printf("%s%.6g", j ? "," : "", e.regtest->tstats[j]);
      std::printf("]");
    }
    std::printf("\n");
  }
  for (const auto& e : out.chain) {
    if (e.outliers) {
      for (const auto& sc : e.outliers->scans)
        std::printf("scan %s pos%zu type%zu t=%.17g beyond=%d\n",
                    e.stage.c_str(), sc.position, sc.type, sc.t_stat,
                    sc.beyond_cv ? 1 : 0);
      if (!e.outliers->stop_reason.empty())
        std::printf("ostop %s %s\n", e.stage.c_str(),
                    e.outliers->stop_reason.c_str());
    }
  }
  for (const auto& e : out.chain) {
    if (e.refine) {
      for (const auto& ev : e.refine->evaluations)
        std::printf("eval %s s%d/%d p%d q%d bp%d bq%d ll=%.17g bic=%.17g\n",
                    e.stage.c_str(), ev.stage, ev.slot, ev.p, ev.q, ev.bp,
                    ev.bq, ev.log_likelihood, ev.bic);
      for (std::size_t i = 0; i < e.refine->preferred_models.size(); ++i) {
        const auto& m = e.refine->preferred_models[i];
        std::printf("pref %s %zu p%d q%d bp%d bq%d bic=%.17g ll=%.17g\n",
                    e.stage.c_str(), i, m.p, m.q, m.bp, m.bq, m.bic,
                    m.log_likelihood);
      }
      for (const auto& st : e.refine->refinement)
        std::printf("step %s cand%d rule%d inc%d moves%d\n", e.stage.c_str(),
                    st.candidate, st.rule, st.incumbent, st.moves);
    }
  }
  for (const auto& e : out.chain) {
    if (!e.atd) continue;
    const auto& a = *e.atd;
    std::printf("evt %s best=%d", e.stage.c_str(), a.best);
    std::printf(" names=[");
    for (std::size_t i = 0; i < a.names.size(); ++i)
      std::printf("%s%s", i ? ", " : "", a.names[i].c_str());
    std::printf("]");
    for (int k = 0; k < 2; ++k)
      std::printf(" ll[%d] n=%.0f eff=%.0f np=%.0f ssq=%.17g aicc=%.17g", k,
                  a.arms[k].observations, a.arms[k].effective,
                  a.arms[k].parameters, a.arms[k].ssq, a.arms[k].aicc);
    std::printf("\n");
  }
  std::printf("fact log=%d mean=%d orders %d %d %d %d %d %d\n", out.log ? 1 : 0,
              out.mean ? 1 : 0, out.p, out.d, out.q, out.bp, out.bd, out.bq);
  std::printf("fact ll=%.17g aicc=%.17g eff=%zu\n", out.log_likelihood,
              out.aicc, out.effective_observations);
  std::printf("landed mean=%d", out.mean ? 1 : 0);
  if (out.auto_calendar)
    for (const auto& name : out.auto_calendar->names)
      std::printf(" %s", name.c_str());
  for (const auto& v : out.variables) std::printf(" %s", v.name.c_str());
  std::printf("\n");
  if (out.auto_calendar) {
    std::printf("auto_calendar preadjust=%d", out.auto_calendar->preadjusted ? 1 : 0);
    for (const double c : out.auto_calendar->coefficients)
      std::printf(" %.17g", c);
    std::printf("\n");
  }
  for (const auto& v : out.variables)
    std::printf("var %s %.17g\n", v.name.c_str(), v.coefficient);
  std::printf("fit_status %s\n", out.fit.status.c_str());
  std::fflush(stdout);
  return 0;
}
