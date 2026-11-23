// Item 24-B: the SEATS decomposition seam as a bit-exact native mirror.
//
// kalos::decompose_ucarima() re-runs the 24-A validation, applies the D6
// envelope (decomposable requests only; approximation_required / Noisy are
// fail-closed refusals with pinned reasons), then mirrors the LIVE jdplus
// 3.9.0 path the probe observed:
//
//   DefaultModelDecomposer.decompose -> selector chain (24-A engine) and
//                                      UcarimaModel.setVarianceMax(-1, false)
//   SeatsModel.compactUcarimaModel(false, false)
//   the probe's dumpModel / dumpUcarima fact grammar
//
// Fact scope: the decomposition seam (d.decompose, the d/cmp/r blocks, and
// the deterministic "status error:...SeatsException:Invalid decomposition"
// line the kernel throws before estimation whenever the decomposition came
// back null). The kernel's k.* lines and the "status ok" / estimation-phase
// statuses (eco.EcoException from DeterminantalTerm, MatrixException during
// variance computation) depend on the estimation machinery beyond this seam
// and are NOT emitted; see docs/contracts/decompose-ucarima-contract.md.
//
// Every routine transcribes the pinned bytecode recorded in
// tmp_opencode/jdprobe (m5/, m524/); the mirrored class is named in the
// comments. Operation order, comparison flavor (dcmpg vs dcmpl -> *_g/*_l)
// and rounding sequence are preserved; the build uses -ffp-contract=off so
// no FMA is introduced.
//
// The numerically interesting part is ArimaModel.getMa()/getIV(), i.e.
// SymmetricFilter.factorize(): the constant shortcut, then the factorizer
// chain EV2 -> EV1 (EigenValuesDecomposer2/EigenValuesDecomposer, installed
// by SymmetricFilter's static{}): UnitRootsSolver -> IVector.sqrt -> exact
// LS divide, the companion-matrix GeneralEigenSystem (balance +
// hessenbergQR), root selection (EV2: real-pair / unit-circle
// classification; EV1: keep |z| > 1), lejaOrder + fromComplexRoots and the
// scaling. A JVM chain probe over every distinct corpus symmetric filter
// (150) settles the mirror's scope: 124 resolve in EV2, 5 in EV1, 1 in the
// fast decomposer, 20 in none (those components carry the corpus maerr
// lines); the single fast-decomposer filter belongs to a validator-refused
// case (v_airline r c1), so EV2 -> EV1 -> null reproduces all 49 cases.

#include "kalos/decompose_ucarima.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "kalos/detail/ucarima_engine.hpp"
#include "kalos/detail/seats_ucarima_model.hpp"

namespace kalos {
namespace {
using namespace detail::seats_model;

// ---------------------------------------------------------------------------
// Fact emission (the probe's dumpModel/dumpUcarima grammar, verbatim).
// ---------------------------------------------------------------------------
class Emitter {
 public:
  void add(const std::string& line) { lines_.push_back(line); }
  const std::vector<std::string>& lines() const { return lines_; }

  void dump_model(const std::string& key, const AModel& m) {
    // dumpModel's grammar (probe 187-222): sar+ur in ONE guarded block, then
    // ma/var/sfd per-getter (failures captured, never silent). sfd only for
    // the ArimaModel variant - every component on the corpus path is one.
    // StringBuilder semantics faithfully: sb.append(key).append(" ma ")
    // already ran when seq(...) throws, so the failure line carries the
    // dangling prefix - corpus lines literally read "r c0 ma r c0 maerr ...".
    std::string pending = key + " sar ";
    try {
      pending += coeffs_row(m.stationary_ar().c);  // throw BEFORE add: the
      add(pending);                                // sar line is a completed
      pending = key + " ur ";                      // JVM line (newline flush)
      pending += coeffs_row(m.non_stationary_ar().c);
      add(pending);
    } catch (const std::exception&) {
      add(pending + key + " arrexcep " + exception_class(std::current_exception()));
      return;
    }
    pending = key + " ma ";
    try {
      add(pending + coeffs_row(m.get_ma().p.c));
    } catch (const std::exception&) {
      add(pending + key + " maerr " + exception_class(std::current_exception()));
    }
    pending = key + " var ";
    try {
      add(pending + fmt_g(m.get_iv()));
    } catch (const std::exception&) {
      add(pending + key + " varr " + exception_class(std::current_exception()));
    }
    pending = key + " sfd ";
    try {
      add(pending + coeffs_row(m.symmetric_ma().a));
    } catch (const std::exception&) {
      add(pending + key + " sfdr " + exception_class(std::current_exception()));
    }
  }

  // The SarimaModel (outer model) variant: sar/ur/ma/var, no sfd line.
  void dump_model(const std::string& key, const UcarimaModel& um) {
    add(key + " sar " + coeffs_row(um.m_sar));
    add(key + " ur " + coeffs_row(um.m_ur));
    add(key + " ma " + coeffs_row(um.m_ma));
    add(key + " var " + fmt_g(um.m_var));
  }

  void dump_ucarima(const std::string& key, const UcarimaModel& um) {
    std::string valid;
    try {
      valid = um_is_valid(um) ? "1" : "0";
    } catch (const std::exception&) {
      add(key + " excep " + exception_class(std::current_exception()));
      return;
    }
    add(key + " count " + std::to_string(um.comps.size()) + " valid " + valid);
    dump_model(key + " model", um);
    for (std::size_t i = 0; i < um.comps.size(); ++i)
      dump_model(key + " c" + std::to_string(i), um.comps[i]);
  }

 private:
  std::vector<std::string> lines_;
};

AModel convert(const detail::Arima& c) {
  AModel a;
  a.ar = c.ar;
  a.delta = c.delta;
  a.symmetric = c.symmetric;
  if (c.symmetric) {
    a.sma = c.sma;
  } else {
    a.has_ma = true;
    a.ma1 = c.ma1;
    a.var = c.var;
  }
  return a;
}

// The 24-B mirror driver: decision envelope -> raw chain -> setVarianceMax
// -> compactUcarimaModel, facts in probe order.
std::vector<std::string> run_seam(const tramo_seats_model_request& request,
                                  const tramo_seats_model_result& dec) {
  detail::Sarima m1;
  m1.period = dec.period;
  m1.d = dec.orders_d;
  m1.bd = dec.orders_bd;
  m1.phi = dec.phi;
  m1.bphi = dec.bphi;
  m1.theta = dec.theta;
  m1.btheta = dec.btheta;
  detail::EngineFacts ef;
  try {
    ef = detail::decompose(m1, request.rmod, request.epsphi, request.smod,
                           request.smod_pi,
                           false);  // Noisy requests were refused above
  } catch (const std::invalid_argument&) {
    throw;  // envelope violation, never swallowed
  } catch (const detail::jvm_decompose_abort&) {
    ef.chain_threw = true;  // JVM: MatrixException -> decompose returns null
  } catch (const std::exception&) {
    // DefaultModelDecomposer.decompose catches every Exception the selector
    // chain throws and returns null (d.decompose null). The probe's raw
    // re-run would print r.excep <class>; the corpus shows r.decompose ok
    // for all 49 cases, so a throwing chain is itself a mirror divergence -
    // emit the r.decompose null line rather than escaping the unit.
    ef.chain_threw = true;
  }

  UcarimaModel raw;
  detail::Arima ma1 = m1.to_arima();
  raw.m_sar = ma1.stationary_ar().c;
  raw.m_ur = ma1.non_stationary_ar().c;
  raw.m_ma = ma1.ma1;
  raw.m_var = 1.0;  // SarimaModel.getInnovationVariance() on every corpus path
  for (const detail::Arima& c : ef.raw) raw.comps.push_back(convert(c));

  std::optional<UcarimaModel> d;
  if (!ef.chain_threw) {
    try {
      d = um_set_variance_max(raw, -1, false);
    } catch (const std::invalid_argument&) {
      throw;  // envelope violation, never swallowed
    } catch (const detail::jvm_decompose_abort&) {
      d = std::nullopt;  // JVM: MatrixException caught in decompose -> null
    } catch (const std::exception&) {
      d = std::nullopt;  // DefaultModelDecomposer: catch (Exception) -> null
    }
  }

  Emitter em;
  em.add(std::string("d.decompose ") + (d ? "ok" : "null"));
  if (d) {
    em.dump_ucarima("d", *d);
    try {
      em.dump_ucarima("cmp", seats_compact_ucarima_model(*d));
    } catch (const std::invalid_argument&) {
      throw;
    } catch (const detail::jvm_decompose_abort&) {
      em.add(std::string("cmp.excep ") + exception_class(std::current_exception()));
    } catch (const std::exception&) {
      // Probe line 470: the compact failure is its own DOTTED key line.
      em.add(std::string("cmp.excep ") + exception_class(std::current_exception()));
    }
  }
  em.add(std::string("r.decompose ") + (ef.chain_threw ? "null" : "ok"));
  if (!ef.chain_threw) {
    em.add("r.count " + std::to_string(raw.comps.size()));
    for (std::size_t i = 0; i < raw.comps.size(); ++i)
      em.dump_model("r c" + std::to_string(i), raw.comps[i]);
  }
  if (!d) {
    // mode None (approximation / Noisy were refused): the kernel throws
    // before estimation, deterministically.
    em.add("status error:jdplus.tramoseats.base.api.seats.SeatsException:Invalid "
           "decomposition");
  }
  return em.lines();
}

}  // namespace

std::string decompose_ucarima_result::fact(const std::string& key) const {
  for (const std::string& line : facts_) {
    if (line.size() >= key.size() && line.compare(0, key.size(), key) == 0 &&
        (line.size() == key.size() || line[key.size()] == ' ')) {
      return line.substr(line.find_first_not_of(' ', key.size()) == std::string::npos
                             ? line.size()
                             : line.find_first_not_of(' ', key.size()));
    }
  }
  return std::string();
}

decompose_ucarima_result decompose_ucarima(const tramo_seats_model_request& request) {
  tramo_seats_model_result dec = validate_tramo_seats_model(request);
  if (request.mode == seats_approx_mode::noisy)
    throw std::invalid_argument(
        "decompose_ucarima: Noisy approximation mode is outside the D6 envelope");
  if (dec.approximation_required)
    throw std::invalid_argument(
        "decompose_ucarima: approximation_required decision is outside the D6 envelope");
  return decompose_ucarima_result(dec, run_seam(request, dec));
}

}  // namespace kalos
