#pragma once

// kalos M5 slice 24-B: the SEATS decomposition seam,
// UcarimaModel.setVarianceMax(-1, false) + compactUcarimaModel(false, false).
//
// Composes on the public 24-A surface (tramo_seats_model_validator.hpp): the
// caller runs validate_tramo_seats_model() on a TramoSeatsRequest, hands the
// request to decompose_ucarima() and receives the canonical UCARIMA split
// exactly as jdplus 3.9.0 prepares it before the SEATS filter construction
// starts - the un-reduced split ("d"), the reduced ucarima ("cmp") and the
// raw selector-chain components ("r") - as JVM-faithful fact strings.
//
// Envelope (pinned D6 policy, docs/contracts/decompose-ucarima-contract.md):
// only admissible requests with decomposition mode != Noisy and an
// approximation_required-free decision are decomposed. Everything else is a
// fail-closed std::invalid_argument with a pinned reason string (the JVM
// DefaultModelApproximator crash path is NOT ported: the validator already
// classifies it, the unit refuses rather than reproducing a crash). A
// request outside the envelope is a programming error at the call site, not
// a data-dependent failure: the validator did the data-dependent work.
//
// Downstream caveat (D1): this lands the UCARIMA seam only. The SEATS
// WK-period Wiener filters still consume the WK-only ucarima, so the cmp/k
// split here is not yet what the seasonal component is built from; the
// 24-C slice wires this seam into the pipeline and re-pins WK-vs-sfd.
//
// Bit-exactness: every emitted numeric is the JVM's double, produced by a
// line-for-line mirror of the pinned arithmetic (EigenValuesDecomposer2 +
// GeneralEigenSystem/hessenbergQR + UnitRootsSolver + LeastSquaresDivision
// + the ArimaModel/UcarimaModel algebra), formatting via the shared
// %.17g/JVM-toString helpers. Fact-line keys and order match
// validation/reference/DecomposeUcarimaProbe.java verbatim.
//
// Error contract: std::invalid_argument (envelope refusal, pinned reason);
// std::overflow_error (numeric abort). No other exception escapes.

#include <string>
#include <vector>

#include "kalos/tramo_seats_model_validator.hpp"

namespace kalos {

// Decomposes `request' through the pinned 24-A validation
// (validate_tramo_seats_model). Throws std::invalid_argument with the pinned
// reason for a Noisy request mode or an approximation_required decision
// (D6 envelope); validator refusals propagate unchanged. Inside the
// envelope the seam itself never throws: JVM abort classes (MatrixException,
// NPE, AIOOBE) are captured as data (null d/r blocks, the per-getter
// maerr/varr/sfdr lines, the SeatsException status row).
struct decompose_ucarima_result {
  // The 24-A decision this decomposition was made on (pass-through).
  const tramo_seats_model_result& decision() const { return decision_; }

  // JVM-faithful fact lines, probe order (see
  // validation/reference/DecomposeUcarimaProbe.java): the decomposition-
  // seam grammar - "d.decompose" (ok|null) and the d block, the reduced
  // "cmp" block (SeatsModel.compactUcarimaModel(false, false)), the raw
  // "r.decompose"/"r.count" and the r cN blocks - ending with the "status"
  // row (the SeatsException line when d came back null under mode None; a
  // null d is a legal, data-dependent outcome, e.g. quarterly D=1/B=1
  // splits). Per-getter failures inside a block appear as "cN maerr
  // <class>" / "cN varr <class>" lines carrying the probe's dangling-key
  // prefix. Record-only rows (toolkit, v.*, model, k.*, the estimation
  // statuses) stay in the csv corpus and are NOT emitted. Each line
  // starts with its key.
  const std::vector<std::string>& facts() const { return facts_; }

  // Convenience: the fact value for `key' (rest of the line after the
  // first whitespace run) or "" when the key is absent.
  std::string fact(const std::string& key) const;

  decompose_ucarima_result() = default;

  decompose_ucarima_result(tramo_seats_model_result decision,
                           std::vector<std::string> facts)
      : decision_(std::move(decision)), facts_(std::move(facts)) {}

 private:
  tramo_seats_model_result decision_;
  std::vector<std::string> facts_;
};

decompose_ucarima_result decompose_ucarima(const tramo_seats_model_request& request);

}  // namespace kalos
