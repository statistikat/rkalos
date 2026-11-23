// Item 24-A: SEATS model validation + canonical-decomposition DECISION mirror.
//
// Every routine below transcribes the pinned jdplus 3.9.0 bytecode recorded in
// tmp_opencode/jdprobe (m5/ and m524/ dumps); the mirrored class is named in
// the comments. Operation order, comparison flavor (dcmpg vs dcmpl -> the *_g /
// *_l helpers) and rounding sequence are preserved; the build uses
// -ffp-contract=off so no FMA is introduced.
//
// Nothing here re-implements DefaultModelApproximator (register D6): where the
// live path would call it, this unit reports approximation_required instead.
// Spectral factorization (SymmetricFilter.factorize) and the post-allocation
// component values are item 24-B scope; this unit exposes only the
// pre-allocation (raw) structure and the decision layer.
//
// The engine itself (everything the driver calls) moved to
// include/kalos/detail/ucarima_engine.hpp at item 24-B so that
// kalos::decompose_ucarima composes on the identical code; this file is now
// the thin 24-A public entry point. Bytecode-order comments live in the
// header.

#include "kalos/tramo_seats_model_validator.hpp"

#include "kalos/detail/ucarima_engine.hpp"

namespace kalos {

tramo_seats_model_result validate_tramo_seats_model(
    const tramo_seats_model_request& request) {
  return detail::engine_run(request);
}

}  // namespace kalos
