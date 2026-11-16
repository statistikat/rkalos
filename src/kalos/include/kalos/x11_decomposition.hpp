#ifndef KALOS_X11_DECOMPOSITION_HPP
#define KALOS_X11_DECOMPOSITION_HPP
// Item 22b.5 M3: the decomposition modes of the pinned X-11 pipeline,
// mirroring jdplus.sa.base.api.DecompositionMode in pinned declaration
// order (Undefined excluded: the pinned steps treat it as additive in every
// decoded branch and no pinned spec default uses it). The fixed passes take
// the mode as an explicit argument; the pinned mode-driven behaviour (the
// ratio remove, the positivity repair of the trend tables, the mode-aware
// corrector, the relative-form MSRE table) is documented with the evidence
// in docs/contracts/x11-multiplicative-contract.md.
//
// Delivered and qualified: additive (items 22b.4 S1-S4), multiplicative
// (item 22b.5 M3a), log_additive (item 22b.5 M3b-i: the pinned kernel
// transforms the input with Math.log before the plain additive steps; the
// delivered tables are those step-level log-domain tables, the results-
// builder exp of D10..D13 is not) and pseudo_additive (item 22b.5 M3b-ii:
// the pinned X11{B,C,D}StepPseudoAdd pipeline — the multiplicative ratio
// context with the PseudoAddUtility forms at the twelve overridden
// composition sites and the removeIrregular global-MSR walk, documented in
// docs/contracts/x11-pseudo-additive-contract.md).
#include <cstddef>
namespace kalos {
enum class X11DecompositionMode {
    additive,
    multiplicative,
    log_additive,
    pseudo_additive
};
// The pinned X11Context.isMultiplicative(): true for BOTH Multiplicative
// and PseudoAdditive — the ratio arithmetic alone does not identify the
// pinned pseudo-additive pipeline, so passes must gate on the full mode.
inline bool x11_mode_is_multiplicative(X11DecompositionMode mode) {
    return mode==X11DecompositionMode::multiplicative
        ||mode==X11DecompositionMode::pseudo_additive;
}
} // namespace kalos
#endif
