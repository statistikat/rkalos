#pragma once

#include "kalos/x11_decomposition.hpp"

namespace kalos::detail {

// X13Kernel.updateSpec preserves PseudoAdditive before inspecting the
// preprocessing transform. Other modes follow the fitted log/level choice.
inline X11DecompositionMode x13_effective_mode(bool log, X11DecompositionMode requested) {
    if (requested == X11DecompositionMode::pseudo_additive) return requested;
    if (!log) return X11DecompositionMode::additive;
    return requested == X11DecompositionMode::log_additive ? requested : X11DecompositionMode::multiplicative;
}

} // namespace kalos::detail
