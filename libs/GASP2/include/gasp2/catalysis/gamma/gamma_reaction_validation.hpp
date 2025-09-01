//===----------------------------------------------------------------------===//
// Gamma reaction validation utilities
//
// Two levels of validation are provided:
//  - static_validate_gamma_reactions: performed once during initialization and
//    checks conditions independent of the wall state (e.g. mass/charge
//    conservation).
//  - verify_gamma_reactions: executed for each call to the flux computation and
//    validates quantities that depend on \f$T_w\f$ or wall densities (gamma
//    values). The latter also populates ReactionInput::computed_gamma so that
//    downstream routines can reuse the values without recomputation.
//===----------------------------------------------------------------------===//

#pragma once
#include "gasp2/types.hpp"
#include <gasp2/catalysis/gamma/types.hpp>
#include <unordered_map>
#include <vector>

namespace gasp2::catalysis::gamma {
namespace details {

using gasp2::SpeciesData;

// Perform initialization-time checks that do not depend on wall state.
void static_validate_gamma_reactions(std::vector<ReactionInput> &reactions,
                                     const std::vector<std::string> &species);

// Compute gamma values using the current wall state and ensure they remain
// physically consistent (0 <= gamma <= 1 with cumulative sums <= 1).
void verify_gamma_reactions(
    std::vector<ReactionInput> &reactions, const SpeciesData &species,
    const std::unordered_map<std::string, std::size_t> &index);

} // namespace details
} // namespace gasp2::catalysis::gamma
