#pragma once
#include "gasp2/catalysis/finite_rate/finite_rate.hpp"
#include <string>
#include <vector>

namespace gasp2::catalysis::finite_rate {
namespace details {

/// Perform static validation checks on the parsed finite-rate reactions.
/// Ensures mass, charge, and site conservation along with basic consistency
/// checks such as species availability, adsorption/desorption pairing, and
/// presence of all required kinetic parameters for each reaction type.
void static_validate_finite_rate_reactions(
    const std::vector<ReactionInput> &reactions,
    const std::vector<std::string> &species_order);

} // namespace details
} // namespace gasp2::catalysis::finite_rate
