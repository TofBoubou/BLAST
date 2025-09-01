#pragma once
#include "gasp2/catalysis/finite_rate/types.hpp"

namespace gasp2::catalysis::finite_rate {
namespace details {

/**
 * @brief Build stoichiometric matrices for all finite-rate reactions.
 *
 * The routine assembles reactant and product stoichiometric matrices ordered by
 * gas species, surface species, and the empty site. It is used by the surface
 * coverage solver and the flux computation to access stoichiometric
 * coefficients quickly.
 *
 * @param props Precomputed reaction properties containing the processed
 *              reaction list.
 * @return Populated StoichiometryData structure.
 */
[[nodiscard]] StoichiometryData
build_stoichiometry(const ReactionProperties &props);

} // namespace details
} // namespace gasp2::catalysis::finite_rate
