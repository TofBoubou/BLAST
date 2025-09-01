#pragma once
#include "gasp2/catalysis/finite_rate/types.hpp"
#include <vector>

namespace gasp2::catalysis::finite_rate {
namespace details {

/**
 * @brief Solve for steady-state surface coverages.
 *
 * The routine enforces zero net production for every adsorbed species and
 * applies a global site-balance constraint to determine the surface coverage
 * of all adsorbed species and the density of empty surface sites.
 *
 * @param props Precomputed reaction properties including rate constants,
 *              gas-phase concentrations and site density.
 * @return Vector of surface coverages (mol/m^2) in the order of
 *         ctx.surface_species followed by the empty site density.
 */
[[nodiscard]] std::vector<double>
solve_surface_coverage(const ReactionProperties &props);

} // namespace details
} // namespace gasp2::catalysis::finite_rate
