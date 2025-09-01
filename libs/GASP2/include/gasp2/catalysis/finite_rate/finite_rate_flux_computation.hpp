#pragma once
#include "gasp2/catalysis/finite_rate/types.hpp"
#include "gasp2/gasp2.hpp"
#include <span>

namespace gasp2::catalysis::finite_rate {
namespace details {

/**
 * @brief Compute mass fluxes for finite-rate surface reactions.
 *
 * Evaluates reaction rates using gas-phase concentrations and solved surface
 * coverages, returning mass fluxes for all gas species.
 *
 * @param props    Precomputed reaction properties (rate constants, etc.).
 * @param coverage Surface coverages (mol/m^2) for surface species followed by
 *                 empty sites.
 * @return CatalysisFluxes containing the finite-rate reaction contributions.
 */
[[nodiscard]] CatalysisFluxes
compute_reaction_fluxes(const ReactionProperties &props,
                        std::span<const double> coverage);

} // namespace details
} // namespace gasp2::catalysis::finite_rate
