#pragma once

#include <gasp2/types.hpp>

namespace gasp2::catalysis::gamma::bose {

using gasp2::CatalysisFluxes;
using gasp2::SpeciesData;

/// Compute mass fluxes for the Bose recombination model.
/// @param species Precomputed species properties including impinging fluxes.
/// @param ns      Number of species (cached from the context).
/// @param idx_O   Index of atomic oxygen.
/// @param idx_CO  Index of carbon monoxide.
/// @param idx_O2  Index of molecular oxygen.
/// @param idx_CO2 Index of carbon dioxide.
/// @param gamma   Base recombination probability.
/// @param p2      Probability of CO oxidation.
/// @return Mass fluxes for all species.
[[nodiscard]] CatalysisFluxes
compute_fluxes(const SpeciesData &species, std::size_t ns, std::size_t idx_O,
               std::size_t idx_CO, std::size_t idx_O2, std::size_t idx_CO2,
               double gamma, double p2);

} // namespace gasp2::catalysis::gamma::bose
