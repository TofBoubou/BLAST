#pragma once

#include "gasp2/types.hpp"
#include <gasp2/catalysis/gamma/types.hpp>
#include <vector>

namespace gasp2::catalysis::gamma {

using gasp2::CatalysisFluxes;
using gasp2::Result;
using gasp2::SpeciesData;

/// Compute the species mass flux at the wall from the recombination probability
/// and precomputed impinging mass flux m*imp_flux.
[[nodiscard]] Result<double> compute_flux(double gamma, double mass_imp_flux,
                                          bool first_order);

/// Combine precomputed recombination probabilities and the nu/mu matrices to
/// obtain the final mass fluxes for all species.
/// @param species        Precomputed species properties.
/// @param ns             Number of species (cached from the context).
/// @param nr             Number of reactions (cached from the context).
/// @param gammas         Per-reaction recombination probabilities.
/// @param nu             Reactant indicator matrix.
/// @param nu_p           Reactant stoichiometric coefficients.
/// @param mu             Production tensor \f$\mu_{l,i,j}\f$.
/// @param first_order    Use first-order kinetics when true.
/// @param limiting_fluxes Adjust impinging fluxes based on limiting species
///                        when true.
/// @param types          Reaction types for each reaction.
/// @param heterogeneous  Flags indicating heterogeneous reactions.
/// @param species_order  Ordering of species names for diagnostics.
/// @param debug          Enable verbose output when true.
[[nodiscard]] Result<CatalysisFluxes>
compute_all_fluxes(const SpeciesData &species, std::size_t ns, std::size_t nr,
                   const std::vector<std::vector<double>> &gammas,
                   const std::vector<std::vector<int>> &nu,
                   const std::vector<std::vector<int>> &nu_p,
                   // mu_l,i,j: for reaction l, production of species i from j
                   const std::vector<std::vector<std::vector<int>>> &mu,
                   bool first_order, bool limiting_fluxes,
                   const std::vector<CatalysisModel> &types,
                   const std::vector<bool> &heterogeneous,
                   const std::vector<std::string> &species_order, bool debug);

} // namespace gasp2::catalysis::gamma
