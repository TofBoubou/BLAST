#pragma once

#include "gasp2/gasp2.hpp"
#include <gasp2/catalysis/gamma/types.hpp>
#include <unordered_map>
#include <vector>

namespace gasp2::catalysis::gamma {

using gasp2::CatalysisFluxes;
using gasp2::Result;
using gasp2::SpeciesData;

/// Dispatch the gamma-model catalysis computation based on the reactions read
/// from the input file. For each reaction the appropriate gamma model is
/// invoked using precomputed stoichiometric data.
///
/// @param species        Precomputed species properties.
/// @param n_species      Total number of species.
/// @param n_reactions    Total number of reactions.
/// @param reactions      Parsed reactions with recombination data.
/// @param first_order    True if the model uses first-order kinetics.
/// @param limiting_fluxes Adjust impinging fluxes when true.
/// @param index          Mapping from species name to its index in
///                       species_order.
/// @param nu            Stoichiometric matrices.
/// @param all_super      True if all reactions are super-catalytic.
/// @param all_rini       True if all reactions use the Rini model.
/// @param species_order  Ordering of gas-phase species.
/// @param debug          Enable verbose diagnostic output.
/// @returns Catalysis fluxes or an error on failure.
[[nodiscard]] Result<CatalysisFluxes>
gamma_dispatcher(const SpeciesData &species, std::size_t n_species,
                 std::size_t n_reactions,
                 const std::vector<ReactionInput> &reactions, bool first_order,
                 bool limiting_fluxes,
                 const std::unordered_map<std::string, std::size_t> &index,
                 const StoichiometricMatrices nu, bool all_super, bool all_rini,
                 const std::vector<std::string> &species_order, bool debug);

} // namespace gasp2::catalysis::gamma
