#pragma once

#include <unordered_map>
#include <vector>

#include <gasp2/catalysis/gamma/types.hpp>
#include <gasp2/types.hpp>

namespace gasp2::catalysis::gamma::rini {

using gasp2::CatalysisFluxes;
using gasp2::SpeciesData;

/// Lookup matrices for the Rini model that do not depend on runtime state.
struct RiniMatrices {
  /// For each species, list of reaction indices where it appears.
  std::vector<std::vector<std::size_t>> species_reactions;
  /// For each reaction, indices of its reactant species.
  std::vector<std::vector<std::size_t>> reaction_reactants;
  /// Flags marking which reactions are heterogeneous.
  std::vector<bool> is_heterogeneous;
  /// Common wall scattering coefficient shared by all reactions.
  double gamma_w{0.0};
};

/// Pre-compute lookup matrices for the Rini model.
/// @param ns         Number of species (cached from the context).
/// @param nr         Number of reactions (cached from the context).
/// @param index      Mapping from species name to its index.
/// @param reactions  List of reactions using the Rini model.
[[nodiscard]] RiniMatrices
precompute_matrices(std::size_t ns, std::size_t nr,
                    const std::unordered_map<std::string, std::size_t> &index,
                    const std::vector<ReactionInput> &reactions);

/// Compute mass fluxes for the Rini recombination model.
/// @param species    Precomputed species properties.
/// @param ns         Number of species (cached from the context).
/// @param nr         Number of reactions (cached from the context).
/// @param nu_p       Stoichiometric coefficients for production terms.
/// @param nu_diff    Differential stoichiometric coefficients.
/// @param matrices   Precomputed lookup tables for the Rini model.
[[nodiscard]] CatalysisFluxes
compute_fluxes(const SpeciesData &species, std::size_t ns, std::size_t nr,
               const std::vector<std::vector<int>> &nu_p,
               const std::vector<std::vector<int>> &nu_diff,
               const RiniMatrices &matrices);

} // namespace gasp2::catalysis::gamma::rini
