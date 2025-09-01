#include "gasp2/catalysis/gamma/rini_flux.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <unordered_set>

namespace gasp2::catalysis::gamma::rini {

// ========== PRECOMPUTATION ==========
/// Build lookup structures for the Rini flux solver.
///
/// Each reaction is mapped to the indices of its reactants while each species
/// stores the list of reactions it participates in. These static mappings are
/// used later during flux evaluation to traverse only relevant reactions for a
/// given species.
///
/// @param ns Number of gas-phase species provided at initialization.
/// @param nr Number of surface reactions.
/// @param index Mapping from species name to its index in the global ordering.
/// @param reactions Parsed reaction list.
///
/// @returns Populated `RiniMatrices` containing species-reaction connectivity.
///
/// @throws std::runtime_error if the closure condition \f$N + M = S\f$ is not
/// satisfied, where \f$N\f$ is the count of homogeneous reactions,
/// \f$M\f$ the count of heterogeneous reactions and \f$S\f$ the number of
/// participating species. The equality ensures that the system of unknown
/// surface fractions is neither under- nor overdetermined.
[[nodiscard]] RiniMatrices
precompute_matrices(std::size_t ns, std::size_t nr,
                    const std::unordered_map<std::string, std::size_t> &index,
                    const std::vector<ReactionInput> &reactions) {
  RiniMatrices matrices;
  matrices.gamma_w = reactions.empty() ? 0.0 : reactions[0].gamma_w.value();
  matrices.species_reactions.assign(ns, {});
  matrices.reaction_reactants.resize(nr);
  matrices.is_heterogeneous.assign(nr, false);
  std::unordered_set<std::size_t> reactant_species;
  std::size_t N = 0; ///< Number of homogeneous reactions.
  std::size_t M = 0; ///< Number of heterogeneous reactions.
  for (std::size_t r = 0; r < nr; ++r) {
    const auto &reaction = reactions[r];
    matrices.is_heterogeneous[r] = reaction.heterogeneous;
    if (reaction.heterogeneous) {
      ++M;
    } else {
      ++N;
    }
    matrices.reaction_reactants[r].reserve(reaction.reactants.size());
    for (const auto &[sp, _] : reaction.reactants) {
      std::size_t idx = index.at(sp);
      matrices.species_reactions[idx].push_back(r);
      matrices.reaction_reactants[r].push_back(idx);
      reactant_species.insert(idx);
    }
  }
  // ========== VALIDATION ==========
  // The Rini model introduces one unknown surface fraction (gamma) for each
  // homogeneous reaction and two for each heterogeneous reaction, yielding
  // \f$N + 2M\f$ unknowns. Species conservation provides \f$S\f$ equations and
  // heterogeneous reactions contribute \f$M\f$ compatibility relations, giving
  // \f$S + M\f$ equations in total. To avoid an under- or overdetermined system
  // we require \f$N + M = S\f$.
  std::size_t S = reactant_species.size();
  if (N + M != S) {
    throw std::runtime_error(
        "Rini flux: N + M must equal number of species (S) for closure");
  }
  return matrices;
}

// ========== FLUX SOLVER ==========
/// Solve the Rini flux model for a given wall state.
///
/// The precomputed connectivity matrices are combined with impinging fluxes
/// and stoichiometric coefficients to determine reaction rates and the
/// resulting catalytic flux for each species.
///
/// @param species Thermodynamic and kinetic properties of gas-phase species.
/// @param ns Number of species.
/// @param nr Number of reactions.
/// @param nu_p Reactant stoichiometric coefficients.
/// @param nu_diff Difference between product and reactant stoichiometry.
/// @param matrices Precomputed lookup tables from `precompute_matrices`.
///
/// @returns Catalytic fluxes for all species.
[[nodiscard]] CatalysisFluxes
compute_fluxes(const SpeciesData &species, std::size_t ns, std::size_t nr,
               const std::vector<std::vector<int>> &nu_p,
               const std::vector<std::vector<int>> &nu_diff,
               const RiniMatrices &matrices) {
  CatalysisFluxes fluxes;
  fluxes.cat_fluxes.assign(ns, 0.0);

  if (nr == 0)
    return fluxes;
  const double gamma_w = matrices.gamma_w;

  // ========== PREPARATION ==========
  // Gamma matrix: rows correspond to species, columns to reactions.
  std::vector<std::vector<double>> gamma_matrix(ns,
                                                std::vector<double>(nr, 0.0));
  // Assign gamma_w to species that participate in only one reaction and track
  // running gamma sums per species.
  std::vector<double> gamma_sums(ns, 0.0);
  for (std::size_t s = 0; s < ns; ++s) {
    if (matrices.species_reactions[s].size() == 1) {
      std::size_t r = matrices.species_reactions[s][0];
      gamma_matrix[s][r] = gamma_w;
      gamma_sums[s] = gamma_w;
    }
  }

  // ========== IMPINGING FLUXES ==========
  // 2/(2-gamma_w) is constant for all species; pre-compute the scale factor.
  const double imp_scale = 2.0 / (2.0 - gamma_w);
  std::vector<double> imp_flux(ns, 0.0);
  for (std::size_t i = 0; i < ns; ++i) {
    imp_flux[i] = imp_scale * species.imp_flux[i];
  }

  // Preallocate the reaction rates chi_r.
  std::vector<double> chi_r(nr, 0.0);

  // ========== HETEROGENEOUS REACTIONS ==========
  for (std::size_t r = 0; r < nr; ++r) {
    if (!matrices.is_heterogeneous[r]) {
      continue;
    }
    // Heterogeneous reactions have two reactants by construction.
    const auto &react_idx = matrices.reaction_reactants[r];
    std::size_t idx1 = react_idx[0];
    std::size_t idx2 = react_idx[1];
    std::size_t single_reactant = 0;
    std::size_t other_reactant = 0;
    // Identify the reactant that only appears in this reaction.
    if (matrices.species_reactions[idx1].size() == 1) {
      single_reactant = idx1;
      other_reactant = idx2;
    }
    if (matrices.species_reactions[idx2].size() == 1) {
      single_reactant = idx2;
      other_reactant = idx1;
    }
    // Compute the reaction rate using the unique reactant.
    // P.S. nu_p and nu_diff are indexed as [reaction][species].
    if (nu_p[r][single_reactant] == 0) {
      throw std::runtime_error("Rini flux: nu_p is zero for a reactant");
    }
    chi_r[r] = gamma_matrix[single_reactant][r] * imp_flux[single_reactant] /
               nu_p[r][single_reactant];
    // Determine the gamma for the second reactant and update its running sum.
    gamma_matrix[other_reactant][r] =
        chi_r[r] * nu_p[r][other_reactant] / imp_flux[other_reactant];
    gamma_sums[other_reactant] += gamma_matrix[other_reactant][r];
  }

  // ========== HOMOGENEOUS REACTIONS ==========
  for (std::size_t r = 0; r < nr; ++r) {
    if (matrices.is_heterogeneous[r]) {
      continue;
    }
    // Only one reactant exists for homogeneous recombination reactions.
    std::size_t species_index = matrices.reaction_reactants[r][0];
    if (nu_p[r][species_index] == 0) {
      throw std::runtime_error("Rini flux: nu_p is zero for a reactant");
    }
    if (matrices.species_reactions[species_index].size() == 1) {
      // Gamma for this species was already assigned to gamma_w during
      // preparation. Compute chi_r directly without modifying the matrix.
      chi_r[r] = gamma_matrix[species_index][r] * imp_flux[species_index] /
                 nu_p[r][species_index];
      continue;
    }
    double gamma = gamma_w - gamma_sums[species_index];
    gamma_matrix[species_index][r] = gamma;
    gamma_sums[species_index] += gamma;
    chi_r[r] = gamma_matrix[species_index][r] * imp_flux[species_index] /
               nu_p[r][species_index];
  }

  // ========== VALIDATION ==========
  for (std::size_t s = 0; s < ns; ++s) {
    double sum_gamma = 0.0;
    for (std::size_t r = 0; r < nr; ++r) {
      double gamma = gamma_matrix[s][r];
      if (gamma < 0.0 || gamma > 1.0) {
        throw std::runtime_error("Rini flux: gamma out of bounds");
      }
      sum_gamma += gamma;
    }
    constexpr double tol = 1e-8;
    if (std::abs(sum_gamma) > tol && std::abs(sum_gamma - gamma_w) > tol) {
      throw std::runtime_error("Rini flux: sum of gammas is not 0 or gamma_w");
    }
  }

  // ========== FLUX ASSEMBLY ==========
  for (std::size_t s = 0; s < ns; ++s) {
    for (std::size_t r = 0; r < nr; ++r) {
      fluxes.cat_fluxes[s] += chi_r[r] * nu_diff[r][s];
    }
    fluxes.cat_fluxes[s] *= species.m[s];
  }

  // The solver expects destruction rates, so flip the sign convention before
  // returning.
  for (auto &f : fluxes.cat_fluxes) {
    f = -f;
  }
  return fluxes;
}

} // namespace gasp2::catalysis::gamma::rini
