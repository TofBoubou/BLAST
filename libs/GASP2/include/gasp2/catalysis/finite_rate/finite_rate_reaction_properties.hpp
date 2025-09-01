#pragma once
#include "gasp2/catalysis/finite_rate/types.hpp"
#include "gasp2/gasp2.hpp"
#include <span>

namespace gasp2::catalysis::finite_rate {
namespace details {

//---------------------------------------------------------------------------
// Reaction property computation
//---------------------------------------------------------------------------
/**
 * @brief Compute rate coefficients and equilibrium constants for all reactions.
 *
 * The routine assembles adsorption/desorption pairs, derives equilibrium
 * constants, and collects auxiliary data needed by the finite-rate solver.
 *
 * @param species           Precomputed thermodynamic data for gas species.
 * @param gibbs_free_energy Gibbs free energy of gas species (J/mol) in the
 *                          ordering supplied during initialization.
 *
 * @return Populated ReactionProperties with placeholders for kinetic rates.
 */
[[nodiscard]] ReactionProperties
compute_reactions_properties(const gasp2::SpeciesData &species,
                             std::span<const double> gibbs_free_energy);

} // namespace details
} // namespace gasp2::catalysis::finite_rate
