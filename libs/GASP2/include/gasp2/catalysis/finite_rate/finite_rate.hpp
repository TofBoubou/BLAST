//===----------------------------------------------------------------------===//
// Finite-rate model data structures and initialization helpers
//
// These utilities mirror the gamma-model infrastructure but are tailored for
// the finite-rate approach. The input parser extracts reaction information once
// during initialization and stores it in compact structures that are later
// consumed by the (yet to be implemented) flux computation routines.
//===----------------------------------------------------------------------===//

#pragma once
#include "gasp2/catalysis/finite_rate/types.hpp"
#include "gasp2/gasp2.hpp"
#include <filesystem>
#include <span>
#include <string>
#include <vector>

namespace gasp2::catalysis::finite_rate {

//--------------------------- Initialization helpers -------------------------//

/// Initialize internal data for the finite-rate model. The routine reads the
/// XML input file, performs static validation (mass/charge/site balance and
/// adsorption/desorption pairing), and caches all data required for subsequent
/// flux computations.
[[nodiscard]] Result<void>
initialize(const std::vector<std::string> &species_order,
           std::span<const double> molar_masses,
           const std::filesystem::path &input_filename, bool debug);

/// Placeholder for the future flux computation routine. The function is already
/// declared so that higher-level code can compile while the actual algorithm is
/// developed. For finite-rate reactions, the Gibbs free energy per mole of each
/// species must also be supplied.
[[nodiscard]] Result<CatalysisFluxes>
compute_fluxes(double T_wall, std::span<const double> rho_wall,
               std::span<const double> gibbs_free_energy);

} // namespace gasp2::catalysis::finite_rate
