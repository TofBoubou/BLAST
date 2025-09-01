#pragma once
#include "gasp2/types.hpp"
#include <filesystem>

namespace gasp2::catalysis::gamma {

using gasp2::CatalysisFluxes;
using gasp2::Result;

// Initialize internal data for the gamma-model catalysis solver. This parses
// the input file, performs static validation, and stores all persistent data.
// Species names, ordering, and molar masses must be validated by the caller
// prior to invocation.
[[nodiscard]] Result<void>
initialize(const std::vector<std::string> &species_order,
           std::span<const double> molar_masses,
           const std::filesystem::path &input_filename, bool debug);

// Compute catalysis fluxes using the data prepared during initialization. Wall
// temperature and species densities must be validated by the caller. The wall
// state may change between calls.
[[nodiscard]] Result<CatalysisFluxes>
compute_fluxes(double T_wall, std::span<const double> rho_wall);

} // namespace gasp2::catalysis::gamma
