#pragma once
#include "types.hpp"

namespace gasp2 {
//---------------------------------------------------------------------------
// This file is the main facade for the GASP2 library.
//---------------------------------------------------------------------------

/**
 * @brief Initialize the catalysis solver.
 *
 * Parses the input file, performs validation independent of the wall state and
 * stores all data required for subsequent flux computations.
 * @param species_order Order of species in the input density array.
 * @param molar_masses Molar masses of the species (kg/mol).
 * @param input_filename Path to the input XML file.
 * @param debug Enable verbose debug printing.
 * @return void on success or an error describing invalid inputs or
 *         initialization failures.
 * @todo Check+implement finite rate
 */
[[nodiscard]] Result<void>
initialize_catalysis(const std::vector<std::string> &species_order,
                     std::span<const double> molar_masses,
                     const std::string &input_filename, bool debug = true);

/**
 * @brief Compute catalysis fluxes using the previously initialized data.
 * An error is returned if the initialization was not successful.
 * The returned fluxes are in [Kg/(m^2*s)] and they are "destruction" rates.
 * For finite-rate chemistry the result also reports the fractional
 * contribution of each reaction to the destruction of every species
 * (\f$w_{ir}/w_i\f$).
 *
 * @param T_wall Wall temperature (K).
 * @param rho_wall Species mass densities at the wall (kg/m^3).
 * @param gibbs_free_energy Gibbs free energies per mole for each species
 *                          (J/mol). Required only when using the finite-rate
 *                          model; ignored otherwise.
 * @return CatalysisFluxes in the order provided during initialization.
 */
[[nodiscard]] Result<CatalysisFluxes>
compute_catalysis_fluxes(double T_wall, std::span<const double> rho_wall,
                         std::span<const double> gibbs_free_energy = {});

/**
 * @brief Retrieve current catalysis library status.
 *
 * Provides a minimal, read-only view of the internal catalysis state so
 * callers can confirm whether initialization succeeded and which modeling
 * options are active.
 *
 * @return CatalysisStatus describing the initialization flag and selected
 *         surface and reactions models.
 */
[[nodiscard]] CatalysisStatus get_catalysis_status();

// (Future) other APIs: ablation, full coupling, etc.

} // namespace gasp2
