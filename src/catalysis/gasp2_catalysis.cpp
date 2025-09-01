#include "blast/catalysis/gasp2_catalysis.hpp"
#include <gasp2/gasp2.hpp>
#include <format>
#include <iostream>

namespace blast::catalysis {

Gasp2Catalysis::Gasp2Catalysis(const std::string& xml_file, const std::vector<std::string>& species_order, 
                               const std::vector<double>& molar_masses)
      : xml_file_(xml_file), species_order_(species_order), molar_masses_(molar_masses), initialized_(false) {
  // Initialize GASP2 using the provided XML definition.
  auto init_result = gasp2::initialize_catalysis(species_order_, molar_masses_, xml_file_);
  if (!init_result) {
    throw CatalysisError(std::format("Failed to initialize GASP2: {}", init_result.error().message));
  }
  
  initialized_ = true;
}

auto Gasp2Catalysis::compute_surface_fluxes(std::span<const double> partial_densities, double wall_temperature) const
    -> std::expected<std::vector<double>, CatalysisError> {
  
  if (!initialized_) {
    return std::unexpected(CatalysisError("GASP2 not initialized"));
  }

  if (partial_densities.size() != species_order_.size()) {
    return std::unexpected(CatalysisError(std::format("Size mismatch: expected {} species, got {}", 
                                                     species_order_.size(), partial_densities.size())));
  }

  // Reinitialize GASP2 every call to ensure gamma values remain intact
  {
    std::cout << "[GASP2] Reinitializing catalysis with XML='" << xml_file_ << "'" << std::endl;
    auto init_result = gasp2::initialize_catalysis(species_order_, molar_masses_, xml_file_);
    if (!init_result) {
      return std::unexpected(CatalysisError(
          std::format("Failed to reinitialize GASP2: {}", init_result.error().message)));
    }
  }


  // Convert partial densities to vector for GASP2
  std::vector<double> rho_wall(partial_densities.begin(), partial_densities.end());

  // Compute fluxes using GASP2 (no Gibbs energies needed for finite rate without reversible reactions)
  auto flux_result = gasp2::compute_catalysis_fluxes(wall_temperature, rho_wall);
  if (!flux_result) {
    return std::unexpected(CatalysisError(std::format("GASP2 flux computation failed: {}", flux_result.error().message)));
  }

  return flux_result.value().cat_fluxes;
}

auto Gasp2Catalysis::species_names() const -> std::vector<std::string> {
  return species_order_;
}

auto Gasp2Catalysis::n_species() const noexcept -> std::size_t {
  return species_order_.size();
}

} // namespace blast::catalysis
