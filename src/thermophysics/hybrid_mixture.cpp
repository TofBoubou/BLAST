#include "blast/thermophysics/hybrid_mixture.hpp"
#include <format>

namespace blast::thermophysics {

HybridMixture::HybridMixture(std::unique_ptr<MixtureInterface> base_mixture, 
                             std::unique_ptr<catalysis::CatalysisInterface> catalysis)
    : base_mixture_(std::move(base_mixture)), catalysis_(std::move(catalysis)) {
  
  if (!base_mixture_) {
    throw ThermophysicsError("Base mixture cannot be null");
  }
  
  if (!catalysis_) {
    throw ThermophysicsError("Catalysis provider cannot be null");
  }

  // Verify species consistency
  if (base_mixture_->n_species() != catalysis_->n_species()) {
    throw ThermophysicsError(std::format("Species count mismatch: base mixture has {} species, catalysis has {}", 
                                        base_mixture_->n_species(), catalysis_->n_species()));
  }
}

// ========== DELEGATE EVERYTHING TO BASE MIXTURE ==========

auto HybridMixture::n_species() const noexcept -> std::size_t {
  return base_mixture_->n_species();
}

auto HybridMixture::has_electrons() const noexcept -> bool {
  return base_mixture_->has_electrons();
}

auto HybridMixture::species_name(std::size_t index) const noexcept -> std::string_view {
  return base_mixture_->species_name(index);
}

auto HybridMixture::mixture_molecular_weight(std::span<const double> mass_fractions) const
    -> std::expected<double, ThermophysicsError> {
  return base_mixture_->mixture_molecular_weight(mass_fractions);
}

auto HybridMixture::species_molecular_weight(std::size_t species_index) const noexcept -> double {
  return base_mixture_->species_molecular_weight(species_index);
}

auto HybridMixture::species_charges() const noexcept -> std::span<const double> {
  return base_mixture_->species_charges();
}

auto HybridMixture::mass_fractions_to_mole_fractions(std::span<const double> mass_fractions) const
    -> std::expected<std::vector<double>, ThermophysicsError> {
  return base_mixture_->mass_fractions_to_mole_fractions(mass_fractions);
}

auto HybridMixture::set_state(std::span<const double> mass_fractions, double temperature, double pressure) const
    -> std::expected<void, ThermophysicsError> {
  return base_mixture_->set_state(mass_fractions, temperature, pressure);
}

auto HybridMixture::mixture_enthalpy(std::span<const double> mass_fractions, double temperature, double pressure) const
    -> std::expected<double, ThermophysicsError> {
  return base_mixture_->mixture_enthalpy(mass_fractions, temperature, pressure);
}

auto HybridMixture::species_enthalpies(double temperature) const -> std::expected<std::vector<double>, ThermophysicsError> {
  return base_mixture_->species_enthalpies(temperature);
}

auto HybridMixture::frozen_cp(std::span<const double> mass_fractions, double temperature, double pressure) const
    -> std::expected<double, ThermophysicsError> {
  return base_mixture_->frozen_cp(mass_fractions, temperature, pressure);
}

auto HybridMixture::viscosity(std::span<const double> mass_fractions, double temperature, double pressure) const
    -> std::expected<double, ThermophysicsError> {
  return base_mixture_->viscosity(mass_fractions, temperature, pressure);
}

auto HybridMixture::frozen_thermal_conductivity(std::span<const double> mass_fractions, double temperature, double pressure) const
    -> std::expected<double, ThermophysicsError> {
  return base_mixture_->frozen_thermal_conductivity(mass_fractions, temperature, pressure);
}

auto HybridMixture::binary_diffusion_coefficients(double temperature, double pressure) const
    -> std::expected<core::Matrix<double>, ThermophysicsError> {
  return base_mixture_->binary_diffusion_coefficients(temperature, pressure);
}

auto HybridMixture::production_rates(std::span<const double> partial_densities, double temperature) const
    -> std::expected<std::vector<double>, ThermophysicsError> {
  return base_mixture_->production_rates(partial_densities, temperature);
}

auto HybridMixture::production_rate_jacobian(std::span<const double> partial_densities, double temperature) const
    -> std::expected<core::Matrix<double>, ThermophysicsError> {
  return base_mixture_->production_rate_jacobian(partial_densities, temperature);
}

auto HybridMixture::thermal_diffusion_ratios(std::span<const double> mass_fractions, double temperature, double pressure) const
    -> std::expected<std::vector<double>, ThermophysicsError> {
  return base_mixture_->thermal_diffusion_ratios(mass_fractions, temperature, pressure);
}

auto HybridMixture::equilibrium_composition(double temperature, double pressure) const
    -> std::expected<std::vector<double>, ThermophysicsError> {
  return base_mixture_->equilibrium_composition(temperature, pressure);
}

auto HybridMixture::reload_gsi() -> std::expected<void, std::string> {
  return base_mixture_->reload_gsi();
}

auto HybridMixture::extract_modal_temperatures(std::span<const double> mass_fractions, double temperature_overall, double pressure) const
    -> std::expected<std::vector<double>, ThermophysicsError> {
  return base_mixture_->extract_modal_temperatures(mass_fractions, temperature_overall, pressure);
}

auto HybridMixture::get_number_energy_modes() const noexcept -> std::size_t {
  return base_mixture_->get_number_energy_modes();
}

// ========== OVERRIDE CATALYSIS TO USE SEPARATE PROVIDER ==========

auto HybridMixture::surface_reaction_rates(std::span<const double> partial_densities, double wall_temperature) const
    -> std::expected<std::vector<double>, ThermophysicsError> {
  
  auto result = catalysis_->compute_surface_fluxes(partial_densities, wall_temperature);
  if (!result) {
    return std::unexpected(ThermophysicsError(std::format("Catalysis computation failed: {}", result.error().message())));
  }
  
  return result.value();
}

} // namespace blast::thermophysics