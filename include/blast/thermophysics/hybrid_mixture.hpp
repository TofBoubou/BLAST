#pragma once
#include "mixture_interface.hpp"
#include "../catalysis/catalysis_interface.hpp"
#include <memory>

namespace blast::thermophysics {

// Hybrid mixture that combines base mixture (thermo/transport) with separate catalysis provider
class HybridMixture : public MixtureInterface {
public:
  // Constructor
  HybridMixture(std::unique_ptr<MixtureInterface> base_mixture, 
                std::unique_ptr<catalysis::CatalysisInterface> catalysis);

  // Destructor
  ~HybridMixture() override = default;

  // Disable copy/move for simplicity  
  HybridMixture(const HybridMixture&) = delete;
  HybridMixture& operator=(const HybridMixture&) = delete;
  HybridMixture(HybridMixture&&) = default;
  HybridMixture& operator=(HybridMixture&&) = default;

  // ========== DELEGATE EVERYTHING TO BASE MIXTURE ==========
  [[nodiscard]] auto n_species() const noexcept -> std::size_t override;
  [[nodiscard]] auto has_electrons() const noexcept -> bool override;
  [[nodiscard]] auto species_name(std::size_t index) const noexcept -> std::string_view override;

  [[nodiscard]] auto mixture_molecular_weight(std::span<const double> mass_fractions) const
      -> std::expected<double, ThermophysicsError> override;

  [[nodiscard]] auto species_molecular_weight(std::size_t species_index) const noexcept -> double override;

  [[nodiscard]] auto species_charges() const noexcept -> std::span<const double> override;

  [[nodiscard]] auto mass_fractions_to_mole_fractions(std::span<const double> mass_fractions) const
      -> std::expected<std::vector<double>, ThermophysicsError> override;

  [[nodiscard]] auto set_state(std::span<const double> mass_fractions, double temperature,
                               double pressure) const -> std::expected<void, ThermophysicsError> override;

  [[nodiscard]] auto mixture_enthalpy(std::span<const double> mass_fractions, double temperature,
                                      double pressure) const -> std::expected<double, ThermophysicsError> override;

  [[nodiscard]] auto species_enthalpies(double temperature) const -> std::expected<std::vector<double>, ThermophysicsError> override;

  [[nodiscard]] auto frozen_cp(std::span<const double> mass_fractions, double temperature,
                               double pressure) const -> std::expected<double, ThermophysicsError> override;

  [[nodiscard]] auto viscosity(std::span<const double> mass_fractions, double temperature,
                               double pressure) const -> std::expected<double, ThermophysicsError> override;

  [[nodiscard]] auto frozen_thermal_conductivity(std::span<const double> mass_fractions, double temperature,
                                                 double pressure) const -> std::expected<double, ThermophysicsError> override;

  [[nodiscard]] auto binary_diffusion_coefficients(double temperature, double pressure) const
      -> std::expected<core::Matrix<double>, ThermophysicsError> override;

  [[nodiscard]] auto production_rates(std::span<const double> partial_densities, double temperature) const
      -> std::expected<std::vector<double>, ThermophysicsError> override;

  [[nodiscard]] auto production_rate_jacobian(std::span<const double> partial_densities, double temperature) const
      -> std::expected<core::Matrix<double>, ThermophysicsError> override;

  [[nodiscard]] auto thermal_diffusion_ratios(std::span<const double> mass_fractions, double temperature,
                                              double pressure) const -> std::expected<std::vector<double>, ThermophysicsError> override;

  [[nodiscard]] auto equilibrium_composition(double temperature, double pressure) const
      -> std::expected<std::vector<double>, ThermophysicsError> override;

  [[nodiscard]] auto reload_gsi() -> std::expected<void, std::string> override;

  [[nodiscard]] auto extract_modal_temperatures(std::span<const double> mass_fractions, double temperature_overall,
                                                double pressure) const -> std::expected<std::vector<double>, ThermophysicsError> override;

  [[nodiscard]] auto get_number_energy_modes() const noexcept -> std::size_t override;

  // ========== OVERRIDE CATALYSIS TO USE SEPARATE PROVIDER ==========
  [[nodiscard]] auto surface_reaction_rates(std::span<const double> partial_densities, double wall_temperature) const
      -> std::expected<std::vector<double>, ThermophysicsError> override;

private:
  std::unique_ptr<MixtureInterface> base_mixture_;
  std::unique_ptr<catalysis::CatalysisInterface> catalysis_;
};

} // namespace blast::thermophysics