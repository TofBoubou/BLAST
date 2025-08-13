#pragma once
#include "../core/containers.hpp"
#include "../core/exceptions.hpp"
#include "../io/config_types.hpp"
#include <expected>
#include <memory>
#include <span>
#include <string_view>
#include <vector>

namespace blast::thermophysics {

namespace constants {
inline constexpr double R_universal = 8.31446261815324;
inline constexpr double k_Boltzmann = 1.380649e-23;
inline constexpr double N_Avogadro = 6.02214076e23;
inline constexpr double Pi = 3.14159265358979323846;
} // namespace constants

// Error type for thermophysics operations
class ThermophysicsError : public core::BlastException {
public:
  explicit ThermophysicsError(std::string_view message, std::source_location location = std::source_location::current())
      : BlastException(std::format("Thermophysics Error: {}", message), location) {}
};

// Abstract interface for mixture properties
class MixtureInterface {
public:
  virtual ~MixtureInterface() = default;

  // Species information
  [[nodiscard]] virtual auto n_species() const noexcept -> std::size_t = 0;
  [[nodiscard]] virtual auto has_electrons() const noexcept -> bool = 0;
  [[nodiscard]] virtual auto species_name(std::size_t index) const noexcept -> std::string_view = 0;

  // Composition and molecular weight
  [[nodiscard]] virtual auto mixture_molecular_weight(std::span<const double> mass_fractions) const
      -> std::expected<double, ThermophysicsError> = 0;

  [[nodiscard]] virtual auto species_molecular_weight(std::size_t species_index) const noexcept -> double = 0;

  [[nodiscard]] virtual auto species_charges() const noexcept -> std::span<const double> = 0;

  // Composition conversions
  [[nodiscard]] virtual auto mass_fractions_to_mole_fractions(std::span<const double> mass_fractions) const
      -> std::expected<std::vector<double>, ThermophysicsError> = 0;

  // State setting
  [[nodiscard]] virtual auto set_state(std::span<const double> mass_fractions, double temperature,
                                       double pressure) const -> std::expected<void, ThermophysicsError> = 0;

  // Thermodynamic properties
  [[nodiscard]] virtual auto mixture_enthalpy(std::span<const double> mass_fractions, double temperature,
                                              double pressure) const -> std::expected<double, ThermophysicsError> = 0;

  [[nodiscard]] virtual auto
  species_enthalpies(double temperature) const -> std::expected<std::vector<double>, ThermophysicsError> = 0;

  [[nodiscard]] virtual auto frozen_cp(std::span<const double> mass_fractions, double temperature,
                                       double pressure) const -> std::expected<double, ThermophysicsError> = 0;

  // Transport properties
  [[nodiscard]] virtual auto viscosity(std::span<const double> mass_fractions, double temperature,
                                       double pressure) const -> std::expected<double, ThermophysicsError> = 0;

  [[nodiscard]] virtual auto
  frozen_thermal_conductivity(std::span<const double> mass_fractions, double temperature,
                              double pressure) const -> std::expected<double, ThermophysicsError> = 0;

  [[nodiscard]] virtual auto binary_diffusion_coefficients(double temperature, double pressure) const
      -> std::expected<core::Matrix<double>, ThermophysicsError> = 0;

  // Chemical properties
  [[nodiscard]] virtual auto production_rates(std::span<const double> partial_densities, double temperature) const
      -> std::expected<std::vector<double>, ThermophysicsError> = 0;

  [[nodiscard]] virtual auto
  production_rate_jacobian(std::span<const double> partial_densities,
                           double temperature) const -> std::expected<core::Matrix<double>, ThermophysicsError> = 0;

  // Thermal diffusion
  [[nodiscard]] virtual auto
  thermal_diffusion_ratios(std::span<const double> mass_fractions, double temperature,
                           double pressure) const -> std::expected<std::vector<double>, ThermophysicsError> = 0;

  // Equilibrium composition
  [[nodiscard]] virtual auto equilibrium_composition(double temperature, double pressure) const
      -> std::expected<std::vector<double>, ThermophysicsError> = 0;

  [[nodiscard]] virtual auto
  surface_reaction_rates(std::span<const double> partial_densities,
                         double wall_temperature) const -> std::expected<std::vector<double>, ThermophysicsError> = 0;

  [[nodiscard]] virtual auto reload_gsi() -> bool = 0;

  // Modal temperature support
  [[nodiscard]] virtual auto
  extract_modal_temperatures(std::span<const double> mass_fractions, double temperature_overall,
                             double pressure) const -> std::expected<std::vector<double>, ThermophysicsError> = 0;

  [[nodiscard]] virtual auto get_number_energy_modes() const noexcept -> std::size_t = 0;
};

// Factory function to create mixture implementation
[[nodiscard]] auto
create_mixture(const io::MixtureConfig& config) -> std::expected<std::unique_ptr<MixtureInterface>, ThermophysicsError>;

} // namespace blast::thermophysics