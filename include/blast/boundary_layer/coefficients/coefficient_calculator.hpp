#pragma once
#include "../../io/config_types.hpp"
#include "../../thermophysics/mixture_interface.hpp"
#include "../conditions/boundary_conditions.hpp"
#include "coefficient_types.hpp"
#include "xi_derivatives.hpp"
#include <concepts>
#include <expected>

namespace blast::boundary_layer::coefficients {

// Derivative computation utilities
namespace derivatives {

template <std::ranges::sized_range Range>
[[nodiscard]] auto compute_eta_derivative(Range&& values,
                                          double d_eta) -> std::expected<std::vector<double>, CoefficientError>;

template <std::ranges::sized_range Range>
[[nodiscard]] auto compute_eta_second_derivative(Range&& values,
                                                 double d_eta) -> std::expected<std::vector<double>, CoefficientError>;

template <typename Matrix>
[[nodiscard]] auto compute_matrix_eta_derivative(const Matrix& values, double d_eta) -> Matrix;

template <typename Matrix>
[[nodiscard]] auto compute_matrix_eta_second_derivative(const Matrix& values,
                                                        double d_eta) -> std::expected<Matrix, CoefficientError>;
} // namespace derivatives

class CoefficientCalculator {
private:
  const thermophysics::MixtureInterface& mixture_;
  const io::SimulationConfig& sim_config_;
  const io::NumericalConfig& num_config_;
  double d_eta_;

  [[nodiscard]] auto calculate_transport_coefficients(
      const CoefficientInputs& inputs, const ThermodynamicCoefficients& thermo,
      const conditions::BoundaryConditions& bc) const -> std::expected<TransportCoefficients, CoefficientError>;

  [[nodiscard]] auto calculate_thermodynamic_coefficients(const CoefficientInputs& inputs,
                                                          const conditions::BoundaryConditions& bc) const
      -> std::expected<ThermodynamicCoefficients, CoefficientError>;

  [[nodiscard]] auto calculate_diffusion_coefficients(
      const CoefficientInputs& inputs, const conditions::BoundaryConditions& bc,
      const XiDerivatives& xi_der) const -> std::expected<DiffusionCoefficients, CoefficientError>;

  [[nodiscard]] auto calculate_chemical_coefficients(
      const CoefficientInputs& inputs, const ThermodynamicCoefficients& thermo,
      const conditions::BoundaryConditions& bc) const -> std::expected<ChemicalCoefficients, CoefficientError>;

  [[nodiscard]] auto calculate_thermal_diffusion(
      const CoefficientInputs& inputs, const ThermodynamicCoefficients& thermo,
      const conditions::BoundaryConditions& bc) const -> std::expected<ThermalDiffusionCoefficients, CoefficientError>;

  [[nodiscard]] auto calculate_wall_properties(
      const CoefficientInputs& inputs, const conditions::BoundaryConditions& bc, const TransportCoefficients& transport,
      const ThermodynamicCoefficients& thermo) const -> std::expected<WallProperties, CoefficientError>;

  [[nodiscard]] auto calculate_species_enthalpies(const CoefficientInputs& inputs) const
      -> std::expected<std::pair<core::Matrix<double>, core::Matrix<double>>, CoefficientError>;

  // Helper functions for chemical coefficients calculation
  [[nodiscard]] auto calculate_station_coefficients(std::size_t station_index, const CoefficientInputs& inputs,
                                                    const ThermodynamicCoefficients& thermo,
                                                    const conditions::BoundaryConditions& bc) const
      -> std::expected<std::pair<std::vector<double>, core::Matrix<double>>, CoefficientError>;

  [[nodiscard]] auto extract_local_composition(std::size_t station_index, const core::Matrix<double>& c,
                                               double rho_total) const
      -> std::expected<std::pair<std::vector<double>, std::vector<double>>, CoefficientError>;

  [[nodiscard]] auto transform_jacobian(
      std::size_t station_index, const CoefficientInputs& inputs, const ThermodynamicCoefficients& thermo,
      const conditions::BoundaryConditions& bc, const std::vector<double>& c_local,
      const std::vector<double>& rho_species, const std::vector<double>& wi_local,
      const core::Matrix<double>& dwi_drho) const -> std::expected<core::Matrix<double>, CoefficientError>;

  [[nodiscard]] auto compute_density_derivatives(const std::vector<double>& rho_species, double rho_total,
                                                 double MW_mixture) const -> core::Matrix<double>;

  [[nodiscard]] auto compute_temperature_derivatives(const std::vector<double>& rho_species, double T,
                                                     const std::vector<double>& wi_base) const
      -> std::expected<std::vector<double>, CoefficientError>;

public:
  CoefficientCalculator(const thermophysics::MixtureInterface& mixture, const io::SimulationConfig& sim_config,
                        const io::NumericalConfig& num_config)
      : mixture_(mixture), sim_config_(sim_config), num_config_(num_config),
        d_eta_(num_config.eta_max / static_cast<double>(num_config.n_eta - 1)) {}

  // Main calculation interface
  [[nodiscard]] auto calculate(const CoefficientInputs& inputs, const conditions::BoundaryConditions& bc,
                               const XiDerivatives& xi_der) const -> std::expected<CoefficientSet, CoefficientError>;
};

} // namespace blast::boundary_layer::coefficients