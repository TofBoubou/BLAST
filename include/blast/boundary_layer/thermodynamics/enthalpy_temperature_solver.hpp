#pragma once
#include "../../thermophysics/mixture_interface.hpp"
#include "../conditions/boundary_conditions.hpp"
#include "thermodynamic_types.hpp"
#include <span>

namespace blast::boundary_layer::thermodynamics {

class EnthalpyTemperatureSolver {
private:
  const thermophysics::MixtureInterface& mixture_;
  const EnthalpyTemperatureSolverConfig config_;

  [[nodiscard]] auto
  solve_single_point(std::span<const double> composition, double target_enthalpy, double pressure,
                     double initial_temperature) const -> std::expected<double, ThermodynamicSolverError>;

  [[nodiscard]] auto brent_method(std::span<const double> composition, double target_enthalpy, double pressure,
                                  double temp_min,
                                  double temp_max) const -> std::expected<double, ThermodynamicSolverError>;

public:
  EnthalpyTemperatureSolver(const thermophysics::MixtureInterface& mixture,
                            const EnthalpyTemperatureSolverConfig& config = {})
      : mixture_(mixture), config_(config) {}

  [[nodiscard]] auto solve(std::span<const double> enthalpy_field, const core::Matrix<double>& composition,
                           const conditions::BoundaryConditions& bc, std::span<const double> initial_temperatures) const
      -> std::expected<TemperatureField, ThermodynamicSolverError>;
};

} // namespace blast::boundary_layer::thermodynamics