#pragma once
#include "../../io/config_types.hpp"
#include "../coefficients/coefficient_types.hpp"
#include "../coefficients/xi_derivatives.hpp"
#include "../conditions/boundary_conditions.hpp"
#include "../../thermophysics/mixture_interface.hpp"
#include "equation_types.hpp"
#include "geometry_factors.hpp"
#include <expected>
#include <span>

namespace blast::boundary_layer::equations {

// Solve energy equation: tridiagonal system for g (dimensionless enthalpy)
// Equation: d/dη[l3 * dg/dη] - V * dg/dη + energy_source_terms = 0
// Boundary conditions: specified temperature or adiabatic wall
[[nodiscard]] auto solve_energy(std::span<const double> g_previous, const coefficients::CoefficientInputs& inputs,
                                const coefficients::CoefficientSet& coeffs, const conditions::BoundaryConditions& bc,
                                const coefficients::XiDerivatives& xi_der, const io::SimulationConfig& sim_config,
                                std::span<const double> F_field, std::span<const double> dF_deta,
                                std::span<const double> V_field, const thermophysics::MixtureInterface& mixture,
                                int station,
                                PhysicalQuantity auto d_eta) -> std::expected<std::vector<double>, EquationError>;

namespace detail {

struct EnergyCoefficients {
  std::vector<double> a; // Lower diagonal
  std::vector<double> b; // Main diagonal
  std::vector<double> c; // Upper diagonal
  std::vector<double> d; // Right-hand side
};

struct EnergyBoundaryConditions {
  double f_bc; // Wall BC coefficient
  double g_bc; // Wall BC coefficient
  double h_bc; // Wall BC RHS
};

[[nodiscard]] auto
build_energy_coefficients(std::span<const double> g_previous, const coefficients::CoefficientInputs& inputs,
                          const coefficients::CoefficientSet& coeffs, const conditions::BoundaryConditions& bc,
                          const coefficients::XiDerivatives& xi_der, const io::SimulationConfig& sim_config,
                          std::span<const double> F_field, std::span<const double> dF_deta,
                          std::span<const double> V_field, const thermophysics::MixtureInterface& mixture,
                          int station,
                          PhysicalQuantity auto d_eta) -> std::expected<EnergyCoefficients, EquationError>;

[[nodiscard]] auto build_energy_boundary_conditions(const coefficients::CoefficientInputs& inputs,
                                                    const coefficients::CoefficientSet& coeffs,
                                                    const conditions::BoundaryConditions& bc,
                                                    const io::SimulationConfig& sim_config, int station,
                                                    PhysicalQuantity auto d_eta) -> EnergyBoundaryConditions;

// Helper to compute species enthalpy terms
[[nodiscard]] auto compute_species_enthalpy_terms(const coefficients::CoefficientInputs& inputs,
                                                  const coefficients::CoefficientSet& coeffs,
                                                  const conditions::BoundaryConditions& bc, double J_fact,
                                                  std::size_t eta_index) -> std::tuple<double, double>;

} // namespace detail

} // namespace blast::boundary_layer::equations