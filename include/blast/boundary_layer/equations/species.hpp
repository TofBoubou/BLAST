#pragma once
#include "../../io/config_types.hpp"
#include "../../thermophysics/mixture_interface.hpp"
#include "../coefficients/coefficient_types.hpp"
#include "../coefficients/xi_derivatives.hpp"
#include "../conditions/boundary_conditions.hpp"
#include "equation_types.hpp"
#include "geometry_factors.hpp"
#include <expected>
#include <span>

namespace blast::boundary_layer::equations {

// Solve species continuity equations: block tridiagonal system for c (mass
// fractions) Equations: d/dη[J_j] + chemistry_terms = 0 for each species j
// Boundary conditions: equilibrium wall
[[nodiscard]] auto solve_species(const core::Matrix<double>& c_previous, const coefficients::CoefficientInputs& inputs,
                                 const coefficients::CoefficientSet& coeffs, const conditions::BoundaryConditions& bc,
                                 const coefficients::XiDerivatives& xi_der,
                                 const thermophysics::MixtureInterface& mixture, const io::SimulationConfig& sim_config,
                                 std::span<const double> F_field, std::span<const double> V_field, int station,
                                 PhysicalQuantity auto d_eta) -> std::expected<core::Matrix<double>, EquationError>;

// Compute equilibrium composition for the entire field
[[nodiscard]] auto compute_equilibrium_composition(std::span<const double> temperature_field, double pressure,
                                                   const thermophysics::MixtureInterface& mixture)
    -> std::expected<core::Matrix<double>, EquationError>;

// Apply charge neutrality constraint for ionized mixtures
auto apply_charge_neutrality(core::Matrix<double>& species_matrix,
                             const thermophysics::MixtureInterface& mixture) -> void;

namespace detail {

struct SpeciesCoefficients {
  core::Matrix<double> a; // Lower diagonal coefficients [n_eta x n_species]
  core::Matrix<double> b; // Main diagonal coefficients [n_eta x n_species]
  core::Matrix<double> c; // Upper diagonal coefficients [n_eta x n_species]
  core::Matrix<double> d; // Right-hand side [n_eta x n_species]
};

struct SpeciesBoundaryConditions {
  std::vector<double> f_bc; // Wall BC coefficients per species
  std::vector<double> g_bc; // Wall BC coefficients per species
  std::vector<double> h_bc; // Wall BC RHS per species
};

[[nodiscard]] auto
build_species_coefficients(const core::Matrix<double>& c_previous, const coefficients::CoefficientInputs& inputs,
                           const coefficients::CoefficientSet& coeffs, const conditions::BoundaryConditions& bc,
                           const coefficients::XiDerivatives& xi_der, const thermophysics::MixtureInterface& mixture,
                           const io::SimulationConfig& sim_config, std::span<const double> F_field,
                           std::span<const double> V_field, int station,
                           PhysicalQuantity auto d_eta) -> std::expected<SpeciesCoefficients, EquationError>;

[[nodiscard]] auto build_species_boundary_conditions(
    const core::Matrix<double>& c_wall, const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc, const thermophysics::MixtureInterface& mixture,
    const io::SimulationConfig& sim_config,
    PhysicalQuantity auto d_eta) -> std::expected<SpeciesBoundaryConditions, EquationError>;

// Compute fake fluxes for Le/Pr approximation
[[nodiscard]] auto
compute_fake_fluxes(const core::Matrix<double>& dc_deta_fixed, const coefficients::CoefficientSet& coeffs,
                    PhysicalQuantity auto d_eta,
                    double Le = 1.2, // Lewis number
                    double Pr = 0.72 // Prandtl number
                    ) -> std::expected<std::pair<core::Matrix<double>, core::Matrix<double>>, EquationError>;

// Fix concentration derivatives for mass conservation
[[nodiscard]] auto fix_concentration_derivatives(const core::Matrix<double>& c_matrix,
                                                 const core::Matrix<double>& dc_deta) -> core::Matrix<double>;

// Compute equilibrium wall composition
[[nodiscard]] auto compute_equilibrium_wall(const conditions::BoundaryConditions& bc,
                                            const thermophysics::MixtureInterface& mixture)
    -> std::expected<std::vector<double>, EquationError>;

[[nodiscard]] auto build_catalytic_boundary_conditions(
    const core::Matrix<double>& c_wall, const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc, const thermophysics::MixtureInterface& mixture,
    const io::SimulationConfig& sim_config,
    PhysicalQuantity auto d_eta) -> std::expected<SpeciesBoundaryConditions, EquationError>;

} // namespace detail

} // namespace blast::boundary_layer::equations