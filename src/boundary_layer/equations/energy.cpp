#include "blast/boundary_layer/equations/energy.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include "blast/boundary_layer/solvers/tridiagonal_solvers.hpp"
#include <algorithm>
#include <cmath>
#include <format>
#include <iomanip>
#include <vector>

namespace blast::boundary_layer::equations {

auto solve_energy(std::span<const double> g_previous, const coefficients::CoefficientInputs& inputs,
                  const coefficients::CoefficientSet& coeffs, const conditions::BoundaryConditions& bc,
                  const coefficients::XiDerivatives& xi_der, const io::SimulationConfig& sim_config,
                  std::span<const double> F_field, std::span<const double> dF_deta, std::span<const double> V_field,
                  const thermophysics::MixtureInterface& mixture, int station,
                  PhysicalQuantity auto d_eta) -> std::expected<std::vector<double>, EquationError> {

  const auto n_eta = g_previous.size();

  // Validation
  if (n_eta != F_field.size() || n_eta != V_field.size()) {
    return std::unexpected(EquationError("Energy: field size mismatch"));
  }

  // Build coefficients
  auto energy_coeffs_result = detail::build_energy_coefficients(g_previous, inputs, coeffs, bc, xi_der, sim_config,
                                                                F_field, dF_deta, V_field, mixture, station, d_eta);
  if (!energy_coeffs_result) {
    return std::unexpected(energy_coeffs_result.error());
  }
  auto energy_coeffs = energy_coeffs_result.value();

  // Build boundary conditions
  auto boundary_conds = detail::build_energy_boundary_conditions(inputs, coeffs, bc, sim_config, station, d_eta);

  // Solve tridiagonal system
  auto solution_result = solvers::solve_momentum_energy_tridiagonal(
      g_previous, boundary_conds.f_bc, boundary_conds.g_bc, boundary_conds.h_bc, energy_coeffs.a, energy_coeffs.b,
      energy_coeffs.c, energy_coeffs.d);

  if (!solution_result) {
    return std::unexpected(EquationError("Energy: tridiagonal solver failed: {}", std::source_location::current(),
                                         solution_result.error().message()));
  }

  return solution_result.value();
}

namespace detail {

[[nodiscard]] auto compute_dufour_term(const coefficients::CoefficientInputs& inputs,
                                       const coefficients::CoefficientSet& coeffs,
                                       const conditions::BoundaryConditions& bc,
                                       const thermophysics::MixtureInterface& mixture, std::size_t eta_index,
                                       PhysicalQuantity auto d_eta) -> std::expected<double, EquationError>;

auto build_energy_coefficients(std::span<const double> g_previous, const coefficients::CoefficientInputs& inputs,
                               const coefficients::CoefficientSet& coeffs, const conditions::BoundaryConditions& bc,
                               const coefficients::XiDerivatives& xi_der, const io::SimulationConfig& sim_config,
                               std::span<const double> F_field, std::span<const double> dF_deta,
                               std::span<const double> V_field, const thermophysics::MixtureInterface& mixture,
                               int station,
                               PhysicalQuantity auto d_eta) -> std::expected<EnergyCoefficients, EquationError> {

  const auto n_eta = g_previous.size();
  const auto n_species = inputs.c.rows();
  const double d_eta_sq = d_eta * d_eta;
  const double xi = inputs.xi;
  const double lambda0 = xi_der.lambda0();
  const auto g_derivatives = xi_der.g_derivative();

  // Compute geometry factor for diffusion fluxes
  const auto J_fact_result = compute_energy_j_factor(station, xi, bc, sim_config);
  if (!J_fact_result) {
    return std::unexpected(J_fact_result.error());
  }
  const double J_fact = J_fact_result.value();

  // Compute Dufour terms if enabled
  std::vector<double> dufour_terms;
  if (sim_config.consider_dufour_effect) {
    dufour_terms.reserve(n_eta);
    for (std::size_t i = 0; i < n_eta; ++i) {
      auto dufour_result = compute_dufour_term(inputs, coeffs, bc, mixture, i, d_eta);
      if (!dufour_result) {
        return std::unexpected(dufour_result.error());
      }
      dufour_terms.push_back(dufour_result.value());
    }

    auto dufour_derivative_result = coefficients::derivatives::compute_eta_derivative(dufour_terms, d_eta);
    if (!dufour_derivative_result) {
      return std::unexpected(EquationError("Failed to compute Dufour term derivative"));
    }
    dufour_terms = std::move(dufour_derivative_result.value());
  }

  EnergyCoefficients energy_coeffs;
  energy_coeffs.a.reserve(n_eta);
  energy_coeffs.b.reserve(n_eta);
  energy_coeffs.c.reserve(n_eta);
  energy_coeffs.d.reserve(n_eta);

  for (std::size_t i = 0; i < n_eta; ++i) {

    // ----- Coefficient a[i] -----
    double l3_i = coeffs.transport.l3[i];
    const double K_bl_sq = coeffs.transport.K_bl * coeffs.transport.K_bl;
    double a_i = l3_i * K_bl_sq / d_eta_sq;
    energy_coeffs.a.push_back(a_i);

    // ----- Coefficient b[i] -----
    double dl3_deta_i = coeffs.transport.dl3_deta[i];
    double V_i = V_field[i];
    double b_i = (dl3_deta_i * K_bl_sq - V_i) / d_eta;
    energy_coeffs.b.push_back(b_i);

    // ----- Coefficient c[i] -----
    double xi_i = xi;
    double F_i = F_field[i];
    double dhe_dxi = bc.d_he_dxi();
    double he = bc.he();
    double c_term = -2.0 * xi_i * F_i * dhe_dxi / he - 2.0 * xi_i * F_i * lambda0;
    energy_coeffs.c.push_back(c_term);

    // Compute species enthalpy terms
    auto [tmp1, tmp2] = compute_species_enthalpy_terms(inputs, coeffs, bc, J_fact, i);

    // Add Dufour effect contribution
    if (sim_config.consider_dufour_effect) {
      const double dufour_contribution = -bc.P_e() / bc.he() * dufour_terms[i];
      tmp2 += dufour_contribution;
    }

    // ----- Coefficient d[i] -----
    const double d_term = -bc.ue() * bc.ue() / bc.he() *
                              (coeffs.transport.l0[i] * dF_deta[i] * dF_deta[i] * K_bl_sq -
                               bc.beta * bc.rho_e() / coeffs.thermodynamic.rho[i] * F_field[i]) +
                          2.0 * xi * F_field[i] * g_derivatives[i] + tmp1 * K_bl_sq + tmp2 * coeffs.transport.K_bl;

    energy_coeffs.d.push_back(d_term);
  }

  return energy_coeffs;
}

[[nodiscard]] auto
build_energy_boundary_conditions(const coefficients::CoefficientSet& coeffs, const conditions::BoundaryConditions& bc,
                                 const io::SimulationConfig& sim_config, const thermophysics::MixtureInterface& mixture,
                                 PhysicalQuantity auto d_eta) -> EnergyBoundaryConditions {
  
  EnergyBoundaryConditions boundary_conds;
  
  if (sim_config.thermal_bc == io::SimulationConfig::ThermalBC::ImposedTemperature) {
    // Given temperature at the wall
    boundary_conds.f_bc = 0.0;
    boundary_conds.g_bc = 1.0;
    boundary_conds.h_bc = coeffs.thermodynamic.h_wall / bc.he();
  } else {
    // Adiabatic wall
    boundary_conds.f_bc = 1.0 / d_eta;
    boundary_conds.g_bc = 0.0;
    boundary_conds.h_bc = 0.0;
    
    // Pour l'instant, on laisse h_bc = 0 pour simplifier
    // Les termes de flux diffusifs seront ajoutés plus tard quand on aura
    // accès aux gradients de concentration et aux flux diffusifs
    // TODO: Ajouter les termes diffusifs comme dans l'ancien code
  }
  
  return boundary_conds;
}

auto compute_species_enthalpy_terms(const coefficients::CoefficientInputs& inputs,
                                    const coefficients::CoefficientSet& coeffs,
                                    const conditions::BoundaryConditions& bc, double J_fact,
                                    std::size_t eta_index) -> std::tuple<double, double> {

  const auto n_species = inputs.c.rows();
  double tmp1 = 0.0;
  double tmp2 = 0.0;

  for (std::size_t j = 0; j < n_species; ++j) {

    // tmp1: concentration and enthalpy derivative terms
    const double dc_deta_j = inputs.dc_deta(j, eta_index);
    const double dc_deta2_j = inputs.dc_deta2(j, eta_index);
    const double h_sp_j = coeffs.h_species(j, eta_index);
    const double dh_sp_deta_j = coeffs.dh_species_deta(j, eta_index);
    const double dl3_deta = coeffs.transport.dl3_deta[eta_index];
    const double l3 = coeffs.transport.l3[eta_index];

    tmp1 += dc_deta_j * h_sp_j / bc.he() * dl3_deta + l3 * h_sp_j / bc.he() * dc_deta2_j +
            l3 * dc_deta_j * dh_sp_deta_j / bc.he();

    // tmp2: diffusion flux terms
    const double J_j = coeffs.diffusion.J(j, eta_index);
    const double dJ_deta_j = coeffs.diffusion.dJ_deta(j, eta_index);

    tmp2 += J_j * dh_sp_deta_j / bc.he() + dJ_deta_j * h_sp_j / bc.he();
  }

  // Apply multiplication factors
  tmp2 *= J_fact;

  return {tmp1, tmp2};
}

[[nodiscard]] auto compute_dufour_term(const coefficients::CoefficientInputs& inputs,
                                       const coefficients::CoefficientSet& coeffs,
                                       const conditions::BoundaryConditions& bc,
                                       const thermophysics::MixtureInterface& mixture, std::size_t eta_index,
                                       PhysicalQuantity auto d_eta) -> std::expected<double, EquationError> {

  const auto n_species = inputs.c.rows();

  // Get mass fractions at this eta point
  std::vector<double> mass_fractions(n_species);
  for (std::size_t j = 0; j < n_species; ++j) {
    mass_fractions[j] = inputs.c(j, eta_index);
  }

  // Convert to mole fractions
  auto mole_fractions_result = mixture.mass_fractions_to_mole_fractions(mass_fractions);
  if (!mole_fractions_result) {
    return std::unexpected(EquationError("Failed to compute mole fractions for Dufour term"));
  }
  const auto& mole_fractions = mole_fractions_result.value();

  // Compute mixture molecular weight
  auto mixture_mw_result = mixture.mixture_molecular_weight(mass_fractions);
  if (!mixture_mw_result) {
    return std::unexpected(EquationError("Failed to compute mixture MW for Dufour term"));
  }
  [[maybe_unused]] const double mixture_mw = mixture_mw_result.value();

  // Calculate Dufour sum: Σ(χᵢ/ρᵢ * Jᵢ)
  double dufour_sum = 0.0;
  const double rho_total = coeffs.thermodynamic.rho[eta_index];

  for (std::size_t j = 0; j < n_species; ++j) {
    const double chi_j = mole_fractions[j];
    const double rho_j = rho_total * mass_fractions[j];
    const double J_j = coeffs.diffusion.J(j, eta_index);

    if (rho_j > 1e-15) {
      dufour_sum += chi_j / rho_j * J_j;
    }
  }

  return dufour_sum;
}

} // namespace detail

// Explicit instantiations for common use cases
template auto solve_energy<double>(std::span<const double> g_previous, const coefficients::CoefficientInputs& inputs,
                                   const coefficients::CoefficientSet& coeffs, const conditions::BoundaryConditions& bc,
                                   const coefficients::XiDerivatives& xi_der, const io::SimulationConfig& sim_config,
                                   std::span<const double> F_field, std::span<const double> dF_deta,
                                   std::span<const double> V_field, const thermophysics::MixtureInterface& mixture,
                                   int station, double d_eta) -> std::expected<std::vector<double>, EquationError>;

} // namespace blast::boundary_layer::equations