#include "blast/boundary_layer/equations/species.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include "blast/boundary_layer/solvers/tridiagonal_solvers.hpp"
#include <algorithm>
#include <cmath>
#include <format>
#include <iomanip>
#include <iostream>
#include <numeric>

namespace blast::boundary_layer::equations {

auto solve_species(const core::Matrix<double>& c_previous, const coefficients::CoefficientInputs& inputs,
                   const coefficients::CoefficientSet& coeffs, const conditions::BoundaryConditions& bc,
                   const coefficients::XiDerivatives& xi_der, const thermophysics::MixtureInterface& mixture,
                   const io::SimulationConfig& sim_config, std::span<const double> F_field,
                   std::span<const double> V_field, int station, PhysicalQuantity auto d_eta)
    -> std::expected<core::Matrix<double>, EquationError> {

  const auto n_eta = c_previous.cols();
  const auto n_species = c_previous.rows();

  if (n_species != mixture.n_species()) {
    return std::unexpected(EquationError("Species: mixture species count mismatch"));
  }

  if (n_species == 1) {
    core::Matrix<double> result(1, n_eta);
    result.setOnes();
    return result;
  }

  if (sim_config.chemical_mode == io::SimulationConfig::ChemicalMode::Equilibrium) {
    return compute_equilibrium_composition(inputs.T, bc.P_e(), mixture);
  }

  if (sim_config.chemical_mode == io::SimulationConfig::ChemicalMode::Frozen) {
  }

  auto species_coeffs_result = detail::build_species_coefficients(c_previous, inputs, coeffs, bc, xi_der, mixture,
                                                                  sim_config, F_field, V_field, station, d_eta);

  if (!species_coeffs_result) {
    return std::unexpected(species_coeffs_result.error());
  }
  auto species_coeffs = species_coeffs_result.value();

  auto boundary_result = detail::build_species_boundary_conditions(c_previous, coeffs, bc, mixture, sim_config, d_eta);

  if (!boundary_result) {
    return std::unexpected(boundary_result.error());
  }
  auto boundary_conds = boundary_result.value();

  auto solution_result = solvers::solve_species_block_tridiagonal(
      c_previous, boundary_conds.f_bc, boundary_conds.g_bc, boundary_conds.h_bc, species_coeffs.a, species_coeffs.b,
      species_coeffs.c, species_coeffs.d, mixture.has_electrons());

  if (!solution_result) {
    return std::unexpected(EquationError("Species: block tridiagonal solver failed: {}",
                                         std::source_location::current(), solution_result.error().message()));
  }

  auto result = solution_result.value();

  if (mixture.has_electrons()) {
    apply_charge_neutrality(result, mixture);
  }

  return result;
}

auto compute_equilibrium_composition(std::span<const double> temperature_field, double pressure,
                                     const thermophysics::MixtureInterface& mixture)
    -> std::expected<core::Matrix<double>, EquationError> {

  const auto n_eta = temperature_field.size();
  const auto n_species = mixture.n_species();

  core::Matrix<double> c_equilibrium(n_species, n_eta);

  for (std::size_t i = 0; i < n_eta; ++i) {
    auto eq_result = mixture.equilibrium_composition(temperature_field[i], pressure);
    if (!eq_result) {
      return std::unexpected(EquationError("Failed to compute equilibrium at eta={}: {}",
                                           std::source_location::current(), i, eq_result.error().message()));
    }

    const auto& eq_composition = eq_result.value();
    for (std::size_t j = 0; j < n_species; ++j) {
      c_equilibrium(j, i) = eq_composition[j];
    }
  }

  return c_equilibrium;
}

auto apply_charge_neutrality(core::Matrix<double>& species_matrix, const thermophysics::MixtureInterface& mixture)
    -> void {

  if (!mixture.has_electrons())
    return;

  const auto n_eta = species_matrix.cols();
  const auto n_species = species_matrix.rows();
  auto charges = mixture.species_charges();

  // Electrons are typically the first species (index 0)
  for (std::size_t i = 0; i < n_eta; ++i) {

    // Get mass fractions for this eta point
    std::vector<double> mass_fractions(n_species);
    for (std::size_t j = 0; j < n_species; ++j) {
      mass_fractions[j] = species_matrix(j, i);
    }

    // Convert mass fractions to mole fractions
    auto mole_fractions_result = mixture.mass_fractions_to_mole_fractions(mass_fractions);
    if (!mole_fractions_result) {
      // If conversion fails, fall back to original (incorrect) method
      double charge_sum = 0.0;
      for (std::size_t j = 1; j < n_species; ++j) {
        charge_sum += species_matrix(j, i) * charges[j];
      }
      if (std::abs(charges[0]) > 1e-15) {
        species_matrix(0, i) = -charge_sum / charges[0];
      }
      continue;
    }
    auto mole_fractions = mole_fractions_result.value();

    // Calculate mixture molecular weight for molar concentration conversion
    auto mixture_mw_result = mixture.mixture_molecular_weight(mass_fractions);
    if (!mixture_mw_result) {
      // Fall back to original method if MW calculation fails
      double charge_sum = 0.0;
      for (std::size_t j = 1; j < n_species; ++j) {
        charge_sum += species_matrix(j, i) * charges[j];
      }
      if (std::abs(charges[0]) > 1e-15) {
        species_matrix(0, i) = -charge_sum / charges[0];
      }
      continue;
    }
    double mixture_mw = mixture_mw_result.value();

    // Convert mole fractions to molar concentrations
    // Assume total molar density = 1.0 mol/m³ (relative calculation)
    const double total_molar_density = 1.0;
    std::vector<double> molar_concentrations(n_species);
    for (std::size_t j = 0; j < n_species; ++j) {
      molar_concentrations[j] = mole_fractions[j] * total_molar_density;
    }

    // Apply charge neutrality: ∑(c_i × z_i) = 0
    // Calculate charge sum from all non-electron species
    double charge_sum = 0.0;
    for (std::size_t j = 1; j < n_species; ++j) {
      // Convert charge from C/kg to elementary charges
      const double species_mw = mixture.species_molecular_weight(j);
      const double elementary_charge_per_mol = charges[j] * species_mw / (1.602176634e-19 * 6.02214076e23);
      charge_sum += molar_concentrations[j] * elementary_charge_per_mol;
    }

    // Set electron molar concentration to ensure neutrality
    const double electron_mw = mixture.species_molecular_weight(0);
    const double electron_elementary_charge_per_mol = charges[0] * electron_mw / (1.602176634e-19 * 6.02214076e23);
    if (std::abs(electron_elementary_charge_per_mol) > 1e-15) {
      molar_concentrations[0] = -charge_sum / electron_elementary_charge_per_mol;
    }

    // Convert back to mole fractions
    double total_moles = 0.0;
    for (std::size_t j = 0; j < n_species; ++j) {
      total_moles += molar_concentrations[j];
    }
    if (total_moles > 1e-15) {
      for (std::size_t j = 0; j < n_species; ++j) {
        mole_fractions[j] = molar_concentrations[j] / total_moles;
      }
    }

    // Convert mole fractions back to mass fractions
    // Y_i = X_i * M_i / MW_mix
    double new_mixture_mw = 0.0;
    for (std::size_t j = 0; j < n_species; ++j) {
      new_mixture_mw += mole_fractions[j] * mixture.species_molecular_weight(j);
    }

    for (std::size_t j = 0; j < n_species; ++j) {
      mass_fractions[j] = mole_fractions[j] * mixture.species_molecular_weight(j) / new_mixture_mw;
    }

    // Renormalize to ensure mass conservation
    double sum = 0.0;
    for (std::size_t j = 0; j < n_species; ++j) {
      sum += mass_fractions[j];
    }
    if (sum > 1e-15) {
      for (std::size_t j = 0; j < n_species; ++j) {
        mass_fractions[j] /= sum;
      }
    }

    // Update species matrix
    for (std::size_t j = 0; j < n_species; ++j) {
      species_matrix(j, i) = mass_fractions[j];
    }
  }
}

namespace detail {

auto build_species_coefficients(const core::Matrix<double>& c_previous, const coefficients::CoefficientInputs& inputs,
                                const coefficients::CoefficientSet& coeffs, const conditions::BoundaryConditions& bc,
                                const coefficients::XiDerivatives& xi_der,
                                const thermophysics::MixtureInterface& mixture, const io::SimulationConfig& sim_config,
                                std::span<const double> F_field, std::span<const double> V_field, int station,
                                PhysicalQuantity auto d_eta) -> std::expected<SpeciesCoefficients, EquationError> {

  const auto n_eta = c_previous.cols();
  const auto n_species = c_previous.rows();
  const double d_eta_sq = d_eta * d_eta;
  const double xi = inputs.xi;
  const double lambda0 = xi_der.lambda0();
  const auto c_derivatives = xi_der.c_derivative();

  // Compute geometry factors
  const auto factors_result = compute_geometry_factors(station, xi, bc, sim_config);
  if (!factors_result) {
    return std::unexpected(factors_result.error());
  }
  const auto factors = factors_result.value();

  // Parameters for fake fluxes
  constexpr double Le = 1.2;  // Lewis number
  constexpr double Pr = 0.72; // Prandtl number

  // Fix concentration derivatives for mass conservation
  auto dc_deta_fixed = fix_concentration_derivatives(c_previous, inputs.dc_deta);

  // Compute fake fluxes and their derivatives
  auto fake_fluxes_result = compute_fake_fluxes(dc_deta_fixed, coeffs, d_eta, Le, Pr);
  if (!fake_fluxes_result) {
    return std::unexpected(fake_fluxes_result.error());
  }

  auto J_fake = std::move(fake_fluxes_result.value().first);
  auto dJ_fake_deta = std::move(fake_fluxes_result.value().second);

  SpeciesCoefficients species_coeffs;
  species_coeffs.a = core::Matrix<double>(n_eta, n_species);
  species_coeffs.b = core::Matrix<double>(n_eta, n_species);
  species_coeffs.c = core::Matrix<double>(n_eta, n_species);
  species_coeffs.d = core::Matrix<double>(n_eta, n_species);

  for (std::size_t i = 0; i < n_eta; ++i) {
    for (std::size_t j = 0; j < n_species; ++j) {

      // ----- Coefficient a[i] -----
      // a[i][j] = -l0[i]*Le/Pr / d_eta²
      species_coeffs.a(i, j) = -coeffs.transport.l0[i] * Le / Pr / d_eta_sq;

      // ----- Coefficient b[i] -----
      species_coeffs.b(i, j) = (V_field[i] - Le / Pr * coeffs.transport.dl0_deta[i]) / d_eta;

      // ----- Coefficient [i] -----
      species_coeffs.c(i, j) = 2.0 * xi * F_field[i] * lambda0;

      // ----- Coefficient d[i] -----
      // d[i][j] = -dJ_fake_deta[j][i] - dJ_deta[j][i]*J_fact -
      // 2*xi*c_der[j][i]*F[i]
      const double d_term = -dJ_fake_deta(j, i) - coeffs.diffusion.dJ_deta(j, i) * factors.J_fact -
                            2.0 * xi * c_derivatives(j, i) * F_field[i];

      const double wi_term = coeffs.chemical.wi(i, j) * factors.W_fact / coeffs.thermodynamic.rho[i];

      species_coeffs.d(i, j) = d_term + wi_term;
    }
  }

  return species_coeffs;
}

auto build_species_boundary_conditions(const core::Matrix<double>& c_wall, const coefficients::CoefficientSet& coeffs,
                                       const conditions::BoundaryConditions& bc,
                                       const thermophysics::MixtureInterface& mixture,
                                       const io::SimulationConfig& sim_config, PhysicalQuantity auto d_eta)
    -> std::expected<SpeciesBoundaryConditions, EquationError> {

  const auto n_species = mixture.n_species();

  SpeciesBoundaryConditions boundary_conds;
  boundary_conds.f_bc.resize(n_species);
  boundary_conds.g_bc.resize(n_species);
  boundary_conds.h_bc.resize(n_species);

  // Check for catalytic wall configuration
  if (sim_config.catalytic_wall) {
    return detail::build_catalytic_boundary_conditions(c_wall, coeffs, bc, mixture, sim_config, d_eta);
  } else {
    // Equilibrium wall boundary conditions (existing behavior)
    for (std::size_t i = 0; i < n_species; ++i) {
      boundary_conds.f_bc[i] = 0.0;
      boundary_conds.g_bc[i] = 1.0;
    }

    auto eq_wall_result = detail::compute_equilibrium_wall(bc, mixture);
    if (!eq_wall_result) {
      return std::unexpected(eq_wall_result.error());
    }

    boundary_conds.h_bc = eq_wall_result.value();
    return boundary_conds;
  }
}

auto compute_fake_fluxes(const core::Matrix<double>& dc_deta_fixed, const coefficients::CoefficientSet& coeffs,
                         PhysicalQuantity auto d_eta, double Le, double Pr)
    -> std::expected<std::pair<core::Matrix<double>, core::Matrix<double>>, EquationError> {

  const auto n_species = dc_deta_fixed.rows();
  const auto n_eta = dc_deta_fixed.cols();

  core::Matrix<double> J_fake(n_species, n_eta);
  core::Matrix<double> dJ_fake_deta(n_species, n_eta);

  // Compute fake fluxes: J_fake[j][i] = Le/Pr * l0[i] * dc_deta_fix[j][i]
  for (std::size_t i = 0; i < n_eta; ++i) {
    for (std::size_t j = 0; j < n_species; ++j) {
      J_fake(j, i) = Le / Pr * coeffs.transport.l0[i] * dc_deta_fixed(j, i);
    }
  }

  // Compute derivatives of fake fluxes
  for (std::size_t j = 0; j < n_species; ++j) {
    std::vector<double> J_row_data(n_eta);
    for (std::size_t i = 0; i < n_eta; ++i) {
      J_row_data[i] = J_fake(j, i);
    }

    auto dJ_result = coefficients::derivatives::compute_eta_derivative(std::span(J_row_data.data(), n_eta), d_eta);

    if (!dJ_result) {
      return std::unexpected(EquationError("Failed to compute diffusion flux derivative"));
    }
    auto dJ = dJ_result.value();

    for (std::size_t i = 0; i < n_eta; ++i) {
      dJ_fake_deta(j, i) = dJ[i];
    }
  }

  return std::make_pair(std::move(J_fake), std::move(dJ_fake_deta));
}

auto fix_concentration_derivatives(const core::Matrix<double>& c_matrix, const core::Matrix<double>& dc_deta)
    -> core::Matrix<double> {

  const auto n_species = c_matrix.rows();
  const auto n_eta = c_matrix.cols();

  core::Matrix<double> dc_deta_fixed(n_species, n_eta);

  for (std::size_t i = 0; i < n_eta; ++i) {
    // Compute sums
    double sum_c = 0.0;
    double sum_dc = 0.0;

    for (std::size_t j = 0; j < n_species; ++j) {
      sum_c += c_matrix(j, i);
      sum_dc += dc_deta(j, i);
    }

    // Apply mass conservation correction
    for (std::size_t j = 0; j < n_species; ++j) {
      dc_deta_fixed(j, i) = dc_deta(j, i) - c_matrix(j, i) / sum_c * sum_dc;
    }
  }

  return dc_deta_fixed;
}

auto compute_equilibrium_wall(const conditions::BoundaryConditions& bc, const thermophysics::MixtureInterface& mixture)
    -> std::expected<std::vector<double>, EquationError> {

  auto eq_result = mixture.equilibrium_composition(bc.Tw(), bc.P_e());
  if (!eq_result) {
    return std::unexpected(EquationError("Failed to compute equilibrium wall composition: {}",
                                         std::source_location::current(), eq_result.error().message()));
  }

  return eq_result.value();
}

// Transport number approximations for catalytic boundary conditions
constexpr double Le = 1.2;  // Lewis number approximation
constexpr double Pr = 0.72; // Prandtl number approximation

auto build_catalytic_boundary_conditions(const core::Matrix<double>& c_wall, const coefficients::CoefficientSet& coeffs,
                                         const conditions::BoundaryConditions& bc,
                                         const thermophysics::MixtureInterface& mixture,
                                         const io::SimulationConfig& sim_config, PhysicalQuantity auto d_eta)
    -> std::expected<SpeciesBoundaryConditions, EquationError> {

  const auto n_species = mixture.n_species();
  const std::size_t start_idx = mixture.has_electrons() ? 1 : 0;

  // Compute partial densities at wall: ρᵢ = cᵢ × ρₜₒₜₐₗ
  std::vector<double> rho_i_wall(n_species);
  const double rho_wall = coeffs.thermodynamic.rho[0]; // Wall density

  for (std::size_t i = 0; i < n_species; ++i) {
    rho_i_wall[i] = c_wall(i, 0) * rho_wall;
  }

  // Compute catalytic fluxes via Mutation++
  auto cat_flux_result = mixture.surface_reaction_rates(rho_i_wall, bc.Tw());
  if (!cat_flux_result) {
    return std::unexpected(EquationError("Failed to compute catalytic flux: {}", std::source_location::current(),
                                         cat_flux_result.error().message()));
  }
  auto cat_flux = cat_flux_result.value();

  // Compute geometry factor for boundary condition scaling
  auto factors_result = compute_geometry_factors(bc.station, bc.xi, bc, sim_config);
  if (!factors_result) {
    return std::unexpected(factors_result.error());
  }
  const double bc_fact = factors_result.value().bc_fact;

  // Build boundary conditions: J_diff|wall = J_cat
  SpeciesBoundaryConditions boundary_conds;
  boundary_conds.f_bc.resize(n_species);
  boundary_conds.g_bc.resize(n_species);
  boundary_conds.h_bc.resize(n_species);

  // Heavy species: catalytic boundary conditions
  for (std::size_t i = start_idx; i < n_species; ++i) {
    // Robin boundary condition: f_bc * dc/dη + g_bc * c = h_bc
    boundary_conds.f_bc[i] = -coeffs.transport.l0[0] * Le / Pr / d_eta;
    boundary_conds.g_bc[i] = 0.0; // No concentration term

    // J_reel
    const double J_reel = coeffs.diffusion.J(i, 0);

    // Fake flux on the wall to compensate f_bc * dc/dη
    // J_fake = Le/Pr * l0 * dc_deta on the wall
    const double dc_deta_wall = (c_wall.cols() > 1) ? (c_wall(i, 1) - c_wall(i, 0)) / d_eta : 0.0;
    const double J_fake_wall = Le / Pr * coeffs.transport.l0[0] * dc_deta_wall;

    // h_bc = J_cat - J_reel + J_fake_wall to ensure the flux balance
    boundary_conds.h_bc[i] = (+cat_flux[i] - J_reel) * 1 + J_fake_wall;
  }

  // Electrons: dummy values (not used in system, determined by charge neutrality)
  if (mixture.has_electrons()) {
    boundary_conds.f_bc[0] = 0.0;
    boundary_conds.g_bc[0] = 0.0;
    boundary_conds.h_bc[0] = 0.0;
  }

  return boundary_conds;
}

} // namespace detail

// Explicit instantiations for common use cases
template auto
solve_species<double>(const core::Matrix<double>& c_previous, const coefficients::CoefficientInputs& inputs,
                      const coefficients::CoefficientSet& coeffs, const conditions::BoundaryConditions& bc,
                      const coefficients::XiDerivatives& xi_der, const thermophysics::MixtureInterface& mixture,
                      const io::SimulationConfig& sim_config, std::span<const double> F_field,
                      std::span<const double> V_field, int station, double d_eta)
    -> std::expected<core::Matrix<double>, EquationError>;

} // namespace blast::boundary_layer::equations