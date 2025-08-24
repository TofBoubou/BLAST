#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include "blast/boundary_layer/coefficients/diffusion.hpp"
#include "blast/boundary_layer/grid/coordinate_transform.hpp"
#include "blast/core/constants.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <ranges>

namespace blast::boundary_layer::coefficients {

using namespace blast::boundary_layer;

void CoefficientCalculator::ensure_workspace_size(std::size_t n_species, std::size_t n_eta) const {
  if (workspace_species_.size() < n_species) {
    workspace_species_.resize(n_species);
  }
  if (workspace_species2_.size() < n_species) {
    workspace_species2_.resize(n_species);
  }
  if (workspace_derivatives_.size() < n_eta) {
    workspace_derivatives_.resize(n_eta);
  }
  if (workspace_integrand_.size() < n_eta) {
    workspace_integrand_.resize(n_eta);
  }
}

auto CoefficientCalculator::calculate(const CoefficientInputs& inputs, const conditions::BoundaryConditions& bc,
                                      const XiDerivatives& xi_der) const
    -> std::expected<CoefficientSet, CoefficientError> {

  // Validate inputs
  if (inputs.T.empty() || inputs.c.rows() == 0 || inputs.c.cols() == 0) {
    return std::unexpected(CoefficientError("Invalid input dimensions"));
  }

  if (inputs.T.size() != static_cast<std::size_t>(inputs.c.cols())) {
    return std::unexpected(CoefficientError("Temperature and concentration dimensions mismatch"));
  }

  CoefficientSet coeffs;

  // 1. Calculate thermodynamic coefficients first (needed by others)
  auto thermo_result = calculate_thermodynamic_coefficients(inputs, bc);
  if (!thermo_result) {
    return std::unexpected(thermo_result.error());
  }
  coeffs.thermodynamic = std::move(thermo_result.value());

  // 3. Transport coefficients
  auto transport_result = calculate_transport_coefficients(inputs, coeffs.thermodynamic, bc);
  if (!transport_result) {
    return std::unexpected(transport_result.error());
  }
  coeffs.transport = std::move(transport_result.value());

  // 4. Diffusion coefficients
  auto diffusion_result = calculate_diffusion_coefficients(inputs, bc, xi_der);
  if (!diffusion_result) {
    return std::unexpected(diffusion_result.error());
  }
  coeffs.diffusion = std::move(diffusion_result.value());

  // 5. Chemical coefficients
  if (sim_config_.chemical_mode == io::SimulationConfig::ChemicalMode::NonEquilibrium) {
    auto chemical_result = calculate_chemical_coefficients(inputs, coeffs.thermodynamic, bc);
    if (!chemical_result) {
      return std::unexpected(chemical_result.error());
    }
    coeffs.chemical = std::move(chemical_result.value());
  } else {
    // Initialize with zeros
    const auto n_eta = inputs.T.size();
    const auto n_species = inputs.c.rows();
    coeffs.chemical.wi = core::Matrix<double>(n_eta, n_species);
    coeffs.chemical.wi.setZero();
    coeffs.chemical.d_wi_dc.resize(n_eta);
    for (auto& mat : coeffs.chemical.d_wi_dc) {
      mat = core::Matrix<double>(n_species, n_species);
      mat.setZero();
    }
  }

  // 6. Thermal diffusion (if enabled)
  if (sim_config_.consider_thermal_diffusion) {
    auto tdr_result = calculate_thermal_diffusion(inputs, coeffs.thermodynamic, bc);
    if (!tdr_result) {
      return std::unexpected(tdr_result.error());
    }
    coeffs.thermal_diffusion = std::move(tdr_result.value());
  }

  // 7. Wall properties
  auto wall_result = calculate_wall_properties(inputs, bc, coeffs.transport, coeffs.thermodynamic);
  if (!wall_result) {
    return std::unexpected(wall_result.error());
  }
  coeffs.wall = std::move(wall_result.value());

  // 8. Species enthalpies
  auto enthalpy_result = calculate_species_enthalpies(inputs);
  if (!enthalpy_result) {
    return std::unexpected(enthalpy_result.error());
  }
  auto [h_sp, dh_sp_deta] = std::move(enthalpy_result.value());
  coeffs.h_species = std::move(h_sp);
  coeffs.dh_species_deta = std::move(dh_sp_deta);

  // 9. Compute diffusion fluxes
  auto flux_result =
      diffusion::compute_stefan_maxwell_fluxes(inputs, coeffs, bc, xi_der, sim_config_, mixture_, d_eta_);
  if (!flux_result) {
    return std::unexpected(flux_result.error());
  }

  return coeffs;
}

auto CoefficientCalculator::calculate_thermodynamic_coefficients(const CoefficientInputs& inputs,
                                                                 const conditions::BoundaryConditions& bc) const
    -> std::expected<ThermodynamicCoefficients, CoefficientError> {

  const auto n_eta = inputs.T.size();
  const auto n_species = inputs.c.rows();
  
  // Ensure workspace is sized correctly
  ensure_workspace_size(n_species, n_eta);

  ThermodynamicCoefficients thermo;
  thermo.rho.reserve(n_eta);
  thermo.MW.reserve(n_eta);

  // Calculate molecular weights and densities
  for (std::size_t i = 0; i < n_eta; ++i) {
    // Validate temperature
    if (!std::isfinite(inputs.T[i]) || inputs.T[i] <= 0.0) {
      return std::unexpected(CoefficientError(std::format("Invalid temperature at eta={}: T={}", i, inputs.T[i])));
    }

    // Extract species concentrations using pre-allocated workspace
    auto& c_local = workspace_species_;  // Reuse workspace
    double sum_c = 0.0;
    for (std::size_t j = 0; j < n_species; ++j) {
      c_local[j] = inputs.c(j, i);
      if (!std::isfinite(c_local[j]) || c_local[j] < -1e-6) {
        return std::unexpected(
            CoefficientError(std::format("Invalid concentration at eta={} species={}: c={}", i, j, c_local[j])));
      }
      sum_c += c_local[j];
    }

    if (sum_c <= 0.0) {
      return std::unexpected(
          CoefficientError(std::format("Zero or negative total concentration at eta={}: sum_c={}", i, sum_c)));
    }

    // Get molecular weight
    auto mw_result = mixture_.mixture_molecular_weight(c_local);
    if (!mw_result) {
      return std::unexpected(
          CoefficientError(std::format("Failed to compute MW at eta={}: {}", i, mw_result.error().message())));
    }
    thermo.MW.push_back(mw_result.value());

    thermo.rho.push_back(bc.P_e() * thermo.MW[i] / (inputs.T[i] * constants::physical::universal_gas_constant));
  }

  // Compute molecular weight derivative
  thermo.d_MW_deta.resize(n_eta);
  for (std::size_t i = 0; i < n_eta; ++i) {
    // Reuse workspace for derivative calculation
    auto& dc_deta_local = workspace_species2_;
    for (std::size_t j = 0; j < n_species; ++j) {
      dc_deta_local[j] = inputs.dc_deta(j, i);
      if (!std::isfinite(dc_deta_local[j])) {
        return std::unexpected(CoefficientError(std::format("Invalid dc_deta at eta={} species={}", i, j)));
      }
    }

    // d(MW)/dη = -MW² · Σ(dc_j/dη / Mw_j)
    double sum = 0.0;
    for (std::size_t j = 0; j < n_species; ++j) {
      double mw_j = mixture_.species_molecular_weight(j);
      if (mw_j <= 0.0) {
        return std::unexpected(CoefficientError(std::format("Invalid species molecular weight for species {}", j)));
      }
      sum += dc_deta_local[j] / mw_j;
    }
    thermo.d_MW_deta[i] = -thermo.MW[i] * thermo.MW[i] * sum;
  }

  // Compute density derivative
  auto d_rho_deta_result = derivatives::compute_eta_derivative(thermo.rho, d_eta_);
  if (!d_rho_deta_result) {
    return std::unexpected(CoefficientError("Failed to compute density derivative"));
  }
  thermo.d_rho_deta = d_rho_deta_result.value();

  // Compute wall enthalpy using workspace
  auto& c_wall = workspace_species_;
  for (std::size_t j = 0; j < n_species; ++j) {
    c_wall[j] = inputs.c(j, 0);
  }

  auto h_wall_result = mixture_.mixture_enthalpy(c_wall, bc.Tw(), bc.P_e());
  if (!h_wall_result) {
    return std::unexpected(CoefficientError("Failed to compute wall enthalpy"));
  }
  thermo.h_wall = h_wall_result.value();

  return thermo;
}

auto CoefficientCalculator::calculate_transport_coefficients(const CoefficientInputs& inputs,
                                                             const ThermodynamicCoefficients& thermo,
                                                             const conditions::BoundaryConditions& bc) const
    -> std::expected<TransportCoefficients, CoefficientError> {

  const auto n_eta = inputs.T.size();
  const auto n_species = inputs.c.rows();

  TransportCoefficients transport;
  transport.l0.reserve(n_eta);
  transport.l3.reserve(n_eta);

  for (std::size_t i = 0; i < n_eta; ++i) {
    // Extract local composition using workspaces
    auto& c_local = workspace_species_;
    auto& c_local_e = workspace_species2_;
    for (std::size_t j = 0; j < n_species; ++j) {
      c_local[j] = inputs.c(j, i);
    }

    for (std::size_t j = 0; j < n_species; ++j) {
      c_local_e[j] = inputs.c(j, n_eta - 1);
    }

    const double P_edge = bc.P_e();

    const double rho_e_calculated =
        bc.P_e() * thermo.MW[n_eta - 1] / (inputs.T[n_eta - 1] * constants::physical::universal_gas_constant);

    // Get transport properties
    auto mu_result = mixture_.viscosity(c_local, inputs.T[i], P_edge);
    auto mu_result_e = mixture_.viscosity(c_local_e, inputs.T[n_eta - 1], P_edge);
    auto cp_result = mixture_.frozen_cp(c_local, inputs.T[i], P_edge);
    auto k_result = mixture_.frozen_thermal_conductivity(c_local, inputs.T[i], P_edge);

    if (!mu_result || !cp_result || !k_result) {
      return std::unexpected(CoefficientError(std::format("Failed to compute transport properties at eta={}", i)));
    }

    const double mu = mu_result.value();
    const double mu_e = mu_result_e.value();
    const double Cp = cp_result.value();
    const double k_fr = k_result.value();

    if (mu <= 0.0 || Cp <= 0.0 || k_fr <= 0.0) {
      return std::unexpected(
          CoefficientError(std::format("Invalid transport properties at eta={}: mu={}, Cp={}, k={}", i, mu, Cp, k_fr)));
    }

    const double Pr = mu * Cp / k_fr;

    if (bc.rho_e() <= 0.0 || bc.mu_e() <= 0.0) {
      return std::unexpected(CoefficientError("Invalid edge conditions: rho_e or mu_e <= 0"));
    }

    double l0_value = thermo.rho[i] * mu / (rho_e_calculated * mu_e);
    transport.l0.push_back(l0_value);
    transport.l3.push_back(transport.l0[i] / Pr);
  }

  // Compute derivatives
  auto dl0_deta_result = derivatives::compute_eta_derivative(transport.l0, d_eta_);
  if (!dl0_deta_result) {
    return std::unexpected(CoefficientError("Failed to compute dl0/deta"));
  }
  transport.dl0_deta = dl0_deta_result.value();

  auto dl3_deta_result = derivatives::compute_eta_derivative(transport.l3, d_eta_);
  if (!dl3_deta_result) {
    return std::unexpected(CoefficientError("Failed to compute dl3/deta"));
  }
  transport.dl3_deta = dl3_deta_result.value();

  // Calculate finite thickness coefficients if enabled
  auto finite_thickness_result = calculate_finite_thickness_coefficients(inputs, thermo, bc);
  if (!finite_thickness_result) {
    return std::unexpected(finite_thickness_result.error());
  }
  auto [K_bl, coeff_finite_thickness] = finite_thickness_result.value();
  transport.K_bl = K_bl;
  transport.coeff_finite_thickness = coeff_finite_thickness;

  return transport;
}

auto CoefficientCalculator::calculate_diffusion_coefficients(const CoefficientInputs& inputs,
                                                             const conditions::BoundaryConditions& bc,
                                                             const XiDerivatives& xi_der) const
    -> std::expected<DiffusionCoefficients, CoefficientError> {

  const auto n_eta = inputs.T.size();
  const auto n_species = inputs.c.rows();

  DiffusionCoefficients diff;

  if (n_species == 1) {
    diff.Dij_bin = core::Matrix<double>(n_eta, 1);
    diff.Dij_bin.setZero();
    diff.J = core::Matrix<double>(1, n_eta);
    diff.J.setZero();
    diff.dJ_deta = core::Matrix<double>(1, n_eta);
    diff.dJ_deta.setZero();
    return diff;
  }

  diff.Dij_bin = core::Matrix<double>(n_eta * n_species, n_species);

  // Use edge pressure consistently
  const double P_edge = bc.P_e();
  if (P_edge <= 0.0) {
    return std::unexpected(CoefficientError("Invalid edge pressure: P_e <= 0"));
  }

  // Calculate binary diffusion coefficients
  for (std::size_t i = 0; i < n_eta; ++i) {
    // Get binary diffusion coefficients using edge pressure
    auto dij_result = mixture_.binary_diffusion_coefficients(inputs.T[i], P_edge);
    if (!dij_result) {
      return std::unexpected(CoefficientError(std::format("Failed to compute Dij at eta={}", i)));
    }

    // Store in matrix (row i*n_species to (i+1)*n_species)
    const auto& dij_local = dij_result.value();

    for (std::size_t j = 0; j < n_species; ++j) {
      for (std::size_t k = 0; k < n_species; ++k) {
        double dij_val = dij_local(j, k);
        if (!std::isfinite(dij_val) || dij_val < 0.0) {
          return std::unexpected(CoefficientError(
              std::format("Invalid diffusion coefficient at eta={}, species {}-{}: Dij={}", i, j, k, dij_val)));
        }
        diff.Dij_bin(i * n_species + j, k) = dij_val;
      }
    }
  }

  return diff;
}

auto CoefficientCalculator::calculate_chemical_coefficients(const CoefficientInputs& inputs,
                                                            const ThermodynamicCoefficients& thermo,
                                                            const conditions::BoundaryConditions& bc) const
    -> std::expected<ChemicalCoefficients, CoefficientError> {

  const auto n_eta = inputs.T.size();
  const auto n_species = inputs.c.rows();

  ChemicalCoefficients chem;
  chem.wi = core::Matrix<double>(n_eta, n_species);
  chem.d_wi_dc.reserve(n_eta);

  // Calculate production rates for each eta station
  for (std::size_t i = 0; i < n_eta; ++i) {
    auto result = calculate_station_coefficients(i, inputs, thermo, bc);
    if (!result) {
      return std::unexpected(result.error());
    }

    // Store results
    const auto& [wi_local, dwi_dc] = result.value();
    for (std::size_t j = 0; j < n_species; ++j) {
      if (!std::isfinite(wi_local[j])) {
        return std::unexpected(CoefficientError(std::format("Invalid production rate at eta={}, species={}", i, j)));
      }
      chem.wi(i, j) = wi_local[j];
    }
    chem.d_wi_dc.push_back(std::move(dwi_dc));
  }

  return chem;
}

auto CoefficientCalculator::calculate_station_coefficients(std::size_t station_index, const CoefficientInputs& inputs,
                                                           const ThermodynamicCoefficients& thermo,
                                                           const conditions::BoundaryConditions& bc) const
    -> std::expected<std::pair<std::vector<double>, core::Matrix<double>>, CoefficientError> {

  const auto n_species = inputs.c.rows();

  if (station_index >= inputs.T.size()) {
    return std::unexpected(CoefficientError("Station index out of bounds"));
  }

  // Extract and normalize composition
  auto composition_result = extract_local_composition(station_index, inputs.c, thermo.rho[station_index]);
  if (!composition_result) {
    return std::unexpected(composition_result.error());
  }
  auto [c_local, rho_species] = composition_result.value();

  // Get production rates
  auto wi_result = mixture_.production_rates(rho_species, inputs.T[station_index]);
  if (!wi_result) {
    return std::unexpected(CoefficientError(
        std::format("Failed to compute production rates at eta={}: {}", station_index, wi_result.error().message())));
  }
  const auto& wi_local = wi_result.value();

  // Get Jacobian d(wi)/d(rho_j)
  auto jac_result = mixture_.production_rate_jacobian(rho_species, inputs.T[station_index]);
  if (!jac_result) {
    return std::unexpected(CoefficientError(
        std::format("Failed to compute Jacobian at eta={}: {}", station_index, jac_result.error().message())));
  }

  // Transform Jacobian from d(wi)/d(rho_j) to d(wi)/d(c_j)
  auto dwi_dc_result =
      transform_jacobian(station_index, inputs, thermo, bc, c_local, rho_species, wi_local, jac_result.value());

  if (!dwi_dc_result) {
    return std::unexpected(dwi_dc_result.error());
  }

  return std::make_pair(wi_local, dwi_dc_result.value());
}

auto CoefficientCalculator::extract_local_composition(std::size_t station_index, const core::Matrix<double>& c,
                                                      double rho_total) const
    -> std::expected<std::pair<std::vector<double>, std::vector<double>>, CoefficientError> {

  const auto n_species = c.rows();
  // Use static thread_local storage to avoid repeated allocations
  thread_local std::vector<double> c_local_tl;
  thread_local std::vector<double> rho_species_tl;
  
  if (c_local_tl.size() < n_species) {
    c_local_tl.resize(n_species);
    rho_species_tl.resize(n_species);
  }
  
  auto& c_local = c_local_tl;
  auto& rho_species = rho_species_tl;

  double sum_c = 0.0;
  for (std::size_t j = 0; j < n_species; ++j) {
    c_local[j] = c(j, station_index);
    sum_c += c_local[j];
  }

  // Normalize and compute partial densities
  if (sum_c > 0.0) {
    for (std::size_t j = 0; j < n_species; ++j) {
      rho_species[j] = rho_total * c_local[j] / sum_c;
    }
  } else {
    // Error: zero composition is invalid
    return std::unexpected(CoefficientError("Zero total composition encountered in extract_local_composition"));
  }

  // Return copies to avoid issues with thread_local storage
  return std::make_pair(
    std::vector<double>(c_local.begin(), c_local.begin() + n_species),
    std::vector<double>(rho_species.begin(), rho_species.begin() + n_species)
  );
}

auto CoefficientCalculator::transform_jacobian(
    std::size_t station_index, const CoefficientInputs& inputs, const ThermodynamicCoefficients& thermo,
    const conditions::BoundaryConditions& bc, const std::vector<double>& c_local,
    const std::vector<double>& rho_species, const std::vector<double>& wi_local,
    const core::Matrix<double>& dwi_drho) const -> std::expected<core::Matrix<double>, CoefficientError> {

  const auto n_species = c_local.size();
  const double T = inputs.T[station_index];
  const double P_edge = bc.P_e();

  // Get species enthalpies and Cp
  auto h_sp_result = mixture_.species_enthalpies(T);
  auto cp_result = mixture_.frozen_cp(c_local, T, P_edge);

  if (!h_sp_result || !cp_result) {
    return std::unexpected(CoefficientError("Failed to get enthalpies or Cp"));
  }

  const auto& h_species = h_sp_result.value();
  const double Cp = cp_result.value();

  if (Cp <= 0.0) {
    return std::unexpected(CoefficientError("Invalid specific heat: Cp <= 0"));
  }

  // Compute dT/dc_j = -h_j/Cp
  std::vector<double> dT_dc(n_species);
  for (std::size_t j = 0; j < n_species; ++j) {
    dT_dc[j] = -h_species[j] / Cp;
  }

  // Compute d(rho_k)/d(c_j)
  auto drho_dc = compute_density_derivatives(rho_species, thermo.rho[station_index], thermo.MW[station_index]);

  // Compute d(wi)/dT
  auto dwi_dT_result = compute_temperature_derivatives(rho_species, T, wi_local);
  if (!dwi_dT_result) {
    return std::unexpected(CoefficientError("Failed to compute temperature derivatives"));
  }
  const auto& dwi_dT = dwi_dT_result.value();

  // Assemble complete Jacobian
  core::Matrix<double> dwi_dc(n_species, n_species);
  for (std::size_t i = 0; i < n_species; ++i) {
    for (std::size_t j = 0; j < n_species; ++j) {
      double sum = 0.0;
      for (std::size_t k = 0; k < n_species; ++k) {
        sum += dwi_drho(i, k) * drho_dc(k, j);
      }
      dwi_dc(i, j) = sum + dwi_dT[i] * dT_dc[j];

      if (!std::isfinite(dwi_dc(i, j))) {
        return std::unexpected(CoefficientError(std::format("Invalid Jacobian element at ({}, {})", i, j)));
      }
    }
  }

  return dwi_dc;
}

auto CoefficientCalculator::compute_density_derivatives(const std::vector<double>& rho_species, double rho_total,
                                                        double MW_mixture) const -> core::Matrix<double> {

  const auto n_species = rho_species.size();
  core::Matrix<double> drho_dc(n_species, n_species);

  for (std::size_t k = 0; k < n_species; ++k) {
    for (std::size_t j = 0; j < n_species; ++j) {
      drho_dc(k, j) = (k == j) ? rho_total : 0.0;
      double mw_j = mixture_.species_molecular_weight(j);
      if (mw_j > 0.0) {
        drho_dc(k, j) -= rho_species[k] * MW_mixture / mw_j;
      }
    }
  }

  return drho_dc;
}

auto CoefficientCalculator::compute_temperature_derivatives(const std::vector<double>& rho_species, double T,
                                                            const std::vector<double>& wi_base) const
    -> std::expected<std::vector<double>, CoefficientError> {

  constexpr double eps = 1e-6;
  const double dT = eps * T;

  if (dT <= 0.0) {
    return std::unexpected(CoefficientError("Invalid temperature for finite difference"));
  }

  auto wi_pert_result = mixture_.production_rates(rho_species, T + dT);
  if (!wi_pert_result) {
    return std::unexpected(CoefficientError("Failed to compute perturbed rates"));
  }

  const auto& wi_pert = wi_pert_result.value();
  
  // Use thread_local storage to avoid allocation
  thread_local std::vector<double> dwi_dT_tl;
  if (dwi_dT_tl.size() < wi_base.size()) {
    dwi_dT_tl.resize(wi_base.size());
  }
  auto& dwi_dT = dwi_dT_tl;

  for (std::size_t j = 0; j < wi_base.size(); ++j) {
    dwi_dT[j] = (wi_pert[j] - wi_base[j]) / dT;
    if (!std::isfinite(dwi_dT[j])) {
      return std::unexpected(CoefficientError(std::format("Invalid temperature derivative for species {}", j)));
    }
  }

  // Return copy to avoid reference to thread_local storage
  return std::vector<double>(dwi_dT.begin(), dwi_dT.begin() + wi_base.size());
}

auto CoefficientCalculator::calculate_thermal_diffusion(const CoefficientInputs& inputs,
                                                        const ThermodynamicCoefficients& thermo,
                                                        const conditions::BoundaryConditions& bc) const
    -> std::expected<ThermalDiffusionCoefficients, CoefficientError> {

  const auto n_eta = inputs.T.size();
  const auto n_species = inputs.c.rows();

  ThermalDiffusionCoefficients tdr;
  tdr.tdr = core::Matrix<double>(n_species, n_eta);

  // Use edge pressure consistently
  const double P_edge = bc.P_e();

  // Calculate thermal diffusion ratios
  for (std::size_t i = 0; i < n_eta; ++i) {
    auto& c_local = workspace_species_;
    for (std::size_t j = 0; j < n_species; ++j) {
      c_local[j] = inputs.c(j, i);
    }

    // Get thermal diffusion ratios using edge pressure
    auto tdr_result = mixture_.thermal_diffusion_ratios(c_local, inputs.T[i], P_edge);
    if (!tdr_result) {
      return std::unexpected(
          CoefficientError(std::format("Failed to compute TDR at eta={}: {}", i, tdr_result.error().message())));
    }

    const auto& tdr_local = tdr_result.value();
    for (std::size_t j = 0; j < n_species; ++j) {
      if (!std::isfinite(tdr_local[j])) {
        return std::unexpected(CoefficientError(std::format("Invalid TDR at eta={}, species={}", i, j)));
      }
      tdr.tdr(j, i) = tdr_local[j];
    }
  }

  // Compute TDR term for diffusion fluxes
  if (sim_config_.consider_thermal_diffusion) {
    // Compute dT/deta
    auto dT_deta_result = derivatives::compute_eta_derivative(inputs.T, d_eta_);
    if (!dT_deta_result) {
      return std::unexpected(CoefficientError("Failed to compute dT/deta"));
    }
    auto dT_deta = dT_deta_result.value();

    // Compute TDR term
    tdr.tdr_term = core::Matrix<double>(n_eta, n_species);
    for (std::size_t i = 0; i < n_eta; ++i) {
      if (inputs.T[i] <= 0.0) {
        return std::unexpected(CoefficientError(std::format("Invalid temperature for TDR calculation at eta={}", i)));
      }
      for (std::size_t j = 0; j < n_species; ++j) {
        tdr.tdr_term(i, j) = tdr.tdr(j, i) / inputs.T[i] * dT_deta[i];
      }
    }
  }

  return tdr;
}

auto CoefficientCalculator::calculate_species_enthalpies(const CoefficientInputs& inputs) const
    -> std::expected<std::pair<core::Matrix<double>, core::Matrix<double>>, CoefficientError> {

  const auto n_eta = inputs.T.size();
  const auto n_species = mixture_.n_species();

  core::Matrix<double> h_species(n_species, n_eta);

  // Get species enthalpies at each temperature
  for (std::size_t i = 0; i < n_eta; ++i) {
    auto h_result = mixture_.species_enthalpies(inputs.T[i]);
    if (!h_result) {
      return std::unexpected(CoefficientError(
          std::format("Failed to get species enthalpies at eta={}: {}", i, h_result.error().message())));
    }

    const auto& h_local = h_result.value();
    for (std::size_t j = 0; j < n_species; ++j) {
      if (!std::isfinite(h_local[j])) {
        return std::unexpected(CoefficientError(std::format("Invalid species enthalpy at eta={}, species={}", i, j)));
      }
      h_species(j, i) = h_local[j];
    }
  }

  // Compute derivatives
  core::Matrix<double> dh_species_deta(n_species, n_eta);

  for (std::size_t j = 0; j < n_species; ++j) {
    auto& h_row = workspace_derivatives_;
    for (std::size_t i = 0; i < n_eta; ++i) {
      h_row[i] = h_species(j, i);
    }

    auto dh_result = derivatives::compute_eta_derivative(h_row, d_eta_);
    if (!dh_result) {
      return std::unexpected(CoefficientError("Failed to compute dh/deta for species"));
    }
    auto dh = dh_result.value();
    for (std::size_t i = 0; i < n_eta; ++i) {
      dh_species_deta(j, i) = dh[i];
    }
  }

  return std::make_pair(std::move(h_species), std::move(dh_species_deta));
}

auto CoefficientCalculator::calculate_wall_properties(const CoefficientInputs& inputs,
                                                      const conditions::BoundaryConditions& bc,
                                                      const TransportCoefficients& transport,
                                                      const ThermodynamicCoefficients& thermo) const
    -> std::expected<WallProperties, CoefficientError> {

  const auto n_species = inputs.c.rows();

  // Extract wall composition using workspace
  auto& c_wall = workspace_species_;
  for (std::size_t j = 0; j < n_species; ++j) {
    c_wall[j] = inputs.c(j, 0);
  }

  // Get wall transport properties using edge pressure
  auto k_result = mixture_.frozen_thermal_conductivity(c_wall, inputs.T[0], bc.P_e());
  if (!k_result) {
    return std::unexpected(CoefficientError("Failed to compute wall thermal conductivity"));
  }

  auto mu_result = mixture_.viscosity(c_wall, inputs.T[0], bc.P_e());
  if (!mu_result) {
    return std::unexpected(CoefficientError("Failed to compute wall viscosity"));
  }

  auto cp_result = mixture_.frozen_cp(c_wall, inputs.T[0], bc.P_e());
  if (!cp_result) {
    return std::unexpected(CoefficientError("Failed to compute wall specific heat"));
  }

  const double k_wall = k_result.value();
  const double mu_wall = mu_result.value();
  const double Cp_wall = cp_result.value();

  if (k_wall <= 0.0 || mu_wall <= 0.0 || Cp_wall <= 0.0) {
    return std::unexpected(
        CoefficientError(std::format("Invalid wall properties: k={}, mu={}, Cp={}", k_wall, mu_wall, Cp_wall)));
  }

  WallProperties wall;
  wall.k_wall = k_wall;
  wall.mu_wall = mu_wall;
  wall.Cp_wall = Cp_wall;
  wall.rho_wall = thermo.rho[0];
  wall.Pr_wall = mu_wall * Cp_wall / k_wall;

  return wall;
}

// Derivative computation implementations
namespace derivatives {

template <std::ranges::sized_range Range>
auto compute_eta_derivative(Range&& values, double d_eta) -> std::expected<std::vector<double>, CoefficientError> {
  const auto n = std::ranges::size(values);
  
  // Use thread_local storage to avoid repeated allocations
  thread_local std::vector<double> derivatives_tl;
  if (derivatives_tl.size() < n) {
    derivatives_tl.resize(n);
  }
  auto& derivatives = derivatives_tl;

  if (d_eta <= 0.0) {
    return std::unexpected(CoefficientError("Invalid grid spacing: d_eta must be positive"));
  }

  if (n < 2) {
    return std::unexpected(CoefficientError("Insufficient points for derivative calculation: need at least 2 points"));
  }

  if (n < 5) {
    // Simple finite differences for small arrays
    derivatives[0] = (values[1] - values[0]) / d_eta;
    for (std::size_t i = 1; i < n - 1; ++i) {
      derivatives[i] = (values[i + 1] - values[i - 1]) / (2.0 * d_eta);
    }
    derivatives[n - 1] = (values[n - 1] - values[n - 2]) / d_eta;
    return std::vector<double>(derivatives.begin(), derivatives.begin() + n);
  }

  const double dx12 = 12.0 * d_eta;

  // Forward 5-point stencil O(h^4) for first 2 points
  derivatives[0] =
      (-25.0 * values[0] + 48.0 * values[1] - 36.0 * values[2] + 16.0 * values[3] - 3.0 * values[4]) / dx12;

  derivatives[1] =
      (-25.0 * values[1] + 48.0 * values[2] - 36.0 * values[3] + 16.0 * values[4] - 3.0 * values[5]) / dx12;

  // Central 5-point stencil O(h^4) for interior points
  for (std::size_t i = 2; i < n - 2; ++i) {
    derivatives[i] = (values[i - 2] - 8.0 * values[i - 1] + 8.0 * values[i + 1] - values[i + 2]) / dx12;
  }

  // Backward 5-point stencil O(h^4) for last 2 points
  derivatives[n - 2] =
      (-1.0 * values[n - 5] + 6.0 * values[n - 4] - 18.0 * values[n - 3] + 10.0 * values[n - 2] + 3.0 * values[n - 1]) /
      dx12;

  derivatives[n - 1] = (3.0 * values[n - 5] - 16.0 * values[n - 4] + 36.0 * values[n - 3] - 48.0 * values[n - 2] +
                        25.0 * values[n - 1]) /
                       dx12;

  // Copy result to avoid returning reference to thread_local
  return std::vector<double>(derivatives.begin(), derivatives.begin() + n);
}

template <std::ranges::sized_range Range>
auto compute_eta_second_derivative(Range&& values, double d_eta)
    -> std::expected<std::vector<double>, CoefficientError> {
  const auto n = std::ranges::size(values);
  
  // Use thread_local storage to avoid repeated allocations
  thread_local std::vector<double> second_derivatives_tl;
  if (second_derivatives_tl.size() < n) {
    second_derivatives_tl.resize(n);
  }
  auto& second_derivatives = second_derivatives_tl;

  if (d_eta <= 0.0) {
    return std::unexpected(CoefficientError("Invalid grid spacing: d_eta must be positive"));
  }

  if (n < 3) {
    return std::unexpected(
        CoefficientError("Insufficient points for second derivative calculation: need at least 3 points"));
  }

  if (n < 5) {
    // Simple 3-point stencil for small arrays
    const double d_eta_sq = d_eta * d_eta;

    second_derivatives[0] = (values[0] - 2.0 * values[1] + values[2]) / d_eta_sq;
    for (std::size_t i = 1; i < n - 1; ++i) {
      second_derivatives[i] = (values[i - 1] - 2.0 * values[i] + values[i + 1]) / d_eta_sq;
    }
    second_derivatives[n - 1] = (values[n - 3] - 2.0 * values[n - 2] + values[n - 1]) / d_eta_sq;

    return std::vector<double>(second_derivatives.begin(), second_derivatives.begin() + n);
  }

  const double dx2_12 = 12.0 * d_eta * d_eta;

  // Forward 5-point stencil O(h^4) for first 2 points
  second_derivatives[0] =
      (35.0 * values[0] - 104.0 * values[1] + 114.0 * values[2] - 56.0 * values[3] + 11.0 * values[4]) / dx2_12;

  second_derivatives[1] =
      (35.0 * values[1] - 104.0 * values[2] + 114.0 * values[3] - 56.0 * values[4] + 11.0 * values[5]) / dx2_12;

  // Central 5-point stencil O(h^4) for interior points
  for (std::size_t i = 2; i < n - 2; ++i) {
    second_derivatives[i] =
        (-values[i - 2] + 16.0 * values[i - 1] - 30.0 * values[i] + 16.0 * values[i + 1] - values[i + 2]) / dx2_12;
  }

  // Backward 5-point stencils O(h^4) for last 2 points
  second_derivatives[n - 2] =
      (-1.0 * values[n - 5] + 4.0 * values[n - 4] + 6.0 * values[n - 3] - 20.0 * values[n - 2] + 11.0 * values[n - 1]) /
      dx2_12;

  second_derivatives[n - 1] = (11.0 * values[n - 5] - 56.0 * values[n - 4] + 114.0 * values[n - 3] -
                               104.0 * values[n - 2] + 35.0 * values[n - 1]) /
                              dx2_12;

  return second_derivatives;
}

template <typename Matrix>
auto compute_matrix_eta_second_derivative(const Matrix& values, double d_eta)
    -> std::expected<Matrix, CoefficientError> {
  const auto n_rows = values.rows();
  const auto n_cols = values.cols();
  Matrix result(n_rows, n_cols);

  if (d_eta <= 0.0) {
    return std::unexpected(CoefficientError("Invalid grid spacing: d_eta must be positive"));
  }

  // Use thread_local storage for row values
  thread_local std::vector<double> row_values_tl;
  if (row_values_tl.size() < n_cols) {
    row_values_tl.resize(n_cols);
  }
  
  for (std::size_t i = 0; i < n_rows; ++i) {
    auto& row_values = row_values_tl;
    for (std::size_t j = 0; j < n_cols; ++j) {
      row_values[j] = values(i, j);
    }
    auto row_derivatives_result = compute_eta_second_derivative(row_values, d_eta);
    if (!row_derivatives_result) {
      return std::unexpected(row_derivatives_result.error());
    }
    const auto& row_derivatives = row_derivatives_result.value();
    for (std::size_t j = 0; j < n_cols; ++j) {
      result(i, j) = row_derivatives[j];
    }
  }

  return result;
}

// Explicit instantiations
template auto compute_eta_derivative(std::span<const double>&&, double)
    -> std::expected<std::vector<double>, CoefficientError>;
template auto compute_eta_derivative(const std::vector<double>&, double)
    -> std::expected<std::vector<double>, CoefficientError>;
template auto compute_eta_derivative(std::vector<double>&&, double)
    -> std::expected<std::vector<double>, CoefficientError>;
template auto compute_eta_derivative(std::span<double>&&, double)
    -> std::expected<std::vector<double>, CoefficientError>;

} // namespace derivatives

auto CoefficientCalculator::calculate_finite_thickness_coefficients(const CoefficientInputs& inputs,
                                                                   const ThermodynamicCoefficients& thermo,
                                                                   const conditions::BoundaryConditions& bc) const
    -> std::expected<std::pair<double, double>, CoefficientError> {
  
  // If finite thickness is not enabled, return default values
  if (!sim_config_.finite_thickness) {
    return std::make_pair(1.0, 0.0);
  }

  // Get finite thickness parameters from configuration
  const auto& ft_params = outer_edge_config_.finite_thickness_params;
  const double v_e = ft_params.v_edge;
  const double d2_ue_dxdy = ft_params.d2_ue_dxdy;
  const double delta_bl = ft_params.delta_bl;

  // Edge properties
  const auto n_eta = inputs.T.size();
  const double rho_e = thermo.rho[n_eta - 1];
  const double d_ue_dx = bc.d_ue_dx();

  // Get edge viscosity using workspace
  const auto n_species = inputs.c.rows();
  ensure_workspace_size(n_species, n_eta);
  auto& c_edge = workspace_species_;
  for (std::size_t j = 0; j < n_species; ++j) {
    c_edge[j] = inputs.c(j, n_eta - 1);
  }
  
  auto mu_e_result = mixture_.viscosity(c_edge, inputs.T[n_eta - 1], bc.P_e());
  if (!mu_e_result) {
    return std::unexpected(CoefficientError("Failed to compute edge viscosity for finite thickness"));
  }
  const double mu_e = mu_e_result.value();

  // Compute finite thickness coefficient
  const double coeff_finite_thickness = v_e * d2_ue_dxdy / (d_ue_dx * d_ue_dx);

  auto& integrand = workspace_integrand_;
  for (std::size_t i = 0; i < n_eta; ++i) {
    integrand[i] = 1.0 / thermo.rho[i];
  }
  
  const double d_eta = num_config_.eta_max / static_cast<double>(num_config_.n_eta - 1);
  
  auto integrated_values = grid::coordinate_transform::simpson_integrate(integrand, d_eta, 0.0);
  
  // Get the final integrated value (integral from 0 to eta_max)
  const double integral = integrated_values.back();
  
  const double K_bl = std::sqrt(rho_e * mu_e / (2.0 * d_ue_dx)) * integral / delta_bl;

  return std::make_pair(K_bl, coeff_finite_thickness);
}

} // namespace blast::boundary_layer::coefficients