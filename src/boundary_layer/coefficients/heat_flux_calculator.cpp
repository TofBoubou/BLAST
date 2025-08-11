#include "blast/boundary_layer/coefficients/heat_flux_calculator.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

namespace blast::boundary_layer::coefficients {

auto HeatFluxCalculator::calculate(const CoefficientInputs& inputs, const CoefficientSet& coeffs,
                                   const conditions::BoundaryConditions& bc, const std::vector<double>& dT_deta,
                                   int station, double xi) const -> std::expected<HeatFluxCoefficients, HeatFluxError> {

  const auto n_eta = inputs.T.size();
  const auto n_species = inputs.c.rows();

  if (n_eta == 0 || dT_deta.size() != n_eta) {
    return std::unexpected(HeatFluxError("Invalid input dimensions for heat flux calculation"));
  }

  HeatFluxCoefficients heat_flux;

  auto geo_factors_result = compute_heat_flux_geometry_factors(station, xi, bc, coeffs.wall);
  if (!geo_factors_result) {
    return std::unexpected(geo_factors_result.error());
  }
  auto geo_factors = geo_factors_result.value();

  if (!geo_factors.valid_geometry) {
    heat_flux.q_conductive_dimensional.resize(n_eta, 0.0);
    heat_flux.q_diffusive_dimensional.resize(n_eta, 0.0);
    heat_flux.q_total_dimensional.resize(n_eta, 0.0);
    heat_flux.q_conductive_nondimensional.resize(n_eta, 0.0);
    heat_flux.q_diffusive_nondimensional.resize(n_eta, 0.0);
    heat_flux.q_total_nondimensional.resize(n_eta, 0.0);

    heat_flux.q_diffusive_species_dimensional = core::Matrix<double>(n_species, n_eta);
    heat_flux.q_diffusive_species_nondimensional = core::Matrix<double>(n_species, n_eta);
    heat_flux.q_diffusive_species_dimensional.setZero();
    heat_flux.q_diffusive_species_nondimensional.setZero();

    heat_flux.q_ref = 1.0;
    return heat_flux;
  }

  heat_flux.q_ref = compute_reference_flux(bc, coeffs);
  // std::cout << "QREF = " << heat_flux.q_ref << std::endl;

  auto k_local_result = compute_local_conductivities(inputs, bc);
  if (!k_local_result) {
    return std::unexpected(k_local_result.error());
  }
  auto k_local = k_local_result.value();

  heat_flux.q_conductive_dimensional = compute_conductive_flux_profile(dT_deta, k_local, geo_factors);

  auto [q_diffusive_total, q_diffusive_species] = compute_diffusive_flux_profile(coeffs);
  heat_flux.q_diffusive_dimensional = std::move(q_diffusive_total);
  heat_flux.q_diffusive_species_dimensional = std::move(q_diffusive_species);

  heat_flux.q_total_dimensional.resize(n_eta);
  for (std::size_t i = 0; i < n_eta; ++i) {
    heat_flux.q_total_dimensional[i] = heat_flux.q_conductive_dimensional[i] + heat_flux.q_diffusive_dimensional[i];
  }

  auto wall_fluxes_result = compute_wall_heat_fluxes(inputs, coeffs, bc, station, xi);
  if (!wall_fluxes_result) {
    return std::unexpected(wall_fluxes_result.error());
  }
  auto [q_wall_cond, q_wall_diff, q_wall_total] = wall_fluxes_result.value();

  heat_flux.q_wall_conductive_dim = q_wall_cond;
  heat_flux.q_wall_diffusive_dim = q_wall_diff;
  heat_flux.q_wall_total_dim = q_wall_total;

  if (heat_flux.q_ref > 0.0) {
    heat_flux.q_conductive_nondimensional.resize(n_eta);
    heat_flux.q_diffusive_nondimensional.resize(n_eta);
    heat_flux.q_total_nondimensional.resize(n_eta);

    for (std::size_t i = 0; i < n_eta; ++i) {
      heat_flux.q_conductive_nondimensional[i] = heat_flux.q_conductive_dimensional[i] / heat_flux.q_ref;
      heat_flux.q_diffusive_nondimensional[i] = heat_flux.q_diffusive_dimensional[i] / heat_flux.q_ref;
      heat_flux.q_total_nondimensional[i] = heat_flux.q_total_dimensional[i] / heat_flux.q_ref;
    }

    heat_flux.q_diffusive_species_nondimensional = core::Matrix<double>(n_species, n_eta);
    for (std::size_t j = 0; j < n_species; ++j) {
      for (std::size_t i = 0; i < n_eta; ++i) {
        heat_flux.q_diffusive_species_nondimensional(j, i) =
            heat_flux.q_diffusive_species_dimensional(j, i) / heat_flux.q_ref;
      }
    }

    heat_flux.q_wall_conductive_nondim = heat_flux.q_wall_conductive_dim / heat_flux.q_ref;
    heat_flux.q_wall_diffusive_nondim = heat_flux.q_wall_diffusive_dim / heat_flux.q_ref;
    heat_flux.q_wall_total_nondim = heat_flux.q_wall_total_dim / heat_flux.q_ref;
  }

  return heat_flux;
}

auto HeatFluxCalculator::compute_heat_flux_geometry_factors(
    int station, double xi, const conditions::BoundaryConditions& bc,
    const WallProperties& wall_props) const -> std::expected<HeatFluxGeometryFactors, HeatFluxError> {

  HeatFluxGeometryFactors factors;
  factors.valid_geometry = true;

  if (station == 0) {
    switch (sim_config_.body_type) {
    case io::SimulationConfig::BodyType::Axisymmetric: {
      if (bc.d_ue_dx() <= 0.0 || bc.rho_e() <= 0.0 || bc.mu_e() <= 0.0) {
        return std::unexpected(HeatFluxError("Invalid boundary conditions for axisymmetric stagnation point"));
      }
      factors.der_fact = std::sqrt(2.0 * bc.d_ue_dx() / (bc.rho_e() * bc.mu_e())) * wall_props.rho_wall;
      factors.dy_deta_factor = std::sqrt(bc.rho_e() * bc.mu_e() / (2.0 * bc.d_ue_dx()));
      break;
    }
    case io::SimulationConfig::BodyType::TwoD: {
      if (bc.d_ue_dx() <= 0.0 || bc.rho_e() <= 0.0 || bc.mu_e() <= 0.0) {
        return std::unexpected(HeatFluxError("Invalid boundary conditions for 2D stagnation point"));
      }
      factors.der_fact = std::sqrt(bc.d_ue_dx() / (bc.rho_e() * bc.mu_e())) * wall_props.rho_wall;
      factors.dy_deta_factor = std::sqrt(bc.rho_e() * bc.mu_e() / bc.d_ue_dx());
      break;
    }
    case io::SimulationConfig::BodyType::Cone:
    case io::SimulationConfig::BodyType::FlatPlate:
      factors.valid_geometry = false;
      factors.der_fact = 0.0;
      factors.dy_deta_factor = 1.0;
      break;
    default:
      return std::unexpected(HeatFluxError("Unknown body type in heat flux calculation"));
    }
  } else {
    if (xi <= 0.0 || bc.ue() <= 0.0 || bc.r_body() <= 0.0) {
      return std::unexpected(HeatFluxError("Invalid conditions for downstream heat flux calculation"));
    }

    factors.der_fact = bc.ue() * bc.r_body() / std::sqrt(2.0 * xi) * wall_props.rho_wall;
    factors.dy_deta_factor = std::sqrt(2.0 * xi) / (bc.ue() * bc.r_body());
  }

  return factors;
}

auto HeatFluxCalculator::compute_local_conductivities(const CoefficientInputs& inputs,
                                                      const conditions::BoundaryConditions& bc) const
    -> std::expected<std::vector<double>, HeatFluxError> {

  const auto n_eta = inputs.T.size();
  const auto n_species = inputs.c.rows();

  std::vector<double> k_local(n_eta);

  for (std::size_t i = 0; i < n_eta; ++i) {
    std::vector<double> c_local(n_species);
    for (std::size_t j = 0; j < n_species; ++j) {
      c_local[j] = inputs.c(j, i);
    }

    auto k_result = mixture_.frozen_thermal_conductivity(c_local, inputs.T[i], bc.P_e());
    if (!k_result) {
      return std::unexpected(HeatFluxError(
          std::format("Failed to compute thermal conductivity at eta point {}: {}", i, k_result.error().message())));
    }

    k_local[i] = k_result.value();

    if (k_local[i] <= 0.0 || !std::isfinite(k_local[i])) {
      return std::unexpected(
          HeatFluxError(std::format("Invalid thermal conductivity at eta point {}: {}", i, k_local[i])));
    }
  }

  return k_local;
}

auto HeatFluxCalculator::compute_reference_flux(const conditions::BoundaryConditions& bc,
                                                const CoefficientSet& coeffs) const -> double {

  const double rho_e = bc.rho_e();
  const double h_wall = coeffs.thermodynamic.h_wall;
  const double delta_h = std::abs(bc.he() + 0.5 * bc.ue() * bc.ue() - h_wall);
  const double q_ref = rho_e * std::sqrt(2.0 * bc.d_ue_dx() * bc.r_body()) * delta_h;

  return q_ref;
}

auto HeatFluxCalculator::compute_conductive_flux_profile(
    const std::vector<double>& dT_deta, const std::vector<double>& k_local,
    const HeatFluxGeometryFactors& geo_factors) const -> std::vector<double> {

  const auto n_eta = dT_deta.size();
  std::vector<double> q_conductive(n_eta);

  for (std::size_t i = 0; i < n_eta; ++i) {
    const double dT_dy = dT_deta[i] / geo_factors.dy_deta_factor;
    q_conductive[i] = std::abs(-k_local[i] * dT_dy);
  }

  return q_conductive;
}

auto HeatFluxCalculator::compute_diffusive_flux_profile(const CoefficientSet& coeffs) const
    -> std::pair<std::vector<double>, core::Matrix<double>> {

  const auto n_eta = coeffs.diffusion.J.cols();
  const auto n_species = coeffs.diffusion.J.rows();

  std::vector<double> q_diffusive_total(n_eta, 0.0);
  core::Matrix<double> q_diffusive_species(n_species, n_eta);

  for (std::size_t i = 0; i < n_eta; ++i) {
    // First compute sum without abs (matching old BLAST behavior)
    for (std::size_t j = 0; j < n_species; ++j) {
      const double q_species_j = coeffs.diffusion.J(j, i) * coeffs.h_species(j, i);
      q_diffusive_total[i] += q_species_j;
    }
    // Apply abs to total sum (matching old BLAST behavior)
    q_diffusive_total[i] = std::abs(q_diffusive_total[i]);
    
    // For individual species, store the absolute values for output
    for (std::size_t j = 0; j < n_species; ++j) {
      q_diffusive_species(j, i) = std::abs(coeffs.diffusion.J(j, i) * coeffs.h_species(j, i));
    }
  }

  return {std::move(q_diffusive_total), std::move(q_diffusive_species)};
}

auto HeatFluxCalculator::compute_wall_heat_fluxes(
    const CoefficientInputs& inputs, const CoefficientSet& coeffs, const conditions::BoundaryConditions& bc,
    int station, double xi) const -> std::expected<std::tuple<double, double, double>, HeatFluxError> {

  auto geo_factors_result = compute_heat_flux_geometry_factors(station, xi, bc, coeffs.wall);
  if (!geo_factors_result) {
    return std::unexpected(geo_factors_result.error());
  }
  auto geo_factors = geo_factors_result.value();

  if (!geo_factors.valid_geometry) {
    return std::make_tuple(0.0, 0.0, 0.0);
  }

  const double dT_deta_wall = inputs.T.size() > 1 ? (inputs.T[1] - inputs.T[0]) / d_eta_ : 0.0;

  // Compute conductive heat flux matching the old BLAST behavior
  double q_wall_conductive = -coeffs.wall.k_wall * dT_deta_wall * geo_factors.der_fact;
  q_wall_conductive = -q_wall_conductive; // To the wall (matching old code)

  double q_wall_diffusive = 0.0;
  const auto n_species = coeffs.diffusion.J.rows();

  // Compute diffusive heat flux matching the old BLAST behavior
  for (std::size_t j = 0; j < n_species; ++j) {
    q_wall_diffusive += coeffs.diffusion.J(j, 0) * coeffs.h_species(j, 0);
  }
  q_wall_diffusive = -q_wall_diffusive; // To the wall (matching old code)

  const double q_wall_total = q_wall_conductive + q_wall_diffusive;

  return std::make_tuple(q_wall_conductive, q_wall_diffusive, q_wall_total);
}

} // namespace blast::boundary_layer::coefficients