#pragma once
#include "../../core/exceptions.hpp"
#include "../../io/config_types.hpp"
#include "../conditions/boundary_conditions.hpp"
#include "equation_types.hpp"
#include <cmath>
#include <expected>
#include <format>
#include <iostream>
#include <source_location>

namespace blast::boundary_layer::equations {

// Compute geometry-dependent factors for equation coefficients
[[nodiscard]] auto
compute_geometry_factors(int station, PhysicalQuantity auto xi, const conditions::BoundaryConditions& bc,
                         const io::SimulationConfig& sim_config) -> std::expected<GeometryFactors, EquationError> {

  if (station < 0) {
    return std::unexpected(EquationError("Invalid negative station number"));
  }

  if (station == 0) {
    // Stagnation point - factors depend on body type
    switch (sim_config.body_type) {
    case io::SimulationConfig::BodyType::Axisymmetric: {
      const double d_ue_dx_val = bc.d_ue_dx();
      if (d_ue_dx_val <= 0.0) {
        return std::unexpected(EquationError(std::format("Invalid d_ue_dx: {} (must be positive)", d_ue_dx_val)));
      }

      const double sqrt_term = std::sqrt(2.0 * bc.rho_e() * bc.mu_e() * d_ue_dx_val);
      return GeometryFactors(1.0 / sqrt_term, 1.0 / (2.0 * d_ue_dx_val),
                             std::sqrt(bc.rho_e() * bc.mu_e() / (2.0 * d_ue_dx_val)));
    }
    case io::SimulationConfig::BodyType::TwoD: {
      const double d_ue_dx_val = bc.d_ue_dx();
      if (d_ue_dx_val <= 0.0) {
        return std::unexpected(EquationError(std::format("Invalid d_ue_dx: {} (must be positive)", d_ue_dx_val)));
      }

      const double sqrt_term = std::sqrt(bc.rho_e() * bc.mu_e() * d_ue_dx_val);
      return GeometryFactors(1.0 / sqrt_term, 1.0 / d_ue_dx_val, std::sqrt(bc.rho_e() * bc.mu_e() / d_ue_dx_val));
    }
    case io::SimulationConfig::BodyType::Cone:
    case io::SimulationConfig::BodyType::FlatPlate:
      return GeometryFactors(1.0, 0.0, 0.0);
    default:
      return std::unexpected(EquationError("Unknown body type"));
    }
  } else {
    // Downstream station - validate xi and other parameters
    if (xi <= 0.0 || !std::isfinite(xi)) {
      return std::unexpected(EquationError(std::format("Invalid xi: {} (must be positive and finite)", xi)));
    }

    if (bc.d_xi_dx() <= 0.0 || bc.ue() <= 0.0 || bc.r_body() <= 0.0) {
      std::cout << " bc.d_xi_dx() " << bc.d_xi_dx() << " bc.ue() " << bc.ue() << " bc.r_body() " << bc.r_body()
                << std::endl;
      return std::unexpected(EquationError("Invalid boundary conditions for downstream station"));
    }

    const double sqrt_2xi = std::sqrt(2.0 * xi);
    return GeometryFactors(bc.r_body() * sqrt_2xi / bc.d_xi_dx(), 2.0 * xi / (bc.ue() * bc.d_xi_dx()),
                           sqrt_2xi / (bc.ue() * bc.r_body()));
  }
}

[[nodiscard]] auto
compute_energy_j_factor(int station, PhysicalQuantity auto xi, const conditions::BoundaryConditions& bc,
                        const io::SimulationConfig& sim_config) -> std::expected<double, EquationError> {

  if (station < 0) {
    return std::unexpected(EquationError("Invalid negative station number"));
  }

  if (station == 0) {
    // Stagnation point
    switch (sim_config.body_type) {
    case io::SimulationConfig::BodyType::Axisymmetric: {
      const double d_ue_dx_val = bc.d_ue_dx();
      if (d_ue_dx_val <= 0.0) {
        return std::unexpected(EquationError(std::format("Invalid d_ue_dx: {} (must be positive)", d_ue_dx_val)));
      }
      return 1.0 / std::sqrt(2.0 * bc.rho_e() * bc.mu_e() * d_ue_dx_val);
    }
    case io::SimulationConfig::BodyType::TwoD: {
      const double d_ue_dx_val = bc.d_ue_dx();
      if (d_ue_dx_val <= 0.0) {
        return std::unexpected(EquationError(std::format("Invalid d_ue_dx: {} (must be positive)", d_ue_dx_val)));
      }
      return 1.0 / std::sqrt(bc.rho_e() * bc.mu_e() * d_ue_dx_val);
    }
    case io::SimulationConfig::BodyType::Cone:
    case io::SimulationConfig::BodyType::FlatPlate:
      return 1.0;
    default:
      return std::unexpected(EquationError("Unknown body type"));
    }
  } else {
    // Downstream station
    if (xi <= 0.0 || !std::isfinite(xi)) {
      return std::unexpected(EquationError(std::format("Invalid xi: {} (must be positive and finite)", xi)));
    }

    if (bc.d_xi_dx() <= 0.0 || bc.r_body() <= 0.0) {
      std::cout << " bc.d_xi_dx() " << bc.d_xi_dx() << " bc.r_body() " << bc.r_body() << std::endl;
      return std::unexpected(EquationError("Invalid boundary conditions for downstream station"));
    }

    return std::sqrt(2.0 * xi) * bc.r_body() / bc.d_xi_dx();
  }
}

} // namespace blast::boundary_layer::equations