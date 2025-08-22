#pragma once

#include "blast/io/config_types.hpp"
#include "blast/io/gsi_manager.hpp"
#include "blast/core/exceptions.hpp"
#include "blast/boundary_layer/conditions/boundary_conditions.hpp"
#include "blast/thermophysics/mixture_interface.hpp"
#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include <expected>
#include <memory>
#include <format>

namespace blast::boundary_layer {

class EdgeTemperatureReconstructor {
public:
  struct ReconstructedEdgeConditions {
    double temperature;           // [K] Edge temperature found
    double pressure;              // [Pa] Edge pressure (input)
    double enthalpy;              // [J/kg] Edge enthalpy
    std::vector<double> mass_fractions;  // Species mass fractions at edge
    double density;               // [kg/m³] Edge density
    double viscosity;             // [Pa·s] Edge viscosity
    double heat_flux_achieved;    // [W/m²] Actual heat flux achieved
    int iterations_used;          // Number of iterations needed
    solver::SolutionResult full_solution;  // Complete boundary layer solution at optimal T_edge
  };

  EdgeTemperatureReconstructor(
      const io::EdgeReconstructionConfig& config,
      const io::Configuration& full_config,
      thermophysics::MixtureInterface& mixture)
      : config_(config),
        full_config_(full_config),
        mixture_(mixture),
        gsi_manager_(full_config) {}

  auto reconstruct() -> std::expected<ReconstructedEdgeConditions, solver::SolverError>;

private:
  auto compute_heat_flux_at_temperature(double T_edge) 
      -> std::expected<double, solver::SolverError>;

  auto setup_boundary_conditions(double T_edge) 
      -> std::expected<io::Configuration, solver::SolverError>;

  const io::EdgeReconstructionConfig& config_;
  const io::Configuration& full_config_;
  thermophysics::MixtureInterface& mixture_;
  io::GsiManager gsi_manager_;
};

} // namespace blast::boundary_layer