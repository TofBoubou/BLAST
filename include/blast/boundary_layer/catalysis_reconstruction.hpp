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

class CatalysisReconstructor {
public:
  struct ReconstructedCatalysisConditions {
    double catalyticity;              // Optimal catalyticity found (0-1)
    double heat_flux_achieved;        // [W/m²] Heat flux achieved
    double edge_temperature;          // [K] Edge temperature (fixed input)
    double pressure;                  // [Pa] Edge pressure (input)
    double enthalpy;                  // [J/kg] Edge enthalpy
    std::vector<double> mass_fractions;  // Species mass fractions at edge
    double density;                   // [kg/m³] Edge density
    double viscosity;                 // [Pa·s] Edge viscosity
    int iterations_used;              // Number of iterations needed
    solver::SolutionResult full_solution;  // Complete boundary layer solution at optimal catalyticity
  };

  CatalysisReconstructor(
      const io::CatalysisReconstructionConfig& config,
      const io::Configuration& full_config,
      thermophysics::MixtureInterface& mixture)
      : config_(config),
        full_config_(full_config),
        mixture_(mixture),
        gsi_manager_(full_config) {}

  auto reconstruct() -> std::expected<ReconstructedCatalysisConditions, solver::SolverError>;

private:
  auto compute_heat_flux_at_catalyticity(double catalyticity) 
      -> std::expected<double, solver::SolverError>;

  auto setup_boundary_conditions(double catalyticity) 
      -> std::expected<io::Configuration, solver::SolverError>;

  const io::CatalysisReconstructionConfig& config_;
  const io::Configuration& full_config_;
  thermophysics::MixtureInterface& mixture_;
  io::GsiManager gsi_manager_;
};

} // namespace blast::boundary_layer