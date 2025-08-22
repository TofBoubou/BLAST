#pragma once
#include "../boundary_layer/solver/boundary_layer_solver.hpp"
#include "../io/abacus_generator.hpp"
#include "../io/config_types.hpp"
#include "../thermophysics/mixture_interface.hpp"
#include "application_types.hpp"
#include <expected>
#include <string>

namespace blast::core {

class SimulationRunner {
public:
  struct SimulationResult {
    boundary_layer::solver::SolutionResult solution;
    bool is_abacus = false;
    bool is_edge_reconstruction = false;
    std::string output_filename;
  };

  // Run simulation (either standard or abacus)
  [[nodiscard]] auto run_simulation(
    boundary_layer::solver::BoundaryLayerSolver& solver,
    const thermophysics::MixtureInterface& mixture,
    const io::Configuration& config,
    const std::string& case_name,
    PerformanceMetrics& metrics) 
    -> std::expected<SimulationResult, ApplicationError>;

  // Display simulation results
  auto display_simulation_results(const SimulationResult& result,
                                  const thermophysics::MixtureInterface& mixture,
                                  const PerformanceMetrics& metrics) const -> void;

private:
  // Run standard boundary layer simulation
  [[nodiscard]] auto run_standard_simulation(
    boundary_layer::solver::BoundaryLayerSolver& solver,
    PerformanceMetrics& metrics) 
    -> std::expected<boundary_layer::solver::SolutionResult, ApplicationError>;
    
  // Run abacus generation
  [[nodiscard]] auto run_abacus_generation(
    boundary_layer::solver::BoundaryLayerSolver& solver,
    const thermophysics::MixtureInterface& mixture,
    const io::Configuration& config,
    const std::string& case_name,
    PerformanceMetrics& metrics) 
    -> std::expected<std::string, ApplicationError>;
    
  // Run edge reconstruction
  [[nodiscard]] auto run_edge_reconstruction(
    const thermophysics::MixtureInterface& mixture,
    const io::Configuration& config,
    PerformanceMetrics& metrics) 
    -> std::expected<void, ApplicationError>;
    
  // Display standard simulation results
  auto display_standard_results(const boundary_layer::solver::SolutionResult& solution,
                                const thermophysics::MixtureInterface& mixture) const -> void;
                                
  // Display abacus generation info
  auto display_abacus_info(const io::Configuration& config) const -> void;
};

} // namespace blast::core