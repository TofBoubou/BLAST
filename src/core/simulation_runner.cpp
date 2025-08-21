#include "blast/core/simulation_runner.hpp"
#include "blast/core/constants.hpp"
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <format>
#include <iomanip>
#include <iostream>

namespace blast::core {

auto SimulationRunner::run_simulation(
  boundary_layer::solver::BoundaryLayerSolver& solver,
  const thermophysics::MixtureInterface& mixture,
  const io::Configuration& config,
  const std::string& case_name,
  PerformanceMetrics& metrics) 
  -> std::expected<SimulationResult, ApplicationError> {
    
  if (config.abacus.enabled) {
    auto abacus_filename = run_abacus_generation(solver, mixture, config, case_name, metrics);
    if (!abacus_filename) {
      return std::unexpected(abacus_filename.error());
    }
    
    return SimulationResult{
      .solution = {},
      .is_abacus = true,
      .output_filename = abacus_filename.value()
    };
  } else {
    auto solution = run_standard_simulation(solver, metrics);
    if (!solution) {
      return std::unexpected(solution.error());
    }
    
    return SimulationResult{
      .solution = std::move(solution.value()),
      .is_abacus = false,
      .output_filename = ""
    };
  }
}

auto SimulationRunner::display_simulation_results(const SimulationResult& result,
                                                  const thermophysics::MixtureInterface& mixture,
                                                  const PerformanceMetrics& metrics) const -> void {
  if (result.is_abacus) {
    std::cout << "✓ Abacus generated successfully!" << std::endl;
    if (!result.output_filename.empty()) {
      auto abacus_file_size = std::filesystem::file_size(result.output_filename);
      std::cout << "  Abacus file: " << std::filesystem::path(result.output_filename).filename().string() 
                << " (" << std::setprecision(constants::string_processing::float_precision_2) << std::fixed
                << (abacus_file_size / constants::io::bytes_to_kb) << " KB)" << std::endl;
    }
  } else {
    display_standard_results(result.solution, mixture);
  }
}

auto SimulationRunner::run_standard_simulation(
  boundary_layer::solver::BoundaryLayerSolver& solver,
  PerformanceMetrics& metrics) 
  -> std::expected<boundary_layer::solver::SolutionResult, ApplicationError> {
    
  std::cout << "\n=== STARTING BOUNDARY LAYER SOLUTION ===" << std::endl;
  std::cout << "Starting solver.solve()..." << std::endl;
  
  auto solve_start = std::chrono::high_resolution_clock::now();
  auto solution_result = solver.solve();
  auto solve_end = std::chrono::high_resolution_clock::now();
  
  metrics.solve_time = std::chrono::duration_cast<std::chrono::milliseconds>(solve_end - solve_start);
  
  if (!solution_result) {
    return std::unexpected(ApplicationError{
      "Solver failed: " + solution_result.error().message(),
      constants::indexing::second
    });
  }
  
  std::cout << "✓ Solution completed successfully!" << std::endl;
  std::cout << "  Solve time: " << metrics.solve_time.count() << " ms" << std::endl;
  
  return std::move(solution_result.value());
}

auto SimulationRunner::run_abacus_generation(
  boundary_layer::solver::BoundaryLayerSolver& solver,
  const thermophysics::MixtureInterface& mixture,
  const io::Configuration& config,
  const std::string& case_name,
  PerformanceMetrics& metrics) 
  -> std::expected<std::string, ApplicationError> {
    
  std::cout << "\n=== GENERATING ABACUS ===" << std::endl;
  std::cout << "Abacus mode enabled - skipping normal simulation" << std::endl;
  
  display_abacus_info(config);
  
  auto abacus_start = std::chrono::high_resolution_clock::now();
  
  io::AbacusGenerator abacus_generator(solver, const_cast<thermophysics::MixtureInterface&>(mixture), config);
  auto abacus_result = abacus_generator.generate();
  
  if (!abacus_result.success) {
    return std::unexpected(ApplicationError{
      "Failed to generate abacus",
      constants::indexing::second
    });
  }
  
  // Save abacus results
  std::string abacus_filename = case_name + "_abacus.h5";
  std::filesystem::path abacus_path = std::filesystem::path(config.output.output_directory) / abacus_filename;
  
  if (!abacus_generator.save_results(abacus_result, abacus_path)) {
    return std::unexpected(ApplicationError{
      "Failed to save abacus results to: " + abacus_path.string(),
      constants::indexing::second
    });
  }
  
  auto abacus_end = std::chrono::high_resolution_clock::now();
  auto abacus_duration = std::chrono::duration_cast<std::chrono::milliseconds>(abacus_end - abacus_start);
  
  std::cout << "  Generation time: " << abacus_duration.count() << " ms" << std::endl;
  
  return abacus_path.string();
}

auto SimulationRunner::display_standard_results(const boundary_layer::solver::SolutionResult& solution,
                                               const thermophysics::MixtureInterface& mixture) const -> void {
  std::cout << "\n=== SOLUTION SUMMARY ===" << std::endl;
  std::cout << "Stations solved: " << solution.stations.size() << std::endl;
  std::cout << "Total iterations: " << solution.total_iterations << std::endl;
  std::cout << "Converged: " << (solution.converged ? "Yes" : "No") << std::endl;

  if (!solution.stations.empty()) {
    const auto& final_station = solution.stations.back();
    const auto n_eta = final_station.F.size();

    std::cout << "\nFinal station results:" << std::endl;
    std::cout << "  ξ = " << solution.xi_solved.back() << std::endl;
    std::cout << "  F at edge: " << final_station.F.back() << std::endl;
    std::cout << "  g at wall: " << final_station.g[constants::indexing::first] << std::endl;
    std::cout << "  g at edge: " << final_station.g.back() << std::endl;

    // Print species concentrations at wall and edge
    std::cout << "\nSpecies concentrations:" << std::endl;
    std::cout << "  Wall → Edge" << std::endl;
    for (std::size_t i = constants::indexing::first; i < mixture.n_species(); ++i) {
      std::cout << "  " << std::setw(constants::string_processing::medium_field_width) << mixture.species_name(i) 
                << ": " << std::setw(constants::string_processing::medium_field_width) << std::fixed
                << std::setprecision(constants::string_processing::float_precision_4) 
                << final_station.c(i, constants::indexing::first) 
                << " → " << std::setw(constants::string_processing::medium_field_width)
                << final_station.c(i, n_eta - constants::indexing::second) << std::endl;
    }

    // Print some profiles
    std::cout << "\nBoundary layer profiles (selected points):" << std::endl;
    std::cout << std::setw(constants::string_processing::wide_field_width) << "η/η_max" 
              << std::setw(constants::string_processing::wide_field_width) << "F" 
              << std::setw(constants::string_processing::wide_field_width) << "g" 
              << std::setw(constants::string_processing::wide_field_width) << "V" << std::endl;
    std::cout << std::string(4 * constants::string_processing::wide_field_width, '-') << std::endl;

    const std::size_t n_print = std::min(static_cast<std::size_t>(constants::indexing::max_profile_print_points), n_eta);
    for (std::size_t i = constants::indexing::first; i < n_print; ++i) {
      const std::size_t idx = (i * (n_eta - constants::indexing::second)) / (n_print - constants::indexing::second);
      const double eta_norm = static_cast<double>(idx) / (n_eta - constants::indexing::second);

      std::cout << std::setw(constants::string_processing::wide_field_width) << std::fixed 
                << std::setprecision(constants::string_processing::float_precision_2) << eta_norm 
                << std::setw(constants::string_processing::wide_field_width) << std::fixed
                << std::setprecision(constants::string_processing::float_precision_4) << final_station.F[idx] 
                << std::setw(constants::string_processing::wide_field_width) << std::fixed 
                << std::setprecision(constants::string_processing::float_precision_4) << final_station.g[idx] 
                << std::setw(constants::string_processing::wide_field_width) << std::fixed 
                << std::setprecision(constants::string_processing::float_precision_4) << final_station.V[idx]
                << std::endl;
    }
  }
}

auto SimulationRunner::display_abacus_info(const io::Configuration& config) const -> void {
  std::cout << "Temperature range: " << config.abacus.temperature_min << " - " 
            << config.abacus.temperature_max << " K" << std::endl;
  std::cout << "Temperature points: " << config.abacus.temperature_points << std::endl;
  std::cout << "Catalyticity values: ";
  
  for (size_t i = constants::indexing::first; i < config.abacus.catalyticity_values.size(); ++i) {
    if (i > constants::indexing::first) {
      std::cout << ", ";
    }
    std::cout << config.abacus.catalyticity_values[i];
  }
  std::cout << std::endl;
}

} // namespace blast::core