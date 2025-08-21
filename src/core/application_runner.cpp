#include "blast/core/application_runner.hpp"
#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/core/constants.hpp"
#include "blast/io/output/hdf5_writer.hpp"
#include <chrono>
#include <filesystem>
#include <iostream>
#include <string>

namespace blast::core {

ApplicationRunner::ApplicationRunner()
  : environment_manager_(std::make_unique<EnvironmentManager>())
  , config_loader_(std::make_unique<ConfigurationLoader>())
  , output_manager_(std::make_unique<OutputManager>())
  , simulation_runner_(std::make_unique<SimulationRunner>()) {
}

ApplicationRunner::~ApplicationRunner() = default;

auto ApplicationRunner::run(int argc, char* argv[]) -> ApplicationResult {
  auto start_time = std::chrono::high_resolution_clock::now();
  PerformanceMetrics metrics;
  
  try {
    // Parse command line arguments
    auto args_result = parse_command_line(argc, argv);
    if (!args_result) {
      return handle_error(args_result.error());
    }
    auto args = args_result.value();
    
    if (args.help_requested) {
      display_usage(argv[constants::indexing::first]);
      return {true, constants::indexing::first, "Help displayed"};
    }
    
    display_header();
    
    // Configure environment
    std::filesystem::path exe_path = std::filesystem::canonical(argv[constants::indexing::first]).parent_path();
    if (auto env_result = environment_manager_->configure_environment(exe_path); !env_result) {
      return handle_error(env_result.error());
    }
    
    // Load configuration and create mixture
    auto config_result = config_loader_->load_configuration(args.config_file);
    if (!config_result) {
      return handle_error(config_result.error());
    }
    auto [config, mixture] = std::move(config_result.value());
    
    // Initialize output system
    if (auto output_init = output_manager_->initialize_output_system(config, args.output_name); !output_init) {
      return handle_error(output_init.error());
    }
    
    // Create solver
    std::cout << "Creating BoundaryLayerSolver..." << std::endl;
    boundary_layer::solver::BoundaryLayerSolver solver(*mixture, config);
    std::cout << "BoundaryLayerSolver created successfully." << std::endl;
    
    // Run simulation
    auto simulation_result = simulation_runner_->run_simulation(solver, *mixture, config, args.output_name, metrics);
    if (!simulation_result) {
      return handle_error(simulation_result.error());
    }
    auto result = std::move(simulation_result.value());
    
    // Write output files (only for standard simulations)
    if (!result.is_abaque) {
      auto output_result = output_manager_->write_simulation_results(
        result.solution, config, *mixture, args.output_name, metrics);
      if (!output_result) {
        return handle_error(output_result.error());
      }
    }
    
    // Calculate total time
    auto end_time = std::chrono::high_resolution_clock::now();
    metrics.total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    // Display results and performance
    simulation_runner_->display_simulation_results(result, *mixture, metrics);
    display_performance_summary(metrics);
    display_completion_message();
    
    cleanup();
    
    return {true, constants::indexing::first, "Success"};
    
  } catch (const std::exception& e) {
    cleanup();
    return handle_error(ApplicationError{"Unexpected error: " + std::string(e.what()), constants::indexing::second});
  }
}

auto ApplicationRunner::parse_command_line(int argc, char* argv[]) 
  -> std::expected<CommandLineArgs, ApplicationError> {
    
  constexpr int min_required_args = 2;
  
  if (argc < min_required_args) {
    return std::unexpected(ApplicationError{
      "Insufficient arguments provided",
      constants::indexing::second
    });
  }
  
  CommandLineArgs args;
  args.config_file = argv[constants::indexing::second];
  
  constexpr int output_name_arg_index = 3;
  if (argc > min_required_args) {
    args.output_name = argv[output_name_arg_index - constants::indexing::second];
  }
  
  return args;
}

auto ApplicationRunner::display_usage(const std::string& program_name) const -> void {
  std::cerr << "Usage: " << program_name << " <config_file.yaml> [output_name]\n";
}

auto ApplicationRunner::display_header() const -> void {
  std::cout << "=== BLAST Boundary Layer Solver ===" << std::endl;
}

auto ApplicationRunner::display_performance_summary(const PerformanceMetrics& metrics) const -> void {
  std::cout << "\n=== PERFORMANCE SUMMARY ===" << std::endl;
  std::cout << "Total runtime: " << metrics.total_time.count() << " ms" << std::endl;
}

auto ApplicationRunner::display_completion_message() const -> void {
  std::cout << "\n=== CALCULATION COMPLETED SUCCESSFULLY ===" << std::endl;
  std::cout << "\nPost-processing recommendations:" << std::endl;
  std::cout << "  • Open .h5 files with HDFView or Python (h5py, pandas)" << std::endl;
}

auto ApplicationRunner::cleanup() -> void {
  io::output::hdf5::finalize();
}

auto ApplicationRunner::handle_error(const ApplicationError& error) -> ApplicationResult {
  std::cerr << "Error: " << error.message << std::endl;
  cleanup();
  return {false, error.exit_code, error.message};
}

} // namespace blast::core