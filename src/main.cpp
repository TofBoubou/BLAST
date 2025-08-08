#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/io/config_manager.hpp"
#include "blast/io/output/hdf5_writer.hpp"
#include "blast/io/output/output_writer.hpp"
#include "blast/thermophysics/mixture_interface.hpp"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <vector>

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <config_file.yaml> [output_name]\n";
    return 1;
  }

  try {
    std::filesystem::path exe_path = std::filesystem::canonical(argv[0]).parent_path();

    std::filesystem::path mpp_data_path = exe_path / "libs" / "mutationpp" / "data";

    if (std::filesystem::exists(mpp_data_path)) {
      std::string mpp_data_str = std::filesystem::canonical(mpp_data_path).string();
      setenv("MPP_DATA_DIRECTORY", mpp_data_str.c_str(), 1);
      std::cout << "MPP_DATA_DIRECTORY auto-set to: " << mpp_data_str << std::endl;
    } else {
      std::cerr << "Warning: Mutation++ data directory not found at: " << mpp_data_path << std::endl;
      std::cerr << "Looking for MPP_DATA_DIRECTORY in environment..." << std::endl;
      if (const char* env_mpp = std::getenv("MPP_DATA_DIRECTORY")) {
        std::cout << "Using existing MPP_DATA_DIRECTORY: " << env_mpp << std::endl;
      } else {
        std::cerr << "Warning: MPP_DATA_DIRECTORY not set. Mutation++ may fail." << std::endl;
      }
    }
  } catch (const std::exception& e) {
    std::cerr << "Warning: Could not auto-configure MPP_DATA_DIRECTORY: " << e.what() << std::endl;
  }

  try {
    std::cout << "=== BLAST Boundary Layer Solver ===" << std::endl;
    std::cout << "Loading configuration from: " << argv[1] << std::endl;

    // Load configuration
    blast::io::ConfigurationManager config_manager;
    auto config_result = config_manager.load(argv[1]);
    if (!config_result) {
      std::cerr << "Failed to load config: " << config_result.error().message() << "\n";
      return 1;
    }
    auto config = config_result.value();
    std::cout << "✓ Configuration loaded successfully" << std::endl;

    // Create mixture
    std::cout << "Creating mixture: " << config.mixture.name << std::endl;
    auto mixture_result = blast::thermophysics::create_mixture(config.mixture);
    if (!mixture_result) {
      std::cerr << "Failed to create mixture: " << mixture_result.error().message() << "\n";
      return 1;
    }
    auto& mixture = *mixture_result.value();
    std::cout << "✓ Mixture created (" << mixture.n_species() << " species)" << std::endl;

    // Print mixture information
    std::cout << "\nMixture species:" << std::endl;
    for (std::size_t i = 0; i < mixture.n_species(); ++i) {
      std::cout << std::format("  [{:2}] {:>8} (MW: {:8.3f} kg/kmol)", i, mixture.species_name(i),
                               mixture.species_molecular_weight(i))
                << std::endl;
    }

    // Print simulation configuration
    std::cout << "\nSimulation setup:" << std::endl;
    std::cout << "  Body type: "
              << (config.simulation.body_type == blast::io::SimulationConfig::BodyType::Axisymmetric ? "Axisymmetric"
                                                                                                     : "Other")
              << std::endl;
    std::cout << "  Stagnation only: " << (config.simulation.only_stagnation_point ? "Yes" : "No") << std::endl;
    std::cout << "  Chemical mode: ";
    switch (config.simulation.chemical_mode) {
    case blast::io::SimulationConfig::ChemicalMode::Equilibrium:
      std::cout << "Equilibrium";
      break;
    case blast::io::SimulationConfig::ChemicalMode::Frozen:
      std::cout << "Frozen";
      break;
    case blast::io::SimulationConfig::ChemicalMode::NonEquilibrium:
      std::cout << "Non-equilibrium";
      break;
    }
    std::cout << std::endl;
    std::cout << "  Thermal diffusion: " << (config.simulation.consider_thermal_diffusion ? "Yes" : "No") << std::endl;
    std::cout << "  Grid points (η): " << config.numerical.n_eta << std::endl;
    std::cout << "  η_max: " << config.numerical.eta_max << std::endl;
    std::cout << "  Convergence tol: " << config.numerical.convergence_tolerance << std::endl;

    // Print edge conditions
    if (!config.outer_edge.edge_points.empty()) {
      auto edge = config.outer_edge.edge_points[0];
      std::cout << "\nEdge conditions:" << std::endl;
      std::cout << "  Pressure: " << edge.pressure << " Pa" << std::endl;
      std::cout << "  Temperature: " << edge.temperature << " K" << std::endl;
      std::cout << "  Enthalpy: " << edge.enthalpy / 1000.0 << " kJ/kg" << std::endl;
      if (!config.wall_parameters.wall_temperatures.empty()) {
        std::cout << "  Wall temp: " << config.wall_parameters.wall_temperatures[0] << " K" << std::endl;
      } else {
        std::cout << "  Wall temp: NOT SET" << std::endl;
      }
    }

    // Initialize HDF5 library
    std::cout << "\nInitializing output system..." << std::endl;
    if (auto hdf5_init = blast::io::output::hdf5::initialize(); !hdf5_init) {
      std::cerr << "Warning: Failed to initialize HDF5: " << hdf5_init.error().message() << std::endl;
    } else {
      if (auto version = blast::io::output::hdf5::check_version()) {
        std::cout << "✓ HDF5 library version: " << version.value() << std::endl;
      }
    }

    // Configure output system
    blast::io::output::OutputConfig output_config;
    output_config.base_directory = config.output.output_directory;
    output_config.primary_format = blast::io::output::OutputFormat::HDF5;
    output_config.additional_formats = {};
    output_config.save_metadata = true;
    output_config.save_derivatives = true;
    output_config.compress_data = true;
    output_config.include_timestamp = true;

    // Select variables to save based on simulation type
    output_config.variables.flow_variables = true;
    output_config.variables.species_concentrations = true;

    blast::io::output::OutputWriter output_writer(output_config);

    // Validate output configuration
    if (auto validation = output_writer.validate_config(); !validation) {
      std::cerr << "Output configuration error: " << validation.error().message() << std::endl;
      return 1;
    }

    std::cout << "✓ Output system configured" << std::endl;

    // Show planned output files
    std::string case_name = (argc > 2) ? argv[2] : "simulation";
    auto output_info = output_writer.get_output_info(case_name);
    std::cout << "\nPlanned output files:" << std::endl;
    for (const auto& [format, path] : output_info) {
      std::string format_name;
      switch (format) {
      case blast::io::output::OutputFormat::HDF5:
        format_name = "HDF5";
        break;
      default:
        format_name = "Unknown";
        break;
      }
      std::cout << "  " << format_name << ": " << path.string() << std::endl;
    }

    // Create and run solver
    std::cout << "\n=== STARTING BOUNDARY LAYER SOLUTION ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    std::cout << "Creating BoundaryLayerSolver..." << std::endl;
    blast::boundary_layer::solver::BoundaryLayerSolver solver(mixture, config);
    std::cout << "BoundaryLayerSolver created successfully." << std::endl;

    std::cout << "Starting solver.solve()..." << std::endl;
    auto solution_result = solver.solve();
    if (!solution_result) {
      std::cerr << "Solver failed: " << solution_result.error().message() << "\n";
      return 1;
    }

    auto& solution = solution_result.value();
    auto solve_time = std::chrono::high_resolution_clock::now();
    auto solve_duration = std::chrono::duration_cast<std::chrono::milliseconds>(solve_time - start_time);

    std::cout << "✓ Solution completed successfully!" << std::endl;
    std::cout << "  Solve time: " << solve_duration.count() << " ms" << std::endl;

    // Print solution summary
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
      std::cout << "  g at wall: " << final_station.g[0] << std::endl;
      std::cout << "  g at edge: " << final_station.g.back() << std::endl;

      // Print species concentrations at wall and edge
      std::cout << "\nSpecies concentrations:" << std::endl;
      std::cout << "  Wall → Edge" << std::endl;
      for (std::size_t i = 0; i < mixture.n_species(); ++i) {
        std::cout << "  " << std::setw(8) << mixture.species_name(i) << ": " << std::setw(8) << std::fixed
                  << std::setprecision(4) << final_station.c(i, 0) << " → " << std::setw(8)
                  << final_station.c(i, n_eta - 1) << std::endl;
      }

      // Print some profiles
      std::cout << "\nBoundary layer profiles (selected points):" << std::endl;
      std::cout << std::setw(6) << "η/η_max" << std::setw(10) << "F" << std::setw(10) << "g" << std::setw(10) << "V"
                << std::endl;
      std::cout << std::string(36, '-') << std::endl;

      const std::size_t n_print = std::min(static_cast<std::size_t>(11), n_eta);
      for (std::size_t i = 0; i < n_print; ++i) {
        const std::size_t idx = (i * (n_eta - 1)) / (n_print - 1);
        const double eta_norm = static_cast<double>(idx) / (n_eta - 1);

        std::cout << std::setw(6) << std::fixed << std::setprecision(2) << eta_norm << std::setw(10) << std::fixed
                  << std::setprecision(4) << final_station.F[idx] << std::setw(10) << std::fixed << std::setprecision(4)
                  << final_station.g[idx] << std::setw(10) << std::fixed << std::setprecision(4) << final_station.V[idx]
                  << std::endl;
      }
    }

    // Write output files
    std::cout << "\n=== WRITING OUTPUT FILES ===" << std::endl;

    // Progress callback for output writing
    auto progress_callback = [](double progress, const std::string& stage) {
      std::cout << "\r" << std::setw(60) << std::left
                << ("  " + stage + " [" + std::to_string(static_cast<int>(progress * 100.0)) + "%]") << std::flush;
    };

    auto output_start = std::chrono::high_resolution_clock::now();
    auto output_result = output_writer.write_solution(solution, config, mixture, case_name, progress_callback);
    auto output_end = std::chrono::high_resolution_clock::now();
    auto output_duration = std::chrono::duration_cast<std::chrono::milliseconds>(output_end - output_start);

    std::cout << std::endl; // New line after progress

    if (!output_result) {
      std::cerr << "Failed to write output: " << output_result.error().message() << std::endl;
      return 1;
    }

    auto& output_files = output_result.value();
    std::cout << "✓ Output written successfully!" << std::endl;
    std::cout << "  Output time: " << output_duration.count() << " ms" << std::endl;
    std::cout << "\nGenerated files:" << std::endl;
    for (const auto& file_path : output_files) {
      auto file_size = std::filesystem::file_size(file_path);
      std::cout << "  " << file_path.filename().string() << " (" << std::setprecision(2) << std::fixed
                << (file_size / 1024.0) << " KB)" << std::endl;
    }

    // Performance summary
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(output_end - start_time);
    std::cout << "\n=== PERFORMANCE SUMMARY ===" << std::endl;
    std::cout << "Total runtime: " << total_time.count() << " ms" << std::endl;
    std::cout << "  Solution: " << solve_duration.count() << " ms (" << std::setprecision(1) << std::fixed
              << (100.0 * solve_duration.count() / total_time.count()) << "%)" << std::endl;
    std::cout << "  Output: " << output_duration.count() << " ms (" << std::setprecision(1) << std::fixed
              << (100.0 * output_duration.count() / total_time.count()) << "%)" << std::endl;

    std::cout << "\n=== CALCULATION COMPLETED SUCCESSFULLY ===" << std::endl;

    // Post-processing recommendations
    std::cout << "\nPost-processing recommendations:" << std::endl;
    std::cout << "  • Open .h5 files with HDFView or Python (h5py, pandas)" << std::endl;

    return 0;

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";

    // Cleanup
    blast::io::output::hdf5::finalize();
    return 1;
  }
}