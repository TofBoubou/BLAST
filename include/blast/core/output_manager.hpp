#pragma once
#include "../io/config_types.hpp"
#include "../io/output/output_writer.hpp"
#include "../boundary_layer/solver/boundary_layer_solver.hpp"
#include "../thermophysics/mixture_interface.hpp"
#include "application_types.hpp"
#include <expected>
#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace blast::core {

class OutputManager {
public:
  using ProgressCallback = std::function<void(double, const std::string&)>;

  // Initialize output system
  [[nodiscard]] auto initialize_output_system(const io::Configuration& config, 
                                              const std::string& case_name) 
    -> std::expected<void, ApplicationError>;

  // Write simulation results
  [[nodiscard]] auto write_simulation_results(
    const boundary_layer::solver::SolutionResult& solution,
    const io::Configuration& config,
    const thermophysics::MixtureInterface& mixture,
    const std::string& case_name,
    PerformanceMetrics& metrics) 
    -> std::expected<std::vector<std::filesystem::path>, ApplicationError>;

  // Display planned output files
  auto display_planned_outputs(const std::string& case_name) const -> void;

private:
  std::unique_ptr<io::output::OutputWriter> output_writer_;
  
  // Initialize HDF5 system
  [[nodiscard]] auto initialize_hdf5() -> std::expected<void, ApplicationError>;
  
  // Create output configuration
  [[nodiscard]] auto create_output_config(const io::Configuration& config) 
    -> io::output::OutputConfig;
    
  // Progress callback for output operations
  auto create_progress_callback() -> ProgressCallback;
};

} // namespace blast::core