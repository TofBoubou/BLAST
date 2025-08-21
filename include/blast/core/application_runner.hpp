#pragma once
#include "application_types.hpp"
#include "configuration_loader.hpp"
#include "environment_manager.hpp"
#include "output_manager.hpp"
#include "simulation_runner.hpp"
#include <expected>
#include <memory>

namespace blast::core {

class ApplicationRunner {
public:
  ApplicationRunner();
  ~ApplicationRunner();

  // Main application entry point
  [[nodiscard]] auto run(int argc, char* argv[]) -> ApplicationResult;

private:
  std::unique_ptr<EnvironmentManager> environment_manager_;
  std::unique_ptr<ConfigurationLoader> config_loader_;
  std::unique_ptr<OutputManager> output_manager_;
  std::unique_ptr<SimulationRunner> simulation_runner_;
  
  // Parse command line arguments
  [[nodiscard]] auto parse_command_line(int argc, char* argv[]) 
    -> std::expected<CommandLineArgs, ApplicationError>;
    
  // Display usage information
  auto display_usage(const std::string& program_name) const -> void;
  
  // Display application header
  auto display_header() const -> void;
  
  // Display performance summary
  auto display_performance_summary(const PerformanceMetrics& metrics) const -> void;
  
  // Display completion message
  auto display_completion_message() const -> void;
  
  // Cleanup resources
  auto cleanup() -> void;
  
  // Convert ApplicationError to ApplicationResult
  auto handle_error(const ApplicationError& error) -> ApplicationResult;
};

} // namespace blast::core