#include "blast/core/output_manager.hpp"
#include "blast/core/constants.hpp"
#include "blast/io/output/hdf5_writer.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>

namespace blast::core {

auto OutputManager::initialize_output_system(const io::Configuration& config, 
                                             const std::string& case_name) 
  -> std::expected<void, ApplicationError> {
    
  std::cout << "\nInitializing output system..." << std::endl;
  
  if (auto hdf5_init = initialize_hdf5(); !hdf5_init) {
    return std::unexpected(hdf5_init.error());
  }
  
  auto output_config = create_output_config(config);
  output_writer_ = std::make_unique<io::output::OutputWriter>(output_config);
  
  if (auto validation = output_writer_->validate_config(); !validation) {
    return std::unexpected(ApplicationError{
      "Output configuration error: " + validation.error().message(),
      constants::indexing::second
    });
  }
  
  std::cout << constants::string_processing::colors::green << "✓ Output system configured" 
            << constants::string_processing::colors::reset << std::endl;
  
  return {};
}

auto OutputManager::write_simulation_results(
  const boundary_layer::solver::SolutionResult& solution,
  const io::Configuration& config,
  const thermophysics::MixtureInterface& mixture,
  const std::string& case_name,
  PerformanceMetrics& metrics) 
  -> std::expected<std::vector<std::filesystem::path>, ApplicationError> {
    
  std::cout << "\n=== WRITING OUTPUT FILES ===" << std::endl;
  
  auto progress_callback = create_progress_callback();
  
  auto output_start = std::chrono::high_resolution_clock::now();
  auto output_result = output_writer_->write_solution(solution, config, mixture, case_name, progress_callback);
  auto output_end = std::chrono::high_resolution_clock::now();
  
  metrics.output_time = std::chrono::duration_cast<std::chrono::milliseconds>(output_end - output_start);
  
  std::cout << std::endl; // New line after progress
  
  if (!output_result) {
    return std::unexpected(ApplicationError{
      "Failed to write output: " + output_result.error().message(),
      constants::indexing::second
    });
  }
  
  auto output_files = std::move(output_result.value());
  metrics.output_files = output_files;
  
  std::cout << constants::string_processing::colors::green << "✓ Output written successfully!" 
            << constants::string_processing::colors::reset << std::endl;
  std::cout << "  Output time: " << metrics.output_time.count() << " ms" << std::endl;
  std::cout << "\nGenerated files:" << std::endl;
  
  for (const auto& file_path : output_files) {
    auto file_size = std::filesystem::file_size(file_path);
    std::cout << "  " << file_path.filename().string() 
              << " (" << std::setprecision(constants::string_processing::float_precision_2) << std::fixed
              << (file_size / constants::io::bytes_to_kb) << " KB)" << std::endl;
  }
  
  return output_files;
}

auto OutputManager::display_planned_outputs(const std::string& case_name) const -> void {
  auto output_info = output_writer_->get_output_info(case_name);
  std::cout << "\nPlanned output files:" << std::endl;
  
  for (const auto& [format, path] : output_info) {
    std::string format_name;
    switch (format) {
    case io::output::OutputFormat::HDF5:
      format_name = "HDF5";
      break;
    default:
      format_name = "Unknown";
      break;
    }
    std::cout << "  " << format_name << ": " << path.filename().string() << " (timestamp will be set at write time)" << std::endl;
  }
}

auto OutputManager::initialize_hdf5() -> std::expected<void, ApplicationError> {
  if (auto hdf5_init = io::output::hdf5::initialize(); !hdf5_init) {
    std::cerr << "Warning: Failed to initialize HDF5: " << hdf5_init.error().message() << std::endl;
  } else {
    if (auto version = io::output::hdf5::check_version()) {
      std::cout << constants::string_processing::colors::green << "✓ HDF5 library version: " 
                << constants::string_processing::colors::reset << version.value() << std::endl;
    }
  }
  return {};
}

auto OutputManager::create_output_config(const io::Configuration& config) -> io::output::OutputConfig {
  io::output::OutputConfig output_config;
  output_config.base_directory = config.output.output_directory;
  output_config.primary_format = io::output::OutputFormat::HDF5;
  output_config.additional_formats = {};
  output_config.save_metadata = true;
  output_config.save_derivatives = true;
  output_config.compress_data = true;
  output_config.include_timestamp = true;
  
  // Select variables to save based on simulation type
  output_config.variables.flow_variables = true;
  output_config.variables.species_concentrations = true;
  
  return output_config;
}

auto OutputManager::create_progress_callback() -> ProgressCallback {
  return [](double progress, const std::string& stage) {
    std::cout << "\r" << std::setw(constants::string_processing::separator_width + 24) << std::left
              << ("  " + stage + " [" + std::to_string(static_cast<int>(progress * constants::conversion::to_percentage)) + "%]") 
              << std::flush;
  };
}

} // namespace blast::core