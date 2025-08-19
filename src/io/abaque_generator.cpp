#include "blast/io/abaque_generator.hpp"
#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/thermophysics/mixture_interface.hpp"
#include <H5Apublic.h>
#include <H5Dpublic.h>
#include <H5Fpublic.h>
#include <H5Gpublic.h>
#include <H5Ppublic.h>
#include <H5Spublic.h>
#include <cstdlib>
#include <fstream>
#include <format>
#include <iostream>
#include <regex>
#include <sstream>

namespace blast::io {

AbaqueGenerator::AbaqueGenerator(boundary_layer::solver::BoundaryLayerSolver& solver,
                                 thermophysics::MixtureInterface& mixture, const Configuration& config)
    : solver_(solver), mixture_(mixture), config_(config) {

  gsi_file_path_ = construct_gsi_path();

  // Backup GSI file on construction
  auto backup_result = backup_gsi_file();
  if (!backup_result) {
    std::cerr << "Warning: Failed to backup GSI file: " << backup_result.error() << std::endl;
  }
}

AbaqueGenerator::~AbaqueGenerator() {
  // Always restore original GSI on destruction
  if (gsi_backed_up_) {
    restore_gsi_file();
  }
}

auto AbaqueGenerator::generate() -> Result {
  Result result;

  if (!config_.abaque.enabled) {
    return result;
  }

  const auto& abaque_config = config_.abaque;
  const int n_temps = abaque_config.temperature_points;
  const int n_gammas = abaque_config.catalyticity_values.size();

  // Setup temperature vector
  result.temperatures.resize(n_temps);
  const double dT = (abaque_config.temperature_max - abaque_config.temperature_min) / (n_temps - 1);
  for (int i = 0; i < n_temps; ++i) {
    result.temperatures[i] = abaque_config.temperature_min + i * dT;
  }

  // Copy catalyticity values
  result.catalyticity_values = abaque_config.catalyticity_values;

  // Initialize heat flux matrix
  result.heat_fluxes = core::Matrix<double>(n_gammas, n_temps);

  // Main computation loop
  std::cout << "\n=== ABAQUE GENERATION ===" << std::endl;

  for (int gamma_idx = 0; gamma_idx < n_gammas; ++gamma_idx) {
    const double gamma = abaque_config.catalyticity_values[gamma_idx];

    std::cout << "\nProcessing gamma = " << gamma << " (" << (gamma_idx + 1) << "/" << n_gammas << ")" << std::endl;

    // Update GSI file with new gamma
    auto update_result = update_gsi_catalyticity(gamma);
    if (!update_result) {
      std::cerr << "Failed to update GSI catalyticity for gamma = " << gamma << ": " << update_result.error() << std::endl;
      return result;
    }

    // Reload mixture with updated GSI
    auto reload_result = mixture_.reload_gsi();
    if (!reload_result) {
      std::cerr << "Failed to reload mixture after GSI update: " << reload_result.error() << std::endl;
      return result;
    }

    // Compute heat flux for each temperature
    for (int temp_idx = 0; temp_idx < n_temps; ++temp_idx) {
      const double Tw = result.temperatures[temp_idx];

      if (temp_idx % 10 == 0) {
        std::cout << "  T_w = " << Tw << " K (" << (temp_idx + 1) << "/" << n_temps << ")" << std::endl;
      }

      double q_wall = solve_for_temperature(Tw);
      result.heat_fluxes(gamma_idx, temp_idx) = q_wall;
    }
  }

  std::cout << "\n=== ABAQUE GENERATION COMPLETE ===" << std::endl;

  result.success = true;
  return result;
}

auto AbaqueGenerator::solve_for_temperature(double wall_temperature) -> double {
  // Update wall temperature in configuration
  solver_.set_wall_temperature(wall_temperature);

  // Solve the boundary layer
  auto solution_result = solver_.solve();

  if (!solution_result) {
    std::cerr << "Solver failed for T_w = " << wall_temperature << std::endl;
    return 0.0;
  }

  // Extract wall heat flux from the first station (stagnation point)
  if (solution_result.value().heat_flux_data.empty()) {
    std::cerr << "No heat flux data available" << std::endl;
    return 0.0;
  }

  return solution_result.value().heat_flux_data[0].q_wall_total_dim;
}

auto AbaqueGenerator::backup_gsi_file() -> std::expected<void, std::string> {
  if (!std::filesystem::exists(gsi_file_path_)) {
    return std::unexpected(std::format("GSI file not found: {}", gsi_file_path_.string()));
  }

  std::ifstream file(gsi_file_path_);
  if (!file.is_open()) {
    return std::unexpected(std::format("Failed to open GSI file: {}", gsi_file_path_.string()));
  }

  std::stringstream buffer;
  buffer << file.rdbuf();
  original_gsi_content_ = buffer.str();
  file.close();

  gsi_backed_up_ = true;
  return {};
}

auto AbaqueGenerator::restore_gsi_file() -> void {
  if (!gsi_backed_up_ || original_gsi_content_.empty()) {
    return;
  }

  std::ofstream file(gsi_file_path_);
  if (file.is_open()) {
    file << original_gsi_content_;
    file.close();
  }
}

auto AbaqueGenerator::update_gsi_catalyticity(double gamma) -> std::expected<void, std::string> {
  if (!gsi_backed_up_) {
    return std::unexpected("GSI file not backed up");
  }

  std::string content = original_gsi_content_;

  // Regex to match gamma_const elements
  std::regex gamma_pattern(R"(<gamma_const>\s*([^<]+)\s*</gamma_const>)");
  std::string result_content;

  auto current_pos = content.cbegin();
  auto content_end = content.cend();

  std::smatch match;
  while (std::regex_search(current_pos, content_end, match, gamma_pattern)) {
    // Append everything before the match
    result_content.append(current_pos, current_pos + match.position());

    // Parse the species list and replace gamma values
    std::string species_list = match[1].str();
    std::stringstream new_gamma_content;
    new_gamma_content << "<gamma_const> ";

    // Match individual species:value pairs
    std::regex species_pattern(R"((\w+):[\d.]+(?:\s+|$))");
    std::sregex_iterator species_begin(species_list.begin(), species_list.end(), species_pattern);
    std::sregex_iterator species_end;

    bool first = true;
    for (auto it = species_begin; it != species_end; ++it) {
      if (!first)
        new_gamma_content << " ";

      std::string species_name = (*it)[1].str();
      new_gamma_content << species_name << ":" << gamma;
      first = false;
    }

    new_gamma_content << " </gamma_const>";
    result_content.append(new_gamma_content.str());

    // Move past this match
    current_pos = current_pos + match.position() + match.length();
  }

  // Append any remaining content
  result_content.append(current_pos, content_end);

  // Write updated content to file
  std::ofstream file(gsi_file_path_);
  if (!file.is_open()) {
    return std::unexpected(std::format("Failed to open GSI file for writing: {}", gsi_file_path_.string()));
  }

  file << result_content;
  file.close();

  return {};
}

auto AbaqueGenerator::construct_gsi_path() const -> std::filesystem::path {
  // Get MPP_DATA_DIRECTORY from environment
  const char* mpp_data = std::getenv("MPP_DATA_DIRECTORY");

  if (!mpp_data) {
    // Try to find it relative to executable
    std::filesystem::path exe_path = std::filesystem::read_symlink("/proc/self/exe");
    std::filesystem::path base_path = exe_path.parent_path();
    std::filesystem::path mpp_path = base_path / "libs" / "mutationpp" / "data";

    if (std::filesystem::exists(mpp_path)) {
      return mpp_path / "gsi" / (config_.mixture.name + "_cata.xml");
    }

    // Fallback to current directory
    return std::filesystem::current_path() / "libs" / "mutationpp" / "data" / "gsi" /
           (config_.mixture.name + "_cata.xml");
  }

  // Use MPP_DATA_DIRECTORY
  std::filesystem::path gsi_path = std::filesystem::path(mpp_data) / "gsi";
  std::string gsi_filename = config_.mixture.name + "_cata.xml";

  return gsi_path / gsi_filename;
}

auto AbaqueGenerator::save_results(const Result& results, const std::filesystem::path& output_path) const -> bool {

  if (!results.success) {
    return false;
  }

  // Create output directory if needed
  std::filesystem::create_directories(output_path.parent_path());

  // Create HDF5 file
  hid_t file_id = H5Fcreate(output_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file_id < 0) {
    std::cerr << "Failed to create HDF5 file: " << output_path << std::endl;
    return false;
  }

  bool success = true;

  // Save temperatures
  {
    hsize_t dims[1] = {results.temperatures.size()};
    hid_t space_id = H5Screate_simple(1, dims, nullptr);
    hid_t dset_id =
        H5Dcreate2(file_id, "temperatures", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (dset_id >= 0) {
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, results.temperatures.data());
      H5Dclose(dset_id);
    } else {
      success = false;
    }
    H5Sclose(space_id);
  }

  // Save catalyticity values
  {
    hsize_t dims[1] = {results.catalyticity_values.size()};
    hid_t space_id = H5Screate_simple(1, dims, nullptr);
    hid_t dset_id =
        H5Dcreate2(file_id, "catalyticity_values", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (dset_id >= 0) {
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, results.catalyticity_values.data());
      H5Dclose(dset_id);
    } else {
      success = false;
    }
    H5Sclose(space_id);
  }

  // Save heat flux matrix
  {
    hsize_t dims[2] = {static_cast<hsize_t>(results.heat_fluxes.rows()),
                       static_cast<hsize_t>(results.heat_fluxes.cols())};
    hid_t space_id = H5Screate_simple(2, dims, nullptr);
    hid_t dset_id =
        H5Dcreate2(file_id, "heat_fluxes", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (dset_id >= 0) {
      // Matrix data is stored row-major
      std::vector<double> flat_data;
      flat_data.reserve(results.heat_fluxes.rows() * results.heat_fluxes.cols());
      for (int i = 0; i < results.heat_fluxes.rows(); ++i) {
        for (int j = 0; j < results.heat_fluxes.cols(); ++j) {
          flat_data.push_back(results.heat_fluxes(i, j));
        }
      }

      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat_data.data());
      H5Dclose(dset_id);
    } else {
      success = false;
    }
    H5Sclose(space_id);
  }

  // Add metadata as attributes
  {
    hid_t root_group = H5Gopen2(file_id, "/", H5P_DEFAULT);

    // Write temperature range
    hsize_t dims[1] = {2};
    hid_t space_id = H5Screate_simple(1, dims, nullptr);
    hid_t attr_id = H5Acreate2(root_group, "temperature_range", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT);
    if (attr_id >= 0) {
      double range[2] = {config_.abaque.temperature_min, config_.abaque.temperature_max};
      H5Awrite(attr_id, H5T_NATIVE_DOUBLE, range);
      H5Aclose(attr_id);
    }
    H5Sclose(space_id);

    // Write number of points
    hid_t scalar_space = H5Screate(H5S_SCALAR);
    attr_id = H5Acreate2(root_group, "n_temperatures", H5T_NATIVE_INT, scalar_space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr_id >= 0) {
      int n_temps = config_.abaque.temperature_points;
      H5Awrite(attr_id, H5T_NATIVE_INT, &n_temps);
      H5Aclose(attr_id);
    }

    attr_id = H5Acreate2(root_group, "n_catalyticity_values", H5T_NATIVE_INT, scalar_space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr_id >= 0) {
      int n_gammas = results.catalyticity_values.size();
      H5Awrite(attr_id, H5T_NATIVE_INT, &n_gammas);
      H5Aclose(attr_id);
    }
    H5Sclose(scalar_space);

    H5Gclose(root_group);
  }

  H5Fclose(file_id);

  if (success) {
    std::cout << "Abaque results saved to: " << output_path << std::endl;
  }

  return success;
}

} // namespace blast::io