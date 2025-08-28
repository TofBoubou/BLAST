#include "blast/io/abacus_generator.hpp"
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

AbacusGenerator::AbacusGenerator(boundary_layer::solver::BoundaryLayerSolver& solver,
                                 thermophysics::MixtureInterface& mixture, const Configuration& config)
    : solver_(solver), mixture_(mixture), config_(config), gsi_manager_(config) {

  // Backup GSI file on construction
  auto backup_result = gsi_manager_.backup_gsi_file();
  if (!backup_result) {
    std::cerr << "Warning: Failed to backup GSI file: " << backup_result.error() << std::endl;
  }
}

AbacusGenerator::~AbacusGenerator() {
  // GSI cleanup is handled by GsiManager destructor
}

auto AbacusGenerator::generate() -> Result {
  Result result;

  if (!config_.abacus.enabled) {
    return result;
  }

  const auto& abacus_config = config_.abacus;
  const int n_temps = abacus_config.temperature_points;
  const int n_gammas = abacus_config.catalyticity_values.size();

  // Setup temperature vector
  result.temperatures.resize(n_temps);
  const double dT = (abacus_config.temperature_max - abacus_config.temperature_min) / (n_temps - 1);
  for (int i = 0; i < n_temps; ++i) {
    result.temperatures[i] = abacus_config.temperature_min + i * dT;
  }

  // Copy catalyticity values
  result.catalyticity_values = abacus_config.catalyticity_values;

  // Initialize heat flux matrix
  result.heat_fluxes = core::Matrix<double>(n_gammas, n_temps);

  // Main computation loop
  std::cout << "\n=== ABACUS GENERATION ===" << std::endl;

  for (int gamma_idx = 0; gamma_idx < n_gammas; ++gamma_idx) {
    const double gamma = abacus_config.catalyticity_values[gamma_idx];

    std::cout << "\nProcessing gamma = " << gamma << " (" << (gamma_idx + 1) << "/" << n_gammas << ")" << std::endl;

    // Update GSI file with new gamma
    auto update_result = gsi_manager_.update_gsi_catalyticity(gamma);
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

  std::cout << "\n=== ABACUS GENERATION COMPLETE ===" << std::endl;

  result.success = true;
  return result;
}

auto AbacusGenerator::solve_for_temperature(double wall_temperature) -> double {
  // Store original configuration
  auto original_config = solver_.get_config();
  
  // Create a copy of the configuration and update wall temperature
  auto modified_config = original_config;
  if (!modified_config.wall_parameters.wall_temperatures.empty()) {
    modified_config.wall_parameters.wall_temperatures[0] = wall_temperature;
  } else {
    modified_config.wall_parameters.wall_temperatures.push_back(wall_temperature);
  }
  
  // Set the modified configuration
  solver_.set_config(modified_config);

  // Solve the boundary layer
  auto solution_result = solver_.solve();

  // Restore original configuration
  solver_.set_config(original_config);

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


auto AbacusGenerator::save_results(const Result& results, const std::filesystem::path& output_path) const -> bool {

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
      double range[2] = {config_.abacus.temperature_min, config_.abacus.temperature_max};
      H5Awrite(attr_id, H5T_NATIVE_DOUBLE, range);
      H5Aclose(attr_id);
    }
    H5Sclose(space_id);

    // Write number of points
    hid_t scalar_space = H5Screate(H5S_SCALAR);
    attr_id = H5Acreate2(root_group, "n_temperatures", H5T_NATIVE_INT, scalar_space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr_id >= 0) {
      int n_temps = config_.abacus.temperature_points;
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
    std::cout << "Abacus results saved to: " << output_path << std::endl;
  }

  return success;
}

} // namespace blast::io