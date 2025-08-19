#pragma once

#include "../core/containers.hpp"
#include "config_types.hpp"
#include <expected>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace blast {

// Forward declarations
namespace thermophysics {
class MixtureInterface;
}

namespace boundary_layer::solver {
class BoundaryLayerSolver;
}

namespace io {

class AbaqueGenerator {
public:
  struct Result {
    std::vector<double> temperatures;
    std::vector<double> catalyticity_values;
    core::Matrix<double> heat_fluxes; // [n_gamma x n_temperatures]
    bool success = false;
  };

  AbaqueGenerator(boundary_layer::solver::BoundaryLayerSolver& solver, thermophysics::MixtureInterface& mixture,
                  const Configuration& config);

  ~AbaqueGenerator();

  // Generate the abaque by varying gamma and Tw
  [[nodiscard]] auto generate() -> Result;

  // Save results to HDF5 file
  [[nodiscard]] auto save_results(const Result& results, const std::filesystem::path& output_path) const -> bool;

private:
  boundary_layer::solver::BoundaryLayerSolver& solver_;
  thermophysics::MixtureInterface& mixture_;
  const Configuration& config_;

  std::filesystem::path gsi_file_path_;
  std::string original_gsi_content_;
  bool gsi_backed_up_ = false;

  // Backup original GSI file content
  [[nodiscard]] auto backup_gsi_file() -> std::expected<void, std::string>;

  // Restore original GSI file
  auto restore_gsi_file() -> void;

  // Update GSI file with new gamma value
  [[nodiscard]] auto update_gsi_catalyticity(double gamma) -> std::expected<void, std::string>;

  // Solve for a single temperature
  [[nodiscard]] auto solve_for_temperature(double wall_temperature) -> double;

  // Build the GSI file path from mixture name
  [[nodiscard]] auto construct_gsi_path() const -> std::filesystem::path;
};

} // namespace io
} // namespace blast