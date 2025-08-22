#pragma once

#include "../core/containers.hpp"
#include "config_types.hpp"
#include "gsi_manager.hpp"
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

class AbacusGenerator {
public:
  struct Result {
    std::vector<double> temperatures;
    std::vector<double> catalyticity_values;
    core::Matrix<double> heat_fluxes; // [n_gamma x n_temperatures]
    bool success = false;
  };

  AbacusGenerator(boundary_layer::solver::BoundaryLayerSolver& solver, thermophysics::MixtureInterface& mixture,
                  const Configuration& config);

  ~AbacusGenerator();

  // Generate the abacus by varying gamma and Tw
  [[nodiscard]] auto generate() -> Result;

  // Save results to HDF5 file
  [[nodiscard]] auto save_results(const Result& results, const std::filesystem::path& output_path) const -> bool;

private:
  boundary_layer::solver::BoundaryLayerSolver& solver_;
  thermophysics::MixtureInterface& mixture_;
  const Configuration& config_;
  GsiManager gsi_manager_;

  // Solve for a single temperature
  [[nodiscard]] auto solve_for_temperature(double wall_temperature) -> double;
};

} // namespace io
} // namespace blast