#pragma once
#include "../../io/config_types.hpp"
#include "../equations/equation_types.hpp"
#include "solver_errors.hpp"
#include <expected>
#include <vector>

namespace blast::boundary_layer::solver {

class BoundaryLayerSolver;

struct ContinuationResult {
  equations::SolutionState solution;
  bool success;
  double final_lambda;
};

class ContinuationMethod {
public:
  ContinuationMethod() = default;

  [[nodiscard]] auto solve_with_continuation(
      BoundaryLayerSolver& solver, int station, double xi, const io::Configuration& target_config,
      const equations::SolutionState& initial_guess) -> std::expected<ContinuationResult, SolverError>;

private:
  // Continuation parameters
  static constexpr double LAMBDA_STEP_INITIAL = 0.01;
  static constexpr double LAMBDA_STEP_MIN = 0.0001;
  static constexpr double LAMBDA_STEP_MAX = 0.2;
  static constexpr double STEP_INCREASE_FACTOR = 1.2;
  static constexpr double STEP_DECREASE_FACTOR = 0.8;
  static constexpr int MAX_STEPS = 10000;

  // Chemical mode switching parameters
  static constexpr int FAILURE_THRESHOLD = 4;  // Switch after 4 consecutive failures
  static constexpr int SUCCESS_THRESHOLD = 2;  // Switch back after 2 consecutive successes

  // Mutable state for chemical mode switching
  mutable int consecutive_failures_ = 0;
  mutable int consecutive_successes_in_equilibrium_ = 0;
  mutable bool using_equilibrium_mode_ = false;
  mutable io::SimulationConfig::ChemicalMode original_chemical_mode_;

  [[nodiscard]] auto interpolate_config(const io::Configuration& target, double lambda) const -> io::Configuration;
  [[nodiscard]] auto create_equilibrium_config(const io::Configuration& config) const -> io::Configuration;
};

} // namespace blast::boundary_layer::solver