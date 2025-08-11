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
  static constexpr double STEP_INCREASE_FACTOR = 2;
  static constexpr double STEP_DECREASE_FACTOR = 0.8;
  static constexpr int MAX_STEPS = 10000;

  [[nodiscard]] auto interpolate_config(const io::Configuration& target, double lambda) const -> io::Configuration;
};

} // namespace blast::boundary_layer::solver