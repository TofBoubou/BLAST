#pragma once
#include "../../io/config_types.hpp"
#include "../equations/equation_types.hpp"
#include "solver_errors.hpp"
#include "../../core/containers.hpp"
#include <expected>
#include <vector>
#include <deque>

namespace blast::boundary_layer::solver {

class BoundaryLayerSolver;

struct ContinuationResult {
  equations::SolutionState solution;
  bool success;
  double final_lambda;
};

struct ContinuationHistoryPoint {
  equations::SolutionState solution;
  double lambda;
};

class ContinuationMethod {
public:
  ContinuationMethod() = default;

  [[nodiscard]] auto solve_with_continuation(
      BoundaryLayerSolver& solver, int station, double xi, const io::Configuration& target_config,
      const equations::SolutionState& initial_guess) -> std::expected<ContinuationResult, SolverError>;

private:
  // Mutable state for chemical mode switching
  mutable int consecutive_failures_ = 0;
  mutable int consecutive_successes_in_equilibrium_ = 0;
  mutable bool using_equilibrium_mode_ = false;
  mutable io::SimulationConfig::ChemicalMode original_chemical_mode_;
  
  // Track if we had success in non-equilibrium since last switch to equilibrium
  mutable bool had_nonequilibrium_success_since_switch_ = true;
  
  // Track last successful non-equilibrium solution above threshold
  mutable std::optional<ContinuationHistoryPoint> last_nonequilibrium_success_;
  
  // Mutable state for polynomial predictor
  mutable std::deque<ContinuationHistoryPoint> history_;
  mutable int step_reductions_for_current_step_ = 0;
  mutable bool predictor_enabled_ = true;

  [[nodiscard]] auto interpolate_config(const io::Configuration& target, double lambda) const -> io::Configuration;
  [[nodiscard]] auto create_equilibrium_config(const io::Configuration& config) const -> io::Configuration;
  
  // Polynomial predictor methods
  [[nodiscard]] auto predict_solution(double target_lambda) const
      -> std::expected<equations::SolutionState, SolverError>;
  void add_to_history(const equations::SolutionState& solution, double lambda) const;
  void reset_predictor_state() const;
};

} // namespace blast::boundary_layer::solver
