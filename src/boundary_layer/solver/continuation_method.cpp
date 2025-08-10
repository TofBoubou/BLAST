#include "blast/boundary_layer/solver/continuation_method.hpp"
#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace blast::boundary_layer::solver {

auto ContinuationMethod::solve_with_continuation(
    BoundaryLayerSolver& solver, int station, double xi, const io::Configuration& target_config,
    const equations::SolutionState& initial_guess) -> std::expected<ContinuationResult, SolverError> {

  double lambda = 0.0;
  double lambda_step = LAMBDA_STEP_INITIAL;
  equations::SolutionState current_solution = initial_guess;
  const auto original_config = solver.config();

  for (int step = 0; step < MAX_STEPS; ++step) {
    if (lambda + lambda_step > 1.0) {
      lambda_step = 1.0 - lambda;
    }

    double lambda_try = lambda + lambda_step;
    auto interp_config = interpolate_config(target_config, lambda_try);

    solver.set_config(interp_config);
    solver.in_continuation_ = true;
    auto result = solver.solve_station(station, xi, current_solution);
    solver.in_continuation_ = false;
    solver.set_config(original_config);

    if (result) {
      current_solution = result.value();
      lambda = lambda_try;

      if (std::abs(lambda - 1.0) < 1e-10) {
        return ContinuationResult{.solution = current_solution, .success = true, .final_lambda = lambda};
      }

      // Increase step after success
      lambda_step = std::min(lambda_step * STEP_INCREASE_FACTOR, LAMBDA_STEP_MAX);

    } else {
      // Check if the failure was due to NaN detection
      const std::string error_msg = result.error().message();
      bool is_nan_error = error_msg.find("NaN detected") != std::string::npos;
      
      if (is_nan_error) {
        std::cout << "[CONTINUATION] NaN detected at lambda=" << std::scientific << std::setprecision(3) << lambda_try 
                  << ", reducing step from " << std::scientific << std::setprecision(3) << lambda_step 
                  << " to " << std::scientific << std::setprecision(3) << lambda_step * STEP_DECREASE_FACTOR << std::endl;
      } else {
        std::cout << "[CONTINUATION] Step failed at lambda=" << std::scientific << std::setprecision(3) << lambda_try 
                  << ", reducing step from " << std::scientific << std::setprecision(3) << lambda_step 
                  << " to " << std::scientific << std::setprecision(3) << lambda_step * STEP_DECREASE_FACTOR << std::endl;
      }
      
      // Decrease step after failure (including NaN errors)
      lambda_step *= STEP_DECREASE_FACTOR;

      if (lambda_step < LAMBDA_STEP_MIN) {
        std::cout << "[CONTINUATION] Step size too small (" << std::scientific << std::setprecision(3) << lambda_step 
                  << " < " << std::scientific << std::setprecision(3) << LAMBDA_STEP_MIN 
                  << "), giving up at lambda=" << std::scientific << std::setprecision(3) << lambda << std::endl;
        return ContinuationResult{.solution = current_solution, .success = false, .final_lambda = lambda};
      }
    }
  }

  return std::unexpected(SolverError("Continuation max steps exceeded", std::source_location::current()));
}

auto ContinuationMethod::interpolate_config(const io::Configuration& target, double lambda) const -> io::Configuration {

  io::Configuration config = target;

  // Interpolate wall temperature
  if (!config.wall_parameters.wall_temperatures.empty()) {
    double Twall_target = target.wall_parameters.wall_temperatures[0];
    config.wall_parameters.wall_temperatures[0] = TWALL_STABLE + lambda * (Twall_target - TWALL_STABLE);
  }

  // Interpolate edge conditions
  if (!config.outer_edge.edge_points.empty()) {
    auto& edge = config.outer_edge.edge_points[0];
    const auto& target_edge = target.outer_edge.edge_points[0];

    edge.temperature = TEDGE_STABLE + lambda * (target_edge.temperature - TEDGE_STABLE);
    edge.pressure = PRESSURE_STABLE + lambda * (target_edge.pressure - PRESSURE_STABLE);
  }

  return config;
}

} // namespace blast::boundary_layer::solver