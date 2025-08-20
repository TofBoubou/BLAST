#include "blast/boundary_layer/solver/continuation_method.hpp"
#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace blast::boundary_layer::solver {

auto ContinuationMethod::solve_with_continuation(
    BoundaryLayerSolver& solver, int station, double xi, const io::Configuration& target_config,
    const equations::SolutionState& initial_guess) -> std::expected<ContinuationResult, SolverError> {

  std::cout << "[DEBUG] $$$$$ ENTERING solve_with_continuation $$$$$" << std::endl;
  
  double lambda = 0.0;
  double lambda_step = LAMBDA_STEP_INITIAL;
  equations::SolutionState current_solution = initial_guess;
  const auto original_config = solver.config();

  std::cout << "[DEBUG] $$$$$ STARTING CONTINUATION LOOP - MAX_STEPS=" << MAX_STEPS << " $$$$$" << std::endl;
  for (int step = 0; step < MAX_STEPS; ++step) {
    std::cout << "[DEBUG] $$$$$ CONTINUATION STEP " << step << " - lambda=" << lambda << " lambda_step=" << lambda_step << " $$$$$" << std::endl;
    if (lambda + lambda_step > 1.0) {
      lambda_step = 1.0 - lambda;
      std::cout << "[DEBUG] $$$$$ Adjusting lambda_step to reach 1.0: lambda_step=" << lambda_step << " $$$$$" << std::endl;
    }

    double lambda_try = lambda + lambda_step;
    std::cout << "[DEBUG] $$$$$ Trying lambda_try=" << lambda_try << " (lambda=" << lambda << " + lambda_step=" << lambda_step << ") $$$$$" << std::endl;
    auto interp_config = interpolate_config(target_config, lambda_try);

    solver.set_config(interp_config);
    solver.in_continuation_ = true;
    std::cout << "[DEBUG] ##### CONTINUATION METHOD CALLING solve_station #####" << std::endl;
    auto result = solver.solve_station(station, xi, current_solution);
    solver.in_continuation_ = false;
    solver.set_config(original_config);

    if (result) {
      std::cout << "[DEBUG] $$$$$ STEP " << step << " SUCCEEDED - lambda_try=" << lambda_try << " $$$$$" << std::endl;
      current_solution = result.value();
      std::cout << "[DEBUG] $$$$$ BEFORE UPDATE: lambda=" << lambda << " $$$$$" << std::endl;
      lambda = lambda_try;
      std::cout << "[DEBUG] $$$$$ AFTER UPDATE: lambda=" << lambda << " $$$$$" << std::endl;

      if (std::abs(lambda - 1.0) < 1e-10) {
        std::cout << "[DEBUG] $$$$$ CONTINUATION COMPLETE - RETURNING SUCCESS $$$$$" << std::endl;
        return ContinuationResult{.solution = current_solution, .success = true, .final_lambda = lambda};
      }
      
      std::cout << "[DEBUG] $$$$$ CONTINUING - lambda=" << lambda << " < 1.0 $$$$$" << std::endl;

      // Increase step after success
      std::cout << "[DEBUG] $$$$$ BEFORE STEP INCREASE: lambda_step=" << lambda_step << " $$$$$" << std::endl;
      lambda_step = std::min(lambda_step * STEP_INCREASE_FACTOR, LAMBDA_STEP_MAX);
      std::cout << "[DEBUG] $$$$$ AFTER STEP INCREASE: lambda_step=" << lambda_step << " (FACTOR=" << STEP_INCREASE_FACTOR << ", MAX=" << LAMBDA_STEP_MAX << ") $$$$$" << std::endl;

    } else {
      // Check if the failure was due to NaN detection
      const std::string error_msg = result.error().message();
      bool is_nan_error = error_msg.find("NaN detected") != std::string::npos;

      if (is_nan_error) {
        std::cout << "[CONTINUATION] NaN detected at lambda=" << std::scientific << std::setprecision(3) << lambda_try
                  << ", reducing step from " << std::scientific << std::setprecision(3) << lambda_step << " to "
                  << std::scientific << std::setprecision(3) << lambda_step * STEP_DECREASE_FACTOR << std::endl;
      } else {
        std::cout << "[CONTINUATION] Step failed at lambda=" << std::scientific << std::setprecision(3) << lambda_try
                  << ", reducing step from " << std::scientific << std::setprecision(3) << lambda_step << " to "
                  << std::scientific << std::setprecision(3) << lambda_step * STEP_DECREASE_FACTOR << std::endl;
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

  return std::unexpected(ConvergenceError("Continuation max steps exceeded"));
}

auto ContinuationMethod::interpolate_config(const io::Configuration& target, double lambda) const -> io::Configuration {

  io::Configuration config = target;

  // Interpolate wall temperature
  if (!config.wall_parameters.wall_temperatures.empty()) {
    double Twall_target = target.wall_parameters.wall_temperatures[0];
    double Twall_stable = target.continuation.wall_temperature_stable;
    config.wall_parameters.wall_temperatures[0] = Twall_stable + lambda * (Twall_target - Twall_stable);
  }

  // Interpolate edge conditions
  if (!config.outer_edge.edge_points.empty()) {
    auto& edge = config.outer_edge.edge_points[0];
    const auto& target_edge = target.outer_edge.edge_points[0];

    double Tedge_stable = target.continuation.edge_temperature_stable;
    double pressure_stable = target.continuation.pressure_stable;
    edge.temperature = Tedge_stable + lambda * (target_edge.temperature - Tedge_stable);
    edge.pressure = pressure_stable + lambda * (target_edge.pressure - pressure_stable);
  }

  return config;
}

} // namespace blast::boundary_layer::solver