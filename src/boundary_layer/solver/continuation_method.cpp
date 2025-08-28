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

  // Initialize switching state
  consecutive_failures_ = 0;
  consecutive_successes_in_equilibrium_ = 0;
  using_equilibrium_mode_ = false;
  original_chemical_mode_ = target_config.simulation.chemical_mode;
  
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

    // Apply chemical mode switching if needed
    if (using_equilibrium_mode_) {
      interp_config = create_equilibrium_config(interp_config);
    }

    solver.set_config(interp_config);
    solver.in_continuation_ = true;
    auto result = solver.solve_station(station, xi, current_solution);
    solver.in_continuation_ = false;
    solver.set_config(original_config);

    if (result) {
      current_solution = result.value();
      lambda = lambda_try;

      // Reset failure counter and handle success in equilibrium mode
      consecutive_failures_ = 0;
      if (using_equilibrium_mode_) {
        consecutive_successes_in_equilibrium_++;
        std::cout << "[CHEMICAL SWITCHING] Success in equilibrium mode (" << consecutive_successes_in_equilibrium_ 
                  << "/" << SUCCESS_THRESHOLD << ")" << std::endl;
        
        // Switch back to original mode after SUCCESS_THRESHOLD successes
        if (consecutive_successes_in_equilibrium_ >= SUCCESS_THRESHOLD) {
          using_equilibrium_mode_ = false;
          consecutive_successes_in_equilibrium_ = 0;
          std::cout << "[CHEMICAL SWITCHING] Switching back to " 
                    << (original_chemical_mode_ == io::SimulationConfig::ChemicalMode::NonEquilibrium ? "non_equilibrium" : "original") 
                    << " mode" << std::endl;
        }
      }

      if (std::abs(lambda - 1.0) < 1e-10) {
        return ContinuationResult{.solution = current_solution, .success = true, .final_lambda = lambda};
      }

      // Increase step after success
      lambda_step = std::min(lambda_step * STEP_INCREASE_FACTOR, LAMBDA_STEP_MAX);

    } else {
      // Handle failure - increment failure counter and check for chemical mode switching
      consecutive_failures_++;
      consecutive_successes_in_equilibrium_ = 0;  // Reset success counter on any failure
      
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

      // Check for chemical mode switching (only for NonEquilibrium -> Equilibrium)
      if (!using_equilibrium_mode_ && 
          original_chemical_mode_ == io::SimulationConfig::ChemicalMode::NonEquilibrium &&
          consecutive_failures_ >= FAILURE_THRESHOLD) {
        using_equilibrium_mode_ = true;
        consecutive_failures_ = 0;
        consecutive_successes_in_equilibrium_ = 0;
        std::cout << "[CHEMICAL SWITCHING] Switching to equilibrium mode after " << FAILURE_THRESHOLD 
                  << " consecutive failures" << std::endl;
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

auto ContinuationMethod::create_equilibrium_config(const io::Configuration& config) const -> io::Configuration {
  io::Configuration equilibrium_config = config;
  equilibrium_config.simulation.chemical_mode = io::SimulationConfig::ChemicalMode::Equilibrium;
  return equilibrium_config;
}

} // namespace blast::boundary_layer::solver