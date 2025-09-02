#include "blast/boundary_layer/solver/continuation_method.hpp"
#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/boundary_layer/solver/continuation_scope.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <Eigen/Dense>

namespace blast::boundary_layer::solver {

auto ContinuationMethod::solve_with_continuation(
    BoundaryLayerSolver& solver, int station, double xi, const io::Configuration& target_config,
    const equations::SolutionState& initial_guess) -> std::expected<ContinuationResult, SolverError> {

  // Initialize switching state
  consecutive_failures_ = 0;
  consecutive_successes_in_equilibrium_ = 0;
  using_equilibrium_mode_ = false;
  original_chemical_mode_ = target_config.simulation.chemical_mode;
  
  // Initialize predictor state
  history_.clear();
  step_reductions_for_current_step_ = 0;
  predictor_enabled_ = target_config.continuation.use_linear_predictor;
  const bool verbose = target_config.verbose || (std::getenv("BLAST_VERBOSE") != nullptr);
  
  double lambda = 0.0;
  double lambda_step = LAMBDA_STEP_INITIAL;
  equations::SolutionState current_solution = initial_guess;
  const auto original_config = solver.config();
  
  // Initialize history with initial solution at lambda=0
  add_to_history(current_solution, lambda);

  for (int step = 0; step < MAX_STEPS; ++step) {
    // Display progress percentage (optional)
    if (verbose) {
      double progress = lambda * 100;
      std::cout << "\r[CONTINUATION] Progress: " << std::fixed << std::setprecision(2) << progress << "%" << std::flush;
    }
    
    if (lambda + lambda_step > 1.0) {
      lambda_step = 1.0 - lambda;
    }

    double lambda_try = lambda + lambda_step;
    
    // Use linear predictor if enabled and we have enough history
    equations::SolutionState predicted_solution = current_solution;
    if (predictor_enabled_ && history_.size() >= 2) {
      try {
        predicted_solution = predict_solution(lambda_try);
        
        // Using linear prediction
      } catch (const std::exception& e) {
        if (verbose) {
          std::cout << "\n[PREDICTOR] Prediction failed: " << e.what() 
                    << ", falling back to previous solution" << std::endl;
        }
        predicted_solution = current_solution;
      }
    }
    
    auto interp_config = interpolate_config(target_config, lambda_try);

    // Apply chemical mode switching if needed
    if (using_equilibrium_mode_) {
      interp_config = create_equilibrium_config(interp_config);
    }

    solver.set_config(interp_config);
    {
      ContinuationScope scope(solver);
      auto result = solver.solve_station(station, xi, predicted_solution);
      if (result) {
        current_solution = result.value();
        lambda = lambda_try;
        add_to_history(current_solution, lambda);
        step_reductions_for_current_step_ = 0;
        consecutive_failures_ = 0;
        if (using_equilibrium_mode_) {
          consecutive_successes_in_equilibrium_++;
          if (verbose) {
            std::cout << "\n[CHEMICAL SWITCHING] Success in equilibrium mode (" << consecutive_successes_in_equilibrium_
                      << "/" << SUCCESS_THRESHOLD << ")" << std::endl;
          }
          if (consecutive_successes_in_equilibrium_ >= SUCCESS_THRESHOLD) {
            using_equilibrium_mode_ = false;
            consecutive_successes_in_equilibrium_ = 0;
            if (verbose) {
              std::cout << "\n[CHEMICAL SWITCHING] Switching back to "
                        << (original_chemical_mode_ == io::SimulationConfig::ChemicalMode::NonEquilibrium ?
                                "non_equilibrium" :
                                "original")
                        << " mode" << std::endl;
            }
          }
        }
        if (std::abs(lambda - 1.0) < 1e-10) {
          if (verbose) {
            std::cout << "\n[CONTINUATION] Completed 100%" << std::endl;
          }
          solver.set_config(original_config);
          return ContinuationResult{.solution = current_solution, .success = true, .final_lambda = lambda};
        }
        lambda_step = std::min(lambda_step * STEP_INCREASE_FACTOR, LAMBDA_STEP_MAX);
        if (!predictor_enabled_ && target_config.continuation.use_linear_predictor) {
          predictor_enabled_ = true;
          if (verbose) {
            std::cout << "\n[PREDICTOR] Re-enabled after successful step" << std::endl;
          }
        }
      } else {
        consecutive_failures_++;
        consecutive_successes_in_equilibrium_ = 0;
        const std::string error_msg = result.error().message();
        bool is_nan_error = error_msg.find("NaN detected") != std::string::npos;
        if (verbose) {
          if (is_nan_error) {
            std::cout << "\n[CONTINUATION] NaN detected at lambda=" << std::scientific << std::setprecision(3) << lambda_try
                      << ", reducing step from " << std::scientific << std::setprecision(3) << lambda_step
                      << " to " << std::scientific << std::setprecision(3) << lambda_step * STEP_DECREASE_FACTOR
                      << std::endl;
          } else {
            std::cout << "\n[CONTINUATION] Step failed at lambda=" << std::scientific << std::setprecision(3) << lambda_try
                      << ", reducing step from " << std::scientific << std::setprecision(3) << lambda_step
                      << " to " << std::scientific << std::setprecision(3) << lambda_step * STEP_DECREASE_FACTOR
                      << std::endl;
          }
        }
        if (!using_equilibrium_mode_ &&
            original_chemical_mode_ == io::SimulationConfig::ChemicalMode::NonEquilibrium &&
            consecutive_failures_ >= FAILURE_THRESHOLD) {
          using_equilibrium_mode_ = true;
          consecutive_failures_ = 0;
          consecutive_successes_in_equilibrium_ = 0;
          if (verbose) {
            std::cout << "\n[CHEMICAL SWITCHING] Switching to equilibrium mode after " << FAILURE_THRESHOLD
                      << " consecutive failures" << std::endl;
          }
        }
        lambda_step *= STEP_DECREASE_FACTOR;
        step_reductions_for_current_step_++;
        if (step_reductions_for_current_step_ >= MAX_STEP_REDUCTIONS && predictor_enabled_) {
          predictor_enabled_ = false;
          if (verbose) {
            std::cout << "\n[PREDICTOR] Disabled after " << MAX_STEP_REDUCTIONS << " step reductions" << std::endl;
          }
        }
        if (lambda_step < LAMBDA_STEP_MIN) {
          if (verbose) {
            std::cout << "\n[CONTINUATION] Step size too small (" << std::scientific << std::setprecision(3) << lambda_step
                      << " < " << std::scientific << std::setprecision(3) << LAMBDA_STEP_MIN
                      << "), giving up at lambda=" << std::scientific << std::setprecision(3) << lambda
                      << std::endl;
          }
          solver.set_config(original_config);
          return ContinuationResult{.solution = current_solution, .success = false, .final_lambda = lambda};
        }
      }
    }
    solver.set_config(original_config);
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

auto ContinuationMethod::predict_solution(double target_lambda) const -> equations::SolutionState {
  if (history_.size() < 2) {
    throw std::runtime_error("Not enough history points for prediction");
  }
  
  // Get the last 2 points from history for linear extrapolation
  const auto& p1 = history_[history_.size() - 2]; // previous
  const auto& p2 = history_[history_.size() - 1]; // current
  
  // Create predicted solution by linear extrapolation
  equations::SolutionState predicted = p2.solution; // Start with last solution as template
  
  // Check for degenerate case
  if (std::abs(p2.lambda - p1.lambda) < 1e-12) {
    return predicted; // Return current solution if lambdas are too close
  }
  
  // Predict V (continuity) using linear extrapolation
  const int n_V_points = p2.solution.V.size();
  for (int i = 0; i < n_V_points; ++i) {
    double slope = (p2.solution.V[i] - p1.solution.V[i]) / (p2.lambda - p1.lambda);
    predicted.V[i] = p2.solution.V[i] + slope * (target_lambda - p2.lambda);
  }
  
  // Predict F (momentum) using linear extrapolation
  const int n_F_points = p2.solution.F.size();
  for (int i = 0; i < n_F_points; ++i) {
    double slope = (p2.solution.F[i] - p1.solution.F[i]) / (p2.lambda - p1.lambda);
    predicted.F[i] = p2.solution.F[i] + slope * (target_lambda - p2.lambda);
  }
  
  // Predict g (energy) using linear extrapolation
  const int n_g_points = p2.solution.g.size();
  for (int i = 0; i < n_g_points; ++i) {
    double slope = (p2.solution.g[i] - p1.solution.g[i]) / (p2.lambda - p1.lambda);
    predicted.g[i] = p2.solution.g[i] + slope * (target_lambda - p2.lambda);
  }
  
  // Predict T (temperature) using linear extrapolation
  const int n_T_points = p2.solution.T.size();
  for (int i = 0; i < n_T_points; ++i) {
    double slope = (p2.solution.T[i] - p1.solution.T[i]) / (p2.lambda - p1.lambda);
    predicted.T[i] = p2.solution.T[i] + slope * (target_lambda - p2.lambda);
  }
  
  // Predict c (species mass fractions) using linear extrapolation
  const int n_species = p2.solution.c.rows();
  const int n_species_points = p2.solution.c.cols();
  for (int i = 0; i < n_species; ++i) {
    for (int j = 0; j < n_species_points; ++j) {
      double slope = (p2.solution.c(i, j) - p1.solution.c(i, j)) / (p2.lambda - p1.lambda);
      predicted.c(i, j) = p2.solution.c(i, j) + slope * (target_lambda - p2.lambda);
    }
  }
  
  return predicted;
}

void ContinuationMethod::add_to_history(const equations::SolutionState& solution, double lambda) const {
  history_.push_back({solution, lambda});
  
  // Keep only the last 2 points for linear extrapolation
  while (history_.size() > 2) {
    history_.pop_front();
  }
}

void ContinuationMethod::reset_predictor_state() const {
  history_.clear();
  step_reductions_for_current_step_ = 0;
  predictor_enabled_ = true;
}

} // namespace blast::boundary_layer::solver
