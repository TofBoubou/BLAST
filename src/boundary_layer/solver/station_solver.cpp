#include "blast/boundary_layer/solver/station_solver.hpp"
#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/boundary_layer/solver/convergence_manager.hpp"
#include "blast/boundary_layer/solver/continuation_method.hpp"
#include "blast/boundary_layer/solver/expected_utils.hpp"
#include "blast/boundary_layer/solver/input_validator.hpp"
#include <format>
#include <cmath>

namespace blast::boundary_layer::solver {

StationSolver::StationSolver(BoundaryLayerSolver& solver,
                            const thermophysics::MixtureInterface& mixture,
                            const io::Configuration& config) noexcept
    : solver_(solver), mixture_(mixture), config_(config) {}

auto StationSolver::solve_station(int station, double xi, const equations::SolutionState& initial_guess)
    -> std::expected<equations::SolutionState, SolverError> {

    // Create validator and validate inputs using unified patterns
    InputValidator validator(solver_.get_grid(), mixture_);
    
    BLAST_TRY_VOID(validator.validate_station_parameters(station, xi));
    BLAST_TRY_VOID(validator.validate_solution_state(initial_guess));
    BLAST_TRY_VOID(validator.validate_xi_coordinate(station, xi));

    // Get boundary conditions
    auto bc_result = get_boundary_conditions(station, xi);
    if (!bc_result) {
        return std::unexpected(BoundaryConditionError(
            std::format("Failed to get boundary conditions for station {}: {}", station, bc_result.error().message())));
    }

    auto bc = bc_result.value();
    auto solution = initial_guess;

    // Get convergence manager from solver (we'll need to add this accessor)
    auto& convergence_manager = solver_.get_convergence_manager();
    
    // Initialize relaxation controller for this station type
    convergence_manager.initialize_relaxation_for_station(station);
    convergence_manager.set_continuation_mode(in_continuation_);

    // Attempt direct solution
    auto convergence_result = convergence_manager.iterate_station_adaptive(station, xi, bc, solution);

    if (!convergence_result) {
        // Try continuation if available and not already in continuation
        if (auto* continuation = solver_.get_continuation(); 
            continuation && !in_continuation_ && should_attempt_continuation(ConvergenceError(convergence_result.error().message()))) {
            
            auto stable_guess = compute_stable_guess(station, xi);
            if (stable_guess) {
                auto cont_result = continuation->solve_with_continuation(solver_, station, xi, 
                                                                       solver_.get_original_config(), 
                                                                       stable_guess.value());

                if (cont_result && cont_result.value().success) {
                    return cont_result.value().solution;
                }
            }
        }

        return std::unexpected(convergence_result.error());
    }

    const auto& conv_info = convergence_result.value();

    if (!conv_info.converged) {
        // Handle non-convergence based on context
        if (in_continuation_) {
            // Return specific error for continuation to handle
            if (std::isnan(conv_info.residual_F) || std::isnan(conv_info.residual_g) || std::isnan(conv_info.residual_c)) {
                return std::unexpected(ConvergenceError(
                    std::format("NaN detected during continuation at station {} iteration {}", station, conv_info.iterations)));
            } else {
                return std::unexpected(ConvergenceError(
                    std::format("Station {} failed to converge during continuation after {} iterations (residual={})", 
                               station, conv_info.iterations, conv_info.max_residual())));
            }
        }

        // Try continuation for non-continuation failures
        if (auto* continuation = solver_.get_continuation(); 
            continuation && !in_continuation_ && conv_info.max_residual() < config_.numerical.residual_guard) {
            
            auto stable_guess = compute_stable_guess(station, xi);
            if (stable_guess) {
                auto cont_result = continuation->solve_with_continuation(solver_, station, xi, 
                                                                       solver_.get_original_config(), 
                                                                       stable_guess.value());

                if (cont_result && cont_result.value().success) {
                    return cont_result.value().solution;
                }
            }
        }

        return std::unexpected(
            ConvergenceError(std::format("Station {} failed to converge after {} iterations (residual={})", 
                                        station, conv_info.iterations, conv_info.max_residual())));
    }

    return solution;
}

// Note: validate_station_inputs method moved to InputValidator

auto StationSolver::get_boundary_conditions(int station, double xi) const
    -> std::expected<conditions::BoundaryConditions, SolverError> {

    if (station == 0) {
        auto result = conditions::create_stagnation_conditions(config_.outer_edge, config_.wall_parameters, 
                                                               config_.simulation, mixture_);
        if (!result) {
            return std::unexpected(BoundaryConditionError(result.error().message()));
        }
        return result.value();
    } else {
        const auto& grid = solver_.get_grid();
        auto result = conditions::interpolate_boundary_conditions(station, xi, grid.xi_coordinates(), 
                                                                 config_.outer_edge, config_.wall_parameters, 
                                                                 config_.simulation, mixture_);
        if (!result) {
            return std::unexpected(BoundaryConditionError(result.error().message()));
        }
        return result.value();
    }
}

auto StationSolver::compute_stable_guess(int station, double xi) const
    -> std::expected<equations::SolutionState, SolverError> {

    // Create stable configuration
    io::Configuration stable_config = solver_.get_original_config();
    
    if (!stable_config.wall_parameters.wall_temperatures.empty()) {
        stable_config.wall_parameters.wall_temperatures[0] = stable_config_.wall_temperature_stable;
    } else {
        // Ensure a valid wall temperature exists for stable boundary conditions
        stable_config.wall_parameters.wall_temperatures.push_back(stable_config_.wall_temperature_stable);
    }
    
    if (!stable_config.outer_edge.edge_points.empty()) {
        stable_config.outer_edge.edge_points[0].temperature = stable_config_.edge_temperature_stable;
        stable_config.outer_edge.edge_points[0].pressure = stable_config_.pressure_stable;
    }

    // Get stable boundary conditions
    auto bc_stable = conditions::create_stagnation_conditions(stable_config.outer_edge, 
                                                             stable_config.wall_parameters,
                                                             stable_config.simulation, mixture_);
    if (!bc_stable) {
        return std::unexpected(BoundaryConditionError(
            std::format("Failed to create stable boundary conditions: {}", bc_stable.error().message())));
    }

    double T_edge_stable = stable_config.outer_edge.edge_points[0].temperature;
    // Use unified initial guess factory
    InitialGuessFactory guess_factory(solver_.get_grid(), mixture_);
    return guess_factory.create_initial_guess(station, xi, bc_stable.value(), T_edge_stable);
}

auto StationSolver::should_attempt_continuation(const ConvergenceError& convergence_error) const noexcept -> bool {
    (void)convergence_error; // Unused; decision based on configured policy
    using CAP = io::NumericalConfig::ContinuationAttemptPolicy;
    switch (config_.numerical.continuation_attempt_policy) {
      case CAP::Always:
        return true;
      case CAP::OnlyIfResidualBelowGuard:
        // When residual is unknown here, be conservative: do not attempt at this stage
        return false;
      case CAP::Never:
      default:
        return false;
    }
}

// Note: create_initial_guess method moved to InitialGuessFactory

} // namespace blast::boundary_layer::solver
