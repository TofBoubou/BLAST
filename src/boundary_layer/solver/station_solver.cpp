#include "blast/boundary_layer/solver/station_solver.hpp"
#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/boundary_layer/solver/convergence_manager.hpp"
#include "blast/boundary_layer/solver/continuation_method.hpp"
#include <format>
#include <cmath>

namespace blast::boundary_layer::solver {

StationSolver::StationSolver(BoundaryLayerSolver& solver,
                            const thermophysics::MixtureInterface& mixture,
                            const io::Configuration& config) noexcept
    : solver_(solver), mixture_(mixture), config_(config) {}

auto StationSolver::solve_station(int station, double xi, const equations::SolutionState& initial_guess)
    -> std::expected<equations::SolutionState, SolverError> {

    // Validate input parameters
    if (auto validation_result = validate_station_inputs(station, xi, initial_guess); !validation_result) {
        return std::unexpected(validation_result.error());
    }

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
            continuation && !in_continuation_ && conv_info.max_residual() < 1e10) {
            
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

auto StationSolver::validate_station_inputs(int station, double xi, const equations::SolutionState& initial_guess) const
    -> std::expected<void, SolverError> {

    // Validate station number
    if (station < 0) {
        return std::unexpected(
            InitializationError(std::format("Invalid station number: {} (must be non-negative)", station)));
    }

    // Validate xi coordinate
    if (!std::isfinite(xi) || xi < 0.0) {
        return std::unexpected(
            InitializationError(std::format("Invalid xi coordinate: {} (must be finite and non-negative)", xi)));
    }

    // Get grid dimensions
    const auto& grid = solver_.get_grid();
    const auto expected_n_eta = grid.n_eta();
    const auto expected_n_species = mixture_.n_species();

    // Validate initial guess field dimensions
    if (initial_guess.F.size() != expected_n_eta || 
        initial_guess.T.size() != expected_n_eta ||
        initial_guess.g.size() != expected_n_eta) {
        return std::unexpected(InitializationError(
            std::format("Initial guess field dimensions mismatch: expected {} eta points", expected_n_eta)));
    }

    // Validate species matrix dimensions
    if (initial_guess.c.rows() != expected_n_species || initial_guess.c.cols() != expected_n_eta) {
        return std::unexpected(InitializationError(std::format(
            "Initial guess species matrix dimensions mismatch: expected {}x{}", expected_n_species, expected_n_eta)));
    }

    // Check consistency between station and xi coordinates for downstream stations
    if (station > 0) {
        const auto& xi_coords = grid.xi_coordinates();
        if (station >= static_cast<int>(xi_coords.size())) {
            return std::unexpected(InitializationError(
                std::format("Station {} exceeds available xi coordinates (max: {})", station, xi_coords.size() - 1)));
        }

        // Allow some tolerance for floating point comparison
        const double expected_xi = xi_coords[station];
        if (std::abs(xi - expected_xi) > 1e-10) {
            return std::unexpected(InitializationError(
                std::format("Xi coordinate mismatch for station {}: provided {}, expected {}", station, xi, expected_xi)));
        }
    }

    return {};
}

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
    return create_initial_guess(station, xi, bc_stable.value(), T_edge_stable);
}

auto StationSolver::should_attempt_continuation(const ConvergenceError& convergence_error) const noexcept -> bool {
    // Simple heuristic: attempt continuation if the error doesn't indicate a fundamental problem
    const std::string& message = convergence_error.message();
    
    // Don't attempt continuation for clearly bad conditions
    if (message.find("NaN") != std::string::npos ||
        message.find("infinite") != std::string::npos ||
        message.find("dimension") != std::string::npos) {
        return false;
    }
    
    return true;
}

auto StationSolver::create_initial_guess(int station, double xi, const conditions::BoundaryConditions& bc,
                                        double T_edge) const
    -> std::expected<equations::SolutionState, SolverError> {

    const auto& grid = solver_.get_grid();
    const auto n_eta = grid.n_eta();
    const auto n_species = mixture_.n_species();
    const double eta_max = grid.eta_max();

    equations::SolutionState guess(n_eta, n_species);

    // Initialize V field to zero
    std::fill(guess.V.begin(), guess.V.end(), 0.0);

    // Handle single species case (simplified initialization)
    if (n_species == 1) {
        guess.c.eigen().setOnes();

        // Compute wall equilibrium enthalpy
        std::array<double, 1> c_wall{{1.0}};
        auto h_wall_eq_result = mixture_.mixture_enthalpy(c_wall, bc.Tw(), bc.P_e());
        if (!h_wall_eq_result) {
            return std::unexpected(NumericError(
                std::format("Failed to compute wall equilibrium enthalpy: {}", h_wall_eq_result.error().message())));
        }
        double g_wall = h_wall_eq_result.value() / bc.he();

        // Create simple profiles
        for (std::size_t i = 0; i < n_eta; ++i) {
            const double eta = static_cast<double>(i) * eta_max / (n_eta - 1);
            const double eta_norm = static_cast<double>(i) / (n_eta - 1);

            guess.F[i] = 1.0 - 1.034 * std::exp(-5.628 * eta / eta_max);
            guess.g[i] = g_wall + eta_norm * (1.0 - g_wall);
            guess.T[i] = bc.Tw() + guess.F[i] * (T_edge - bc.Tw());
        }

        return guess;
    }

    // Multi-species case - compute equilibrium composition at wall
    auto equilibrium_result = mixture_.equilibrium_composition(bc.Tw(), bc.P_e());
    if (!equilibrium_result) {
        return std::unexpected(NumericError(std::format("Failed to compute equilibrium composition at wall conditions: {}",
                                                        equilibrium_result.error().message())));
    }
    auto c_wall_equilibrium = equilibrium_result.value();

    // Extract edge composition
    const auto& c_edge = bc.c_e();
    if (c_edge.size() != n_species) {
        return std::unexpected(InitializationError(
            std::format("Edge composition size mismatch: expected {}, got {}", n_species, c_edge.size())));
    }

    // Compute wall equilibrium enthalpy
    auto h_wall_eq_result = mixture_.mixture_enthalpy(c_wall_equilibrium, bc.Tw(), bc.P_e());
    if (!h_wall_eq_result) {
        return std::unexpected(NumericError(
            std::format("Failed to compute wall equilibrium enthalpy: {}", h_wall_eq_result.error().message())));
    }
    double g_wall = h_wall_eq_result.value() / bc.he();

    // Create profiles for multi-species case
    for (std::size_t i = 0; i < n_eta; ++i) {
        const double eta = static_cast<double>(i) * eta_max / (n_eta - 1);
        const double eta_norm = static_cast<double>(i) / (n_eta - 1);

        // Momentum profile
        guess.F[i] = 1.0 - 1.034 * std::exp(-5.628 * eta / eta_max);
        
        // Energy profile
        guess.g[i] = g_wall + eta_norm * (1.0 - g_wall);
        
        // Temperature profile
        guess.T[i] = bc.Tw() + guess.F[i] * (T_edge - bc.Tw());

        // Species profiles - linear interpolation between wall equilibrium and edge
        for (std::size_t s = 0; s < n_species; ++s) {
            guess.c(s, i) = c_wall_equilibrium[s] + eta_norm * (c_edge[s] - c_wall_equilibrium[s]);
        }
    }

    return guess;
}

} // namespace blast::boundary_layer::solver