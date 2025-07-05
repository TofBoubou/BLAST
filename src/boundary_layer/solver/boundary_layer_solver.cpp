#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/boundary_layer/equations/continuity.hpp"
#include "blast/boundary_layer/equations/momentum.hpp"
#include "blast/boundary_layer/equations/energy.hpp"
#include "blast/boundary_layer/equations/species.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace blast::boundary_layer::solver {

BoundaryLayerSolver::BoundaryLayerSolver(
    const thermophysics::MixtureInterface& mixture,
    const io::Configuration& config
) : mixture_(mixture), config_(config) {
    
    // Create grid
    if (config.simulation.only_stagnation_point) {
        auto grid_result = grid::BoundaryLayerGrid::create_stagnation_grid(
            config.numerical, config.outer_edge
        );
        if (!grid_result) {
            throw SolverError("Failed to create stagnation grid: {}", 
                            std::source_location::current(), grid_result.error().message());
        }
        grid_ = std::make_unique<grid::BoundaryLayerGrid>(std::move(grid_result.value()));
    } else {
        auto grid_result = grid::BoundaryLayerGrid::create_downstream_grid(
            config.numerical, config.outer_edge, config.output
        );
        if (!grid_result) {
            throw SolverError("Failed to create downstream grid: {}", 
                            std::source_location::current(), grid_result.error().message());
        }
        grid_ = std::make_unique<grid::BoundaryLayerGrid>(std::move(grid_result.value()));
    }
    
    // Create coefficient calculator
    coeff_calculator_ = std::make_unique<coefficients::CoefficientCalculator>(
        mixture_, config_.simulation, config_.numerical
    );
    
    // Create enthalpy-temperature solver
    thermodynamics::EnthalpyTemperatureSolverConfig h2t_config{
        .tolerance = config_.numerical.solvers.h2t_tolerance,
        .max_iterations = config_.numerical.solvers.h2t_max_iterations,
        .max_bracket_expansions = config_.numerical.solvers.h2t_max_bracket_expansions
    };
    h2t_solver_ = std::make_unique<thermodynamics::EnthalpyTemperatureSolver>(
        mixture_, h2t_config
    );
    
    // Create xi derivatives manager
    xi_derivatives_ = std::make_unique<coefficients::XiDerivatives>(
        grid_->n_eta(), mixture_.n_species()
    );
}

auto BoundaryLayerSolver::solve() -> std::expected<SolutionResult, SolverError> {
    
    SolutionResult result;
    result.xi_solved.reserve(grid_->xi_coordinates().size());
    result.stations.reserve(grid_->xi_coordinates().size());
    result.temperature_fields.reserve(grid_->xi_coordinates().size());
    
    const auto xi_stations = grid_->xi_coordinates();
    
    // Solve each xi station
    for (std::size_t station_idx = 0; station_idx < xi_stations.size(); ++station_idx) {
        const int station = static_cast<int>(station_idx);
        const double xi = xi_stations[station_idx];
        
        // Create initial guess for this station
        equations::SolutionState initial_guess;
        if (station == 0 || result.stations.empty()) {
            // Create boundary conditions
            auto bc_result = conditions::create_stagnation_conditions(
                config_.outer_edge, config_.wall_parameters, config_.simulation
            );
            if (!bc_result) {
                return std::unexpected(SolverError(
                    "Failed to create boundary conditions for station {}: {}", 
                    std::source_location::current(), station, bc_result.error().message()
                ));
            }
            
            auto guess_result = create_initial_guess(station, xi, bc_result.value());
            if (!guess_result) {
                return std::unexpected(guess_result.error());
            }
            initial_guess = std::move(guess_result.value());
        } else {
            // Extrapolate from previous station
            initial_guess = extrapolate_from_previous(
                result.stations.back(), 
                result.xi_solved.back(), 
                xi
            );
        }
        
        // Solve this station
        auto station_result = solve_station(station, xi, initial_guess);
        if (!station_result) {
            return std::unexpected(SolverError(
                "Failed to solve station {} (xi={}): {}", 
                std::source_location::current(), station, xi, station_result.error().message()
            ));
        }
        
        // Store results
        result.xi_solved.push_back(xi);
        result.stations.push_back(std::move(station_result.value()));
        result.total_iterations++; // This would be accumulated from all stations
        
        // Update xi derivatives for next station
        if (station_idx + 1 < xi_stations.size()) {
            const auto& solution = result.stations.back();
            xi_derivatives_->update_station(station, xi, solution.F, solution.g, solution.c);
        }
    }
    
    result.converged = true; // All stations solved successfully
    return result;
}

auto BoundaryLayerSolver::solve_station(
    int station,
    double xi,
    const equations::SolutionState& initial_guess
) -> std::expected<equations::SolutionState, SolverError> {
    
    // Get boundary conditions for this station
    auto bc_result = (station == 0) ?
        conditions::create_stagnation_conditions(
            config_.outer_edge, config_.wall_parameters, config_.simulation
        ) :
        conditions::interpolate_boundary_conditions(
            station, xi, grid_->xi_coordinates(),
            config_.outer_edge, config_.wall_parameters, config_.simulation
        );
    
    if (!bc_result) {
        return std::unexpected(SolverError(
            "Failed to get boundary conditions for station {}: {}", 
            std::source_location::current(), station, bc_result.error().message()
        ));
    }
    
    auto bc = bc_result.value();
    auto solution = initial_guess;
    
    // Iterative solution at this station
    auto convergence_result = iterate_station(station, xi, bc, solution);
    if (!convergence_result) {
        return std::unexpected(convergence_result.error());
    }
    
    const auto& conv_info = convergence_result.value();
    if (!conv_info.converged) {
        return std::unexpected(SolverError(
            "Station {} failed to converge after {} iterations (residual={})", 
            std::source_location::current(), station, conv_info.iterations, conv_info.max_residual()
        ));
    }
    
    return solution;
}

auto BoundaryLayerSolver::iterate_station(
    int station,
    double xi,
    const conditions::BoundaryConditions& bc,
    equations::SolutionState& solution
) -> std::expected<ConvergenceInfo, SolverError> {
    
    const auto n_eta = grid_->n_eta();
    const auto n_species = mixture_.n_species();
    std::vector<double> temperature_field(n_eta);
    
    // Initialize temperature field from wall temperature
    std::ranges::fill(temperature_field, bc.Tw());
    
    ConvergenceInfo conv_info;
    
    for (int iter = 0; iter < config_.numerical.max_iterations; ++iter) {
        const auto solution_old = solution;
        
        // 1. Solve continuity equation: dV/dη = -y
        auto V_result = solve_continuity_equation(solution);
        if (!V_result) {
            return std::unexpected(SolverError(
                "Continuity solver failed at station {} iteration {}: {}", 
                std::source_location::current(), station, iter, V_result.error().message()
            ));
        }
        solution.V = std::move(V_result.value());
        
        // 2. Compute eta derivatives
        auto derivatives_result = compute_eta_derivatives(solution);
        if (!derivatives_result) {
            return std::unexpected(derivatives_result.error());
        }
        auto solution_with_derivatives = derivatives_result.value();
        
        // 3. Create coefficient inputs
        coefficients::CoefficientInputs inputs{
            .xi = xi,
            .F = solution.F,
            .c = solution.c,
            .dc_deta = solution_with_derivatives.c, // derivatives stored here temporarily
            .T = temperature_field
        };
        
        // 4. Calculate all coefficients
        auto coeffs_result = coeff_calculator_->calculate(inputs, bc, *xi_derivatives_);
        if (!coeffs_result) {
            return std::unexpected(SolverError(
                "Coefficient calculation failed at station {} iteration {}: {}", 
                std::source_location::current(), station, iter, coeffs_result.error().message()
            ));
        }
        auto coeffs = coeffs_result.value();
        
        // 5. Solve momentum equation
        auto F_result = solve_momentum_equation(solution, coeffs, bc, xi);
        if (!F_result) {
            return std::unexpected(F_result.error());
        }
        auto F_new = F_result.value();
        
        // 6. Solve energy equation  
        auto g_result = solve_energy_equation(solution, inputs, coeffs, bc, station);
        if (!g_result) {
            return std::unexpected(g_result.error());
        }
        auto g_new = g_result.value();
        
        // 7. Update temperature from enthalpy
        auto T_result = update_temperature_field(g_new, solution.c, bc, temperature_field);
        if (!T_result) {
            return std::unexpected(T_result.error());
        }
        temperature_field = T_result.value();
        
        // 8. Update inputs with new temperature
        inputs.T = temperature_field;
        
        // 9. Recalculate coefficients with updated temperature
        auto coeffs_updated_result = coeff_calculator_->calculate(inputs, bc, *xi_derivatives_);
        if (!coeffs_updated_result) {
            return std::unexpected(SolverError("Coefficient recalculation failed"));
        }
        coeffs = coeffs_updated_result.value();
        
        // 10. Solve species equations
        auto c_result = solve_species_equations(solution, inputs, coeffs, bc, station);
        if (!c_result) {
            return std::unexpected(c_result.error());
        }
        auto c_new = c_result.value();
        
        // 11. Apply relaxation
        equations::SolutionState solution_new(n_eta, n_species);
        solution_new.V = solution.V;
        solution_new.F = std::move(F_new);
        solution_new.g = std::move(g_new);
        solution_new.c = std::move(c_new);
        
        solution = apply_relaxation(solution_old, solution_new, config_.numerical.under_relaxation);
        
        // 12. Check convergence
        conv_info = check_convergence(solution_old, solution);
        conv_info.iterations = iter + 1;
        
        if (conv_info.converged) {
            break;
        }
    }
    
    return conv_info;
}

auto BoundaryLayerSolver::solve_continuity_equation(
    const equations::SolutionState& solution
) -> std::expected<std::vector<double>, SolverError> {
    
    const auto n_eta = grid_->n_eta();
    const double d_eta = grid_->d_eta();
    
    // Get current xi value and lambda0 from xi_derivatives
    const double xi = (xi_derivatives_->station() == 0) ? 0.0 : 
                      grid_->xi_coordinates()[xi_derivatives_->station()];
    const double lambda0 = xi_derivatives_->lambda0();
    const auto F_derivatives = xi_derivatives_->F_derivative();
    
    // Compute y field according to the formula:
    // y[i] = -(2ξλ₀ + 1)F[i] - 2ξ∂F/∂ξ
    std::vector<double> y_field(n_eta);
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        y_field[i] = -(2.0 * xi * lambda0 + 1.0) * solution.F[i] - 
                      2.0 * xi * F_derivatives[i];
    }
    
    // Integrate dV/dη = -y to get V
    // Note: The negative sign is handled in the continuity equation solver
    return equations::solve_continuity(y_field, d_eta, 0.0);
}

auto BoundaryLayerSolver::solve_momentum_equation(
    const equations::SolutionState& solution,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    double xi
) -> std::expected<std::vector<double>, SolverError> {
    
    auto result = equations::solve_momentum(
        solution.F, coeffs, bc, *xi_derivatives_, solution.V, xi, grid_->d_eta()
    );
    
    if (!result) {
        return std::unexpected(SolverError(
            "Momentum equation failed: {}", 
            std::source_location::current(), result.error().message()
        ));
    }
    
    return result.value();
}

auto BoundaryLayerSolver::solve_energy_equation(
    const equations::SolutionState& solution,
    const coefficients::CoefficientInputs& inputs,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    int station
) -> std::expected<std::vector<double>, SolverError> {
    
    // Compute dF/deta for energy equation  
    std::vector<double> dF_deta(solution.F.size());
    for (std::size_t i = 1; i < solution.F.size() - 1; ++i) {
        dF_deta[i] = (solution.F[i+1] - solution.F[i-1]) / (2.0 * grid_->d_eta());
    }
    dF_deta[0] = (solution.F[1] - solution.F[0]) / grid_->d_eta();
    dF_deta.back() = (solution.F.back() - solution.F[solution.F.size()-2]) / grid_->d_eta();
    
    auto result = equations::solve_energy(
        solution.g, inputs, coeffs, bc, *xi_derivatives_,
        config_.simulation, solution.F, dF_deta, solution.V, station, grid_->d_eta()
    );
    
    if (!result) {
        return std::unexpected(SolverError(
            "Energy equation failed: {}", 
            std::source_location::current(), result.error().message()
        ));
    }
    
    return result.value();
}

auto BoundaryLayerSolver::solve_species_equations(
    const equations::SolutionState& solution,
    const coefficients::CoefficientInputs& inputs,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    int station
) -> std::expected<core::Matrix<double>, SolverError> {
    
    auto result = equations::solve_species(
        solution.c, inputs, coeffs, bc, *xi_derivatives_,
        mixture_, config_.simulation, solution.F, solution.V, station, grid_->d_eta()
    );
    
    if (!result) {
        return std::unexpected(SolverError(
            "Species equations failed: {}", 
            std::source_location::current(), result.error().message()
        ));
    }
    
    return result.value();
}

auto BoundaryLayerSolver::update_temperature_field(
    std::span<const double> g_field,
    const core::Matrix<double>& composition,
    const conditions::BoundaryConditions& bc,
    std::span<const double> current_temperatures
) -> std::expected<std::vector<double>, SolverError> {
    
    const auto n_eta = g_field.size();
    std::vector<double> enthalpy_field(n_eta);
    
    // Convert g (dimensionless enthalpy) to dimensional enthalpy
    for (std::size_t i = 0; i < n_eta; ++i) {
        enthalpy_field[i] = g_field[i] * bc.he();
    }
    
    auto result = h2t_solver_->solve(enthalpy_field, composition, bc, current_temperatures);
    if (!result) {
        return std::unexpected(SolverError(
            "Temperature solve failed: {}", 
            std::source_location::current(), result.error().message()
        ));
    }
    
    return result.value().temperatures;
}

auto BoundaryLayerSolver::check_convergence(
    const equations::SolutionState& old_solution,
    const equations::SolutionState& new_solution
) const noexcept -> ConvergenceInfo {
    
    ConvergenceInfo info;
    const double tol = config_.numerical.convergence_tolerance;
    
    // Compute L2 norms of residuals
    auto compute_residual = [](const auto& old_field, const auto& new_field) {
        double sum = 0.0;
        for (std::size_t i = 0; i < old_field.size(); ++i) {
            const double diff = new_field[i] - old_field[i];
            sum += diff * diff;
        }
        return std::sqrt(sum / old_field.size());
    };
    
    info.residual_F = compute_residual(old_solution.F, new_solution.F);
    info.residual_g = compute_residual(old_solution.g, new_solution.g);
    
    // Species residual (matrix)
    double c_sum = 0.0;
    std::size_t c_count = 0;
    for (std::size_t i = 0; i < old_solution.c.rows(); ++i) {
        for (std::size_t j = 0; j < old_solution.c.cols(); ++j) {
            const double diff = new_solution.c(i, j) - old_solution.c(i, j);
            c_sum += diff * diff;
            ++c_count;
        }
    }
    info.residual_c = std::sqrt(c_sum / c_count);
    
    info.converged = (info.residual_F < tol) && (info.residual_g < tol) && (info.residual_c < tol);
    
    return info;
}

auto BoundaryLayerSolver::create_initial_guess(
    int station,
    double xi,
    const conditions::BoundaryConditions& bc
) const -> std::expected<equations::SolutionState, SolverError> {
    
    const auto n_eta = grid_->n_eta();
    const auto n_species = mixture_.n_species();
    const double eta_max = grid_->eta_max();
    
    equations::SolutionState guess(n_eta, n_species);
    
    // Get wall composition from boundary conditions
    std::vector<double> c_wall(n_species);
    for (std::size_t j = 0; j < n_species; ++j) {
        c_wall[j] = (j < bc.c_e().size()) ? bc.c_e()[j] : 0.0;
    }
    
    // Compute wall enthalpy
    auto h_wall_result = mixture_.mixture_enthalpy(c_wall, bc.Tw(), bc.P_e());
    if (!h_wall_result) {
        return std::unexpected(SolverError(
            "Failed to compute wall enthalpy for initial guess: {}", 
            std::source_location::current(), h_wall_result.error().message()
        ));
    }
    const double h_wall = h_wall_result.value();
    const double g_wall = h_wall / bc.he();  // Dimensionless wall enthalpy
    const double g_edge = 1.0;               // Dimensionless edge enthalpy (h_e/h_e = 1)
    
    // Initialize with reasonable profiles
    for (std::size_t i = 0; i < n_eta; ++i) {
        const double eta = static_cast<double>(i) * eta_max / (n_eta - 1);
        const double eta_norm = static_cast<double>(i) / (n_eta - 1);
        
        // Boundary layer-like profiles
        guess.V[i] = 0.0; // Will be computed from continuity
        guess.F[i] = eta_norm * (2.0 - eta_norm); // Smooth profile
        guess.g[i] = g_wall + eta_norm * (g_edge - g_wall); // Linear dimensionless enthalpy
        
        // Species from boundary conditions
        for (std::size_t j = 0; j < n_species; ++j) {
            guess.c(j, i) = c_wall[j];  // Use consistent composition
        }
    }
    
    return guess;
}

auto BoundaryLayerSolver::extrapolate_from_previous(
    const equations::SolutionState& previous_solution,
    double xi_prev,
    double xi_current
) const -> equations::SolutionState {
    
    // Simple extrapolation - could be improved with higher-order methods
    const double factor = (xi_current > xi_prev) ? 1.1 : 0.9; // Modest extrapolation
    
    auto extrapolated = previous_solution;
    
    // Extrapolate F and g fields
    for (std::size_t i = 0; i < extrapolated.F.size(); ++i) {
        extrapolated.F[i] *= factor;
        extrapolated.g[i] = std::max(0.1, extrapolated.g[i] * factor); // Keep positive
    }
    
    // Species concentrations stay the same (conservative)
    
    return extrapolated;
}

auto BoundaryLayerSolver::apply_relaxation(
    const equations::SolutionState& old_solution,
    const equations::SolutionState& new_solution,
    double relaxation_factor
) const -> equations::SolutionState {
    
    auto relaxed = new_solution; // Start with new solution
    
    // Apply relaxation: relaxed = (1-α)*old + α*new
    const double alpha = relaxation_factor;
    const double one_minus_alpha = 1.0 - alpha;
    
    for (std::size_t i = 0; i < relaxed.F.size(); ++i) {
        relaxed.F[i] = one_minus_alpha * old_solution.F[i] + alpha * new_solution.F[i];
        relaxed.g[i] = one_minus_alpha * old_solution.g[i] + alpha * new_solution.g[i];
    }
    
    for (std::size_t i = 0; i < relaxed.c.rows(); ++i) {
        for (std::size_t j = 0; j < relaxed.c.cols(); ++j) {
            relaxed.c(i, j) = one_minus_alpha * old_solution.c(i, j) + alpha * new_solution.c(i, j);
        }
    }
    
    return relaxed;
}

auto BoundaryLayerSolver::compute_eta_derivatives(
    const equations::SolutionState& solution
) const -> std::expected<equations::SolutionState, SolverError> {
    
    const auto n_eta = grid_->n_eta();
    const auto n_species = mixture_.n_species();
    const double d_eta = grid_->d_eta();
    
    equations::SolutionState derivatives(n_eta, n_species);
    
    // Simple finite difference derivatives
    auto compute_derivative = [d_eta](const std::vector<double>& field) {
        std::vector<double> deriv(field.size());
        
        // Forward difference at start
        deriv[0] = (field[1] - field[0]) / d_eta;
        
        // Central difference in interior
        for (std::size_t i = 1; i < field.size() - 1; ++i) {
            deriv[i] = (field[i+1] - field[i-1]) / (2.0 * d_eta);
        }
        
        // Backward difference at end
        deriv.back() = (field.back() - field[field.size()-2]) / d_eta;
        
        return deriv;
    };
    
    derivatives.F = compute_derivative(solution.F);
    derivatives.g = compute_derivative(solution.g);
    derivatives.V = compute_derivative(solution.V);
    
    // Species derivatives
    for (std::size_t j = 0; j < n_species; ++j) {
        std::vector<double> c_row(n_eta);
        for (std::size_t i = 0; i < n_eta; ++i) {
            c_row[i] = solution.c(j, i);
        }
        
        auto dc_deta = compute_derivative(c_row);
        for (std::size_t i = 0; i < n_eta; ++i) {
            derivatives.c(j, i) = dc_deta[i];
        }
    }
    
    return derivatives;
}

} // namespace blast::boundary_layer::solver