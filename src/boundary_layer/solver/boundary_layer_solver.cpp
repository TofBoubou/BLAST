#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/boundary_layer/equations/continuity.hpp"
#include "blast/boundary_layer/equations/momentum.hpp"
#include "blast/boundary_layer/equations/energy.hpp"
#include "blast/boundary_layer/equations/species.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>
#include <iomanip>

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
        .max_iterations = config_.numerical.solvers.h2t_max_iterations
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
    
    // Variables to store previous station results for xi derivatives
    double prev_xi = 0.0;
    std::vector<double> prev_F, prev_g;
    core::Matrix<double> prev_c;
    
    // Solve each xi station
    for (std::size_t station_idx = 0; station_idx < xi_stations.size(); ++station_idx) {
        const int station = static_cast<int>(station_idx);
        const double xi = xi_stations[station_idx];
        
        // CRITICAL: Update xi derivatives BEFORE solving (except for station 0)
        if (station_idx > 0) {
            xi_derivatives_->update_station(station - 1, xi, prev_F, prev_g, prev_c);
        }
        
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
        
        // Store current results for next iteration
        prev_xi = xi;
        prev_F = result.stations.back().F;
        prev_g = result.stations.back().g;
        prev_c = result.stations.back().c;
        
        result.total_iterations++;
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
    
    // Create a mutable copy of boundary conditions for dynamic updates
    auto bc_dynamic = bc;
    
    // Initialize temperature field from wall temperature
    std::ranges::fill(temperature_field, bc_dynamic.Tw());
    
    ConvergenceInfo conv_info;
    
    for (int iter = 0; iter < config_.numerical.max_iterations; ++iter) {
        std::cout << "=== ITERATION " << iter << " at station " << station << " ===" << std::endl;
        const auto solution_old = solution;
        
        // 1. Solve continuity equation: dV/dη = -solve_continuity
        auto V_result = solve_continuity_equation(solution, xi);
        if (!V_result) {
            return std::unexpected(SolverError(
                "Continuity solver failed at station {} iteration {}: {}", 
                std::source_location::current(), station, iter, V_result.error().message()
            ));
        }
        solution.V = std::move(V_result.value());
        // std::cout << solution.V[5] << std::endl;
        
        // 2. Compute eta derivatives
        auto derivatives_result = compute_eta_derivatives(solution);
        if (!derivatives_result) {
            return std::unexpected(derivatives_result.error());
        }
        auto solution_with_derivatives = derivatives_result.value();
        
        // 2b. Compute concentration derivatives (first and second)
        auto conc_derivatives_result = compute_concentration_derivatives(solution);
        if (!conc_derivatives_result) {
            return std::unexpected(conc_derivatives_result.error());
        }
        auto concentration_derivatives = conc_derivatives_result.value();
        
        // 3. Create coefficient inputs
        coefficients::CoefficientInputs inputs{
            .xi = xi,
            .F = solution.F,
            .c = solution.c,
            .dc_deta = concentration_derivatives.dc_deta,
            .dc_deta2 = concentration_derivatives.dc_deta2,
            .T = temperature_field
        };
        
        // 4. Calculate all coefficients
        auto coeffs_result = coeff_calculator_->calculate(inputs, bc_dynamic, *xi_derivatives_);
        if (!coeffs_result) {
            return std::unexpected(SolverError(
                "Coefficient calculation failed at station {} iteration {}: {}", 
                std::source_location::current(), station, iter, coeffs_result.error().message()
            ));
        }
        auto coeffs = coeffs_result.value();
        
        // 5. Solve momentum equation
        auto F_result = solve_momentum_equation(solution, coeffs, bc_dynamic, xi);
        if (!F_result) {
            return std::unexpected(F_result.error());
        }
        auto F_new = F_result.value();
        
        // 5.5. CRITICAL: Force F[edge] = 1.0 immediately after momentum solve
        if (!F_new.empty()) {
            F_new.back() = 1.0;  // Force dimensionless velocity = 1 at edge
        }
        
        // 6. Solve energy equation  
        auto g_result = solve_energy_equation(solution, inputs, coeffs, bc_dynamic, station);
        if (!g_result) {
            return std::unexpected(g_result.error());
        }
        auto g_new = g_result.value();
        
        // 6.5. CRITICAL: Force boundary conditions BEFORE temperature update
        if (!g_new.empty()) {
            g_new.back() = 1.0;  // Force dimensionless enthalpy = 1 at edge
        }
        
        // 7. Update temperature from enthalpy
        auto T_result = update_temperature_field(g_new, solution.c, bc_dynamic, temperature_field);
        if (!T_result) {
            return std::unexpected(T_result.error());
        }
        temperature_field = T_result.value();
        
        // 8. Update inputs with new temperature
        inputs.T = temperature_field;
        
        // 8.5. Update edge properties dynamically for thermodynamic consistency
        update_edge_properties(bc_dynamic, inputs, solution.c);
        
        // 9. Recalculate coefficients with updated temperature
        auto coeffs_updated_result = coeff_calculator_->calculate(inputs, bc_dynamic, *xi_derivatives_);
        if (!coeffs_updated_result) {
            return std::unexpected(SolverError("Coefficient recalculation failed"));
        }
        coeffs = coeffs_updated_result.value();
        
        // 10. Solve species equations
        auto c_result = solve_species_equations(solution, inputs, coeffs, bc_dynamic, station);
        if (!c_result) {
            return std::unexpected(c_result.error());
        }
        auto c_new = c_result.value();
        
        // 11. Build new solution state
        equations::SolutionState solution_new(n_eta, n_species);
        solution_new.V = solution.V;
        solution_new.F = std::move(F_new);
        solution_new.g = std::move(g_new);
        solution_new.c = std::move(c_new);
        
        // 11.5. CRITICAL: Enforce boundary conditions BEFORE relaxation
        // This ensures the forced values are not corrupted by relaxation
        enforce_edge_boundary_conditions(solution_new, bc_dynamic);
        
        // 12. Apply relaxation AFTER enforcing boundary conditions
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
    const equations::SolutionState& solution,
    double xi  // <- Paramètre ajouté
) -> std::expected<std::vector<double>, SolverError> {
    
    const auto n_eta = grid_->n_eta();
    const double d_eta = grid_->d_eta();
    
    // Get current xi value and lambda0 from xi_derivatives
    const double lambda0 = xi_derivatives_->lambda0();
    const auto F_derivatives = xi_derivatives_->F_derivative();
    
    // Compute y field according to the formula:
    // y[i] = -(2ξλ₀ + 1)F[i] - 2ξ∂F/∂ξ
    std::vector<double> y_field(n_eta);
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        y_field[i] = -(2.0 * xi * lambda0 + 1.0) * solution.F[i] - 
                      2.0 * xi * F_derivatives[i];
    }

    std::cout << std::scientific << "From continuity lambda0 = " << lambda0 << "   -------   " << "xi = " << xi << std::endl;
    
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
    
    // Compute dF/deta for energy equation using 4th-order scheme
    auto dF_deta_result = coefficients::derivatives::compute_eta_derivative(solution.F, grid_->d_eta());
    if (!dF_deta_result) {
        return std::unexpected(SolverError("Failed to compute dF/deta: {}", 
                                          std::source_location::current(), dF_deta_result.error().message()));
    }
    auto dF_deta = dF_deta_result.value();
    
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
    
    auto g_solution = result.value();
/*     std::cout << "DEBUG: Energy equation result (station " << station << "): ";
    for (size_t i = 0; i < std::min(g_solution.size(), size_t(5)); ++i) {
        std::cout << g_solution[i] << " ";
    }
    std::cout << std::endl; */
    
    return g_solution;
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
    // std::cout << "DEBUG: bc.he() = " << bc.he() << std::endl;
    for (std::size_t i = 0; i < n_eta; ++i) {
        enthalpy_field[i] = g_field[i] * bc.he();
        if (!std::isfinite(enthalpy_field[i])) {
            std::cout << "DEBUG: NaN detected at i=" << i << ", g_field[i]=" << g_field[i] << ", bc.he()=" << bc.he() << std::endl;
        }
    }

    // std::cout << "h_w = " << enthalpy_field[0] << " ---------- " << "h_e = " << enthalpy_field[19] << std::endl;
    
    // DEBUG: Print enthalpy field values
/*     std::cout << "[DEBUG] update_temperature_field: Starting temperature solve..." << std::endl;
    std::cout << "[DEBUG] Enthalpy field values:" << std::endl;
    for (std::size_t i = 0; i < std::min(enthalpy_field.size(), size_t(20)); ++i) {
        std::cout << "[DEBUG] h[" << i << "] = " << enthalpy_field[i] << std::endl;
    } */
    
    auto result = h2t_solver_->solve(enthalpy_field, composition, bc, current_temperatures);
    if (!result) {
        // std::cout << "[DEBUG] Temperature solve failed!" << std::endl;
        return std::unexpected(SolverError(
            "Temperature solve failed: {}", 
            std::source_location::current(), result.error().message()
        ));
    }
    
    // std::cout << "[DEBUG] Temperature solve succeeded!" << std::endl;
    
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
    
    // Get equilibrium composition at wall conditions
    auto equilibrium_result = mixture_.equilibrium_composition(bc.Tw(), bc.P_e());
    if (!equilibrium_result) {
        return std::unexpected(SolverError(
            "Failed to compute equilibrium composition at wall conditions: {}",
            std::source_location::current(), equilibrium_result.error().message()
        ));
    }
    auto c_wall_equilibrium = equilibrium_result.value();
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        const double eta = static_cast<double>(i) * eta_max / (n_eta - 1);
        const double eta_norm = static_cast<double>(i) / (n_eta - 1);
        
        // Boundary layer-like profiles  
        guess.F[i] = eta_norm * (2.0 - eta_norm);
        
        guess.g[i] = 1.0; 
        
        for (std::size_t j = 0; j < n_species; ++j) {
            guess.c(j, i) = c_wall_equilibrium[j]; 
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
    const double factor = (xi_current > xi_prev) ? 1.1 : 0.9;
    
    auto extrapolated = previous_solution;
    
    // Extrapolate F and g fields
    for (std::size_t i = 0; i < extrapolated.F.size(); ++i) {
        extrapolated.F[i] *= factor;
        extrapolated.g[i] = std::max(0.1, extrapolated.g[i] * factor); // Keep positive
    }
    
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
    
    // High-order finite difference derivatives using unified function
    using namespace coefficients::derivatives;
    auto dF_result = compute_eta_derivative(solution.F, d_eta);
    if (!dF_result) {
        return std::unexpected(SolverError("Failed to compute dF/deta: {}", 
                                          std::source_location::current(), dF_result.error().message()));
    }
    derivatives.F = dF_result.value();
    
    auto dg_result = compute_eta_derivative(solution.g, d_eta);
    if (!dg_result) {
        return std::unexpected(SolverError("Failed to compute dg/deta: {}", 
                                          std::source_location::current(), dg_result.error().message()));
    }
    derivatives.g = dg_result.value();
    
    auto dV_result = compute_eta_derivative(solution.V, d_eta);
    if (!dV_result) {
        return std::unexpected(SolverError("Failed to compute dV/deta: {}", 
                                          std::source_location::current(), dV_result.error().message()));
    }
    derivatives.V = dV_result.value();
    
    // Species derivatives - use the high-order derivative functions
    for (std::size_t j = 0; j < n_species; ++j) {
        std::vector<double> c_row(n_eta);
        for (std::size_t i = 0; i < n_eta; ++i) {
            c_row[i] = solution.c(j, i);
        }
        
        auto dc_deta_result = compute_eta_derivative(c_row, d_eta);
        if (!dc_deta_result) {
            return std::unexpected(SolverError("Failed to compute dc/deta for species {}: {}", 
                                              std::source_location::current(), j, dc_deta_result.error().message()));
        }
        auto dc_deta = dc_deta_result.value();
        for (std::size_t i = 0; i < n_eta; ++i) {
            derivatives.c(j, i) = dc_deta[i];
        }
    }
    
    return derivatives;
}

auto BoundaryLayerSolver::compute_concentration_derivatives(
    const equations::SolutionState& solution
) const -> std::expected<DerivativeState, SolverError> {
    
    const auto n_eta = grid_->n_eta();
    const auto n_species = mixture_.n_species();
    const double d_eta = grid_->d_eta();
    
    DerivativeState derivatives(n_species, n_eta);
    
    // Use high-order finite difference schemes for both first and second derivatives
    using namespace coefficients::derivatives;
    
    for (std::size_t j = 0; j < n_species; ++j) {
        std::vector<double> c_row(n_eta);
        for (std::size_t i = 0; i < n_eta; ++i) {
            c_row[i] = solution.c(j, i);
        }
        
        // Compute first derivatives
        auto dc_deta_result = compute_eta_derivative(c_row, d_eta);
        if (!dc_deta_result) {
            return std::unexpected(SolverError("Failed to compute dc/deta for species {}: {}", 
                                              std::source_location::current(), j, dc_deta_result.error().message()));
        }
        auto dc_deta = dc_deta_result.value();
        for (std::size_t i = 0; i < n_eta; ++i) {
            derivatives.dc_deta(j, i) = dc_deta[i];
        }
        
        // Compute second derivatives
        auto dc_deta2 = compute_eta_second_derivative(c_row, d_eta);
        for (std::size_t i = 0; i < n_eta; ++i) {
            derivatives.dc_deta2(j, i) = dc_deta2[i];
        }
    }
    
    return derivatives;
}

auto BoundaryLayerSolver::enforce_edge_boundary_conditions(
    equations::SolutionState& solution,
    const conditions::BoundaryConditions& bc
) const -> void {
    
    const auto n_eta = grid_->n_eta();
    const auto n_species = mixture_.n_species();
    
    if (n_eta == 0) return;
    
    const std::size_t edge_idx = n_eta - 1;  // Last eta point is the edge
    
    // Force species composition at edge only
    
    // Force species composition at edge
    const auto& edge_composition = bc.c_e();
    for (std::size_t j = 0; j < n_species && j < edge_composition.size(); ++j) {
        solution.c(j, edge_idx) = edge_composition[j];
    }
    
    // Ensure species mass conservation (normalize if needed)
    double total_mass_fraction = 0.0;
    for (std::size_t j = 0; j < n_species; ++j) {
        total_mass_fraction += solution.c(j, edge_idx);
    }
    
    // DEBUG: Print species mass conservation details
/*     std::cout << "[DEBUG] boundary_layer_solver.cpp:718 - Species mass fraction sum: " << total_mass_fraction 
              << " (difference from 1.0: " << std::abs(total_mass_fraction - 1.0) << ")" << std::endl; */
    
    // Normalize if total is significantly different from 1.0
    if (std::abs(total_mass_fraction - 1.0) > 1e-12 && total_mass_fraction > 1e-12) {
        // std::cout << "[DEBUG] Normalizing species mass fractions..." << std::endl;
        for (std::size_t j = 0; j < n_species; ++j) {
            solution.c(j, edge_idx) /= total_mass_fraction;
        }
        // std::cout << "[DEBUG] Normalization complete." << std::endl;
    }
}

auto BoundaryLayerSolver::update_edge_properties(
    conditions::BoundaryConditions& bc,
    const coefficients::CoefficientInputs& inputs,
    const core::Matrix<double>& species_matrix
) const -> void {
    
    const auto n_eta = grid_->n_eta();
    const auto n_species = mixture_.n_species();
    
    if (n_eta == 0 || inputs.T.empty()) {
        return; // No data to update
    }
    
    // Get edge conditions (last point in eta grid)
    const auto edge_idx = n_eta - 1;
    const double T_edge = inputs.T[edge_idx];
    const double P_edge = bc.P_e(); // Pressure remains constant
    
    // Get edge composition
    std::vector<double> edge_composition(n_species);
    for (std::size_t j = 0; j < n_species; ++j) {
        edge_composition[j] = species_matrix(j, edge_idx);
    }
    
    // Calculate new edge density using equation of state
    auto MW_result = mixture_.mixture_molecular_weight(edge_composition);
    if (!MW_result) {
        throw SolverError("Failed to compute edge molecular weight: {}", 
                         std::source_location::current(), MW_result.error().message());
    }
    const double MW_edge = MW_result.value();
    const double rho_e_new = P_edge * MW_edge / 
                            (T_edge * thermophysics::constants::R_universal);
    
    // Calculate new edge viscosity using mixture properties
    auto mu_result = mixture_.viscosity(edge_composition, T_edge, P_edge);
    if (!mu_result) {
        throw SolverError("Failed to compute edge viscosity: {}", 
                         std::source_location::current(), mu_result.error().message());
    }
    const double mu_e_new = mu_result.value();
    
    // Update boundary conditions with new values
    bc.update_edge_density(rho_e_new);
    bc.update_edge_viscosity(mu_e_new);
}

} // namespace blast::boundary_layer::solver