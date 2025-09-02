#include "blast/boundary_layer/solver/convergence_manager.hpp"
#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/boundary_layer/solver/radiative_equilibrium_solver.hpp"
#include <algorithm>
#include <cmath>
#include <format>
#include <iomanip>
#include <iostream>

namespace blast::boundary_layer::solver {

ConvergenceManager::ConvergenceManager(BoundaryLayerSolver& solver,
                                      const thermophysics::MixtureInterface& mixture,
                                      const io::Configuration& config) noexcept
    : solver_(solver), mixture_(mixture), config_(config) {}

auto ConvergenceManager::iterate_station_adaptive(int station, double xi, 
                                                 conditions::BoundaryConditions& bc,
                                                 equations::SolutionState& solution)
    -> std::expected<ConvergenceInfo, SolverError> {

    const auto n_eta = solver_.get_grid().n_eta();
    const auto n_species = mixture_.n_species();

    auto bc_dynamic = bc;
    ConvergenceInfo conv_info;

    for (int iter = 0; iter < config_.numerical.max_iterations; ++iter) {
        const auto solution_old = solution;

        // Execute solver pipeline for this iteration
        auto pipeline_result = execute_solver_pipeline(station, xi, bc_dynamic, solution, iter);
        if (!pipeline_result) {
            return std::unexpected(pipeline_result.error());
        }

        // New solution after pipeline execution
        equations::SolutionState solution_new = solution;

        // Convergence check
        conv_info = check_convergence(solution_old, solution_new);
        conv_info.iterations = iter + 1;

        // Handle NaN detection
        if (auto nan_result = handle_nan_detection(conv_info, station, iter); !nan_result) {
            if (in_continuation_) {
                conv_info.converged = false;
                conv_info.iterations = iter + 1;
                return conv_info;
            } else {
                return std::unexpected(nan_result.error());
            }
        }

        // Adaptive relaxation
        auto adaptive_factor_result = relaxation_controller_->adapt_relaxation_factor(conv_info, iter);
        if (!adaptive_factor_result) {
            return std::unexpected<SolverError>(adaptive_factor_result.error());
        }
        double adaptive_factor = adaptive_factor_result.value();

        // Check convergence first
        if (conv_info.converged) {
            // Already at new solution state from pipeline
            return conv_info;
        }

        // Apply relaxation only if not converged
        solution = apply_relaxation_differential(solution_old, solution_new, adaptive_factor);
        
        // Update wall temperature for radiative equilibrium AFTER checking convergence
        if (radiative_solver_ && radiative_solver_->is_radiative_equilibrium_enabled()) {
            auto radiative_result = update_wall_temperature_radiative(station, xi, bc_dynamic, solution, iter);
            if (!radiative_result) {
                return std::unexpected(radiative_result.error());
            }
        }

        // Check for divergence
        if (auto divergence_result = check_divergence(conv_info, station, iter); !divergence_result) {
            return std::unexpected(divergence_result.error());
        }
    }

    // Print final heat flux summary if radiative equilibrium was used
    if (radiative_solver_ && radiative_solver_->is_radiative_equilibrium_enabled()) {
        print_final_heat_flux_summary(station, xi, bc_dynamic, solution);
    }

    // Only reach here if max iterations reached without convergence
    return conv_info;
}

auto ConvergenceManager::check_convergence(const equations::SolutionState& old_solution,
                                          const equations::SolutionState& new_solution) const noexcept
    -> ConvergenceInfo {

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

/*     std::cout << "CONVERGENCE : " << std::setprecision(12) << info.residual_F << " " << info.residual_g << " " << info.residual_c
            << " | tol=" << tol << " | converged=" << info.converged
            << " | g_wall=" << std::format("{:.12g}", new_solution.g[0])
            << " | g_mid="  << std::format("{:.12g}", new_solution.g[new_solution.g.size()/2])
            << std::endl; */

    
    // Print mass fractions at wall (eta=0)
/*     std::cout << "Mass fractions at wall: ";
    for (std::size_t s = 0; s < new_solution.c.rows(); ++s) {
        std::cout << mixture_.species_name(s) << "=" << std::format("{:.6e}", new_solution.c(s, 0)) << " ";
    }
    std::cout << std::endl; */

    return info;
}

auto ConvergenceManager::apply_relaxation_differential(const equations::SolutionState& old_solution,
                                                      const equations::SolutionState& new_solution,
                                                      double base_factor) const -> equations::SolutionState {
    
    const auto n_eta = old_solution.F.size();
    const auto n_species = old_solution.c.rows();
    
    equations::SolutionState relaxed_solution(n_eta, n_species);
    
    // Apply relaxation to each field
    for (std::size_t i = 0; i < n_eta; ++i) {
        relaxed_solution.F[i] = old_solution.F[i] + base_factor * (new_solution.F[i] - old_solution.F[i]);
        relaxed_solution.g[i] = old_solution.g[i] + base_factor * (new_solution.g[i] - old_solution.g[i]);
        relaxed_solution.T[i] = old_solution.T[i] + base_factor * (new_solution.T[i] - old_solution.T[i]);
        relaxed_solution.V[i] = old_solution.V[i] + base_factor * (new_solution.V[i] - old_solution.V[i]);
        
        for (std::size_t s = 0; s < n_species; ++s) {
            relaxed_solution.c(s, i) = old_solution.c(s, i) + base_factor * (new_solution.c(s, i) - old_solution.c(s, i));
        }
    }
    
    return relaxed_solution;
}

auto ConvergenceManager::initialize_relaxation_for_station(int station) -> void {
    if (station == 0) {
        relaxation_controller_ = 
            std::make_unique<AdaptiveRelaxationController>(AdaptiveRelaxationController::Config::for_stagnation_point());
    } else {
        relaxation_controller_ = 
            std::make_unique<AdaptiveRelaxationController>(AdaptiveRelaxationController::Config::for_downstream_station());
    }
}

auto ConvergenceManager::execute_solver_pipeline(int station, double xi,
                                                conditions::BoundaryConditions& bc,
                                                equations::SolutionState& solution,
                                                int iteration) -> std::expected<void, SolverError> {

    // Initialize coefficients container (filled by steps)
    coefficients::CoefficientSet coeffs;

    // Create context (no const_cast; pipeline mutates state explicitly)
    SolverContext ctx{.solution = solution,
                     .solution_old = solution,
                     .bc = bc,
                     .coeffs = coeffs,
                     .mixture = mixture_,
                     .station = station,
                     .xi = xi,
                     .iteration = iteration,
                     .solver = solver_,
                     .grid = solver_.get_grid(),
                     .xi_derivatives = solver_.get_xi_derivatives()};

    // Create and execute pipeline
    auto pipeline = SolverPipeline::create_for_solver(solver_);
    auto pipeline_result = pipeline.execute_all(ctx);
    if (!pipeline_result) {
        return std::unexpected(NumericError(std::format("Pipeline execution failed at station {} iteration {}: {}",
                                                        station, iteration, pipeline_result.error().what())));
    }

    return {};
}

auto ConvergenceManager::update_wall_temperature_radiative(int station, double xi,
                                                          conditions::BoundaryConditions& bc,
                                                          const equations::SolutionState& solution,
                                                          int iteration) -> std::expected<void, SolverError> {
    if (!radiative_solver_) {
        return std::unexpected(NumericError("Radiative solver not available"));
    }

    return radiative_solver_->update_wall_temperature_iteration(station, xi, bc, solution, iteration);
}

auto ConvergenceManager::check_divergence(const ConvergenceInfo& conv_info, int station, int iteration) const
    -> std::expected<void, SolverError> {
    
    if (conv_info.max_residual() > config_.numerical.divergence_threshold 1e6) {
        return std::unexpected(NumericError(std::format("Solution diverged at station {} iteration {} (residual={})",
                                                        station, iteration, conv_info.max_residual())));
    }
    return {};
}

auto ConvergenceManager::handle_nan_detection(const ConvergenceInfo& conv_info, int station, int iteration) const
    -> std::expected<void, SolverError> {
    
    if (std::isnan(conv_info.residual_F) || std::isnan(conv_info.residual_g) || std::isnan(conv_info.residual_c)) {
        return std::unexpected(NumericError(std::format("NaN detected in residuals at station {} iteration {}", 
                                                        station, iteration)));
    }
    return {};
}

auto ConvergenceManager::print_final_heat_flux_summary(int station, double xi,
                                                      const conditions::BoundaryConditions& bc,
                                                      const equations::SolutionState& solution) const -> void {
    if (!radiative_solver_) return;

    auto heat_flux_result = radiative_solver_->compute_heat_flux_analysis(station, xi, bc, solution);
    if (heat_flux_result) {
        RadiativeEquilibriumSolver::print_heat_flux_summary(station, heat_flux_result.value());
    }
}

} // namespace blast::boundary_layer::solver
