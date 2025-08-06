#include "blast/boundary_layer/solver/continuation_method.hpp"
#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include <algorithm>
#include <cmath>

namespace blast::boundary_layer::solver {

auto ContinuationMethod::solve_with_continuation(
    BoundaryLayerSolver& solver,
    int station,
    double xi,
    const io::Configuration& target_config,
    const equations::SolutionState& initial_guess
) -> std::expected<ContinuationResult, SolverError> {
    
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
        auto result = solver.solve_station(station, xi, current_solution);
        solver.set_config(original_config);
        
        if (result) {
            current_solution = result.value();
            lambda = lambda_try;
            
            if (std::abs(lambda - 1.0) < 1e-10) {
                return ContinuationResult{
                    .solution = current_solution,
                    .success = true,
                    .final_lambda = lambda
                };
            }
            
            // Increase step after success
            lambda_step = std::min(lambda_step * STEP_INCREASE_FACTOR, LAMBDA_STEP_MAX);
            
        } else {
            // Decrease step after failure
            lambda_step *= STEP_DECREASE_FACTOR;
            
            if (lambda_step < LAMBDA_STEP_MIN) {
                return ContinuationResult{
                    .solution = current_solution,
                    .success = false,
                    .final_lambda = lambda
                };
            }
        }
    }
    
    return std::unexpected(SolverError(
        "Continuation max steps exceeded", 
        std::source_location::current()
    ));
}

auto ContinuationMethod::interpolate_config(
    const io::Configuration& target,
    double lambda
) const -> io::Configuration {
    
    io::Configuration config = target;
    
    // Interpolate wall temperature
    if (!config.wall_parameters.wall_temperatures.empty()) {
        double Twall_target = target.wall_parameters.wall_temperatures[0];
        config.wall_parameters.wall_temperatures[0] = 
            TWALL_STABLE + lambda * (Twall_target - TWALL_STABLE);
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