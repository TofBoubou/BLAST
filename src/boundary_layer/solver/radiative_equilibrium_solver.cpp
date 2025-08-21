#include "blast/boundary_layer/solver/radiative_equilibrium_solver.hpp"
#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/boundary_layer/solver/expected_utils.hpp"
#include "blast/boundary_layer/solver/input_validator.hpp"
#include <cmath>
#include <format>
#include <iomanip>
#include <iostream>

namespace blast::boundary_layer::solver {

RadiativeEquilibriumSolver::RadiativeEquilibriumSolver(BoundaryLayerSolver& solver,
                                                      const io::Configuration& config) noexcept
    : solver_(solver), config_(config) {}

auto RadiativeEquilibriumSolver::solve_radiative_equilibrium(double q_wall, double emissivity, double T_infinity)
    -> std::expected<double, std::string> {
    
    // Validate inputs using unified validator
    if (auto validation_result = InputValidator::validate_radiative_equilibrium_inputs(emissivity, T_infinity, q_wall); !validation_result) {
        return std::unexpected(validation_result.error().message());
    }
    
    const double T_inf_4 = std::pow(T_infinity, 4);
    const double T_wall_4 = q_wall / (emissivity * STEFAN_BOLTZMANN) + T_inf_4;
    
    if (T_wall_4 <= 0.0) {
        return std::unexpected(std::format("Radiative equilibrium solution invalid: T_wall^4 = {} <= 0 (q_wall={}, ε={}, T_∞={})", 
                                          T_wall_4, q_wall, emissivity, T_infinity));
    }
    
    const double T_wall = std::pow(T_wall_4, 0.25);
    
    if (!std::isfinite(T_wall) || T_wall <= 0.0) {
        return std::unexpected(std::format("Invalid computed wall temperature: {} K", T_wall));
    }
    
    return T_wall;
}

auto RadiativeEquilibriumSolver::update_wall_temperature_iteration(int station, double xi,
                                                                  conditions::BoundaryConditions& bc,
                                                                  const equations::SolutionState& solution,
                                                                  int iteration) -> std::expected<void, SolverError> {
    
    // Use unified heat flux computer to eliminate duplication
    auto& heat_flux_computer = solver_.get_heat_flux_computer();
    auto heat_flux_data = BLAST_TRY_WITH_CONTEXT(
        heat_flux_computer.compute_heat_flux_only(solution, bc, xi, station),
        "Failed to compute heat flux for radiative equilibrium"
    );
    
    double q_wall = heat_flux_data.q_wall_total_dim;
    
    // Calculate current radiative flux for debugging
    const double T_wall_current = bc.wall.temperature;
    const double T_inf = config_.wall_parameters.environment_temperature;
    const double q_rad = calculate_radiative_flux(T_wall_current, T_inf, config_.wall_parameters.emissivity);
    
    // Solve for new wall temperature
    auto T_wall_result = solve_radiative_equilibrium(
        q_wall, 
        config_.wall_parameters.emissivity,
        config_.wall_parameters.environment_temperature
    );

    if (!T_wall_result) {
        return std::unexpected(NumericError(std::format("Radiative equilibrium failed: {}", T_wall_result.error())));
    }

    // Update wall temperature
    double T_wall_old = bc.wall.temperature;
    bc.wall.temperature = T_wall_result.value();
    
    // Optional debug output (commented out for production)
/*     std::cout << std::format("[RADIATIVE] Station {} Iter {}: q_wall={:.2e} W/m², q_rad={:.2e} W/m², "
                           "q_wall-q_rad={:.2e} W/m², T_wall={:.1f} K", 
                           station, iteration, q_wall, q_rad, q_wall - q_rad, T_wall_current) << std::endl;
    
    std::cout << std::format("[DEBUG] T_wall update: {:.1f} K -> {:.1f} K (delta={:.1f} K)", 
                           T_wall_old, bc.wall.temperature, 
                           bc.wall.temperature - T_wall_old) << std::endl; */
    
    return {};
}

auto RadiativeEquilibriumSolver::compute_heat_flux_analysis(int station, double xi,
                                                           const conditions::BoundaryConditions& bc,
                                                           const equations::SolutionState& solution) const
    -> std::expected<HeatFluxComponents, SolverError> {
    
    // Use unified heat flux computer for analysis
    auto& heat_flux_computer = solver_.get_heat_flux_computer();
    auto heat_flux = BLAST_TRY_WITH_CONTEXT(
        heat_flux_computer.compute_heat_flux_only(solution, bc, xi, station),
        "Failed to compute heat flux for analysis"
    );
    
    // Calculate radiative flux
    const double T_wall = bc.wall.temperature;
    const double T_inf = config_.wall_parameters.environment_temperature;
    const double q_rad = calculate_radiative_flux(T_wall, T_inf, config_.wall_parameters.emissivity);
    
    return HeatFluxComponents{
        .q_wall_conductive = heat_flux.q_wall_conductive_dim,
        .q_wall_diffusive = heat_flux.q_wall_diffusive_dim,
        .q_wall_total = heat_flux.q_wall_total_dim,
        .q_radiative = q_rad,
        .q_balance = heat_flux.q_wall_total_dim - q_rad,
        .wall_temperature = T_wall
    };
}

auto RadiativeEquilibriumSolver::print_heat_flux_summary(int station, const HeatFluxComponents& components) -> void {
    std::cout << std::format("\n=== FINAL HEAT FLUX SUMMARY - Station {} ===", station) << std::endl;
    std::cout << std::format("  Wall Temperature: {:.1f} K", components.wall_temperature) << std::endl;
    std::cout << std::format("  Conductive flux:  {:.2e} W/m²", components.q_wall_conductive) << std::endl;
    std::cout << std::format("  Diffusive flux:   {:.2e} W/m²", components.q_wall_diffusive) << std::endl;
    std::cout << std::format("  TOTAL flux:       {:.2e} W/m²", components.q_wall_total) << std::endl;
    std::cout << std::format("  Radiative flux:   {:.2e} W/m²", components.q_radiative) << std::endl;
    std::cout << std::format("  Balance (q-qrad): {:.2e} W/m²", components.q_balance) << std::endl;
    std::cout << "============================================\n" << std::endl;
}

auto RadiativeEquilibriumSolver::calculate_radiative_flux(double T_wall, double T_infinity, double emissivity) noexcept
    -> double {
    return emissivity * STEFAN_BOLTZMANN * (std::pow(T_wall, 4) - std::pow(T_infinity, 4));
}

// Note: validate_inputs method moved to InputValidator::validate_radiative_equilibrium_inputs

} // namespace blast::boundary_layer::solver