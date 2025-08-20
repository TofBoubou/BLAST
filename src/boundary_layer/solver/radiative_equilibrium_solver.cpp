#include "blast/boundary_layer/solver/radiative_equilibrium_solver.hpp"
#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
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
    
    // Validate inputs
    if (auto validation_result = validate_inputs(q_wall, emissivity, T_infinity); !validation_result) {
        return std::unexpected(validation_result.error());
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
    
    // Compute all derivatives
    auto derivatives_result = solver_.compute_all_derivatives(solution);
    if (!derivatives_result) {
        return std::unexpected(NumericError(std::format("Failed to compute derivatives for radiative equilibrium: {}",
                                                        derivatives_result.error().message())));
    }
    auto derivatives = derivatives_result.value();
    
    // Create coefficient inputs
    auto final_inputs = coefficients::CoefficientInputs{
        .xi = xi,
        .F = solution.F,
        .c = solution.c,
        .dc_deta = derivatives.dc_deta,
        .dc_deta2 = derivatives.dc_deta2,
        .T = solution.T
    };
    
    // Calculate coefficients
    auto& coeff_calculator = solver_.get_coeff_calculator();
    auto& xi_derivatives = solver_.get_xi_derivatives();
    auto coeffs_result = coeff_calculator.calculate(final_inputs, bc, xi_derivatives);
    if (!coeffs_result) {
        return std::unexpected(NumericError(std::format("Failed to calculate coefficients for radiative equilibrium: {}",
                                                        coeffs_result.error().message())));
    }
    auto coeffs = coeffs_result.value();
    
    // Calculate heat flux
    auto& heat_flux_calculator = solver_.get_heat_flux_calculator();
    auto heat_flux_result = heat_flux_calculator.calculate(
        final_inputs, coeffs, bc, derivatives.dT_deta, station, xi);
    
    if (!heat_flux_result) {
        return std::unexpected(NumericError(std::format("Failed to calculate heat flux for radiative equilibrium: {}",
                                                        heat_flux_result.error().message())));
    }
    
    double q_wall = heat_flux_result.value().q_wall_total_dim;
    
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
    std::cout << std::format("[RADIATIVE] Station {} Iter {}: q_wall={:.2e} W/m², q_rad={:.2e} W/m², "
                           "q_wall-q_rad={:.2e} W/m², T_wall={:.1f} K", 
                           station, iteration, q_wall, q_rad, q_wall - q_rad, T_wall_current) << std::endl;
    
    std::cout << std::format("[DEBUG] T_wall update: {:.1f} K -> {:.1f} K (delta={:.1f} K)", 
                           T_wall_old, bc.wall.temperature, 
                           bc.wall.temperature - T_wall_old) << std::endl;
    
    return {};
}

auto RadiativeEquilibriumSolver::compute_heat_flux_analysis(int station, double xi,
                                                           const conditions::BoundaryConditions& bc,
                                                           const equations::SolutionState& solution) const
    -> std::expected<HeatFluxComponents, SolverError> {
    
    // Compute derivatives
    auto derivatives_result = solver_.compute_all_derivatives(solution);
    if (!derivatives_result) {
        return std::unexpected(NumericError(std::format("Failed to compute derivatives for heat flux analysis: {}",
                                                        derivatives_result.error().message())));
    }
    auto derivatives = derivatives_result.value();
    
    // Create coefficient inputs
    auto final_inputs = coefficients::CoefficientInputs{
        .xi = xi,
        .F = solution.F,
        .c = solution.c,
        .dc_deta = derivatives.dc_deta,
        .dc_deta2 = derivatives.dc_deta2,
        .T = solution.T
    };
    
    // Calculate coefficients
    auto& coeff_calculator = solver_.get_coeff_calculator();
    auto& xi_derivatives = solver_.get_xi_derivatives();
    auto coeffs_result = coeff_calculator.calculate(final_inputs, bc, xi_derivatives);
    if (!coeffs_result) {
        return std::unexpected(NumericError(std::format("Failed to calculate coefficients for heat flux analysis: {}",
                                                        coeffs_result.error().message())));
    }
    auto coeffs = coeffs_result.value();
    
    // Calculate heat flux
    auto& heat_flux_calculator = solver_.get_heat_flux_calculator();
    auto heat_flux_result = heat_flux_calculator.calculate(
        final_inputs, coeffs, bc, derivatives.dT_deta, station, xi);
    
    if (!heat_flux_result) {
        return std::unexpected(NumericError(std::format("Failed to calculate heat flux for analysis: {}",
                                                        heat_flux_result.error().message())));
    }
    
    auto heat_flux = heat_flux_result.value();
    
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

auto RadiativeEquilibriumSolver::validate_inputs(double q_wall, double emissivity, double T_infinity)
    -> std::expected<void, std::string> {
    
    if (emissivity <= 0.0) {
        return std::unexpected(std::format("Invalid emissivity: {} (must be > 0)", emissivity));
    }
    
    if (T_infinity <= 0.0) {
        return std::unexpected(std::format("Invalid environment temperature: {} K (must be > 0)", T_infinity));
    }
    
    if (!std::isfinite(q_wall)) {
        return std::unexpected(std::format("Invalid wall heat flux: {} (not finite)", q_wall));
    }
    
    return {};
}

} // namespace blast::boundary_layer::solver