#include "blast/boundary_layer/solver/heat_flux_computer.hpp"
#include <format>

namespace blast::boundary_layer::solver {

HeatFluxComputer::HeatFluxComputer(Config config) noexcept
    : config_(config) {}

auto HeatFluxComputer::compute_complete_heat_flux(
    const equations::SolutionState& solution,
    const conditions::BoundaryConditions& bc,
    double xi,
    int station
) -> std::expected<Result, SolverError> {

    // Create result with proper dimensions
    const auto n_eta = solution.F.size();
    const auto n_species = solution.c.rows();
    Result result(n_eta, n_species);

    // Step 1: Compute all derivatives
    auto derivatives_result = config_.derivative_calculator.compute_all_derivatives(solution);
    if (!derivatives_result) {
        return std::unexpected(NumericError(
            std::format("Failed to compute derivatives for heat flux at station {}: {}", 
                       station, derivatives_result.error().message())));
    }
    result.derivatives = derivatives_result.value();

    // Step 2: Create coefficient inputs (unified pattern)
    auto inputs = coefficients::CoefficientInputs{
        .xi = xi,
        .F = solution.F,
        .c = solution.c,
        .dc_deta = result.derivatives.dc_deta,
        .dc_deta2 = result.derivatives.dc_deta2,
        .T = solution.T
    };

    // Step 3: Calculate coefficients
    auto coeffs_result = config_.coeff_calculator.calculate(inputs, bc, config_.xi_derivatives);
    if (!coeffs_result) {
        return std::unexpected(NumericError(
            std::format("Failed to compute coefficients for heat flux at station {}: {}", 
                       station, coeffs_result.error().message())));
    }
    result.coefficients = coeffs_result.value();

    // Step 4: Calculate heat flux
    auto heat_flux_result = config_.heat_flux_calculator.calculate(
        inputs, result.coefficients, bc, result.derivatives.dT_deta, station, xi);
    
    if (!heat_flux_result) {
        return std::unexpected(NumericError(
            std::format("Failed to compute heat flux at station {}: {}", 
                       station, heat_flux_result.error().message())));
    }
    result.heat_flux = heat_flux_result.value();

    return result;
}

auto HeatFluxComputer::compute_heat_flux_only(
    const equations::SolutionState& solution,
    const conditions::BoundaryConditions& bc,
    double xi,
    int station
) -> std::expected<coefficients::HeatFluxCoefficients, SolverError> {

    auto complete_result = compute_complete_heat_flux(solution, bc, xi, station);
    if (!complete_result) {
        return std::unexpected(complete_result.error());
    }

    return complete_result.value().heat_flux;
}

} // namespace blast::boundary_layer::solver