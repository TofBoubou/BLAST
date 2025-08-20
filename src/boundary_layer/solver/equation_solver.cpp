#include "blast/boundary_layer/solver/equation_solver.hpp"
#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/boundary_layer/equations/continuity.hpp"
#include "blast/boundary_layer/equations/energy.hpp"
#include "blast/boundary_layer/equations/momentum.hpp"
#include "blast/boundary_layer/equations/species.hpp"
#include <format>

namespace blast::boundary_layer::solver {

EquationSolver::EquationSolver(BoundaryLayerSolver& solver,
                              const thermophysics::MixtureInterface& mixture,
                              const io::Configuration& config) noexcept
    : solver_(solver), mixture_(mixture), config_(config) {}

auto EquationSolver::solve_momentum_equation(const equations::SolutionState& solution,
                                            const coefficients::CoefficientSet& coeffs,
                                            const conditions::BoundaryConditions& bc,
                                            double xi) -> std::expected<std::vector<double>, SolverError> {

    const auto& grid = solver_.get_grid();
    auto& xi_derivatives = solver_.get_xi_derivatives();
    
    auto result = equations::solve_momentum(solution.F, coeffs, bc, xi_derivatives, solution.V, xi, grid.d_eta());

    if (!result) {
        return std::unexpected(NumericError(std::format("Momentum equation failed: {}", result.error().message())));
    }

    return result.value();
}

auto EquationSolver::solve_energy_equation(const equations::SolutionState& solution,
                                          const coefficients::CoefficientInputs& inputs,
                                          const coefficients::CoefficientSet& coeffs,
                                          const conditions::BoundaryConditions& bc,
                                          const thermophysics::MixtureInterface& mixture,
                                          int station) -> std::expected<std::vector<double>, SolverError> {

    // Get dF/deta from unified derivative calculation
    auto all_derivatives_result = solver_.compute_all_derivatives(solution);
    if (!all_derivatives_result) {
        return std::unexpected(NumericError(std::format("Failed to compute derivatives for energy equation: {}",
                                                        all_derivatives_result.error().message())));
    }
    const auto& dF_deta = all_derivatives_result.value().dF_deta;

    const auto& grid = solver_.get_grid();
    auto& xi_derivatives = solver_.get_xi_derivatives();

    auto result = equations::solve_energy(solution.g, inputs, coeffs, bc, xi_derivatives, config_.simulation,
                                         solution.F, dF_deta, solution.V, mixture, station, grid.d_eta());

    if (!result) {
        return std::unexpected(NumericError(std::format("Energy equation failed: {}", result.error().message())));
    }

    return result.value();
}

auto EquationSolver::solve_species_equations(const equations::SolutionState& solution,
                                            const coefficients::CoefficientInputs& inputs,
                                            const coefficients::CoefficientSet& coeffs,
                                            const conditions::BoundaryConditions& bc,
                                            int station) -> std::expected<core::Matrix<double>, SolverError> {

    const auto& grid = solver_.get_grid();
    auto& xi_derivatives = solver_.get_xi_derivatives();

    auto result = equations::solve_species(solution.c, inputs, coeffs, bc, xi_derivatives, mixture_, config_.simulation,
                                          solution.F, solution.V, station, grid.d_eta());

    if (!result) {
        return std::unexpected(NumericError(std::format("Species equations failed: {}", result.error().message())));
    }

    return result.value();
}

auto EquationSolver::update_temperature_field(std::span<const double> g_field,
                                             const core::Matrix<double>& composition,
                                             const conditions::BoundaryConditions& bc,
                                             std::span<const double> current_temperatures)
    -> std::expected<std::vector<double>, SolverError> {

    const auto n_eta = g_field.size();
    std::vector<double> enthalpy_field(n_eta);

    // Convert g (dimensionless enthalpy) to dimensional enthalpy
    for (std::size_t i = 0; i < n_eta; ++i) {
        enthalpy_field[i] = g_field[i] * bc.he();
    }

    auto& h2t_solver = solver_.get_h2t_solver();
    auto result = h2t_solver.solve(enthalpy_field, composition, bc, current_temperatures);
    if (!result) {
        return std::unexpected(NumericError(std::format("Temperature solve failed: {}", result.error().message())));
    }

    return result.value().temperatures;
}

} // namespace blast::boundary_layer::solver