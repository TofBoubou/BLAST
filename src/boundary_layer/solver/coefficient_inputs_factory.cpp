#include "blast/boundary_layer/solver/coefficient_inputs_factory.hpp"
#include <format>

namespace blast::boundary_layer::solver {

CoefficientInputsFactory::CoefficientInputsFactory(coefficients::DerivativeCalculator& derivative_calculator) noexcept
    : derivative_calculator_(derivative_calculator) {}

auto CoefficientInputsFactory::create_with_derivatives(
    const equations::SolutionState& solution,
    double xi
) -> std::expected<coefficients::CoefficientInputs, SolverError> {

    // Compute derivatives first
    auto derivatives_result = derivative_calculator_.compute_all_derivatives(solution);
    if (!derivatives_result) {
        return std::unexpected(NumericError(
            std::format("Failed to compute derivatives for coefficient inputs at xi={}: {}", 
                       xi, derivatives_result.error().message())));
    }

    auto derivatives = derivatives_result.value();

    // Create complete coefficient inputs
    return coefficients::CoefficientInputs{
        .xi = xi,
        .F = solution.F,
        .c = solution.c,
        .dc_deta = derivatives.dc_deta,
        .dc_deta2 = derivatives.dc_deta2,
        .T = solution.T
    };
}

auto CoefficientInputsFactory::create_basic(
    const equations::SolutionState& solution,
    double xi
) noexcept -> coefficients::CoefficientInputs {

    return coefficients::CoefficientInputs{
        .xi = xi,
        .F = solution.F,
        .c = solution.c,
        .dc_deta = core::Matrix<double>(),  // Empty derivatives
        .dc_deta2 = core::Matrix<double>(), // Empty derivatives
        .T = solution.T
    };
}

auto CoefficientInputsFactory::create_with_provided_derivatives(
    const equations::SolutionState& solution,
    double xi,
    const coefficients::UnifiedDerivativeState& derivatives
) noexcept -> coefficients::CoefficientInputs {

    return coefficients::CoefficientInputs{
        .xi = xi,
        .F = solution.F,
        .c = solution.c,
        .dc_deta = derivatives.dc_deta,
        .dc_deta2 = derivatives.dc_deta2,
        .T = solution.T
    };
}

} // namespace blast::boundary_layer::solver