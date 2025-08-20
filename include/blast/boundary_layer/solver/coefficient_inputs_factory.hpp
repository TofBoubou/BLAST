#pragma once

#include "blast/boundary_layer/coefficients/coefficient_types.hpp"
#include "blast/boundary_layer/coefficients/derivative_calculator.hpp"
#include "blast/boundary_layer/equations/equation_types.hpp"
#include "blast/boundary_layer/solver/solver_errors.hpp"

#include <expected>

namespace blast::boundary_layer::solver {

/**
 * @brief Factory for creating CoefficientInputs objects with computed derivatives
 * 
 * This class eliminates the duplication of creating CoefficientInputs objects
 * which appeared 9+ times throughout the codebase with identical patterns.
 */
class CoefficientInputsFactory {
public:
    explicit CoefficientInputsFactory(coefficients::DerivativeCalculator& derivative_calculator) noexcept;

    /**
     * @brief Create coefficient inputs with computed derivatives
     * 
     * This replaces the common pattern:
     * auto derivatives_result = compute_all_derivatives(solution);
     * auto inputs = CoefficientInputs{.xi = xi, .F = solution.F, .c = solution.c, 
     *                                .dc_deta = derivatives.dc_deta, ...};
     * 
     * @param solution Current solution state
     * @param xi Streamwise coordinate
     * @return Complete coefficient inputs or error
     */
    [[nodiscard]] auto create_with_derivatives(
        const equations::SolutionState& solution,
        double xi
    ) -> std::expected<coefficients::CoefficientInputs, SolverError>;

    /**
     * @brief Create coefficient inputs without derivatives (for initial steps)
     * 
     * Used when derivatives are not yet available or needed.
     * 
     * @param solution Current solution state
     * @param xi Streamwise coordinate
     * @return Basic coefficient inputs
     */
    [[nodiscard]] auto create_basic(
        const equations::SolutionState& solution,
        double xi
    ) noexcept -> coefficients::CoefficientInputs;

    /**
     * @brief Create coefficient inputs with provided derivatives
     * 
     * For cases where derivatives are already computed externally.
     * 
     * @param solution Current solution state
     * @param xi Streamwise coordinate
     * @param derivatives Pre-computed derivatives
     * @return Coefficient inputs with provided derivatives
     */
    [[nodiscard]] auto create_with_provided_derivatives(
        const equations::SolutionState& solution,
        double xi,
        const coefficients::UnifiedDerivativeState& derivatives
    ) noexcept -> coefficients::CoefficientInputs;

private:
    coefficients::DerivativeCalculator& derivative_calculator_;
};

} // namespace blast::boundary_layer::solver