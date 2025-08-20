#pragma once

#include "blast/boundary_layer/conditions/boundary_conditions.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include "blast/boundary_layer/coefficients/coefficient_types.hpp"
#include "blast/boundary_layer/coefficients/derivative_calculator.hpp"
#include "blast/boundary_layer/coefficients/heat_flux_calculator.hpp"
#include "blast/boundary_layer/coefficients/xi_derivatives.hpp"
#include "blast/boundary_layer/equations/equation_types.hpp"
#include "blast/boundary_layer/solver/solver_errors.hpp"

#include <expected>
#include <optional>

namespace blast::boundary_layer::solver {

/**
 * @brief Unified heat flux computer that eliminates duplication in heat flux calculations
 * 
 * This class encapsulates the common pattern of:
 * 1. Computing all derivatives
 * 2. Creating coefficient inputs
 * 3. Computing coefficients 
 * 4. Computing heat flux
 * 
 * This pattern was duplicated in 3+ locations throughout the codebase.
 */
class HeatFluxComputer {
public:
    /**
     * @brief Configuration for the heat flux computer
     */
    struct Config {
        coefficients::CoefficientCalculator& coeff_calculator;
        coefficients::HeatFluxCalculator& heat_flux_calculator;
        coefficients::DerivativeCalculator& derivative_calculator;
        coefficients::XiDerivatives& xi_derivatives;
    };

    /**
     * @brief Result of heat flux computation with all intermediate values
     */
    struct Result {
        coefficients::UnifiedDerivativeState derivatives;
        coefficients::CoefficientSet coefficients;
        coefficients::HeatFluxCoefficients heat_flux;
        
        // Constructor that properly initializes derivatives
        Result(std::size_t n_eta, std::size_t n_species) 
            : derivatives(n_eta, n_species) {}
    };

    explicit HeatFluxComputer(Config config) noexcept;

    /**
     * @brief Compute complete heat flux with all intermediate steps
     * 
     * This replaces the common pattern found throughout the codebase:
     * - auto derivatives_result = compute_all_derivatives(solution);
     * - auto final_inputs = CoefficientInputs{...};
     * - auto coeffs_result = coeff_calculator_->calculate(...);
     * - auto heat_flux_result = heat_flux_calculator_->calculate(...);
     * 
     * @param solution Current solution state
     * @param bc Boundary conditions
     * @param xi Streamwise coordinate
     * @param station Station number
     * @return Complete heat flux result or error
     */
    [[nodiscard]] auto compute_complete_heat_flux(
        const equations::SolutionState& solution,
        const conditions::BoundaryConditions& bc,
        double xi,
        int station
    ) -> std::expected<Result, SolverError>;

    /**
     * @brief Compute only the heat flux value (simplified interface)
     * 
     * @param solution Current solution state
     * @param bc Boundary conditions
     * @param xi Streamwise coordinate
     * @param station Station number
     * @return Heat flux data or error
     */
    [[nodiscard]] auto compute_heat_flux_only(
        const equations::SolutionState& solution,
        const conditions::BoundaryConditions& bc,
        double xi,
        int station
    ) -> std::expected<coefficients::HeatFluxCoefficients, SolverError>;

private:
    Config config_;
};

} // namespace blast::boundary_layer::solver