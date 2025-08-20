#pragma once

#include "../../core/exceptions.hpp"
#include "../../core/containers.hpp"
#include "../../thermophysics/mixture_interface.hpp"
#include "../coefficients/coefficient_calculator.hpp"
#include "../conditions/boundary_interpolator.hpp"
#include "../equations/equation_types.hpp"
#include "../thermodynamics/enthalpy_temperature_solver.hpp"
#include "solver_errors.hpp"
#include <expected>
#include <span>

namespace blast::boundary_layer::solver {

// Forward declarations
class BoundaryLayerSolver;

/**
 * @brief Handles the solution of individual governing equations
 * 
 * This class encapsulates the responsibility of solving:
 * - Momentum equation for F field
 * - Energy equation for g field  
 * - Species transport equations for composition matrix
 * - Temperature field updates from enthalpy
 */
class EquationSolver {
public:
    /**
     * @brief Constructor
     * @param solver Reference to the main boundary layer solver
     * @param mixture Reference to the thermophysics mixture
     * @param config Solver configuration
     */
    explicit EquationSolver(BoundaryLayerSolver& solver,
                           const thermophysics::MixtureInterface& mixture,
                           const io::Configuration& config) noexcept;

    /**
     * @brief Solve momentum equation for F field
     * @param solution Current solution state
     * @param coeffs Coefficient set
     * @param bc Boundary conditions
     * @param xi Streamwise coordinate
     * @return New F field or error
     */
    [[nodiscard]] auto solve_momentum_equation(const equations::SolutionState& solution,
                                              const coefficients::CoefficientSet& coeffs,
                                              const conditions::BoundaryConditions& bc,
                                              double xi) -> std::expected<std::vector<double>, SolverError>;

    /**
     * @brief Solve energy equation for g field
     * @param solution Current solution state
     * @param inputs Coefficient inputs
     * @param coeffs Coefficient set
     * @param bc Boundary conditions
     * @param mixture Thermophysics mixture
     * @param station Station number
     * @return New g field or error
     */
    [[nodiscard]] auto solve_energy_equation(const equations::SolutionState& solution,
                                            const coefficients::CoefficientInputs& inputs,
                                            const coefficients::CoefficientSet& coeffs,
                                            const conditions::BoundaryConditions& bc,
                                            const thermophysics::MixtureInterface& mixture,
                                            int station) -> std::expected<std::vector<double>, SolverError>;

    /**
     * @brief Solve species transport equations
     * @param solution Current solution state
     * @param inputs Coefficient inputs
     * @param coeffs Coefficient set
     * @param bc Boundary conditions
     * @param station Station number
     * @return New composition matrix or error
     */
    [[nodiscard]] auto solve_species_equations(const equations::SolutionState& solution,
                                              const coefficients::CoefficientInputs& inputs,
                                              const coefficients::CoefficientSet& coeffs,
                                              const conditions::BoundaryConditions& bc,
                                              int station) -> std::expected<core::Matrix<double>, SolverError>;

    /**
     * @brief Update temperature field from enthalpy
     * @param g_field Dimensionless enthalpy field
     * @param composition Species composition matrix
     * @param bc Boundary conditions
     * @param current_temperatures Current temperature field for initial guess
     * @return New temperature field or error
     */
    [[nodiscard]] auto update_temperature_field(std::span<const double> g_field,
                                               const core::Matrix<double>& composition,
                                               const conditions::BoundaryConditions& bc,
                                               std::span<const double> current_temperatures)
        -> std::expected<std::vector<double>, SolverError>;

private:
    BoundaryLayerSolver& solver_;
    const thermophysics::MixtureInterface& mixture_;
    const io::Configuration& config_;
};

} // namespace blast::boundary_layer::solver