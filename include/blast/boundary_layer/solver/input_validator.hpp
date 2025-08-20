#pragma once

#include "blast/boundary_layer/equations/equation_types.hpp"
#include "blast/boundary_layer/grid/grid.hpp"
#include "blast/boundary_layer/solver/solver_errors.hpp"
#include "blast/thermophysics/mixture_interface.hpp"

#include <expected>

namespace blast::boundary_layer::solver {

/**
 * @brief Utility class for validating solver inputs
 * 
 * This class eliminates duplication in input validation patterns
 * found throughout the solver components.
 */
class InputValidator {
public:
    InputValidator(const grid::BoundaryLayerGrid& grid,
                  const thermophysics::MixtureInterface& mixture) noexcept;

    /**
     * @brief Validate station parameters
     * 
     * @param station Station number
     * @param xi Streamwise coordinate
     * @return Success or validation error
     */
    [[nodiscard]] auto validate_station_parameters(int station, double xi) const noexcept
        -> std::expected<void, SolverError>;

    /**
     * @brief Validate solution state dimensions
     * 
     * @param solution Solution state to validate
     * @return Success or validation error
     */
    [[nodiscard]] auto validate_solution_state(const equations::SolutionState& solution) const noexcept
        -> std::expected<void, SolverError>;

    /**
     * @brief Validate xi coordinate against grid
     * 
     * @param station Station number
     * @param xi Streamwise coordinate
     * @return Success or validation error
     */
    [[nodiscard]] auto validate_xi_coordinate(int station, double xi) const noexcept
        -> std::expected<void, SolverError>;

    /**
     * @brief Validate input parameters for radiative equilibrium
     * 
     * @param emissivity Surface emissivity
     * @param T_infinity Environment temperature
     * @param q_wall Wall heat flux
     * @return Success or validation error
     */
    [[nodiscard]] static auto validate_radiative_equilibrium_inputs(
        double emissivity,
        double T_infinity,
        double q_wall
    ) noexcept -> std::expected<void, SolverError>;

    /**
     * @brief Validate solution for NaN or infinite values
     * 
     * @param solution Solution state to check
     * @return Success or validation error
     */
    [[nodiscard]] static auto validate_solution_finite(const equations::SolutionState& solution) noexcept
        -> std::expected<void, SolverError>;

private:
    const grid::BoundaryLayerGrid& grid_;
    const thermophysics::MixtureInterface& mixture_;
};

} // namespace blast::boundary_layer::solver