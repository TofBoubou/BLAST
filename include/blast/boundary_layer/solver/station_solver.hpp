#pragma once

#include "../../core/exceptions.hpp"
#include "../../io/config_types.hpp"
#include "../../thermophysics/mixture_interface.hpp"
#include "../conditions/boundary_interpolator.hpp"
#include "../equations/equation_types.hpp"
#include "../grid/grid.hpp"
#include "solver_errors.hpp"
#include <expected>
#include <functional>

namespace blast::boundary_layer::solver {

// Forward declarations
class BoundaryLayerSolver;
class ConvergenceManager;

/**
 * @brief Handles the resolution of individual stations in the boundary layer solution
 * 
 * This class encapsulates the responsibility of solving a single station, including:
 * - Input validation for station parameters
 * - Boundary condition creation/interpolation
 * - Initial guess generation
 * - Coordination with continuation methods when direct solution fails
 */
class StationSolver {
public:
    /**
     * @brief Configuration for stable guess computation used in continuation
     */
    struct StableGuessConfig {
        double wall_temperature_stable;
        double edge_temperature_stable;
        double pressure_stable;
    };

    /**
     * @brief Constructor
     * @param solver Reference to the main boundary layer solver
     * @param mixture Reference to the thermophysics mixture
     * @param config Solver configuration
     */
    explicit StationSolver(BoundaryLayerSolver& solver,
                          const thermophysics::MixtureInterface& mixture,
                          const io::Configuration& config) noexcept;

    /**
     * @brief Solve a single station
     * @param station Station number (0-based)
     * @param xi Streamwise coordinate
     * @param initial_guess Initial solution guess
     * @return Solution state or error
     */
    [[nodiscard]] auto solve_station(int station, double xi, const equations::SolutionState& initial_guess)
        -> std::expected<equations::SolutionState, SolverError>;

    // Note: create_initial_guess method moved to InitialGuessFactory

    /**
     * @brief Set stable configuration for continuation method
     * @param stable_config Configuration for stable initial guess
     */
    auto set_stable_config(const StableGuessConfig& stable_config) noexcept -> void {
        stable_config_ = stable_config;
    }

    /**
     * @brief Set continuation flag
     * @param in_continuation True if currently in continuation mode
     */
    auto set_continuation_mode(bool in_continuation) noexcept -> void {
        in_continuation_ = in_continuation;
    }

private:
    /**
     * @brief Validate input parameters for station solving
     * @param station Station number
     * @param xi Streamwise coordinate
     * @param initial_guess Initial solution guess
     * @return Error if validation fails
     */
    // Note: validate_station_inputs method moved to InputValidator

    /**
     * @brief Get boundary conditions for a station
     * @param station Station number
     * @param xi Streamwise coordinate
     * @return Boundary conditions or error
     */
    [[nodiscard]] auto get_boundary_conditions(int station, double xi) const
        -> std::expected<conditions::BoundaryConditions, SolverError>;

    /**
     * @brief Compute stable guess for continuation method
     * @param station Station number
     * @param xi Streamwise coordinate
     * @return Stable initial guess or error
     */
    [[nodiscard]] auto compute_stable_guess(int station, double xi) const
        -> std::expected<equations::SolutionState, SolverError>;

    /**
     * @brief Check if continuation should be attempted
     * @param convergence_error The convergence error from direct solution
     * @return True if continuation should be tried
     */
    [[nodiscard]] auto should_attempt_continuation(const ConvergenceError& convergence_error) const noexcept -> bool;

private:
    BoundaryLayerSolver& solver_;
    const thermophysics::MixtureInterface& mixture_;
    const io::Configuration& config_;
    StableGuessConfig stable_config_{};
    bool in_continuation_{false};
};

} // namespace blast::boundary_layer::solver