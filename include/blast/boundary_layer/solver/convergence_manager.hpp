#pragma once

#include "../../core/exceptions.hpp"
#include "../../io/config_types.hpp"
#include "../../thermophysics/mixture_interface.hpp"
#include "../conditions/boundary_interpolator.hpp"
#include "../equations/equation_types.hpp"
#include "adaptive_relaxation_controller.hpp"
#include "solver_errors.hpp"
#include "solver_steps.hpp"
#include <expected>

namespace blast::boundary_layer::solver {

// Forward declarations
class BoundaryLayerSolver;
class RadiativeEquilibriumSolver;

/**
 * @brief Manages the iterative convergence process for boundary layer stations
 * 
 * This class encapsulates:
 * - Adaptive relaxation control
 * - Convergence checking with customizable tolerances
 * - Iteration management and divergence detection
 * - Integration with radiative equilibrium updates
 */
class ConvergenceManager {
public:
    /**
     * @brief Constructor
     * @param solver Reference to the main boundary layer solver
     * @param mixture Reference to the thermophysics mixture
     * @param config Solver configuration
     */
    explicit ConvergenceManager(BoundaryLayerSolver& solver,
                               const thermophysics::MixtureInterface& mixture,
                               const io::Configuration& config) noexcept;

    /**
     * @brief Perform adaptive iteration for a station until convergence
     * @param station Station number
     * @param xi Streamwise coordinate
     * @param bc Boundary conditions (may be modified for radiative equilibrium)
     * @param solution Solution state (modified in-place)
     * @return Convergence information or error
     */
    [[nodiscard]] auto iterate_station_adaptive(int station, double xi, 
                                               conditions::BoundaryConditions& bc,
                                               equations::SolutionState& solution)
        -> std::expected<ConvergenceInfo, SolverError>;

    /**
     * @brief Check convergence between old and new solutions
     * @param old_solution Previous iteration solution
     * @param new_solution Current iteration solution
     * @return Convergence information
     */
    [[nodiscard]] auto check_convergence(const equations::SolutionState& old_solution,
                                        const equations::SolutionState& new_solution) const noexcept 
        -> ConvergenceInfo;

    /**
     * @brief Apply differential relaxation between solutions
     * @param old_solution Previous solution
     * @param new_solution New solution
     * @param base_factor Base relaxation factor
     * @return Relaxed solution
     */
    [[nodiscard]] auto apply_relaxation_differential(const equations::SolutionState& old_solution,
                                                    const equations::SolutionState& new_solution,
                                                    double base_factor) const -> equations::SolutionState;

    // Continuation mode is queried directly from BoundaryLayerSolver

    /**
     * @brief Initialize relaxation controller for station type
     * @param station Station number (0 = stagnation point, >0 = downstream)
     */
    auto initialize_relaxation_for_station(int station) -> void;

    /**
     * @brief Set radiative equilibrium solver for wall temperature updates
     * @param radiative_solver Pointer to radiative equilibrium solver
     */
    auto set_radiative_solver(RadiativeEquilibriumSolver* radiative_solver) noexcept -> void {
        radiative_solver_ = radiative_solver;
    }

private:
    /**
     * @brief Execute solver pipeline for one iteration
     * @param station Station number
     * @param xi Streamwise coordinate
     * @param bc Boundary conditions
     * @param solution Current solution
     * @param iteration Iteration number
     * @return New solution state or error
     */
    [[nodiscard]] auto execute_solver_pipeline(int station, double xi,
                                              conditions::BoundaryConditions& bc,
                                              equations::SolutionState& solution,
                                              int iteration) -> std::expected<void, SolverError>;

    /**
     * @brief Update wall temperature for radiative equilibrium
     * @param station Station number
     * @param xi Streamwise coordinate
     * @param bc Boundary conditions (modified in-place)
     * @param solution Current solution
     * @param iteration Iteration number
     * @return Success or error
     */
    [[nodiscard]] auto update_wall_temperature_radiative(int station, double xi,
                                                        conditions::BoundaryConditions& bc,
                                                        const equations::SolutionState& solution,
                                                        int iteration) -> std::expected<void, SolverError>;

    /**
     * @brief Check for solution divergence
     * @param conv_info Convergence information
     * @param station Station number
     * @param iteration Iteration number
     * @return Error if solution has diverged
     */
    [[nodiscard]] auto check_divergence(const ConvergenceInfo& conv_info, int station, int iteration) const
        -> std::expected<void, SolverError>;

    /**
     * @brief Print final heat flux summary for radiative equilibrium
     * @param station Station number
     * @param xi Streamwise coordinate
     * @param bc Final boundary conditions
     * @param solution Final solution
     */
    auto print_final_heat_flux_summary(int station, double xi,
                                      const conditions::BoundaryConditions& bc,
                                      const equations::SolutionState& solution) const -> void;

private:
    BoundaryLayerSolver& solver_;
    const thermophysics::MixtureInterface& mixture_;
    const io::Configuration& config_;
    std::unique_ptr<AdaptiveRelaxationController> relaxation_controller_;
    RadiativeEquilibriumSolver* radiative_solver_{nullptr};
    // No local continuation state; rely on solver_.is_in_continuation()
};

} // namespace blast::boundary_layer::solver
