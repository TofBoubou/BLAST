#pragma once

#include "../../core/exceptions.hpp"
#include "../../io/config_types.hpp"
#include "../coefficients/coefficient_calculator.hpp"
#include "../coefficients/heat_flux_calculator.hpp"
#include "../conditions/boundary_interpolator.hpp"
#include "../equations/equation_types.hpp"
#include "solver_errors.hpp"
#include <expected>
#include <string>

namespace blast::boundary_layer::solver {

// Forward declarations
class BoundaryLayerSolver;

/**
 * @brief Handles radiative equilibrium calculations for wall temperature
 * 
 * This class encapsulates:
 * - Radiative heat flux calculations
 * - Wall temperature computation from heat balance
 * - Stefan-Boltzmann law application
 * - Heat flux balance validation
 */
class RadiativeEquilibriumSolver {
public:
    /**
     * @brief Heat flux components for analysis
     */
    struct HeatFluxComponents {
        double q_wall_conductive;    ///< Conductive heat flux [W/m²]
        double q_wall_diffusive;     ///< Diffusive heat flux [W/m²]
        double q_wall_total;         ///< Total wall heat flux [W/m²]
        double q_radiative;          ///< Radiative heat flux [W/m²]
        double q_balance;            ///< Heat balance (total - radiative) [W/m²]
        double wall_temperature;     ///< Wall temperature [K]
    };

    /**
     * @brief Constructor
     * @param solver Reference to the main boundary layer solver
     * @param config Solver configuration
     */
    explicit RadiativeEquilibriumSolver(BoundaryLayerSolver& solver,
                                       const io::Configuration& config) noexcept;

    /**
     * @brief Solve radiative equilibrium for wall temperature
     * @param q_wall Wall heat flux from boundary layer [W/m²]
     * @param emissivity Surface emissivity [-]
     * @param T_infinity Environment temperature [K]
     * @return Wall temperature in radiative equilibrium [K] or error
     */
    [[nodiscard]] static auto solve_radiative_equilibrium(double q_wall, double emissivity, double T_infinity)
        -> std::expected<double, std::string>;

    /**
     * @brief Update wall temperature for radiative equilibrium during iteration
     * @param station Station number
     * @param xi Streamwise coordinate
     * @param bc Boundary conditions (modified in-place)
     * @param solution Current solution state
     * @param iteration Iteration number
     * @return Success or error
     */
    [[nodiscard]] auto update_wall_temperature_iteration(int station, double xi,
                                                        conditions::BoundaryConditions& bc,
                                                        const equations::SolutionState& solution,
                                                        int iteration) -> std::expected<void, SolverError>;

    /**
     * @brief Compute complete heat flux analysis for final summary
     * @param station Station number
     * @param xi Streamwise coordinate
     * @param bc Boundary conditions
     * @param solution Final solution state
     * @return Heat flux components or error
     */
    [[nodiscard]] auto compute_heat_flux_analysis(int station, double xi,
                                                 const conditions::BoundaryConditions& bc,
                                                 const equations::SolutionState& solution) const
        -> std::expected<HeatFluxComponents, SolverError>;

    /**
     * @brief Print detailed heat flux summary
     * @param station Station number
     * @param components Heat flux components
     */
    static auto print_heat_flux_summary(int station, const HeatFluxComponents& components) -> void;

    /**
     * @brief Check if radiative equilibrium is enabled
     * @return True if wall mode is Radiative
     */
    [[nodiscard]] auto is_radiative_equilibrium_enabled() const noexcept -> bool {
        return config_.simulation.wall_mode == io::SimulationConfig::WallMode::Radiative;
    }

    /**
     * @brief Get environment temperature
     * @return Environment temperature [K]
     */
    [[nodiscard]] auto get_environment_temperature() const noexcept -> double {
        return config_.wall_parameters.environment_temperature;
    }

    /**
     * @brief Get surface emissivity
     * @return Surface emissivity [-]
     */
    [[nodiscard]] auto get_emissivity() const noexcept -> double {
        return config_.wall_parameters.emissivity;
    }

private:
    /**
     * @brief Calculate radiative heat flux using Stefan-Boltzmann law
     * @param T_wall Wall temperature [K]
     * @param T_infinity Environment temperature [K]
     * @param emissivity Surface emissivity [-]
     * @return Radiative heat flux [W/m²]
     */
    [[nodiscard]] static auto calculate_radiative_flux(double T_wall, double T_infinity, double emissivity) noexcept
        -> double;

    /**
     * @brief Validate radiative equilibrium inputs
     * @param q_wall Wall heat flux [W/m²]
     * @param emissivity Surface emissivity [-]
     * @param T_infinity Environment temperature [K]
     * @return Error message if invalid
     */
    // Note: validate_inputs method moved to InputValidator::validate_radiative_equilibrium_inputs

private:
    BoundaryLayerSolver& solver_;
    const io::Configuration& config_;
    
    static constexpr double STEFAN_BOLTZMANN = 5.670374419e-8; ///< Stefan-Boltzmann constant [W/(m²·K⁴)]
};

} // namespace blast::boundary_layer::solver