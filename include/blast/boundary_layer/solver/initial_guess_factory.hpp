#pragma once

#include "blast/boundary_layer/conditions/boundary_conditions.hpp"
#include "blast/boundary_layer/equations/equation_types.hpp"
#include "blast/boundary_layer/grid/grid.hpp"
#include "blast/boundary_layer/solver/solver_errors.hpp"
#include "blast/thermophysics/mixture_interface.hpp"

#include <expected>

namespace blast::boundary_layer::solver {

/**
 * @brief Factory for creating initial guess solution states
 *
 * This class eliminates the duplication of create_initial_guess methods
 * that appeared in both BoundaryLayerSolver and StationSolver with identical
 * logic.
 */
class InitialGuessFactory {
public:
  InitialGuessFactory(const grid::BoundaryLayerGrid &grid,
                      const thermophysics::MixtureInterface &mixture) noexcept;

  /**
   * @brief Create initial guess for a station
   *
   * This replaces the duplicated create_initial_guess methods with unified
   * logic. Handles both single-species and multi-species cases.
   *
   * @param station Station number
   * @param xi Streamwise coordinate
   * @param bc Boundary conditions
   * @param T_edge Edge temperature
   * @return Initial guess solution state or error
   */
  [[nodiscard]] auto
  create_initial_guess(int station, double xi,
                       const conditions::BoundaryConditions &bc,
                       double T_edge) const
      -> std::expected<equations::SolutionState, SolverError>;

  auto create_initial_guess_from_previous(
      const std::vector<double> &prev_F, const std::vector<double> &prev_g,
      const std::vector<double> &prev_T, const core::Matrix<double> &prev_c,
      const conditions::BoundaryConditions &bc, double xi) const
      -> equations::SolutionState;

private:
  /**
   * @brief Create simple single-species initial guess
   */
  [[nodiscard]] auto
  create_single_species_guess(const conditions::BoundaryConditions &bc,
                              double T_edge) const
      -> std::expected<equations::SolutionState, SolverError>;

  /**
   * @brief Create multi-species initial guess with equilibrium composition
   */
  [[nodiscard]] auto
  create_multi_species_guess(const conditions::BoundaryConditions &bc,
                             double T_edge) const
      -> std::expected<equations::SolutionState, SolverError>;

  /**
   * @brief Smooth transition function for species interpolation
   */
  [[nodiscard]] static auto
  transition_function(double eta_norm, double eta_center = 0.35,
                      double sharpness = 10.0) noexcept -> double;

  const grid::BoundaryLayerGrid &grid_;
  const thermophysics::MixtureInterface &mixture_;
};

} // namespace blast::boundary_layer::solver