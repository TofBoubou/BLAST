#pragma once
#include "../../core/exceptions.hpp"
#include "../../io/config_types.hpp"
#include "../../thermophysics/mixture_interface.hpp"
#include "../coefficients/coefficient_calculator.hpp"
#include "../coefficients/derivative_calculator.hpp"
#include "../coefficients/heat_flux_calculator.hpp"
#include "../coefficients/xi_derivatives.hpp"
#include "../conditions/boundary_interpolator.hpp"
#include "../equations/equation_types.hpp"
#include "../grid/grid.hpp"
#include "../thermodynamics/enthalpy_temperature_solver.hpp"
#include "adaptive_relaxation_controller.hpp"
#include "continuation_method.hpp"
#include "solver_errors.hpp"
#include "solver_steps.hpp"
#include "station_solver.hpp"
#include "convergence_manager.hpp"
#include "equation_solver.hpp"
#include "radiative_equilibrium_solver.hpp"
#include "heat_flux_computer.hpp"
#include "initial_guess_factory.hpp"
#include <expected>
#include <memory>
#include <string>
#include <vector>

namespace blast::boundary_layer::solver {

using DerivativeState = coefficients::UnifiedDerivativeState;

// Complete solution data for all stations
struct SolutionResult {
  std::vector<equations::SolutionState> stations;
  std::vector<double> xi_solved;
  std::vector<std::vector<double>> temperature_fields;
  bool converged = false;
  int total_iterations = 0;
  std::vector<coefficients::HeatFluxCoefficients> heat_flux_data;
  std::vector<std::vector<std::vector<double>>> modal_temperature_fields; // [n_stations][n_modes][n_eta]
  std::vector<std::string> temperature_mode_names;
};

// ConvergenceInfo is now defined in adaptive_relaxation_controller.hpp

class BoundaryLayerSolver {
private:
  const thermophysics::MixtureInterface& mixture_;
  io::Configuration config_;
  std::unique_ptr<grid::BoundaryLayerGrid> grid_;
  std::unique_ptr<coefficients::CoefficientCalculator> coeff_calculator_;
  std::unique_ptr<thermodynamics::EnthalpyTemperatureSolver> h2t_solver_;
  std::unique_ptr<coefficients::XiDerivatives> xi_derivatives_;
  std::unique_ptr<coefficients::DerivativeCalculator> derivative_calculator_;
  std::unique_ptr<coefficients::HeatFluxCalculator> heat_flux_calculator_;
  std::unique_ptr<ContinuationMethod> continuation_;
  io::Configuration original_config_;
  bool in_continuation_ = false;
  
  // Specialized solver components
  std::unique_ptr<StationSolver> station_solver_;
  std::unique_ptr<ConvergenceManager> convergence_manager_;
  std::unique_ptr<EquationSolver> equation_solver_;
  std::unique_ptr<RadiativeEquilibriumSolver> radiative_solver_;
  
  // Utility components for eliminating duplication
  std::unique_ptr<HeatFluxComputer> heat_flux_computer_;
  std::unique_ptr<InitialGuessFactory> initial_guess_factory_;

public:
  friend class ContinuationMethod;
  friend class StationSolver;
  friend class ConvergenceManager;
  friend class EquationSolver;
  friend class RadiativeEquilibriumSolver;

  explicit BoundaryLayerSolver(const thermophysics::MixtureInterface& mixture, const io::Configuration& config);

  [[nodiscard]] auto solve() -> std::expected<SolutionResult, SolverError>;

  // =============================================================================
  // PUBLIC ACCESSORS FOR PIPELINE STEPS
  // =============================================================================

  [[nodiscard]] auto get_h2t_solver() -> thermodynamics::EnthalpyTemperatureSolver& { return *h2t_solver_; }

  [[nodiscard]] auto get_grid() -> const grid::BoundaryLayerGrid& { return *grid_; }

  [[nodiscard]] auto get_coeff_calculator() -> coefficients::CoefficientCalculator& { return *coeff_calculator_; }
  
  [[nodiscard]] auto get_heat_flux_calculator() -> coefficients::HeatFluxCalculator& { return *heat_flux_calculator_; }
  
  [[nodiscard]] auto get_heat_flux_computer() -> HeatFluxComputer& { return *heat_flux_computer_; }
  
  [[nodiscard]] auto get_xi_derivatives() -> coefficients::XiDerivatives& { return *xi_derivatives_; }
  
  [[nodiscard]] auto get_convergence_manager() -> ConvergenceManager& { return *convergence_manager_; }
  
  [[nodiscard]] auto get_continuation() -> ContinuationMethod* { return continuation_.get(); }
  
  [[nodiscard]] auto get_original_config() const -> const io::Configuration& { return original_config_; }

  // =============================================================================
  // PUBLIC METHODS FOR STEP ACCESS
  // =============================================================================

  auto update_edge_properties(conditions::BoundaryConditions& bc, const coefficients::CoefficientInputs& inputs,
                              const core::Matrix<double>& species_matrix) const -> std::expected<void, SolverError>;

  auto enforce_edge_boundary_conditions(equations::SolutionState& solution,
                                        const conditions::BoundaryConditions& bc) const -> void;

  [[nodiscard]] auto compute_all_derivatives(const equations::SolutionState& solution) const
      -> std::expected<coefficients::UnifiedDerivativeState, SolverError>;

  [[nodiscard]] auto solve_momentum_equation(const equations::SolutionState& solution,
                                             const coefficients::CoefficientSet& coeffs,
                                             const conditions::BoundaryConditions& bc,
                                             double xi) -> std::expected<std::vector<double>, SolverError>;

  [[nodiscard]] auto solve_energy_equation(const equations::SolutionState& solution,
                                           const coefficients::CoefficientInputs& inputs,
                                           const coefficients::CoefficientSet& coeffs,
                                           const conditions::BoundaryConditions& bc,
                                           const thermophysics::MixtureInterface& mixture,
                                           int station) -> std::expected<std::vector<double>, SolverError>;

  [[nodiscard]] auto solve_species_equations(const equations::SolutionState& solution,
                                             const coefficients::CoefficientInputs& inputs,
                                             const coefficients::CoefficientSet& coeffs,
                                             const conditions::BoundaryConditions& bc,
                                             int station) -> std::expected<core::Matrix<double>, SolverError>;

  [[nodiscard]] auto update_temperature_field(std::span<const double> g_field, const core::Matrix<double>& composition,
                                              const conditions::BoundaryConditions& bc,
                                              std::span<const double> current_temperatures)
      -> std::expected<std::vector<double>, SolverError>;

  auto set_wall_temperature(double Tw) -> void {
    if (!config_.wall_parameters.wall_temperatures.empty()) {
      config_.wall_parameters.wall_temperatures[0] = Tw;
    } else {
      config_.wall_parameters.wall_temperatures.push_back(Tw);
    }
  }

  [[nodiscard]] auto get_wall_temperature() const -> double {
    if (!config_.wall_parameters.wall_temperatures.empty()) {
      return config_.wall_parameters.wall_temperatures[0];
    }
    return 300.0; // Default
  }

  [[nodiscard]] auto get_config() -> io::Configuration& { return config_; }

  [[nodiscard]] auto get_mixture() -> thermophysics::MixtureInterface& {
    return const_cast<thermophysics::MixtureInterface&>(mixture_);
  }

  [[nodiscard]] auto config() const noexcept -> const io::Configuration& { return config_; }

  auto set_config(const io::Configuration& config) noexcept -> void { config_ = config; }

private:
  // =============================================================================
  // PRIVATE ORCHESTRATION METHODS
  // =============================================================================

  [[nodiscard]] auto solve_station(int station, double xi, const equations::SolutionState& initial_guess)
      -> std::expected<equations::SolutionState, SolverError>;

  // Note: create_initial_guess method moved to InitialGuessFactory
                                          
  // Initialize specialized solver components
  auto initialize_solver_components() -> void;
};

} // namespace blast::boundary_layer::solver