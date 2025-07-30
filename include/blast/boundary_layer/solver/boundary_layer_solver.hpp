#pragma once
#include "../../core/exceptions.hpp"
#include "../../io/config_types.hpp"
#include "../../thermophysics/mixture_interface.hpp"
#include "../coefficients/coefficient_calculator.hpp"
#include "../coefficients/xi_derivatives.hpp"
#include "../conditions/boundary_interpolator.hpp"
#include "../equations/equation_types.hpp"
#include "../grid/grid.hpp"
#include "../thermodynamics/enthalpy_temperature_solver.hpp"
#include "adaptive_relaxation_controller.hpp"
#include "solver_steps.hpp"
#include <expected>
#include <memory>

namespace blast::boundary_layer::solver {

// Structure to hold both first and second derivatives
struct DerivativeState {
  core::Matrix<double> dc_deta;  // First derivatives [n_species x n_eta]
  core::Matrix<double> dc_deta2; // Second derivatives [n_species x n_eta]

  DerivativeState(std::size_t n_species, std::size_t n_eta) : dc_deta(n_species, n_eta), dc_deta2(n_species, n_eta) {}
};

// Complete solution data for all stations
struct SolutionResult {
  std::vector<equations::SolutionState> stations;
  std::vector<double> xi_solved;
  std::vector<std::vector<double>> temperature_fields;
  bool converged = false;
  int total_iterations = 0;
};

// ConvergenceInfo is now defined in adaptive_relaxation_controller.hpp

// Error type for solver operations
class SolverError : public core::BlastException {
public:
  explicit SolverError(std::string_view message, std::source_location location = std::source_location::current())
      : BlastException(std::format("Solver Error: {}", message), location) {}

  template <typename... Args>
  explicit SolverError(std::string_view format_str, std::source_location location, Args&&... args)
      : BlastException(std::format("Solver Error: {}", std::vformat(format_str, std::make_format_args(args...))),
                       location) {}
};

class BoundaryLayerSolver {
private:
  const thermophysics::MixtureInterface& mixture_;
  const io::Configuration& config_;
  std::unique_ptr<grid::BoundaryLayerGrid> grid_;
  std::unique_ptr<coefficients::CoefficientCalculator> coeff_calculator_;
  std::unique_ptr<thermodynamics::EnthalpyTemperatureSolver> h2t_solver_;
  std::unique_ptr<coefficients::XiDerivatives> xi_derivatives_;
  std::unique_ptr<AdaptiveRelaxationController> relaxation_controller_;

public:
  explicit BoundaryLayerSolver(const thermophysics::MixtureInterface& mixture, const io::Configuration& config);

  [[nodiscard]] auto solve() -> std::expected<SolutionResult, SolverError>;

  // =============================================================================
  // PUBLIC ACCESSORS FOR PIPELINE STEPS
  // =============================================================================
  
  [[nodiscard]] auto get_h2t_solver() -> thermodynamics::EnthalpyTemperatureSolver& { 
      return *h2t_solver_; 
  }
  
  [[nodiscard]] auto get_grid() -> const grid::BoundaryLayerGrid& { 
      return *grid_; 
  }
  
  [[nodiscard]] auto get_coeff_calculator() -> coefficients::CoefficientCalculator& { 
      return *coeff_calculator_; 
  }

  // =============================================================================
  // PUBLIC METHODS FOR STEP ACCESS
  // =============================================================================
  
  auto update_edge_properties(conditions::BoundaryConditions& bc, 
                             const coefficients::CoefficientInputs& inputs,
                             const core::Matrix<double>& species_matrix) const -> void;
  
  auto enforce_edge_boundary_conditions(equations::SolutionState& solution,
                                      const conditions::BoundaryConditions& bc) const -> void;
  
  [[nodiscard]] auto compute_eta_derivatives(const equations::SolutionState& solution) const
      -> std::expected<equations::SolutionState, SolverError>;
  
  [[nodiscard]] auto compute_concentration_derivatives(const equations::SolutionState& solution) const
      -> std::expected<DerivativeState, SolverError>;
  
  [[nodiscard]] auto solve_momentum_equation(const equations::SolutionState& solution,
                                           const coefficients::CoefficientSet& coeffs,
                                           const conditions::BoundaryConditions& bc,
                                           double xi) -> std::expected<std::vector<double>, SolverError>;
  
  [[nodiscard]] auto solve_energy_equation(const equations::SolutionState& solution,
                                         const coefficients::CoefficientInputs& inputs,
                                         const coefficients::CoefficientSet& coeffs,
                                         const conditions::BoundaryConditions& bc,
                                         int station) -> std::expected<std::vector<double>, SolverError>;
  
  [[nodiscard]] auto solve_species_equations(const equations::SolutionState& solution,
                                           const coefficients::CoefficientInputs& inputs,
                                           const coefficients::CoefficientSet& coeffs,
                                           const conditions::BoundaryConditions& bc,
                                           int station) -> std::expected<core::Matrix<double>, SolverError>;

private:
  // =============================================================================
  // PRIVATE ORCHESTRATION METHODS
  // =============================================================================
  
  [[nodiscard]] auto solve_station(int station, double xi, const equations::SolutionState& initial_guess)
      -> std::expected<equations::SolutionState, SolverError>;

  [[nodiscard]] auto iterate_station_adaptive(int station, double xi, 
                                             const conditions::BoundaryConditions& bc,
                                             equations::SolutionState& solution)
      -> std::expected<ConvergenceInfo, SolverError>;

  [[nodiscard]] auto update_temperature_field(std::span<const double> g_field, 
                                             const core::Matrix<double>& composition,
                                             const conditions::BoundaryConditions& bc,
                                             std::span<const double> current_temperatures)
      -> std::expected<std::vector<double>, SolverError>;

  [[nodiscard]] auto check_convergence(const equations::SolutionState& old_solution,
                                     const equations::SolutionState& new_solution) const noexcept -> ConvergenceInfo;

  [[nodiscard]] auto create_initial_guess(int station, double xi, const conditions::BoundaryConditions& bc,
                                        double T_edge) const -> std::expected<equations::SolutionState, SolverError>;

  [[nodiscard]] auto apply_relaxation_differential(const equations::SolutionState& old_solution,
                                                 const equations::SolutionState& new_solution,
                                                 double base_factor) const -> equations::SolutionState;
};

} // namespace blast::boundary_layer::solver