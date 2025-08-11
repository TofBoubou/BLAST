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
#include "solver_errors.hpp"
#include <expected>
#include <memory>

namespace blast::boundary_layer::solver {

// Forward declarations
class BoundaryLayerSolver;

// =============================================================================
// CORE STEP INTERFACES
// =============================================================================

struct SolverContext {
  equations::SolutionState& solution;
  equations::SolutionState& solution_old;
  conditions::BoundaryConditions& bc;
  coefficients::CoefficientSet& coeffs;
  const thermophysics::MixtureInterface& mixture;

  const int station;
  const double xi;
  const int iteration;

  // Services (references to avoid copying)
  BoundaryLayerSolver& solver;
  const grid::BoundaryLayerGrid& grid;
  coefficients::XiDerivatives& xi_derivatives;
};

class SolverStep {
public:
  virtual ~SolverStep() = default;
  virtual auto execute(SolverContext& ctx) -> std::expected<void, StepExecutionError> = 0;
  virtual auto name() const -> std::string_view = 0;
};

// =============================================================================
// CONCRETE STEP IMPLEMENTATIONS
// =============================================================================

class ThermodynamicConsistencyStep : public SolverStep {
private:
  thermodynamics::EnthalpyTemperatureSolver& h2t_solver_;

public:
  explicit ThermodynamicConsistencyStep(thermodynamics::EnthalpyTemperatureSolver& solver);
  auto execute(SolverContext& ctx) -> std::expected<void, StepExecutionError> override;
  auto name() const -> std::string_view override { return "ThermodynamicConsistency"; }
};

class MechanicalResolutionStep : public SolverStep {
private:
  const grid::BoundaryLayerGrid& grid_;
  coefficients::CoefficientCalculator& coeff_calculator_;

  auto solve_continuity(SolverContext& ctx) -> std::expected<std::vector<double>, SolverError>;

public:
  MechanicalResolutionStep(const grid::BoundaryLayerGrid& grid, coefficients::CoefficientCalculator& calc);
  auto execute(SolverContext& ctx) -> std::expected<void, StepExecutionError> override;
  auto name() const -> std::string_view override { return "MechanicalResolution"; }
};

class ThermalResolutionStep : public SolverStep {
public:
  auto execute(SolverContext& ctx) -> std::expected<void, StepExecutionError> override;
  auto name() const -> std::string_view override { return "ThermalResolution"; }
};

class ChemicalResolutionStep : public SolverStep {
public:
  auto execute(SolverContext& ctx) -> std::expected<void, StepExecutionError> override;
  auto name() const -> std::string_view override { return "ChemicalResolution"; }
};

class BoundaryEnforcementStep : public SolverStep {
public:
  auto execute(SolverContext& ctx) -> std::expected<void, StepExecutionError> override;
  auto name() const -> std::string_view override { return "BoundaryEnforcement"; }
};

class InputCalculationStep : public SolverStep {
private:
  std::unique_ptr<coefficients::CoefficientInputs> current_inputs_;

public:
  auto execute(SolverContext& ctx) -> std::expected<void, StepExecutionError> override;
  auto get_inputs() const -> const coefficients::CoefficientInputs&;
  auto name() const -> std::string_view override { return "InputCalculation"; }
};

// =============================================================================
// PIPELINE ORCHESTRATION
// =============================================================================

class SolverPipeline {
private:
  std::vector<std::unique_ptr<SolverStep>> steps_;

public:
  static auto create_for_solver(BoundaryLayerSolver& solver) -> SolverPipeline;
  auto execute_all(SolverContext& ctx) -> std::expected<void, StepExecutionError>;
};

} // namespace blast::boundary_layer::solver