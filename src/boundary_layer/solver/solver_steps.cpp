#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/boundary_layer/equations/continuity.hpp"
#include "blast/boundary_layer/solver/solver_steps.hpp"
#include <iomanip>
#include <iostream>

namespace blast::boundary_layer::solver {

// =============================================================================
// THERMODYNAMIC CONSISTENCY STEP
// =============================================================================

ThermodynamicConsistencyStep::ThermodynamicConsistencyStep(thermodynamics::EnthalpyTemperatureSolver& solver) 
    : h2t_solver_(solver) {}

auto ThermodynamicConsistencyStep::execute(SolverContext& ctx) -> StepResult {
    // 1. Temperature update from enthalpy
    std::vector<double> enthalpy_field(ctx.solution.g.size());
    for (std::size_t i = 0; i < enthalpy_field.size(); ++i) {
        enthalpy_field[i] = ctx.solution.g[i] * ctx.bc.he();
    }

    auto T_result = h2t_solver_.solve(enthalpy_field, ctx.solution.c, ctx.bc, ctx.solution.T);
    if (!T_result) {
        std::cerr << "Temperature update failed: " << T_result.error().message() << std::endl;
        return StepResult::Failed;
    }
    
    ctx.solution.T = std::move(T_result.value().temperatures);
    
    // 2. Update edge properties for thermodynamic consistency
    auto dummy_inputs = coefficients::CoefficientInputs{
        .xi = ctx.xi,
        .F = ctx.solution.F,
        .c = ctx.solution.c,
        .dc_deta = core::Matrix<double>(),
        .dc_deta2 = core::Matrix<double>(),
        .T = ctx.solution.T
    };
    ctx.solver.update_edge_properties(ctx.bc, dummy_inputs, ctx.solution.c);
    
    return StepResult::Success;
}

// =============================================================================
// MECHANICAL RESOLUTION STEP
// =============================================================================

MechanicalResolutionStep::MechanicalResolutionStep(const grid::BoundaryLayerGrid& grid, 
                                                 coefficients::CoefficientCalculator& calc)
    : grid_(grid), coeff_calculator_(calc) {}

auto MechanicalResolutionStep::execute(SolverContext& ctx) -> StepResult {
    // 1. Continuity equation
    auto V_result = solve_continuity(ctx);
    if (!V_result) {
        std::cerr << "Continuity equation failed: " << V_result.error().message() << std::endl;
        return StepResult::Failed;
    }
    ctx.solution.V = std::move(V_result.value());
    
    // 2. Compute all derivatives after V update
    auto all_derivatives_result = ctx.solver.compute_all_derivatives(ctx.solution);
    if (!all_derivatives_result) {
        std::cerr << "All derivatives computation failed: " << all_derivatives_result.error().message() << std::endl;
        return StepResult::Failed;
    }
    auto all_derivatives = all_derivatives_result.value();
    
    // 3. Create updated inputs with consistent state
    auto updated_inputs = coefficients::CoefficientInputs{
        .xi = ctx.xi,
        .F = ctx.solution.F,
        .c = ctx.solution.c,
        .dc_deta = all_derivatives.dc_deta,
        .dc_deta2 = all_derivatives.dc_deta2,
        .T = ctx.solution.T
    };
    
    // 4. Calculate coefficients
    auto coeffs_result = coeff_calculator_.calculate(updated_inputs, ctx.bc, ctx.xi_derivatives);
    if (!coeffs_result) {
        std::cerr << "Coefficient calculation failed: " << coeffs_result.error().message() << std::endl;
        return StepResult::Failed;
    }
    ctx.coeffs = std::move(coeffs_result.value());
    
    // 5. Solve momentum equation
    auto F_result = ctx.solver.solve_momentum_equation(ctx.solution, ctx.coeffs, ctx.bc, ctx.xi);
    if (!F_result) {
        std::cerr << "Momentum equation failed: " << F_result.error().message() << std::endl;
        return StepResult::Failed;
    }
    auto F_new = F_result.value();
    
    // Apply momentum boundary conditions
    if (!F_new.empty()) {
        F_new.back() = 1.0;
    }
    ctx.solution.F = std::move(F_new);
    
    return StepResult::Success;
}

auto MechanicalResolutionStep::solve_continuity(SolverContext& ctx) -> std::expected<std::vector<double>, SolverError> {
    const auto n_eta = grid_.n_eta();
    const double d_eta = grid_.d_eta();
    const double lambda0 = ctx.xi_derivatives.lambda0();
    const auto F_derivatives = ctx.xi_derivatives.F_derivative();
    
    std::vector<double> y_field(n_eta);
    for (std::size_t i = 0; i < n_eta; ++i) {
        y_field[i] = -(2.0 * ctx.xi * lambda0 + 1.0) * ctx.solution.F[i] - 2.0 * ctx.xi * F_derivatives[i];
    }
    
    return equations::solve_continuity(y_field, d_eta, 0.0);
}

// =============================================================================
// THERMAL RESOLUTION STEP
// =============================================================================

auto ThermalResolutionStep::execute(SolverContext& ctx) -> StepResult {
    auto all_derivatives_result = ctx.solver.compute_all_derivatives(ctx.solution);
    if (!all_derivatives_result) {
        std::cerr << "All derivatives computation failed in ThermalResolution" << std::endl;
        return StepResult::Failed;
    }
    auto derivatives = all_derivatives_result.value();

    auto thermal_inputs = coefficients::CoefficientInputs{
        .xi = ctx.xi,
        .F = ctx.solution.F,
        .c = ctx.solution.c,
        .dc_deta = derivatives.dc_deta,
        .dc_deta2 = derivatives.dc_deta2,
        .T = ctx.solution.T
    };

    auto g_result = ctx.solver.solve_energy_equation(ctx.solution, thermal_inputs, ctx.coeffs, ctx.bc, ctx.station);
    if (!g_result) {
        std::cerr << "Energy equation failed: " << g_result.error().message() << std::endl;
        return StepResult::Failed;
    }
    
    auto g_new = g_result.value();
    
    // Apply energy boundary conditions
    if (!g_new.empty()) {
        g_new.back() = 1.0;
    }
    ctx.solution.g = std::move(g_new);
    
    return StepResult::Success;
}

// =============================================================================
// CHEMICAL RESOLUTION STEP
// =============================================================================

auto ChemicalResolutionStep::execute(SolverContext& ctx) -> StepResult {
    auto all_derivatives_result = ctx.solver.compute_all_derivatives(ctx.solution);
    if (!all_derivatives_result) {
        std::cerr << "All derivatives computation failed in ThermalResolution" << std::endl;
        return StepResult::Failed;
    }
    auto derivatives = all_derivatives_result.value();

    auto chemical_inputs = coefficients::CoefficientInputs{
        .xi = ctx.xi,
        .F = ctx.solution.F,
        .c = ctx.solution.c,
        .dc_deta = derivatives.dc_deta,
        .dc_deta2 = derivatives.dc_deta2,
        .T = ctx.solution.T
    };

    auto c_result = ctx.solver.solve_species_equations(ctx.solution, chemical_inputs, ctx.coeffs, ctx.bc, ctx.station);
    if (!c_result) {
        std::cerr << "Species equations failed: " << c_result.error().message() << std::endl;
        return StepResult::Failed;
    }
    
    ctx.solution.c = std::move(c_result.value());
    
    return StepResult::Success;
}

// =============================================================================
// INPUT CALCULATION STEP  
// =============================================================================

auto InputCalculationStep::execute(SolverContext& ctx) -> StepResult {
    auto all_derivatives_result = ctx.solver.compute_all_derivatives(ctx.solution);
    if (!all_derivatives_result) {
        std::cerr << "All derivatives computation failed in ThermalResolution" << std::endl;
        return StepResult::Failed;
    }
    auto derivatives = all_derivatives_result.value();
    
    current_inputs_ = std::make_unique<coefficients::CoefficientInputs>(
        coefficients::CoefficientInputs{
            .xi = ctx.xi,
            .F = ctx.solution.F,
            .c = ctx.solution.c,
            .dc_deta = derivatives.dc_deta,
            .dc_deta2 = derivatives.dc_deta2,
            .T = ctx.solution.T
        }
    );
    return StepResult::Success;
}

auto InputCalculationStep::get_inputs() const -> const coefficients::CoefficientInputs& {
    return *current_inputs_;
}

// =============================================================================
// BOUNDARY ENFORCEMENT STEP
// =============================================================================

auto BoundaryEnforcementStep::execute(SolverContext& ctx) -> StepResult {
    ctx.solver.enforce_edge_boundary_conditions(ctx.solution, ctx.bc);
    return StepResult::Success;
}

// =============================================================================
// SOLVER PIPELINE
// =============================================================================

auto SolverPipeline::create_for_solver(BoundaryLayerSolver& solver) -> SolverPipeline {
    SolverPipeline pipeline;
    
    pipeline.steps_.push_back(
        std::make_unique<ThermodynamicConsistencyStep>(solver.get_h2t_solver())
    );
    pipeline.steps_.push_back(
        std::make_unique<MechanicalResolutionStep>(solver.get_grid(), solver.get_coeff_calculator())
    );
    pipeline.steps_.push_back(
        std::make_unique<ThermalResolutionStep>()
    );
    pipeline.steps_.push_back(
        std::make_unique<ChemicalResolutionStep>()
    );
    pipeline.steps_.push_back(
        std::make_unique<BoundaryEnforcementStep>()
    );
    
    return pipeline;
}

auto SolverPipeline::execute_all(SolverContext& ctx) -> StepResult {
    last_failed_step_.clear(); // Reset previous failure state
    
    for (auto& step : steps_) {
        // std::cout << "Executing step: " << step->name() << std::endl;
        
        auto result = step->execute(ctx);
        if (result == StepResult::Failed) {
            last_failed_step_ = step->name();
            std::cerr << "Step failed: " << step->name() << std::endl;
            return StepResult::Failed;
        }
    }
    return StepResult::Success;
}

} // namespace blast::boundary_layer::solver