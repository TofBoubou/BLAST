#pragma once
#include "../equations/equation_types.hpp"
#include "../coefficients/coefficient_calculator.hpp"
#include "../coefficients/xi_derivatives.hpp"
#include "../conditions/boundary_interpolator.hpp"
#include "../thermodynamics/enthalpy_temperature_solver.hpp"
#include "../grid/grid.hpp"
#include "../../thermophysics/mixture_interface.hpp"
#include "../../io/config_types.hpp"
#include "../../core/exceptions.hpp"
#include "adaptive_relaxation_controller.hpp"
#include <expected>
#include <memory>

namespace blast::boundary_layer::solver {

// Configuration for adaptive subdivision strategy
struct SubdivisionConfig {
    int max_depth = 10;                      // Maximum subdivision depth
    double min_step_fraction = 0.000000000000001;        // Minimum step as fraction of total domain
    int max_retry_attempts = 15;             // Maximum retry attempts per station
    bool enable_physics_simplification = true;  // Allow temporary physics simplification
    double subdivision_factor = 0.3;        // Factor for first subdivision attempt (xi_new = xi_start + factor*(xi_end - xi_start))
    
    // Factory methods for different strategies
    [[nodiscard]] static auto conservative() -> SubdivisionConfig {
        SubdivisionConfig config;
        config.max_depth = 2;
        config.subdivision_factor = 0.2;
        config.max_retry_attempts = 3;
        return config;
    }
    
    [[nodiscard]] static auto aggressive() -> SubdivisionConfig {
        SubdivisionConfig config;
        config.max_depth = 4;
        config.subdivision_factor = 0.5;
        config.max_retry_attempts = 8;
        return config;
    }
};

// Structure to hold both first and second derivatives
struct DerivativeState {
    core::Matrix<double> dc_deta;   // First derivatives [n_species x n_eta]
    core::Matrix<double> dc_deta2;  // Second derivatives [n_species x n_eta]
    
    DerivativeState(std::size_t n_species, std::size_t n_eta) 
        : dc_deta(n_species, n_eta), dc_deta2(n_species, n_eta) {}
};

// Complete solution data for all stations
struct SolutionResult {
    std::vector<equations::SolutionState> stations;
    std::vector<double> xi_solved;
    std::vector<std::vector<double>> temperature_fields;
    bool converged = false;
    int total_iterations = 0;
    
    // Subdivision statistics
    int total_subdivisions = 0;
    int max_subdivision_depth_used = 0;
    std::vector<std::pair<double, int>> subdivision_points; // (xi, depth) pairs
};

// Recovery strategy information
struct RecoveryInfo {
    enum class Strategy {
        None,
        Subdivision,
        RelaxationReduction,
        PhysicsSimplification,
        Extrapolation
    };
    
    Strategy strategy_used = Strategy::None;
    int attempts_made = 0;
    double final_xi_step = 0.0;
    std::string details;
};

// Error type for solver operations
class SolverError : public core::BlastException {
public:
    explicit SolverError(std::string_view message,
                        std::source_location location = std::source_location::current())
        : BlastException(std::format("Solver Error: {}", message), location) {}
        
    template<typename... Args>
    explicit SolverError(std::string_view format_str, 
                        std::source_location location,
                        Args&&... args) 
        : BlastException(std::format("Solver Error: {}", 
                        std::vformat(format_str, std::make_format_args(args...))), location) {}
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
    
    // Subdivision configuration
    SubdivisionConfig subdivision_config_;

public:
    explicit BoundaryLayerSolver(
        const thermophysics::MixtureInterface& mixture,
        const io::Configuration& config,
        const SubdivisionConfig& subdivision_config = SubdivisionConfig{}
    );
    
    // Main solving interface
    [[nodiscard]] auto solve() -> std::expected<SolutionResult, SolverError>;

private:
    // Adaptive station-level solving with subdivision support
    [[nodiscard]] auto solve_stations_adaptive() -> std::expected<SolutionResult, SolverError>;
    
    // Attempt to solve a specific xi interval with subdivision if needed
    [[nodiscard]] auto solve_xi_interval(
        std::size_t start_station_idx,
        std::size_t end_station_idx,
        const equations::SolutionState& initial_guess,
        SolutionResult& result,
        int depth = 0
    ) -> std::expected<RecoveryInfo, SolverError>;
    
    // Single station solving (existing logic)
    [[nodiscard]] auto solve_station(
        int station,
        double xi,
        const equations::SolutionState& initial_guess
    ) -> std::expected<equations::SolutionState, SolverError>;
    
    // Iterative solving at one station
    [[nodiscard]] auto iterate_station_adaptive(
        int station,
        double xi,
        const conditions::BoundaryConditions& bc,
        equations::SolutionState& solution
    ) -> std::expected<ConvergenceInfo, SolverError>;
    
    // Individual equation solving with proper sequencing
    auto solve_continuity_equation(
        const equations::SolutionState& solution,
        double xi
    ) -> std::expected<std::vector<double>, SolverError>;
    
    [[nodiscard]] auto solve_momentum_equation(
        const equations::SolutionState& solution,
        const coefficients::CoefficientSet& coeffs,
        const conditions::BoundaryConditions& bc,
        double xi
    ) -> std::expected<std::vector<double>, SolverError>;
    
    [[nodiscard]] auto solve_energy_equation(
        const equations::SolutionState& solution,
        const coefficients::CoefficientInputs& inputs,
        const coefficients::CoefficientSet& coeffs,
        const conditions::BoundaryConditions& bc,
        int station
    ) -> std::expected<std::vector<double>, SolverError>;
    
    [[nodiscard]] auto solve_species_equations(
        const equations::SolutionState& solution,
        const coefficients::CoefficientInputs& inputs,
        const coefficients::CoefficientSet& coeffs,
        const conditions::BoundaryConditions& bc,
        int station
    ) -> std::expected<core::Matrix<double>, SolverError>;
    
    // Temperature update from enthalpy
    [[nodiscard]] auto update_temperature_field(
        std::span<const double> g_field,
        const core::Matrix<double>& composition,
        const conditions::BoundaryConditions& bc,
        std::span<const double> current_temperatures
    ) -> std::expected<std::vector<double>, SolverError>;
    
    // Convergence checking
    [[nodiscard]] auto check_convergence(
        const equations::SolutionState& old_solution,
        const equations::SolutionState& new_solution
    ) const noexcept -> ConvergenceInfo;
    
    // Solution initialization and extrapolation
    [[nodiscard]] auto create_initial_guess(
        int station,
        double xi,
        const conditions::BoundaryConditions& bc
    ) const -> std::expected<equations::SolutionState, SolverError>;
    
    [[nodiscard]] auto extrapolate_from_previous(
        const equations::SolutionState& previous_solution,
        double xi_prev,
        double xi_current
    ) const -> equations::SolutionState;
    
    // Recovery strategies
    [[nodiscard]] auto attempt_subdivision_recovery(
        std::size_t start_station_idx,
        std::size_t end_station_idx,
        const equations::SolutionState& initial_guess,
        SolutionResult& result,
        int depth
    ) -> std::expected<RecoveryInfo, SolverError>;
    
    [[nodiscard]] auto attempt_relaxation_recovery(
        int station,
        double xi,
        const equations::SolutionState& initial_guess
    ) -> std::expected<equations::SolutionState, SolverError>;
    
    [[nodiscard]] auto attempt_physics_simplification_recovery(
        int station,
        double xi,
        const equations::SolutionState& initial_guess
    ) -> std::expected<equations::SolutionState, SolverError>;
    
    // Utility functions
    [[nodiscard]] auto apply_relaxation_differential(
        const equations::SolutionState& old_solution,
        const equations::SolutionState& new_solution,
        double base_factor
    ) const -> equations::SolutionState;
    
    [[nodiscard]] auto compute_eta_derivatives(
        const equations::SolutionState& solution
    ) const -> std::expected<equations::SolutionState, SolverError>;
    
    [[nodiscard]] auto compute_concentration_derivatives(
        const equations::SolutionState& solution
    ) const -> std::expected<DerivativeState, SolverError>;
    
    // Enforce edge boundary conditions (critical for convergence)
    auto enforce_edge_boundary_conditions(
        equations::SolutionState& solution,
        const conditions::BoundaryConditions& bc
    ) const -> void;
    
    // Dynamic edge properties update for thermodynamic consistency
    auto update_edge_properties(
        conditions::BoundaryConditions& bc,
        const coefficients::CoefficientInputs& inputs,
        const core::Matrix<double>& species_matrix
    ) const -> void;
    
    // Validation and diagnostics
    [[nodiscard]] auto validate_subdivision_feasibility(
        double xi_start, 
        double xi_end, 
        int depth
    ) const noexcept -> bool;
    
    auto log_subdivision_attempt(
        double xi_start, 
        double xi_end, 
        double xi_new, 
        int depth
    ) const -> void;
};

} // namespace blast::boundary_layer::solver