#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include "blast/boundary_layer/equations/continuity.hpp"
#include "blast/boundary_layer/equations/energy.hpp"
#include "blast/boundary_layer/equations/momentum.hpp"
#include "blast/boundary_layer/equations/species.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

namespace blast::boundary_layer::solver {

BoundaryLayerSolver::BoundaryLayerSolver(const thermophysics::MixtureInterface& mixture,
                                         const io::Configuration& config)
    : mixture_(mixture), config_(config) {

  // Create grid
  if (config.simulation.only_stagnation_point) {
    auto grid_result = grid::BoundaryLayerGrid::create_stagnation_grid(config.numerical, config.outer_edge);
    if (!grid_result) {
      throw SolverError("Failed to create stagnation grid: {}", std::source_location::current(),
                        grid_result.error().message());
    }
    grid_ = std::make_unique<grid::BoundaryLayerGrid>(std::move(grid_result.value()));
  } else {
    auto grid_result =
        grid::BoundaryLayerGrid::create_downstream_grid(config.numerical, config.outer_edge, config.output);
    if (!grid_result) {
      throw SolverError("Failed to create downstream grid: {}", std::source_location::current(),
                        grid_result.error().message());
    }
    grid_ = std::make_unique<grid::BoundaryLayerGrid>(std::move(grid_result.value()));
  }

  // Create coefficient calculator
  coeff_calculator_ =
      std::make_unique<coefficients::CoefficientCalculator>(mixture_, config_.simulation, config_.numerical);

  // Create enthalpy-temperature solver
  thermodynamics::EnthalpyTemperatureSolverConfig h2t_config{.tolerance = config_.numerical.solvers.h2t_tolerance,
                                                             .max_iterations =
                                                                 config_.numerical.solvers.h2t_max_iterations};
  h2t_solver_ = std::make_unique<thermodynamics::EnthalpyTemperatureSolver>(mixture_, h2t_config);

  // Create xi derivatives manager
  xi_derivatives_ = std::make_unique<coefficients::XiDerivatives>(grid_->n_eta(), mixture_.n_species());

  relaxation_controller_ =
      std::make_unique<AdaptiveRelaxationController>(AdaptiveRelaxationController::Config::for_stagnation_point());
}

auto BoundaryLayerSolver::solve() -> std::expected<SolutionResult, SolverError> {

  SolutionResult result;
  result.xi_solved.reserve(grid_->xi_coordinates().size());
  result.stations.reserve(grid_->xi_coordinates().size());
  result.temperature_fields.reserve(grid_->xi_coordinates().size());

  const auto xi_stations = grid_->xi_coordinates();

  // Variables to store previous station results for xi derivatives
  double prev_xi = 0.0;
  std::vector<double> prev_F, prev_g;
  core::Matrix<double> prev_c;

  // Solve each xi station
  for (std::size_t station_idx = 0; station_idx < xi_stations.size(); ++station_idx) {
    const int station = static_cast<int>(station_idx);
    const double xi = xi_stations[station_idx];

    // CRITICAL: Update xi derivatives BEFORE solving (except for station 0)
    if (station_idx > 0) {
      xi_derivatives_->update_station(station, xi, prev_F, prev_g, prev_c);
    }

    // Create initial guess for this station
    equations::SolutionState initial_guess;
    if (station == 0 || result.stations.empty()) {
      // Create boundary conditions
      auto bc_result =
          conditions::create_stagnation_conditions(config_.outer_edge, config_.wall_parameters, config_.simulation);
      if (!bc_result) {
        return std::unexpected(SolverError("Failed to create boundary conditions for station {}: {}",
                                           std::source_location::current(), station, bc_result.error().message()));
      }

      auto guess_result = create_initial_guess(station, xi, bc_result.value());
      if (!guess_result) {
        return std::unexpected(guess_result.error());
      }
      initial_guess = std::move(guess_result.value());
    } else {
      // Extrapolate from previous station
      initial_guess = extrapolate_from_previous(result.stations.back(), result.xi_solved.back(), xi);
    }

    // Solve this station
    auto station_result = solve_station(station, xi, initial_guess);
    if (!station_result) {
      return std::unexpected(SolverError("Failed to solve station {} (xi={}): {}", std::source_location::current(),
                                         station, xi, station_result.error().message()));
    }

    // Store results
    result.xi_solved.push_back(xi);
    result.stations.push_back(std::move(station_result.value()));
    result.temperature_fields.push_back(result.stations.back().T);

    // Store current results for next iteration
    prev_xi = xi;
    prev_F = result.stations.back().F;
    prev_g = result.stations.back().g;
    prev_c = result.stations.back().c;

    result.total_iterations++;
  }

  result.converged = true; // All stations solved successfully
  return result;
}

auto BoundaryLayerSolver::solve_station(int station, double xi, const equations::SolutionState& initial_guess)
    -> std::expected<equations::SolutionState, SolverError> {

  // Get boundary conditions for this station
  auto bc_result =
      (station == 0)
          ? conditions::create_stagnation_conditions(config_.outer_edge, config_.wall_parameters, config_.simulation)
          : conditions::interpolate_boundary_conditions(station, xi, grid_->xi_coordinates(), config_.outer_edge,
                                                        config_.wall_parameters, config_.simulation);

  if (!bc_result) {
    return std::unexpected(SolverError("Failed to get boundary conditions for station {}: {}",
                                       std::source_location::current(), station, bc_result.error().message()));
  }

  auto bc = bc_result.value();
  auto solution = initial_guess;

  if (station == 0) {
    relaxation_controller_ =
        std::make_unique<AdaptiveRelaxationController>(AdaptiveRelaxationController::Config::for_stagnation_point());
  } else {
    relaxation_controller_ =
        std::make_unique<AdaptiveRelaxationController>(AdaptiveRelaxationController::Config::for_downstream_station());
  }

  auto convergence_result = iterate_station_adaptive(station, xi, bc, solution);
  if (!convergence_result) {
    return std::unexpected(convergence_result.error());
  }

  const auto& conv_info = convergence_result.value();
  if (!conv_info.converged) {
    return std::unexpected(SolverError("Station {} failed to converge after {} iterations (residual={})",
                                       std::source_location::current(), station, conv_info.iterations,
                                       conv_info.max_residual()));
  }

  return solution;
}

/* auto BoundaryLayerSolver::iterate_station_adaptive(int station, double xi, const conditions::BoundaryConditions& bc,
                                                   equations::SolutionState& solution)
    -> std::expected<ConvergenceInfo, SolverError> {

  const auto n_eta = grid_->n_eta();
  const auto n_species = mixture_.n_species();

  auto bc_dynamic = bc;
  std::ranges::fill(solution.T, bc_dynamic.Tw());

  ConvergenceInfo conv_info;

  for (int iter = 0; iter < config_.numerical.max_iterations; ++iter) {
    std::cout << "=== ITERATION " << iter << " at station " << station << " ===" << std::endl;
    const auto solution_old = solution;

    auto T_result = update_temperature_field(solution.g, solution.c, bc_dynamic, solution.T);
    if (!T_result) {
      return std::unexpected(SolverError("Temperature update failed at station {} iteration {}: {}",
                                         std::source_location::current(), station, iter, T_result.error().message()));
    }
    solution.T = T_result.value();

    // 2. Compute eta derivatives
    auto derivatives_result = compute_eta_derivatives(solution);
    if (!derivatives_result) {
      return std::unexpected(derivatives_result.error());
    }
    auto solution_with_derivatives = derivatives_result.value();

    // 2b. Compute concentration derivatives
    auto conc_derivatives_result = compute_concentration_derivatives(solution);
    if (!conc_derivatives_result) {
      return std::unexpected(conc_derivatives_result.error());
    }
    auto concentration_derivatives = conc_derivatives_result.value();

    // 3. Create coefficient inputs
    coefficients::CoefficientInputs inputs{.xi = xi,
                                           .F = solution.F,
                                           .c = solution.c,
                                           .dc_deta = concentration_derivatives.dc_deta,
                                           .dc_deta2 = concentration_derivatives.dc_deta2,
                                           .T = solution.T};

    // 4. Calculate coefficients
    auto coeffs_result = coeff_calculator_->calculate(inputs, bc_dynamic, *xi_derivatives_);
    if (!coeffs_result) {
      return std::unexpected(SolverError("Coefficient calculation failed at station {} iteration {}: {}",
                                         std::source_location::current(), station, iter,
                                         coeffs_result.error().message()));
    }
    auto coeffs = coeffs_result.value();

    // 1. Solve continuity equation
    auto V_result = solve_continuity_equation(solution, xi);
    if (!V_result) {
      return std::unexpected(SolverError("Continuity solver failed at station {} iteration {}: {}",
                                         std::source_location::current(), station, iter, V_result.error().message()));
    }
    solution.V = std::move(V_result.value());

    // 5. Solve momentum equation
    auto F_result = solve_momentum_equation(solution, coeffs, bc_dynamic, xi);
    if (!F_result)
      return std::unexpected(F_result.error());
    auto F_new = F_result.value();
    if (!F_new.empty())
      F_new.back() = 1.0;
    std::cout << "Verification f[10] juste après le calcul : " << F_new[10] << std::endl;

    // 6. Solve energy equation
    auto g_result = solve_energy_equation(solution, inputs, coeffs, bc_dynamic, station);
    if (!g_result)
      return std::unexpected(g_result.error());
    auto g_new = g_result.value();
    std::cout << "Print de g après résolution : " << g_new[0] << "  " << g_new[10] << "  " << g_new[19] << std::endl;
    if (!g_new.empty())
      g_new.back() = 1.0;

    // 10. Solve species equations
    auto c_result = solve_species_equations(solution, inputs, coeffs, bc_dynamic, station);
    if (!c_result) {
      return std::unexpected(c_result.error());
    }
    auto c_new = c_result.value();

    // 11. Build new solution state
    equations::SolutionState solution_new(n_eta, n_species);
    solution_new.V = solution.V;
    solution_new.F = std::move(F_new);
    solution_new.g = std::move(g_new);
    solution_new.c = std::move(c_new);
    solution_new.T = solution.T;

    // 11.5. Enforce BC
    enforce_edge_boundary_conditions(solution_new, bc_dynamic);

    // 12. Check convergence
    conv_info = check_convergence(solution_old, solution_new);
    conv_info.iterations = iter + 1;

    // 13. Adapt relaxation
    double adaptive_factor = relaxation_controller_->adapt_relaxation_factor(conv_info, iter);
    std::cout << "Adaptive relaxation factor: " << std::scientific << std::setprecision(3) << adaptive_factor
              << " (max residual: " << conv_info.max_residual() << ")" << std::endl;

    // 14. Apply relaxation
    solution = apply_relaxation_differential(solution_old, solution_new, adaptive_factor);

    // 15. Exit if converged
    if (conv_info.converged) {
      std::cout << "✓ Converged with adaptive factor: " << adaptive_factor << std::endl;
      break;
    }
  }

  return conv_info;
} */

auto BoundaryLayerSolver::iterate_station_adaptive(
    int station, 
    double xi, 
    const conditions::BoundaryConditions& bc,
    equations::SolutionState& solution
) -> std::expected<ConvergenceInfo, SolverError> {

    const auto n_eta = grid_->n_eta();
    const auto n_species = mixture_.n_species();

    auto bc_dynamic = bc;
    ConvergenceInfo conv_info;

    for (int iter = 0; iter < config_.numerical.max_iterations; ++iter) {
        std::cout << "=== ITERATION " << iter << " at station " << station << " ===" << std::endl;
        const auto solution_old = solution;

        // ================================================================
        // PHASE 1: ESTABLISH BASIC THERMODYNAMIC CONSISTENCY
        // ================================================================
        
        // 1. Temperature field update for T-g-c consistency
        auto T_result = update_temperature_field(solution.g, solution.c, bc_dynamic, solution.T);
        if (!T_result) {
            return std::unexpected(SolverError(
                "Temperature update failed at station {} iteration {}: {}",
                std::source_location::current(), station, iter, T_result.error().message()
            ));
        }
        solution.T = T_result.value();

        // 2. Dynamic update of boundary properties for thermodynamic consistency
        update_edge_properties(bc_dynamic, 
            coefficients::CoefficientInputs{.xi = xi, .F = solution.F, .c = solution.c, 
                                          .dc_deta = core::Matrix<double>(), .dc_deta2 = core::Matrix<double>(), 
                                          .T = solution.T}, 
            solution.c);

        // ================================================================
        // PHASE 2: MECHANICAL SYSTEM RESOLUTION (KINEMATICS + DYNAMICS)
        // ================================================================
        
        // 3. Continuity resolution BEFORE derivative calculation (V depends on F and ∂F/∂ξ)
        auto V_result = solve_continuity_equation(solution, xi);
        if (!V_result) {
            return std::unexpected(SolverError(
                "Continuity solver failed at station {} iteration {}: {}",
                std::source_location::current(), station, iter, V_result.error().message()
            ));
        }
        solution.V = std::move(V_result.value());

        // 4. Spatial derivatives calculation AFTER V update
        auto derivatives_result = compute_eta_derivatives(solution);
        if (!derivatives_result) {
            return std::unexpected(derivatives_result.error());
        }

        auto conc_derivatives_result = compute_concentration_derivatives(solution);
        if (!conc_derivatives_result) {
            return std::unexpected(conc_derivatives_result.error());
        }
        auto concentration_derivatives = conc_derivatives_result.value();

        // 5. Input construction with consistent state
        coefficients::CoefficientInputs inputs{
            .xi = xi,
            .F = solution.F,
            .c = solution.c,
            .dc_deta = concentration_derivatives.dc_deta,
            .dc_deta2 = concentration_derivatives.dc_deta2,
            .T = solution.T
        };

        // 6. Coefficient calculation with consistent derivatives and state
        auto coeffs_result = coeff_calculator_->calculate(inputs, bc_dynamic, *xi_derivatives_);
        if (!coeffs_result) {
            return std::unexpected(SolverError(
                "Coefficient calculation failed at station {} iteration {}: {}",
                std::source_location::current(), station, iter, coeffs_result.error().message()
            ));
        }
        auto coeffs = coeffs_result.value();

        // 7. Momentum resolution with accurate coefficients
        auto F_result = solve_momentum_equation(solution, coeffs, bc_dynamic, xi);
        if (!F_result) {
            return std::unexpected(F_result.error());
        }
        auto F_new = F_result.value();
        
        // Immediate application of momentum boundary conditions
        if (!F_new.empty()) {
            F_new.back() = 1.0;  // Dimensionless velocity = 1 at edge
        }

        // ================================================================
        // PHASE 3: THERMAL SYSTEM RESOLUTION
        // ================================================================
        
        // 8. Solution update for energy equation
        solution.F = std::move(F_new);

        // 9. Energy resolution with updated mechanical state
        auto g_result = solve_energy_equation(solution, inputs, coeffs, bc_dynamic, station);
        if (!g_result) {
            return std::unexpected(g_result.error());
        }
        auto g_new = g_result.value();
        
        // Immediate application of energy boundary conditions
        if (!g_new.empty()) {
            g_new.back() = 1.0;  // Dimensionless enthalpy = 1 at edge
        }

        // ================================================================
        // PHASE 4: CHEMICAL SYSTEM RESOLUTION
        // ================================================================
        
        // 10. Solution update for species equations
        solution.g = std::move(g_new);

        // 11. Species resolution with updated thermomechanical state
        auto c_result = solve_species_equations(solution, inputs, coeffs, bc_dynamic, station);
        if (!c_result) {
            return std::unexpected(c_result.error());
        }
        auto c_new = c_result.value();

        // ================================================================
        // PHASE 5: ASSEMBLY AND CONVERGENCE CONTROL
        // ================================================================
        
        // 12. Complete new solution construction
        equations::SolutionState solution_new(n_eta, n_species);
        solution_new.V = solution.V;        // Already updated
        solution_new.F = solution.F;        // Already updated
        solution_new.g = solution.g;        // Already updated
        solution_new.c = std::move(c_new);  // New species solution
        solution_new.T = solution.T;        // Consistent temperature

        // 13. Final application of all boundary conditions
        enforce_edge_boundary_conditions(solution_new, bc_dynamic);

        // 14. Convergence check with fully consistent solution
        conv_info = check_convergence(solution_old, solution_new);
        conv_info.iterations = iter + 1;

        // 15. Relaxation factor adaptation based on convergence history
        double adaptive_factor = relaxation_controller_->adapt_relaxation_factor(conv_info, iter);
        std::cout << "Adaptive relaxation factor: " << std::scientific << std::setprecision(3) 
                  << adaptive_factor << " (max residual: " << conv_info.max_residual() << ")" << std::endl;

        // 16. Differential relaxation application
        solution = apply_relaxation_differential(solution_old, solution_new, adaptive_factor);

        // 17. Convergence test and exit
        if (conv_info.converged) {
            std::cout << "✓ Converged with adaptive factor: " << adaptive_factor << std::endl;
            break;
        }

        // 18. Divergence detection for early termination
        if (conv_info.max_residual() > 1e6) {
            return std::unexpected(SolverError(
                "Solution diverged at station {} iteration {} (residual={})",
                std::source_location::current(), station, iter, conv_info.max_residual()
            ));
        }
    }

    return conv_info;
}

auto BoundaryLayerSolver::solve_continuity_equation(const equations::SolutionState& solution,
                                                    double xi) -> std::expected<std::vector<double>, SolverError> {

  const auto n_eta = grid_->n_eta();
  const double d_eta = grid_->d_eta();

  // Get current xi value and lambda0 from xi_derivatives
  const double lambda0 = xi_derivatives_->lambda0();
  const auto F_derivatives = xi_derivatives_->F_derivative();

  std::vector<double> y_field(n_eta);

  for (std::size_t i = 0; i < n_eta; ++i) {
    y_field[i] = -(2.0 * xi * lambda0 + 1.0) * solution.F[i] - 2.0 * xi * F_derivatives[i];
  }

  // Integrate dV/dη = -y to get V
  // Note: The negative sign is handled in the continuity equation solver
  return equations::solve_continuity(y_field, d_eta, 0.0);
}

auto BoundaryLayerSolver::solve_momentum_equation(const equations::SolutionState& solution,
                                                  const coefficients::CoefficientSet& coeffs,
                                                  const conditions::BoundaryConditions& bc,
                                                  double xi) -> std::expected<std::vector<double>, SolverError> {

  auto result = equations::solve_momentum(solution.F, coeffs, bc, *xi_derivatives_, solution.V, xi, grid_->d_eta());

  if (!result) {
    return std::unexpected(
        SolverError("Momentum equation failed: {}", std::source_location::current(), result.error().message()));
  }

  return result.value();
}

auto BoundaryLayerSolver::solve_energy_equation(const equations::SolutionState& solution,
                                                const coefficients::CoefficientInputs& inputs,
                                                const coefficients::CoefficientSet& coeffs,
                                                const conditions::BoundaryConditions& bc,
                                                int station) -> std::expected<std::vector<double>, SolverError> {

  // Compute dF/deta for energy equation using 4th-order scheme
  auto dF_deta_result = coefficients::derivatives::compute_eta_derivative(solution.F, grid_->d_eta());
  if (!dF_deta_result) {
    return std::unexpected(SolverError("Failed to compute dF/deta: {}", std::source_location::current(),
                                       dF_deta_result.error().message()));
  }
  auto dF_deta = dF_deta_result.value();

  auto result = equations::solve_energy(solution.g, inputs, coeffs, bc, *xi_derivatives_, config_.simulation,
                                        solution.F, dF_deta, solution.V, station, grid_->d_eta());

  if (!result) {
    return std::unexpected(
        SolverError("Energy equation failed: {}", std::source_location::current(), result.error().message()));
  }

  auto g_solution = result.value();

  return g_solution;
}

auto BoundaryLayerSolver::solve_species_equations(const equations::SolutionState& solution,
                                                  const coefficients::CoefficientInputs& inputs,
                                                  const coefficients::CoefficientSet& coeffs,
                                                  const conditions::BoundaryConditions& bc,
                                                  int station) -> std::expected<core::Matrix<double>, SolverError> {

  auto result = equations::solve_species(solution.c, inputs, coeffs, bc, *xi_derivatives_, mixture_, config_.simulation,
                                         solution.F, solution.V, station, grid_->d_eta());

  if (!result) {
    return std::unexpected(
        SolverError("Species equations failed: {}", std::source_location::current(), result.error().message()));
  }

  return result.value();
}

auto BoundaryLayerSolver::update_temperature_field(
    std::span<const double> g_field, const core::Matrix<double>& composition, const conditions::BoundaryConditions& bc,
    std::span<const double> current_temperatures) -> std::expected<std::vector<double>, SolverError> {

  const auto n_eta = g_field.size();
  std::vector<double> enthalpy_field(n_eta);

  // Convert g (dimensionless enthalpy) to dimensional enthalpy
  for (std::size_t i = 0; i < n_eta; ++i) {
    enthalpy_field[i] = g_field[i] * bc.he();
  }

  auto result = h2t_solver_->solve(enthalpy_field, composition, bc, current_temperatures);
  if (!result) {
    return std::unexpected(
        SolverError("Temperature solve failed: {}", std::source_location::current(), result.error().message()));
  }

  return result.value().temperatures;
}

auto BoundaryLayerSolver::check_convergence(const equations::SolutionState& old_solution,
                                            const equations::SolutionState& new_solution) const noexcept
    -> ConvergenceInfo {

  ConvergenceInfo info;
  const double tol = config_.numerical.convergence_tolerance;

  // Compute L2 norms of residuals
  auto compute_residual = [](const auto& old_field, const auto& new_field) {
    double sum = 0.0;
    for (std::size_t i = 0; i < old_field.size(); ++i) {
      const double diff = new_field[i] - old_field[i];
      sum += diff * diff;
    }
    return std::sqrt(sum / old_field.size());
  };

  info.residual_F = compute_residual(old_solution.F, new_solution.F);
  info.residual_g = compute_residual(old_solution.g, new_solution.g);

  // Species residual (matrix)
  double c_sum = 0.0;
  std::size_t c_count = 0;
  for (std::size_t i = 0; i < old_solution.c.rows(); ++i) {
    for (std::size_t j = 0; j < old_solution.c.cols(); ++j) {
      const double diff = new_solution.c(i, j) - old_solution.c(i, j);
      c_sum += diff * diff;
      ++c_count;
    }
  }
  info.residual_c = std::sqrt(c_sum / c_count);

  info.converged = (info.residual_F < tol) && (info.residual_g < tol) && (info.residual_c < tol);

  std::cout << "CONVERGENCE : " << info.residual_F << " " << info.residual_g << " " << info.residual_c << std::endl;

  return info;
}

/* auto BoundaryLayerSolver::create_initial_guess(
    int station,
    double xi,
    const conditions::BoundaryConditions& bc
) const -> std::expected<equations::SolutionState, SolverError> {

    const auto n_eta = grid_->n_eta();
    const auto n_species = mixture_.n_species();
    const double eta_max = grid_->eta_max();

    equations::SolutionState guess(n_eta, n_species);

    // Get equilibrium composition at wall conditions
    auto equilibrium_result = mixture_.equilibrium_composition(bc.Tw(),
bc.P_e()); if (!equilibrium_result) { return std::unexpected(SolverError(
            "Failed to compute equilibrium composition at wall conditions: {}",
            std::source_location::current(),
equilibrium_result.error().message()
        ));
    }
    auto c_wall_equilibrium = equilibrium_result.value();

    for (std::size_t i = 0; i < n_eta; ++i) {
        const double eta = static_cast<double>(i) * eta_max / (n_eta - 1);
        const double eta_norm = static_cast<double>(i) / (n_eta - 1);

        // Boundary layer-like profiles
        guess.F[i] = eta_norm * (2.0 - eta_norm);

        guess.g[i] = 1.0;

        for (std::size_t j = 0; j < n_species; ++j) {
            guess.c(j, i) = c_wall_equilibrium[j];
        }
    }

    return guess;
} */

auto BoundaryLayerSolver::create_initial_guess(int station, double xi, const conditions::BoundaryConditions& bc) const
    -> std::expected<equations::SolutionState, SolverError> {

  const auto n_eta = grid_->n_eta();
  const auto n_species = mixture_.n_species();
  const double eta_max = grid_->eta_max();

  equations::SolutionState guess(n_eta, n_species);

  // ===== STEP 1: RETRIEVE BOUNDARY COMPOSITIONS =====

  // Compute equilibrium composition at the wall (Tw, P_e)
  auto equilibrium_result = mixture_.equilibrium_composition(bc.Tw(), bc.P_e());
  if (!equilibrium_result) {
    return std::unexpected(SolverError("Failed to compute equilibrium composition at wall conditions: {}",
                                       std::source_location::current(), equilibrium_result.error().message()));
  }
  auto c_wall_equilibrium = equilibrium_result.value();

  // Extract edge composition from input
  const auto& c_edge = bc.c_e();
  if (c_edge.size() != n_species) {
    return std::unexpected(SolverError("Edge composition size mismatch: expected {}, got {}",
                                       std::source_location::current(), n_species, c_edge.size()));
  }

  // ===== STEP 1.5: COMPUTE WALL EQUILIBRIUM ENTHALPY =====

  // Calculate enthalpy of equilibrium mixture at wall conditions
  auto h_wall_eq_result = mixture_.mixture_enthalpy(c_wall_equilibrium, bc.Tw(), bc.P_e());
  if (!h_wall_eq_result) {
    return std::unexpected(SolverError("Failed to compute wall equilibrium enthalpy: {}",
                                       std::source_location::current(), h_wall_eq_result.error().message()));
  }
  double h_wall_equilibrium = h_wall_eq_result.value();

  std::cout << "ENTHALPIE AU MUR = " << h_wall_equilibrium << std::endl;

  // Compute dimensionless enthalpy at wall
  double g_wall = h_wall_equilibrium / bc.he();

  // ===== STEP 2: TRANSITION FUNCTION (TANH INTERPOLATION) =====

  // Smooth transition function from wall (0) to edge (1)
  // Centered at 35% of the boundary layer thickness, with steep slope
  auto transition_function = [](double eta_norm, double eta_center = 0.35, double sharpness = 10.0) -> double {
    return 0.5 * (1.0 + std::tanh(sharpness * (eta_norm - eta_center)));
  };

  // ===== STEP 3: COMPUTE SPECIES PROFILES =====

  for (std::size_t i = 0; i < n_eta; ++i) {
    const double eta = static_cast<double>(i) * eta_max / (n_eta - 1);
    const double eta_norm = static_cast<double>(i) / (n_eta - 1); // Normalize to [0,1]

    // Analytical profiles for velocity F (unchanged)
    guess.F[i] = 1.0 - 1.034 * std::exp(-5.628 * eta / eta_max);

    // Physical enthalpy profile: linear interpolation from wall equilibrium to
    // edge
    guess.g[i] = g_wall + eta_norm * (1.0 - g_wall);

    // ===== STEP 4: INTERPOLATE SPECIES CONCENTRATIONS =====

    // Compute transition factor from wall to edge
    const double f = transition_function(eta_norm);

    // Linear interpolation between wall and edge compositions
    double sum_interpolated = 0.0;
    for (std::size_t j = 0; j < n_species; ++j) {
      guess.c(j, i) = c_wall_equilibrium[j] * (1.0 - f) + c_edge[j] * f;
      sum_interpolated += guess.c(j, i);
    }

    // Enforce mass conservation: normalize so that sum_j c_j = 1
    if (sum_interpolated > 1e-15) {
      for (std::size_t j = 0; j < n_species; ++j) {
        guess.c(j, i) /= sum_interpolated;
      }
    } else {
      std::cout << "GROS BUG" << std::endl;
      abort();
    }
  }
  return guess;
}

auto BoundaryLayerSolver::extrapolate_from_previous(const equations::SolutionState& previous_solution, double xi_prev,
                                                    double xi_current) const -> equations::SolutionState {

  // Simple extrapolation - could be improved with higher-order methods
  const double factor = (xi_current > xi_prev) ? 1 : 1;

  auto extrapolated = previous_solution;

  // Extrapolate F and g fields
  for (std::size_t i = 0; i < extrapolated.F.size(); ++i) {
    extrapolated.F[i] *= factor;
    extrapolated.g[i] = std::max(0.1, extrapolated.g[i] * factor); // Keep positive
  }

  return extrapolated;
}

auto BoundaryLayerSolver::apply_relaxation_differential(const equations::SolutionState& old_solution,
                                                        const equations::SolutionState& new_solution,
                                                        double base_factor) const -> equations::SolutionState {

  auto relaxed = new_solution;

  const double alpha_F = base_factor;
  const double alpha_g = base_factor * 1;
  const double alpha_c = base_factor * 1;
  const double alpha_T = base_factor * 1;

  for (std::size_t i = 0; i < relaxed.F.size() - 1; ++i) {
    relaxed.F[i] = (1.0 - alpha_F) * old_solution.F[i] + alpha_F * new_solution.F[i];
    relaxed.g[i] = (1.0 - alpha_g) * old_solution.g[i] + alpha_g * new_solution.g[i];
    // relaxed.T[i] = (1.0 - alpha_T) * old_solution.T[i] + alpha_T *
    // new_solution.T[i];
  }

  for (std::size_t i = 0; i < relaxed.c.rows(); ++i) {
    for (std::size_t j = 0; j < relaxed.c.cols(); ++j) {
      relaxed.c(i, j) = (1.0 - alpha_c) * old_solution.c(i, j) + alpha_c * new_solution.c(i, j);
    }
  }

  return relaxed;
}

auto BoundaryLayerSolver::compute_eta_derivatives(const equations::SolutionState& solution) const
    -> std::expected<equations::SolutionState, SolverError> {

  const auto n_eta = grid_->n_eta();
  const auto n_species = mixture_.n_species();
  const double d_eta = grid_->d_eta();

  equations::SolutionState derivatives(n_eta, n_species);

  // High-order finite difference derivatives using unified function
  using namespace coefficients::derivatives;
  auto dF_result = compute_eta_derivative(solution.F, d_eta);
  if (!dF_result) {
    return std::unexpected(
        SolverError("Failed to compute dF/deta: {}", std::source_location::current(), dF_result.error().message()));
  }
  derivatives.F = dF_result.value();

  auto dg_result = compute_eta_derivative(solution.g, d_eta);
  if (!dg_result) {
    return std::unexpected(
        SolverError("Failed to compute dg/deta: {}", std::source_location::current(), dg_result.error().message()));
  }
  derivatives.g = dg_result.value();

  auto dV_result = compute_eta_derivative(solution.V, d_eta);
  if (!dV_result) {
    return std::unexpected(
        SolverError("Failed to compute dV/deta: {}", std::source_location::current(), dV_result.error().message()));
  }
  derivatives.V = dV_result.value();

  // Species derivatives - use the high-order derivative functions
  for (std::size_t j = 0; j < n_species; ++j) {
    std::vector<double> c_row(n_eta);
    for (std::size_t i = 0; i < n_eta; ++i) {
      c_row[i] = solution.c(j, i);
    }

    auto dc_deta_result = compute_eta_derivative(c_row, d_eta);
    if (!dc_deta_result) {
      return std::unexpected(SolverError("Failed to compute dc/deta for species {}: {}",
                                         std::source_location::current(), j, dc_deta_result.error().message()));
    }
    auto dc_deta = dc_deta_result.value();
    for (std::size_t i = 0; i < n_eta; ++i) {
      derivatives.c(j, i) = dc_deta[i];
    }
  }

  const std::array<std::size_t, 3> indices_to_print{0, 10, 19};

  return derivatives;
}

/* auto BoundaryLayerSolver::compute_concentration_derivatives(
    const equations::SolutionState& solution
) const -> std::expected<DerivativeState, SolverError> {

    const auto n_eta = grid_->n_eta();
    const auto n_species = mixture_.n_species();
    const double d_eta = grid_->d_eta();

    DerivativeState derivatives(n_species, n_eta);

    // Use high-order finite difference schemes for both first and second
derivatives using namespace coefficients::derivatives;

    for (std::size_t j = 0; j < n_species; ++j) {
        std::vector<double> c_row(n_eta);
        for (std::size_t i = 0; i < n_eta; ++i) {
            c_row[i] = solution.c(j, i);
        }

        // Compute first derivatives
        auto dc_deta_result = compute_eta_derivative(c_row, d_eta);
        if (!dc_deta_result) {
            return std::unexpected(SolverError("Failed to compute dc/deta for
species {}: {}", std::source_location::current(), j,
dc_deta_result.error().message()));
        }
        auto dc_deta = dc_deta_result.value();
        for (std::size_t i = 0; i < n_eta; ++i) {
            derivatives.dc_deta(j, i) = dc_deta[i];
        }

        // Compute second derivatives
        auto dc_deta2 = compute_eta_second_derivative(c_row, d_eta);
        for (std::size_t i = 0; i < n_eta; ++i) {
            derivatives.dc_deta2(j, i) = dc_deta2[i];
        }
    }

    return derivatives;
} */

auto BoundaryLayerSolver::compute_concentration_derivatives(const equations::SolutionState& solution) const
    -> std::expected<DerivativeState, SolverError> {

  const auto n_eta = grid_->n_eta();
  const auto n_species = mixture_.n_species();
  const double d_eta = grid_->d_eta();

  DerivativeState derivatives(n_species, n_eta);

  // Use cascaded approach: apply the same high-order scheme twice
  using namespace coefficients::derivatives;

  for (std::size_t j = 0; j < n_species; ++j) {
    std::vector<double> c_row(n_eta);
    for (std::size_t i = 0; i < n_eta; ++i) {
      c_row[i] = solution.c(j, i);
    }

    // Compute first derivatives
    auto dc_deta_result = compute_eta_derivative(c_row, d_eta);
    if (!dc_deta_result) {
      return std::unexpected(SolverError("Failed to compute dc/deta for species {}: {}",
                                         std::source_location::current(), j, dc_deta_result.error().message()));
    }
    auto dc_deta = dc_deta_result.value();
    for (std::size_t i = 0; i < n_eta; ++i) {
      derivatives.dc_deta(j, i) = dc_deta[i];
    }

    // Compute second derivatives by applying the same scheme to first
    // derivatives
    auto dc_deta2_result = compute_eta_derivative(dc_deta, d_eta);
    if (!dc_deta2_result) {
      return std::unexpected(SolverError("Failed to compute dc/deta2 for species {}: {}",
                                         std::source_location::current(), j, dc_deta2_result.error().message()));
    }
    auto dc_deta2 = dc_deta2_result.value();
    for (std::size_t i = 0; i < n_eta; ++i) {
      derivatives.dc_deta2(j, i) = dc_deta2[i];
    }
  }

  return derivatives;
}

auto BoundaryLayerSolver::enforce_edge_boundary_conditions(equations::SolutionState& solution,
                                                           const conditions::BoundaryConditions& bc) const -> void {

  const auto n_eta = grid_->n_eta();
  const auto n_species = mixture_.n_species();

  if (n_eta == 0)
    return;

  const std::size_t edge_idx = n_eta - 1; // Last eta point is the edge

  // Use dynamic edge composition (updated by update_edge_properties)
  // The edge composition is now calculated from equilibrium conditions
  const auto& edge_composition = bc.c_e();
  for (std::size_t j = 0; j < n_species && j < edge_composition.size(); ++j) {
    solution.c(j, edge_idx) = edge_composition[j];
  }

  // Ensure species mass conservation (normalize if needed)
  double total_mass_fraction = 0.0;
  for (std::size_t j = 0; j < n_species; ++j) {
    total_mass_fraction += solution.c(j, edge_idx);
  }

  // Normalize if total is significantly different from 1.0
  if (std::abs(total_mass_fraction - 1.0) > 1e-12 && total_mass_fraction > 1e-12) {
    for (std::size_t j = 0; j < n_species; ++j) {
      solution.c(j, edge_idx) /= total_mass_fraction;
    }
  }
}

auto BoundaryLayerSolver::update_edge_properties(conditions::BoundaryConditions& bc,
                                                 const coefficients::CoefficientInputs& inputs,
                                                 const core::Matrix<double>& species_matrix) const -> void {

  const auto n_eta = grid_->n_eta();
  const auto n_species = mixture_.n_species();

  if (n_eta == 0 || inputs.T.empty()) {
    return; // No data to update
  }

  // Get edge conditions (last point in eta grid)
  const auto edge_idx = n_eta - 1;
  const double T_edge = inputs.T[edge_idx];
  const double P_edge = bc.P_e();

  // Get edge composition
  std::vector<double> edge_composition(n_species);
  for (std::size_t j = 0; j < n_species; ++j) {
    edge_composition[j] = species_matrix(j, edge_idx);
  }

  // Calculate new edge density using equation of state
  auto MW_result = mixture_.mixture_molecular_weight(edge_composition);
  if (!MW_result) {
    throw SolverError("Failed to compute edge molecular weight: {}", std::source_location::current(),
                      MW_result.error().message());
  }
  const double MW_edge = MW_result.value();
  const double rho_e_new = P_edge * MW_edge / (T_edge * thermophysics::constants::R_universal);

  // Calculate equilibrium composition at edge conditions
  auto eq_result = mixture_.equilibrium_composition(T_edge, P_edge);
  if (!eq_result) {
    throw SolverError("Failed to compute edge equilibrium composition: {}", std::source_location::current(),
                      eq_result.error().message());
  }
  auto edge_composition_eq = eq_result.value();

  // Calculate new edge viscosity using equilibrium composition
  auto mu_result = mixture_.viscosity(edge_composition, T_edge, P_edge);
  if (!mu_result) {
    throw SolverError("Failed to compute edge viscosity: {}", std::source_location::current(),
                      mu_result.error().message());
  }
  const double mu_e_new = mu_result.value();

  // Update boundary conditions with new equilibrium values
  bc.update_edge_density(rho_e_new);
  bc.update_edge_viscosity(mu_e_new);
  bc.edge.species_fractions = edge_composition_eq;
}

} // namespace blast::boundary_layer::solver