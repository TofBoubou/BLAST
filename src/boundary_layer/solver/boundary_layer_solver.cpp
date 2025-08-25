#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include "blast/boundary_layer/coefficients/heat_flux_calculator.hpp"
#include "blast/boundary_layer/equations/continuity.hpp"
#include "blast/boundary_layer/equations/energy.hpp"
#include "blast/boundary_layer/equations/momentum.hpp"
#include "blast/boundary_layer/equations/species.hpp"
#include "blast/boundary_layer/solver/expected_utils.hpp"
#include "blast/boundary_layer/solver/heat_flux_computer.hpp"
#include "blast/boundary_layer/solver/initial_guess_factory.hpp"
#include "blast/core/constants.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

// Note: solve_radiative_equilibrium function has been moved to
// RadiativeEquilibriumSolver

namespace blast::boundary_layer::solver {

BoundaryLayerSolver::BoundaryLayerSolver(
    const thermophysics::MixtureInterface &mixture,
    const io::Configuration &config)
    : mixture_(mixture), config_(config), original_config_(config) {

  // Create grid
  if (config.simulation.only_stagnation_point) {
    auto grid_result = grid::BoundaryLayerGrid::create_stagnation_grid(
        config.numerical, config.outer_edge, mixture_);
    if (!grid_result) {
      throw GridError(std::format("Failed to create stagnation grid: {}",
                                  grid_result.error().message()));
    }
    grid_ = std::make_unique<grid::BoundaryLayerGrid>(
        std::move(grid_result.value()));
  } else {
    auto grid_result = grid::BoundaryLayerGrid::create_downstream_grid(
        config.numerical, config.outer_edge, config.output, mixture_);
    if (!grid_result) {
      throw GridError(std::format("Failed to create downstream grid: {}",
                                  grid_result.error().message()));
    }
    grid_ = std::make_unique<grid::BoundaryLayerGrid>(
        std::move(grid_result.value()));
  }

  // Create coefficient calculator
  coeff_calculator_ = std::make_unique<coefficients::CoefficientCalculator>(
      mixture_, config_.simulation, config_.numerical, config_.outer_edge);

  heat_flux_calculator_ = std::make_unique<coefficients::HeatFluxCalculator>(
      mixture_, config_.simulation, config_.numerical);

  // Create enthalpy-temperature solver
  thermodynamics::EnthalpyTemperatureSolverConfig h2t_config{
      .tolerance = config_.numerical.solvers.h2t_tolerance,
      .max_iterations = config_.numerical.solvers.h2t_max_iterations};
  h2t_solver_ = std::make_unique<thermodynamics::EnthalpyTemperatureSolver>(
      mixture_, h2t_config);

  // Create xi derivatives manager
  xi_derivatives_ = std::make_unique<coefficients::XiDerivatives>(
      grid_->n_eta(), mixture_.n_species());

  // Create derivative calculator
  derivative_calculator_ =
      std::make_unique<coefficients::DerivativeCalculator>(grid_->d_eta());

  continuation_ = std::make_unique<ContinuationMethod>();

  // Create utility objects for eliminating duplication
  heat_flux_computer_ = std::make_unique<HeatFluxComputer>(
      HeatFluxComputer::Config{.coeff_calculator = *coeff_calculator_,
                               .heat_flux_calculator = *heat_flux_calculator_,
                               .derivative_calculator = *derivative_calculator_,
                               .xi_derivatives = *xi_derivatives_});

  initial_guess_factory_ =
      std::make_unique<InitialGuessFactory>(*grid_, mixture_);

  // Initialize specialized solver components
  initialize_solver_components();
}

auto BoundaryLayerSolver::initialize_solver_components() -> void {
  // Create specialized solver components
  station_solver_ = std::make_unique<StationSolver>(*this, mixture_, config_);
  convergence_manager_ =
      std::make_unique<ConvergenceManager>(*this, mixture_, config_);
  equation_solver_ = std::make_unique<EquationSolver>(*this, mixture_, config_);
  radiative_solver_ =
      std::make_unique<RadiativeEquilibriumSolver>(*this, config_);

  // Set up cross-dependencies
  convergence_manager_->set_radiative_solver(radiative_solver_.get());

  // Set stable configuration for station solver if continuation is configured
  if (original_config_.continuation.wall_temperature_stable > 0) {
    StationSolver::StableGuessConfig stable_config{
        .wall_temperature_stable =
            original_config_.continuation.wall_temperature_stable,
        .edge_temperature_stable =
            original_config_.continuation.edge_temperature_stable,
        .pressure_stable = original_config_.continuation.pressure_stable};
    station_solver_->set_stable_config(stable_config);
  }
}

auto BoundaryLayerSolver::solve()
    -> std::expected<SolutionResult, SolverError> {

  SolutionResult result;
  const auto xi_stations_original = grid_->xi_coordinates();

  double prev_xi = 0.0;
  std::vector<double> prev_F, prev_g, prev_T;
  core::Matrix<double> prev_c;

  // Solve stagnation point first
  {
    const double xi = 0.0;
    auto bc_result = conditions::create_stagnation_conditions(
        config_.outer_edge, config_.wall_parameters, config_.simulation,
        mixture_);
    if (!bc_result) {
      return std::unexpected(bc_result.error());
    }

    const auto &edge_point = config_.outer_edge.edge_points[0];
    const double T_edge = edge_point.temperature;

    auto guess_result = initial_guess_factory_->create_initial_guess(
        0, xi, bc_result.value(), T_edge);
    if (!guess_result) {
      return std::unexpected(guess_result.error());
    }

    auto station_result = solve_station(0, xi, guess_result.value());
    if (!station_result) {
      return std::unexpected(station_result.error());
    }

    result.xi_solved.push_back(xi);
    result.stations.push_back(station_result.value());

    auto bc = bc_result.value();
    if (config_.simulation.wall_mode ==
        io::SimulationConfig::WallMode::Radiative) {
      bc.wall.temperature = result.stations.back().T[0];
    }

    auto heat_flux = heat_flux_computer_->compute_heat_flux_only(
        result.stations.back(), bc, xi, 0);
    if (!heat_flux) {
      return std::unexpected(heat_flux.error());
    }
    result.heat_flux_data.push_back(heat_flux.value());

    prev_xi = xi;
    prev_F = result.stations.back().F;
    prev_g = result.stations.back().g;
    prev_T = result.stations.back().T;
    prev_c = result.stations.back().c;
  }

  // Solve other stations adaptively
  for (std::size_t i = 1; i < xi_stations_original.size(); ++i) {
    StationInterval interval{xi_stations_original[i - 1],
                             xi_stations_original[i], 0};
    auto interval_result = solve_interval_adaptive(
        interval, result, prev_F, prev_g, prev_T, prev_c, prev_xi);
    if (!interval_result) {
      return std::unexpected(interval_result.error());
    }
  }

  result.converged = true;
  result.temperature_fields.reserve(result.stations.size());
  for (const auto &station : result.stations) {
    result.temperature_fields.push_back(station.T);
  }
  return result;
}

auto BoundaryLayerSolver::solve_interval_adaptive(
    const StationInterval &interval, SolutionResult &result,
    std::vector<double> &prev_F, std::vector<double> &prev_g,
    std::vector<double> &prev_T, core::Matrix<double> &prev_c, double prev_xi)
    -> std::expected<void, SolverError> {

  (void)prev_xi;
  const double interval_size = interval.xi_end - interval.xi_start;

  int station_idx = static_cast<int>(result.xi_solved.size());
  xi_derivatives_->update_station(station_idx, interval.xi_end, prev_F, prev_g,
                                  prev_c);

  auto bc = interpolate_edge_conditions(interval.xi_end);

  auto initial_guess =
      initial_guess_factory_->create_initial_guess_from_previous(
          prev_F, prev_g, prev_T, prev_c, bc, interval.xi_end);

  auto conv_result = convergence_manager_->iterate_station_adaptive(
      station_idx, interval.xi_end, bc, initial_guess);
  if (conv_result && conv_result->converged) {
    result.xi_solved.push_back(interval.xi_end);
    result.stations.push_back(initial_guess);

    auto heat_flux = heat_flux_computer_->compute_heat_flux_only(
        initial_guess, bc, interval.xi_end, station_idx);
    if (!heat_flux) {
      return std::unexpected(heat_flux.error());
    }
    result.heat_flux_data.push_back(heat_flux.value());
    result.total_iterations += conv_result->iterations;

    prev_F = initial_guess.F;
    prev_g = initial_guess.g;
    prev_T = initial_guess.T;
    prev_c = initial_guess.c;
    return {};
  }

  if (interval.depth >= MAX_SUBDIVISION_DEPTH ||
      interval_size < MIN_INTERVAL_SIZE) {
    return std::unexpected(ConvergenceError(std::format(
        "Cannot subdivide further: interval [{:.4f}, {:.4f}], depth {}, size "
        "{:.6f}",
        interval.xi_start, interval.xi_end, interval.depth, interval_size)));
  }

  std::cout << std::format("Station at xi={:.4f} failed to converge. "
                           "Subdividing interval [{:.4f}, {:.4f}]...\n",
                           interval.xi_end, interval.xi_start, interval.xi_end);

  const double xi_mid = interval.midpoint();

  StationInterval first_half{interval.xi_start, xi_mid, interval.depth + 1};
  auto first_result = solve_interval_adaptive(first_half, result, prev_F,
                                              prev_g, prev_T, prev_c, prev_xi);
  if (!first_result) {
    return first_result;
  }

  StationInterval second_half{xi_mid, interval.xi_end, interval.depth + 1};
  return solve_interval_adaptive(second_half, result, prev_F, prev_g, prev_T,
                                 prev_c, prev_xi);
}

auto BoundaryLayerSolver::interpolate_edge_conditions(double xi) const
    -> conditions::BoundaryConditions {
  auto bc_result = conditions::interpolate_boundary_conditions(
      1, xi, grid_->xi_coordinates(), config_.outer_edge,
      config_.wall_parameters, config_.simulation, mixture_);
  if (!bc_result) {
    throw BoundaryConditionError(
        std::format("Failed to interpolate boundary conditions at xi={}", xi));
  }
  return bc_result.value();
}

auto BoundaryLayerSolver::solve_station(
    int station, double xi, const equations::SolutionState &initial_guess)
    -> std::expected<equations::SolutionState, SolverError> {

  // Set continuation mode for station solver
  station_solver_->set_continuation_mode(in_continuation_);

  // Delegate to specialized StationSolver
  return station_solver_->solve_station(station, xi, initial_guess);
}

// Note: iterate_station_adaptive method has been moved to ConvergenceManager

auto BoundaryLayerSolver::solve_momentum_equation(
    const equations::SolutionState &solution,
    const coefficients::CoefficientSet &coeffs,
    const conditions::BoundaryConditions &bc, double xi)
    -> std::expected<std::vector<double>, SolverError> {

  // Delegate to specialized EquationSolver
  return equation_solver_->solve_momentum_equation(solution, coeffs, bc, xi);
}

auto BoundaryLayerSolver::solve_energy_equation(
    const equations::SolutionState &solution,
    const coefficients::CoefficientInputs &inputs,
    const coefficients::CoefficientSet &coeffs,
    const conditions::BoundaryConditions &bc,
    const thermophysics::MixtureInterface &mixture, int station)
    -> std::expected<std::vector<double>, SolverError> {

  // Delegate to specialized EquationSolver
  return equation_solver_->solve_energy_equation(solution, inputs, coeffs, bc,
                                                 mixture, station);
}

auto BoundaryLayerSolver::solve_species_equations(
    const equations::SolutionState &solution,
    const coefficients::CoefficientInputs &inputs,
    const coefficients::CoefficientSet &coeffs,
    const conditions::BoundaryConditions &bc, int station)
    -> std::expected<core::Matrix<double>, SolverError> {

  // Delegate to specialized EquationSolver
  return equation_solver_->solve_species_equations(solution, inputs, coeffs, bc,
                                                   station);
}

auto BoundaryLayerSolver::update_temperature_field(
    std::span<const double> g_field, const core::Matrix<double> &composition,
    const conditions::BoundaryConditions &bc,
    std::span<const double> current_temperatures)
    -> std::expected<std::vector<double>, SolverError> {

  // Delegate to specialized EquationSolver
  return equation_solver_->update_temperature_field(g_field, composition, bc,
                                                    current_temperatures);
}

// Note: check_convergence method has been moved to ConvergenceManager

// Note: create_initial_guess method moved to InitialGuessFactory

// Note: apply_relaxation_differential method has been moved to
// ConvergenceManager

auto BoundaryLayerSolver::compute_all_derivatives(
    const equations::SolutionState &solution) const
    -> std::expected<coefficients::UnifiedDerivativeState, SolverError> {

  auto result = derivative_calculator_->compute_all_derivatives(solution);
  if (!result) {
    return std::unexpected(NumericError(std::format(
        "Failed to compute derivatives: {}", result.error().message())));
  }
  return result.value();
}

auto BoundaryLayerSolver::enforce_edge_boundary_conditions(
    equations::SolutionState &solution,
    const conditions::BoundaryConditions &bc) const -> void {

  const auto n_eta = grid_->n_eta();
  const auto n_species = mixture_.n_species();

  if (n_eta == 0)
    return;

  const std::size_t edge_idx = n_eta - 1; // Last eta point is the edge

  // Use dynamic edge composition (updated by update_edge_properties)
  // The edge composition is now calculated from equilibrium conditions
  const auto &edge_composition = bc.c_e();
  for (std::size_t j = 0; j < n_species && j < edge_composition.size(); ++j) {
    solution.c(j, edge_idx) = edge_composition[j];
  }

  // Ensure species mass conservation (normalize if needed)
  double total_mass_fraction = 0.0;
  for (std::size_t j = 0; j < n_species; ++j) {
    total_mass_fraction += solution.c(j, edge_idx);
  }

  // Normalize if total is significantly different from 1.0
  if (std::abs(total_mass_fraction - 1.0) > 1e-12 &&
      total_mass_fraction > 1e-12) {
    for (std::size_t j = 0; j < n_species; ++j) {
      solution.c(j, edge_idx) /= total_mass_fraction;
    }
  }
}

auto BoundaryLayerSolver::update_edge_properties(
    conditions::BoundaryConditions &bc,
    const coefficients::CoefficientInputs &inputs,
    const core::Matrix<double> &species_matrix) const
    -> std::expected<void, SolverError> {

  const auto n_eta = grid_->n_eta();
  const auto n_species = mixture_.n_species();

  if (n_eta == 0 || inputs.T.empty()) {
    return {}; // No data to update
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
    return std::unexpected(
        NumericError(std::format("Failed to compute edge molecular weight: {}",
                                 MW_result.error().message())));
  }
  const double MW_edge = MW_result.value();
  const double rho_e_new =
      P_edge * MW_edge / (T_edge * constants::physical::universal_gas_constant);

  // Calculate equilibrium composition at edge conditions
  auto eq_result = mixture_.equilibrium_composition(T_edge, P_edge);
  if (!eq_result) {
    return std::unexpected(NumericError(
        std::format("Failed to compute edge equilibrium composition: {}",
                    eq_result.error().message())));
  }
  auto edge_composition_eq = eq_result.value();

  // Calculate new edge viscosity using equilibrium composition
  auto mu_result = mixture_.viscosity(edge_composition, T_edge, P_edge);
  if (!mu_result) {
    return std::unexpected(NumericError(std::format(
        "Failed to compute edge viscosity: {}", mu_result.error().message())));
  }
  const double mu_e_new = mu_result.value();

  // Update boundary conditions with new values
  bc.update_edge_density(rho_e_new);
  bc.update_edge_viscosity(mu_e_new);

  // Only update species fractions if boundary_override is false
  if (!bc.edge.boundary_override) {
    bc.edge.species_fractions = edge_composition_eq;
  }

  return {};
}

} // namespace blast::boundary_layer::solver