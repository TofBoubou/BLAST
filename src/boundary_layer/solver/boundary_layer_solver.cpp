#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include "blast/boundary_layer/coefficients/heat_flux_calculator.hpp"
#include "blast/boundary_layer/equations/continuity.hpp"
#include "blast/boundary_layer/equations/energy.hpp"
#include "blast/boundary_layer/equations/momentum.hpp"
#include "blast/boundary_layer/equations/species.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

// Note: solve_radiative_equilibrium function has been moved to RadiativeEquilibriumSolver

namespace blast::boundary_layer::solver {

BoundaryLayerSolver::BoundaryLayerSolver(const thermophysics::MixtureInterface& mixture,
                                         const io::Configuration& config)
    : mixture_(mixture), config_(config), original_config_(config) {

  // Create grid
  if (config.simulation.only_stagnation_point) {
    auto grid_result = grid::BoundaryLayerGrid::create_stagnation_grid(config.numerical, config.outer_edge, mixture_);
    if (!grid_result) {
      throw GridError(std::format("Failed to create stagnation grid: {}", grid_result.error().message()));
    }
    grid_ = std::make_unique<grid::BoundaryLayerGrid>(std::move(grid_result.value()));
  } else {
    auto grid_result =
        grid::BoundaryLayerGrid::create_downstream_grid(config.numerical, config.outer_edge, config.output, mixture_);
    if (!grid_result) {
      throw GridError(std::format("Failed to create downstream grid: {}", grid_result.error().message()));
    }
    grid_ = std::make_unique<grid::BoundaryLayerGrid>(std::move(grid_result.value()));
  }

  // Create coefficient calculator
  coeff_calculator_ =
      std::make_unique<coefficients::CoefficientCalculator>(mixture_, config_.simulation, config_.numerical, config_.outer_edge);

  heat_flux_calculator_ =
      std::make_unique<coefficients::HeatFluxCalculator>(mixture_, config_.simulation, config_.numerical);

  // Create enthalpy-temperature solver
  thermodynamics::EnthalpyTemperatureSolverConfig h2t_config{.tolerance = config_.numerical.solvers.h2t_tolerance,
                                                             .max_iterations =
                                                                 config_.numerical.solvers.h2t_max_iterations};
  h2t_solver_ = std::make_unique<thermodynamics::EnthalpyTemperatureSolver>(mixture_, h2t_config);

  // Create xi derivatives manager
  xi_derivatives_ = std::make_unique<coefficients::XiDerivatives>(grid_->n_eta(), mixture_.n_species());

  // Create derivative calculator
  derivative_calculator_ = std::make_unique<coefficients::DerivativeCalculator>(grid_->d_eta());

  continuation_ = std::make_unique<ContinuationMethod>();
  
  // Initialize specialized solver components
  initialize_solver_components();
}

auto BoundaryLayerSolver::initialize_solver_components() -> void {
  // Create specialized solver components
  station_solver_ = std::make_unique<StationSolver>(*this, mixture_, config_);
  convergence_manager_ = std::make_unique<ConvergenceManager>(*this, mixture_, config_);
  equation_solver_ = std::make_unique<EquationSolver>(*this, mixture_, config_);
  radiative_solver_ = std::make_unique<RadiativeEquilibriumSolver>(*this, config_);
  
  // Set up cross-dependencies
  convergence_manager_->set_radiative_solver(radiative_solver_.get());
  
  // Set stable configuration for station solver if continuation is configured
  if (original_config_.continuation.wall_temperature_stable > 0) {
    StationSolver::StableGuessConfig stable_config{
      .wall_temperature_stable = original_config_.continuation.wall_temperature_stable,
      .edge_temperature_stable = original_config_.continuation.edge_temperature_stable,
      .pressure_stable = original_config_.continuation.pressure_stable
    };
    station_solver_->set_stable_config(stable_config);
  }
}

auto BoundaryLayerSolver::solve() -> std::expected<SolutionResult, SolverError> {

  SolutionResult result;
  result.xi_solved.reserve(grid_->xi_coordinates().size());
  result.stations.reserve(grid_->xi_coordinates().size());
  result.temperature_fields.reserve(grid_->xi_coordinates().size());
  result.heat_flux_data.reserve(grid_->xi_coordinates().size());

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
      auto bc_result = conditions::create_stagnation_conditions(config_.outer_edge, config_.wall_parameters,
                                                                config_.simulation, mixture_);
      if (!bc_result) {
        return std::unexpected(BoundaryConditionError(std::format(
            "Failed to create boundary conditions for station {}: {}", station, bc_result.error().message())));
      }

      // Get T_edge from configuration
      const auto& edge_point = config_.outer_edge.edge_points[0];
      const double T_edge = edge_point.temperature;

      auto guess_result = create_initial_guess(station, xi, bc_result.value(), T_edge);
      if (!guess_result) {
        return std::unexpected(guess_result.error());
      }
      initial_guess = std::move(guess_result.value());
    } else {
      // For downstream stations, use previous station as initial guess
      initial_guess = result.stations.back();
    }

    // Solve this station
    auto station_result = solve_station(station, xi, initial_guess);
    if (!station_result) {
      return std::unexpected(NumericError(
          std::format("Failed to solve station {} (xi={}): {}", station, xi, station_result.error().message())));
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

    // Calculate heat flux for this converged station
    auto derivatives_result = compute_all_derivatives(result.stations.back());
    if (!derivatives_result) {
      return std::unexpected(NumericError(std::format("Failed to compute derivatives for heat flux at station {}: {}",
                                                      station, derivatives_result.error().message())));
    }
    auto derivatives = derivatives_result.value();

    auto final_inputs = coefficients::CoefficientInputs{.xi = xi,
                                                        .F = result.stations.back().F,
                                                        .c = result.stations.back().c,
                                                        .dc_deta = derivatives.dc_deta,
                                                        .dc_deta2 = derivatives.dc_deta2,
                                                        .T = result.stations.back().T};

    auto bc_result =
        (station == 0)
            ? conditions::create_stagnation_conditions(config_.outer_edge, config_.wall_parameters, config_.simulation,
                                                       mixture_)
            : conditions::interpolate_boundary_conditions(station, xi, grid_->xi_coordinates(), config_.outer_edge,
                                                          config_.wall_parameters, config_.simulation, mixture_);

    if (!bc_result) {
      return std::unexpected(BoundaryConditionError(std::format(
          "Failed to get boundary conditions for heat flux at station {}: {}", station, bc_result.error().message())));
    }
    auto bc = bc_result.value();
    
    // Update wall temperature if radiative equilibrium was used
    if (config_.wall_parameters.emissivity > 0.0) {
      // Get the final temperature from the converged solution
      bc.wall.temperature = result.stations.back().T[0];  // T at wall (eta=0)
    }
    
    // DEBUG: Print wall temperature used for final heat flux calculation
    /* std::cout << std::format("[FINAL_BC] Station {}: T_wall={:.1f} K (updated for radiative equilibrium)", station, bc.wall.temperature) << std::endl; */

    auto coeffs_result = coeff_calculator_->calculate(final_inputs, bc, *xi_derivatives_);
    if (!coeffs_result) {
      return std::unexpected(NumericError(std::format("Failed to compute coefficients for heat flux at station {}: {}",
                                                      station, coeffs_result.error().message())));
    }
    auto coeffs = coeffs_result.value();

    // DEBUG: Print temperatures being used for profile calculation
/*     std::cout << std::format("[PROFILE_CALC] Station {}: T_wall_in_inputs={:.1f} K, T_wall_in_bc={:.1f} K", 
                             station, final_inputs.T[0], bc.wall.temperature) << std::endl; */
    
    auto heat_flux_result =
        heat_flux_calculator_->calculate(final_inputs, coeffs, bc, derivatives.dT_deta, station, xi);
    if (!heat_flux_result) {
      return std::unexpected(NumericError(
          std::format("Failed to compute heat flux at station {}: {}", station, heat_flux_result.error().message())));
    }
    
    // DEBUG: Print heat flux calculated for final output
    auto heat_flux_final = heat_flux_result.value();
/*     std::cout << std::format("[FINAL_CALC] Station {}: q_cond={:.2e} W/m², q_diff={:.2e} W/m², q_total={:.2e} W/m²", 
                             station, heat_flux_final.q_wall_conductive_dim, heat_flux_final.q_wall_diffusive_dim, heat_flux_final.q_wall_total_dim) << std::endl; */

    result.heat_flux_data.push_back(std::move(heat_flux_result.value()));

    // Extract modal temperatures if multiple energy modes are present
    if (mixture_.get_number_energy_modes() > 1) {
      auto extract_modal_temperatures_for_station =
          [&](const equations::SolutionState& sol) -> std::expected<std::vector<std::vector<double>>, SolverError> {
        const auto n_eta = grid_->n_eta();
        const auto n_species = mixture_.n_species();
        const auto n_modes = mixture_.get_number_energy_modes();
        std::vector<std::vector<double>> modal_temps(n_modes, std::vector<double>(n_eta));

        for (std::size_t i = 0; i < n_eta; ++i) {
          std::vector<double> local_composition(n_species);
          for (std::size_t j = 0; j < n_species; ++j) {
            local_composition[j] = sol.c(j, i);
          }

          auto modal_result = mixture_.extract_modal_temperatures(local_composition, sol.T[i], bc.P_e());
          if (!modal_result) {
            return std::unexpected(NumericError(
                std::format("Failed to extract modal temperatures at eta={}: {}", i, modal_result.error().message())));
          }

          const auto& temps = modal_result.value();
          for (std::size_t mode = 0; mode < n_modes; ++mode) {
            modal_temps[mode][i] = temps[mode];
          }
        }

        return modal_temps;
      };

      auto modal_temps_result = extract_modal_temperatures_for_station(result.stations.back());
      if (!modal_temps_result) {
        return std::unexpected(modal_temps_result.error());
      }
      result.modal_temperature_fields.push_back(std::move(modal_temps_result.value()));

      if (result.temperature_mode_names.empty()) {
        const auto n_modes = mixture_.get_number_energy_modes();
        if (n_modes == 2) {
          result.temperature_mode_names = {"T_translation_rotation", "T_vibrational_electronic"};
        } else {
          result.temperature_mode_names.resize(n_modes);
          for (std::size_t m = 0; m < n_modes; ++m) {
            result.temperature_mode_names[m] = std::format("T_mode_{}", m);
          }
        }
      }
    }
  }

  result.converged = true; // All stations solved successfully
  return result;
}

auto BoundaryLayerSolver::solve_station(int station, double xi, const equations::SolutionState& initial_guess)
    -> std::expected<equations::SolutionState, SolverError> {
  
  // Set continuation mode for station solver
  station_solver_->set_continuation_mode(in_continuation_);
  
  // Delegate to specialized StationSolver
  return station_solver_->solve_station(station, xi, initial_guess);
}

// Note: iterate_station_adaptive method has been moved to ConvergenceManager

auto BoundaryLayerSolver::solve_momentum_equation(const equations::SolutionState& solution,
                                                  const coefficients::CoefficientSet& coeffs,
                                                  const conditions::BoundaryConditions& bc, double xi)
    -> std::expected<std::vector<double>, SolverError> {
  
  // Delegate to specialized EquationSolver
  return equation_solver_->solve_momentum_equation(solution, coeffs, bc, xi);
}

auto BoundaryLayerSolver::solve_energy_equation(const equations::SolutionState& solution,
                                                const coefficients::CoefficientInputs& inputs,
                                                const coefficients::CoefficientSet& coeffs,
                                                const conditions::BoundaryConditions& bc,
                                                const thermophysics::MixtureInterface& mixture, int station)
    -> std::expected<std::vector<double>, SolverError> {
  
  // Delegate to specialized EquationSolver
  return equation_solver_->solve_energy_equation(solution, inputs, coeffs, bc, mixture, station);
}

auto BoundaryLayerSolver::solve_species_equations(const equations::SolutionState& solution,
                                                  const coefficients::CoefficientInputs& inputs,
                                                  const coefficients::CoefficientSet& coeffs,
                                                  const conditions::BoundaryConditions& bc, int station)
    -> std::expected<core::Matrix<double>, SolverError> {
  
  // Delegate to specialized EquationSolver
  return equation_solver_->solve_species_equations(solution, inputs, coeffs, bc, station);
}

auto BoundaryLayerSolver::update_temperature_field(std::span<const double> g_field,
                                                   const core::Matrix<double>& composition,
                                                   const conditions::BoundaryConditions& bc,
                                                   std::span<const double> current_temperatures)
    -> std::expected<std::vector<double>, SolverError> {
  
  // Delegate to specialized EquationSolver
  return equation_solver_->update_temperature_field(g_field, composition, bc, current_temperatures);
}

// Note: check_convergence method has been moved to ConvergenceManager

auto BoundaryLayerSolver::create_initial_guess(int station, double xi, const conditions::BoundaryConditions& bc,
                                               double T_edge) const
    -> std::expected<equations::SolutionState, SolverError> {

  const auto n_eta = grid_->n_eta();
  const auto n_species = mixture_.n_species();
  const double eta_max = grid_->eta_max();

  equations::SolutionState guess(n_eta, n_species);

  std::fill(guess.V.begin(), guess.V.end(), 0.0);

  if (n_species == 1) {
    guess.c.eigen().setOnes();

    std::array<double, 1> c_wall{{1.0}};
    auto h_wall_eq_result = mixture_.mixture_enthalpy(c_wall, bc.Tw(), bc.P_e());
    if (!h_wall_eq_result) {
      return std::unexpected(NumericError(
          std::format("Failed to compute wall equilibrium enthalpy: {}", h_wall_eq_result.error().message())));
    }
    double g_wall = h_wall_eq_result.value() / bc.he();

    for (std::size_t i = 0; i < n_eta; ++i) {
      const double eta = static_cast<double>(i) * eta_max / (n_eta - 1);
      const double eta_norm = static_cast<double>(i) / (n_eta - 1);

      guess.F[i] = 1.0 - 1.034 * std::exp(-5.628 * eta / eta_max);
      guess.g[i] = g_wall + eta_norm * (1.0 - g_wall);
      guess.T[i] = bc.Tw() + guess.F[i] * (T_edge - bc.Tw());
    }

    return guess;
  }

  // ===== STEP 1: RETRIEVE BOUNDARY COMPOSITIONS =====

  // Compute equilibrium composition at the wall (Tw, P_e)
  auto equilibrium_result = mixture_.equilibrium_composition(bc.Tw(), bc.P_e());
  if (!equilibrium_result) {
    return std::unexpected(NumericError(std::format("Failed to compute equilibrium composition at wall conditions: {}",
                                                    equilibrium_result.error().message())));
  }
  auto c_wall_equilibrium = equilibrium_result.value();

  // Extract edge composition from input
  const auto& c_edge = bc.c_e();
  if (c_edge.size() != n_species) {
    return std::unexpected(InitializationError(
        std::format("Edge composition size mismatch: expected {}, got {}", n_species, c_edge.size())));
  }

  // ===== STEP 1.5: COMPUTE WALL EQUILIBRIUM ENTHALPY =====

  // Calculate enthalpy of equilibrium mixture at wall conditions
  auto h_wall_eq_result = mixture_.mixture_enthalpy(c_wall_equilibrium, bc.Tw(), bc.P_e());
  if (!h_wall_eq_result) {
    return std::unexpected(NumericError(
        std::format("Failed to compute wall equilibrium enthalpy: {}", h_wall_eq_result.error().message())));
  }
  double h_wall_equilibrium = h_wall_eq_result.value();

  // std::cout << "ENTHALPIE AU MUR = " << h_wall_equilibrium << std::endl;

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
      return std::unexpected(InitializationError(std::format(
          "Mass conservation violation: sum of species concentrations is too small ({})", sum_interpolated)));
    }
  }

  // ===== STEP 5: INITIALIZE TEMPERATURE FIELD =====
  for (std::size_t i = 0; i < n_eta; ++i) {
    const double eta = static_cast<double>(i) * eta_max / (n_eta - 1);
    double F_normalized = 1.0 - 1.034 * std::exp(-5.628 * eta / eta_max);
    guess.T[i] = bc.Tw() + F_normalized * (T_edge - bc.Tw());
  }

  return guess;
}

// Note: apply_relaxation_differential method has been moved to ConvergenceManager

auto BoundaryLayerSolver::compute_all_derivatives(const equations::SolutionState& solution) const
    -> std::expected<coefficients::UnifiedDerivativeState, SolverError> {

  auto result = derivative_calculator_->compute_all_derivatives(solution);
  if (!result) {
    return std::unexpected(NumericError(std::format("Failed to compute derivatives: {}", result.error().message())));
  }
  return result.value();
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
                                                 const core::Matrix<double>& species_matrix) const
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
        NumericError(std::format("Failed to compute edge molecular weight: {}", MW_result.error().message())));
  }
  const double MW_edge = MW_result.value();
  const double rho_e_new = P_edge * MW_edge / (T_edge * thermophysics::constants::R_universal);

  // Calculate equilibrium composition at edge conditions
  auto eq_result = mixture_.equilibrium_composition(T_edge, P_edge);
  if (!eq_result) {
    return std::unexpected(
        NumericError(std::format("Failed to compute edge equilibrium composition: {}", eq_result.error().message())));
  }
  auto edge_composition_eq = eq_result.value();

  // Calculate new edge viscosity using equilibrium composition
  auto mu_result = mixture_.viscosity(edge_composition, T_edge, P_edge);
  if (!mu_result) {
    return std::unexpected(
        NumericError(std::format("Failed to compute edge viscosity: {}", mu_result.error().message())));
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