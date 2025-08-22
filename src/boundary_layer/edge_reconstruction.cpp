#include "blast/boundary_layer/edge_reconstruction.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

namespace blast::boundary_layer {

auto EdgeTemperatureReconstructor::reconstruct() 
    -> std::expected<ReconstructedEdgeConditions, solver::SolverError> {
  
  // Backup and update GSI file with the catalyticity value
  auto backup_result = gsi_manager_.backup_gsi_file();
  if (!backup_result) {
    std::cerr << "Warning: Failed to backup GSI file: " << backup_result.error() << std::endl;
  }
  
  // Update GSI file with the specified catalyticity
  if (config_.boundary_conditions.catalyticity > 0) {
    auto update_result = gsi_manager_.update_gsi_catalyticity(config_.boundary_conditions.catalyticity);
    if (!update_result) {
      std::cerr << "Warning: Failed to update GSI catalyticity: " << update_result.error() << std::endl;
    } else {
      std::cout << "Updated GSI file with catalyticity = " << config_.boundary_conditions.catalyticity << std::endl;
      
      // Reload mixture with updated GSI
      auto reload_result = mixture_.reload_gsi();
      if (!reload_result) {
        return std::unexpected(solver::SolverError(
            std::format("Failed to reload mixture after GSI update: {}", reload_result.error())));
      }
      std::cout << "Mixture reloaded with updated GSI file" << std::endl;
    }
  }
  
  // Brent's method parameters
  bool mflag = true;
  double a = config_.solver.temperature_min;
  double b = config_.solver.temperature_max;
  double d = 0;
  double c_brent = a;
  
  std::cout << "\n=== Edge Temperature Reconstruction ===" << std::endl;
  std::cout << "Target heat flux: " << config_.target_heat_flux << " W/m²" << std::endl;
  std::cout << "Temperature search range: [" << a << ", " << b << "] K" << std::endl;
  
  // Compute flux at bounds
  auto qa_result = compute_heat_flux_at_temperature(a);
  if (!qa_result) {
    return std::unexpected(qa_result.error());
  }
  
  auto qb_result = compute_heat_flux_at_temperature(b);
  if (!qb_result) {
    return std::unexpected(qb_result.error());
  }
  
  double fa = *qa_result - config_.target_heat_flux;
  double fb = *qb_result - config_.target_heat_flux;
  double fc = fb;
  double fs;
  
  std::cout << "\n>>> INITIAL BOUNDS EVALUATION <<<" << std::endl;
  std::cout << "Lower bound: T_edge = " << a << " K → q = " << *qa_result << " W/m² (error: " << fa << ")" << std::endl;
  std::cout << "Upper bound: T_edge = " << b << " K → q = " << *qb_result << " W/m² (error: " << fb << ")" << std::endl;
  
  // Check if target is within bounds
  if (fa * fb > 0) {
    return std::unexpected(solver::SolverError(
        std::format("Target heat flux {} W/m² is outside the range [{}, {}] W/m²",
                   config_.target_heat_flux, *qa_result, *qb_result)));
  }
  
  // Brent's method variables
  double R, S, T, U, Q;
  double s;
  
  int count = 0;
  int max_iter = config_.solver.max_iterations;
  double conv = std::pow(config_.solver.tolerance, 1.0/3.0);
  double step_tol = std::pow(config_.solver.tolerance, 2.0/3.0);
  
  while (count < max_iter && std::abs(fb) > conv) {
    count++;
    
    if (std::abs(fb) < conv || std::abs(b - a) < step_tol) {
      break;  // Converged
    }
    
    // Inverse quadratic interpolation
    if (std::abs(fa - fc) > 1e-10 && std::abs(fb - fc) > 1e-10) {
      R = fb / fc;
      S = fb / fa;
      T = fa / fc;
      
      U = S * (T * (R - T) * (c_brent - b) - (1 - R) * (b - a));
      Q = (T - 1) * (R - 1) * (S - 1);
      s = b + U / Q;
      std::cout << "  Using inverse quadratic interpolation → T_candidate=" << s << " K" << std::endl;
    } else {
      // Secant method
      s = b - fb * (b - a) / (fb - fa);
      std::cout << "  Using secant method → T_candidate=" << s << " K" << std::endl;
    }
    
    // Check conditions for bisection
    bool condition1 = (s < (3 * a + b) / 4) || (s > b);
    bool condition2 = mflag && (std::abs(s - b) >= std::abs(b - c_brent) / 2);
    bool condition3 = !mflag && (std::abs(s - b) >= std::abs(c_brent - d) / 2);
    bool condition4 = mflag && (std::abs(b - c_brent) < std::abs(step_tol));
    bool condition5 = !mflag && (std::abs(c_brent - d) < std::abs(step_tol));
    
    if (condition1 || condition2 || condition3 || condition4 || condition5) {
      s = (a + b) / 2;  // Bisection
      std::cout << "  Switching to bisection → T_candidate=" << s << " K" << std::endl;
      mflag = true;
    } else {
      mflag = false;
    }
    
    // Compute flux at new temperature
    std::cout << "\n>>> INVERSION STEP " << count << " <<<" << std::endl;
    std::cout << "Testing edge temperature: T_edge = " << s << " K" << std::endl;
    std::cout << "Computing boundary layer solution..." << std::flush;
    auto qs_result = compute_heat_flux_at_temperature(s);
    if (!qs_result) {
      return std::unexpected(qs_result.error());
    }
    fs = *qs_result - config_.target_heat_flux;
    
    std::cout << " DONE!" << std::endl;
    std::cout << "Heat flux calculated: q = " << *qs_result << " W/m²" << std::endl;
    std::cout << "Target flux:          q = " << config_.target_heat_flux << " W/m²" << std::endl;
    std::cout << "Error:                Δq = " << fs << " W/m²" << std::endl;
    std::cout << "Search interval:      [" << a << ", " << b << "] K" << std::endl;
    
    // Update intervals
    d = c_brent;
    c_brent = b;
    fc = fb;
    
    if (fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }
    
    // Swap to maintain |f(b)| <= |f(a)|
    if (std::abs(fa) < std::abs(fb)) {
      std::swap(a, b);
      std::swap(fa, fb);
    }
  }
  
  // Final temperature found
  double T_edge = b;
  
  std::cout << "\n=== Reconstruction Complete ===" << std::endl;
  std::cout << "Edge temperature found: " << T_edge << " K" << std::endl;
  std::cout << "Final heat flux: " << (*qb_result + fb) << " W/m²" << std::endl;
  std::cout << "Iterations used: " << count << std::endl;
  
  // Compute final edge conditions
  auto mass_fractions_result = mixture_.equilibrium_composition(
      T_edge, config_.boundary_conditions.pressure);
  if (!mass_fractions_result) {
    return std::unexpected(solver::SolverError(
        "Failed to compute equilibrium composition at reconstructed edge temperature"));
  }
  
  auto enthalpy_result = mixture_.mixture_enthalpy(
      *mass_fractions_result, T_edge, config_.boundary_conditions.pressure);
  if (!enthalpy_result) {
    return std::unexpected(solver::SolverError(
        "Failed to compute enthalpy at reconstructed edge temperature"));
  }
  
  // Compute density using ideal gas law: rho = P * Mmix / (R * T)
  auto mol_weight_result = mixture_.mixture_molecular_weight(*mass_fractions_result);
  if (!mol_weight_result) {
    return std::unexpected(solver::SolverError(
        "Failed to compute molecular weight at reconstructed edge temperature"));
  }
  
  const double R_universal = 8314.4621;  // J/(kmol·K)
  double density_edge = (config_.boundary_conditions.pressure * (*mol_weight_result)) / 
                        (R_universal * T_edge);
  
  auto viscosity_result = mixture_.viscosity(
      *mass_fractions_result, T_edge, config_.boundary_conditions.pressure);
  if (!viscosity_result) {
    return std::unexpected(solver::SolverError(
        "Failed to compute viscosity at reconstructed edge temperature"));
  }
  
  // Compute final heat flux at converged temperature
  auto final_heat_flux_result = compute_heat_flux_at_temperature(T_edge);
  if (!final_heat_flux_result) {
    return std::unexpected(solver::SolverError(
        "Failed to compute final heat flux at reconstructed edge temperature"));
  }
  
  return ReconstructedEdgeConditions{
      .temperature = T_edge,
      .pressure = config_.boundary_conditions.pressure,
      .enthalpy = *enthalpy_result,
      .mass_fractions = *mass_fractions_result,
      .density = density_edge,
      .viscosity = *viscosity_result,
      .heat_flux_achieved = *final_heat_flux_result,
      .iterations_used = count
  };
}

auto EdgeTemperatureReconstructor::compute_heat_flux_at_temperature(double T_edge)
    -> std::expected<double, solver::SolverError> {
  
  // Setup boundary conditions with current T_edge
  auto temp_config_result = setup_boundary_conditions(T_edge);
  if (!temp_config_result) {
    return std::unexpected(temp_config_result.error());
  }
  
  // Create solver  
  solver::BoundaryLayerSolver solver(mixture_, temp_config_result.value());
  
  // Solve boundary layer
  auto solution = solver.solve();
  if (!solution) {
    return std::unexpected(solver::SolverError(
        std::format("Failed to solve boundary layer at T_edge={} K: {}", 
                   T_edge, solution.error().message())));
  }
  
  // Return total heat flux at wall
  if (solution->heat_flux_data.empty()) {
    return std::unexpected(solver::SolverError(
        "No heat flux data available in solution"));
  }
  
  return solution->heat_flux_data[0].q_wall_total_dim;
}

auto EdgeTemperatureReconstructor::setup_boundary_conditions(double T_edge)
    -> std::expected<io::Configuration, solver::SolverError> {
  
  // Get equilibrium composition at edge
  auto mass_fractions_result = mixture_.equilibrium_composition(
      T_edge, config_.boundary_conditions.pressure);
  if (!mass_fractions_result) {
    return std::unexpected(solver::SolverError(
        std::format("Failed to compute equilibrium composition at T_edge={} K", T_edge)));
  }
  
  // Compute enthalpy at edge
  auto enthalpy_result = mixture_.mixture_enthalpy(
      *mass_fractions_result, T_edge, config_.boundary_conditions.pressure);
  if (!enthalpy_result) {
    return std::unexpected(solver::SolverError(
        std::format("Failed to compute enthalpy at T_edge={} K", T_edge)));
  }
  
  // Create modified configuration for this iteration
  io::Configuration temp_config = full_config_;
  
  // Override with reconstruction parameters
  temp_config.simulation.finite_thickness = true;
  temp_config.simulation.catalytic_wall = (config_.boundary_conditions.catalyticity > 0);
  temp_config.wall_parameters.catalytic_wall = (config_.boundary_conditions.catalyticity > 0);
  
  // Set edge conditions
  temp_config.outer_edge.edge_points[0].temperature = T_edge;
  temp_config.outer_edge.edge_points[0].pressure = config_.boundary_conditions.pressure;
  temp_config.outer_edge.edge_points[0].enthalpy = *enthalpy_result;
  temp_config.outer_edge.edge_points[0].velocity = 0.0;  // Stagnation point
  
  // Set flow parameters
  temp_config.outer_edge.velocity_gradient_stagnation = 
      config_.flow_parameters.velocity_gradient_stagnation;
  temp_config.outer_edge.freestream_density = 
      config_.flow_parameters.freestream_density;
  temp_config.outer_edge.freestream_velocity = 
      config_.flow_parameters.freestream_velocity;
  
  // Set finite thickness parameters
  temp_config.outer_edge.finite_thickness_params.v_edge = 
      config_.finite_thickness_params.v_edge;
  temp_config.outer_edge.finite_thickness_params.d2_ue_dxdy = 
      config_.finite_thickness_params.d2_ue_dxdy;
  temp_config.outer_edge.finite_thickness_params.delta_bl = 
      config_.finite_thickness_params.delta_bl;
  
  // Set wall conditions
  temp_config.wall_parameters.wall_temperatures[0] = 
      config_.boundary_conditions.wall_temperature;
  
  return temp_config;
}

} // namespace blast::boundary_layer