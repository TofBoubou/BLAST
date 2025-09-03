#include "blast/boundary_layer/catalysis_reconstruction.hpp"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>

namespace blast::boundary_layer {

auto CatalysisReconstructor::reconstruct() 
    -> std::expected<ReconstructedCatalysisConditions, solver::SolverError> {
  
  // Backup GSI file (will be modified with different catalyticity values)
  auto backup_result = gsi_manager_.backup_gsi_file();
  if (!backup_result) {
    std::cerr << "Warning: Failed to backup GSI file: " << backup_result.error() << std::endl;
  }
  
  // Brent's method parameters
  bool mflag = true;
  double a = config_.solver.catalyticity_min;
  double b = config_.solver.catalyticity_max;
  double d = 0;
  double c_brent = a;
  
  if (full_config_.verbose || std::getenv("BLAST_VERBOSE")) {
    std::cout << "\n=== Catalysis Reconstruction ===" << std::endl;
    std::cout << "Target heat flux: " << config_.target_heat_flux << " W/m²" << std::endl;
    std::cout << "Fixed edge temperature: " << config_.boundary_conditions.edge_temperature << " K" << std::endl;
    std::cout << "Catalyticity search range: [" << a << ", " << b << "]" << std::endl;
  }
  
  // Compute heat flux at bounds
  auto qa_result = compute_heat_flux_at_catalyticity(a);
  if (!qa_result) {
    return std::unexpected(qa_result.error());
  }
  
  auto qb_result = compute_heat_flux_at_catalyticity(b);
  if (!qb_result) {
    return std::unexpected(qb_result.error());
  }
  
  double fa = *qa_result - config_.target_heat_flux;
  double fb = *qb_result - config_.target_heat_flux;
  double fc = fb;
  double fs;
  
  if (full_config_.verbose || std::getenv("BLAST_VERBOSE")) {
    std::cout << "\n>>> INITIAL BOUNDS EVALUATION <<<" << std::endl;
    std::cout << "Lower bound: γ = " << a << " → q = " << *qa_result << " W/m² (error: " << fa << ")" << std::endl;
    std::cout << "Upper bound: γ = " << b << " → q = " << *qb_result << " W/m² (error: " << fb << ")" << std::endl;
  }
  
  // Check if target is within bounds
  if (fa * fb > 0) {
    std::string error_msg = std::format(
        "\n"
        "═══════════════════════════════════════════════════════════════════════════════════\n"
        "                        CATALYTICITY RANGE ERROR\n"
        "═══════════════════════════════════════════════════════════════════════════════════\n"
        "\n"
        "Target heat flux cannot be achieved within the specified catalyticity range!\n"
        "\n"
        "  Target heat flux:            {} W/m²\n"
        "  Flux at γ_min ({}):         {} W/m²\n"
        "  Flux at γ_max ({}):         {} W/m²\n"
        "  Current flux range:          [{}, {}] W/m²\n"
        "\n"
        "SOLUTION: Update your YAML configuration file\n"
        "\n",
        config_.target_heat_flux,
        config_.solver.catalyticity_min, *qa_result,
        config_.solver.catalyticity_max, *qb_result,
        std::min(*qa_result, *qb_result), std::max(*qa_result, *qb_result));
    
    error_msg += "\n═══════════════════════════════════════════════════════════════════════════════════\n";
    
    return std::unexpected(solver::SolverError(error_msg));
  }
  
  // Brent's method variables
  double R, S, T, U, Q;
  double s;
  
  int count = 0;
  int max_iter = config_.solver.max_iterations;
  double conv = std::pow(config_.solver.tolerance, 1.0/3.0);
  double step_tol = std::max(0.001, std::pow(config_.solver.tolerance, 2.0/3.0)); // Minimum step for catalyticity
  
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
      if (full_config_.verbose || std::getenv("BLAST_VERBOSE")) { std::cout << "  Using inverse quadratic interpolation → γ_candidate=" << s << std::endl; }
    } else {
      // Secant method
      s = b - fb * (b - a) / (fb - fa);
      if (full_config_.verbose || std::getenv("BLAST_VERBOSE")) { std::cout << "  Using secant method → γ_candidate=" << s << std::endl; }
    }
    
    bool condition2 = mflag && (std::abs(s - b) >= std::abs(b - c_brent) * 0.75);
    bool condition3 = !mflag && (std::abs(s - b) >= std::abs(c_brent - d) * 0.75);
    bool condition4 = mflag && (std::abs(b - c_brent) < std::abs(step_tol));
    bool condition5 = !mflag && (std::abs(c_brent - d) < std::abs(step_tol));
    
    if (condition2 || condition3 || condition4 || condition5) {
      s = (a + b) / 2;
      if (full_config_.verbose || std::getenv("BLAST_VERBOSE")) { std::cout << "  Switching to bisection → γ_candidate=" << s << std::endl; }
      mflag = true;
    } else {
      mflag = false;
    }
    
    // Compute heat flux at new catalyticity
    if (full_config_.verbose || std::getenv("BLAST_VERBOSE")) {
      std::cout << "\n>>> INVERSION STEP " << count << " <<<" << std::endl;
      std::cout << "Testing catalyticity: γ = " << s << std::endl;
      std::cout << "Computing boundary layer solution..." << std::flush;
    }
    auto qs_result = compute_heat_flux_at_catalyticity(s);
    if (!qs_result) {
      return std::unexpected(qs_result.error());
    }
    fs = *qs_result - config_.target_heat_flux;
    
    if (full_config_.verbose || std::getenv("BLAST_VERBOSE")) {
      std::cout << " DONE!" << std::endl;
      std::cout << "Heat flux calculated: q = " << *qs_result << " W/m²" << std::endl;
      std::cout << "Target flux:          q = " << config_.target_heat_flux << " W/m²" << std::endl;
      std::cout << "Error:                Δq = " << fs << " W/m²" << std::endl;
      std::cout << "Search interval:      [" << a << ", " << b << "]" << std::endl;
    }
    
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
  
  // Final catalyticity found
  double optimal_catalyticity = b;
  
  if (full_config_.verbose || std::getenv("BLAST_VERBOSE")) {
    std::cout << "\n=== Reconstruction Complete ===" << std::endl;
    std::cout << "Optimal catalyticity found: " << optimal_catalyticity << std::endl;
    std::cout << "Final heat flux: " << (*qb_result + fb) << " W/m²" << std::endl;
    std::cout << "Iterations used: " << count << std::endl;
  }
  
  // Generate complete boundary layer solution at optimal catalyticity
  if (full_config_.verbose || std::getenv("BLAST_VERBOSE")) {
    std::cout << "\n=== Generating Complete Solution ===" << std::endl;
    std::cout << "Computing full boundary layer solution at γ = " << optimal_catalyticity << "..." << std::endl;
  }
  
  auto final_config_result = setup_boundary_conditions(optimal_catalyticity);
  if (!final_config_result) {
    return std::unexpected(final_config_result.error());
  }
  
  solver::BoundaryLayerSolver final_solver(mixture_, final_config_result.value());
  auto final_solution_result = final_solver.solve();
  if (!final_solution_result) {
    return std::unexpected(solver::SolverError(
        std::format("Failed to generate complete solution at γ={}: {}", 
                   optimal_catalyticity, final_solution_result.error().message())));
  }
  
  // Get final heat flux from solution
  double final_heat_flux = 0.0;
  if (!final_solution_result->heat_flux_data.empty()) {
    final_heat_flux = final_solution_result->heat_flux_data[0].q_wall_total_dim;
  }
  
  // Use fixed edge temperature for final conditions
  double final_edge_temperature = config_.boundary_conditions.edge_temperature;
  
  // Compute edge conditions at fixed edge temperature
  auto mass_fractions_result = mixture_.equilibrium_composition(
      final_edge_temperature, config_.boundary_conditions.pressure);
  if (!mass_fractions_result) {
    return std::unexpected(solver::SolverError(
        "Failed to compute equilibrium composition at reconstructed catalyticity"));
  }
  
  auto enthalpy_result = mixture_.mixture_enthalpy(
      *mass_fractions_result, final_edge_temperature, config_.boundary_conditions.pressure);
  if (!enthalpy_result) {
    return std::unexpected(solver::SolverError(
        "Failed to compute enthalpy at reconstructed catalyticity"));
  }
  
  // Compute density using ideal gas law: rho = P * Mmix / (R * T)
  auto mol_weight_result = mixture_.mixture_molecular_weight(*mass_fractions_result);
  if (!mol_weight_result) {
    return std::unexpected(solver::SolverError(
        "Failed to compute molecular weight at reconstructed catalyticity"));
  }
  
  const double R_universal = 8314.4621;  // J/(kmol·K)
  double density_edge = (config_.boundary_conditions.pressure * (*mol_weight_result)) / 
                        (R_universal * final_edge_temperature);
  
  auto viscosity_result = mixture_.viscosity(
      *mass_fractions_result, final_edge_temperature, config_.boundary_conditions.pressure);
  if (!viscosity_result) {
    return std::unexpected(solver::SolverError(
        "Failed to compute viscosity at reconstructed catalyticity"));
  }
  
  if (full_config_.verbose || std::getenv("BLAST_VERBOSE")) {
    std::cout << "Complete solution generated successfully!" << std::endl;
    std::cout << "Final heat flux verification: " << final_heat_flux << " W/m²" << std::endl;
  }
  
  // Restore GSI file before returning
  gsi_manager_.restore_gsi_file();
  
  return ReconstructedCatalysisConditions{
      .catalyticity = optimal_catalyticity,
      .heat_flux_achieved = final_heat_flux,
      .edge_temperature = final_edge_temperature,
      .pressure = config_.boundary_conditions.pressure,
      .enthalpy = *enthalpy_result,
      .mass_fractions = *mass_fractions_result,
      .density = density_edge,
      .viscosity = *viscosity_result,
      .iterations_used = count,
      .full_solution = std::move(final_solution_result.value())
  };
}

auto CatalysisReconstructor::compute_heat_flux_at_catalyticity(double catalyticity)
    -> std::expected<double, solver::SolverError> {
  
  // Setup boundary conditions with current catalyticity
  auto temp_config_result = setup_boundary_conditions(catalyticity);
  if (!temp_config_result) {
    return std::unexpected(temp_config_result.error());
  }
  
  // Create solver  
  solver::BoundaryLayerSolver solver(mixture_, temp_config_result.value());
  
  // Solve boundary layer
  auto solution = solver.solve();
  if (!solution) {
    return std::unexpected(solver::SolverError(
        std::format("Failed to solve boundary layer at γ={}: {}", 
                   catalyticity, solution.error().message())));
  }
  
  // Return total heat flux at wall
  if (solution->heat_flux_data.empty()) {
    return std::unexpected(solver::SolverError(
        "No heat flux data available in solution"));
  }
  
  return solution->heat_flux_data[0].q_wall_total_dim;
}

auto CatalysisReconstructor::setup_boundary_conditions(double catalyticity)
    -> std::expected<io::Configuration, solver::SolverError> {
  
  // Update GSI file with the current catalyticity value
  if (catalyticity > 0) {
    auto update_result = gsi_manager_.update_gsi_catalyticity(catalyticity);
    if (!update_result) {
      return std::unexpected(solver::SolverError(
          std::format("Failed to update GSI catalyticity to {}: {}", catalyticity, update_result.error())));
    }
    
    // Reload mixture with updated GSI
    auto reload_result = mixture_.reload_gsi();
    if (!reload_result) {
      return std::unexpected(solver::SolverError(
          std::format("Failed to reload mixture after GSI update: {}", reload_result.error())));
    }
  }
  
  // Get equilibrium composition at fixed edge temperature
  auto mass_fractions_result = mixture_.equilibrium_composition(
      config_.boundary_conditions.edge_temperature, config_.boundary_conditions.pressure);
  if (!mass_fractions_result) {
    return std::unexpected(solver::SolverError(
        std::format("Failed to compute equilibrium composition at T_edge={} K", config_.boundary_conditions.edge_temperature)));
  }
  
  // Compute enthalpy at edge temperature
  auto enthalpy_result = mixture_.mixture_enthalpy(
      *mass_fractions_result, config_.boundary_conditions.edge_temperature, config_.boundary_conditions.pressure);
  if (!enthalpy_result) {
    return std::unexpected(solver::SolverError(
        std::format("Failed to compute enthalpy at T_edge={} K", config_.boundary_conditions.edge_temperature)));
  }
  
  // Create modified configuration for this iteration
  io::Configuration temp_config = full_config_;
  
  // Override with reconstruction parameters
  temp_config.simulation.finite_thickness = true;
  temp_config.simulation.catalytic_wall = true;  // Always catalytic for catalysis reconstruction
  temp_config.simulation.catalysis_provider = io::SimulationConfig::CatalysisProvider::MutationPP;
  
  // Set edge conditions with fixed edge temperature
  if (temp_config.outer_edge.edge_points.empty()) {
    io::OuterEdgeConfig::EdgePoint edge_point;
    edge_point.x = 0.0;  // Stagnation point
    edge_point.temperature = config_.boundary_conditions.edge_temperature;
    edge_point.pressure = config_.boundary_conditions.pressure;
    edge_point.enthalpy = *enthalpy_result;
    edge_point.velocity = 0.0;  // Always 0 for catalysis_reconstruction
    edge_point.radius = config_.boundary_conditions.radius;
    temp_config.outer_edge.edge_points.push_back(edge_point);
  } else {
    temp_config.outer_edge.edge_points[0].temperature = config_.boundary_conditions.edge_temperature;
    temp_config.outer_edge.edge_points[0].pressure = config_.boundary_conditions.pressure;
    temp_config.outer_edge.edge_points[0].enthalpy = *enthalpy_result;
    temp_config.outer_edge.edge_points[0].velocity = 0.0;  // Stagnation point
  }
  
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
  
  // Set wall conditions - fixed wall temperature
  if (temp_config.wall_parameters.wall_temperatures.empty()) {
    temp_config.wall_parameters.wall_temperatures.push_back(
        config_.boundary_conditions.wall_temperature);
  } else {
    temp_config.wall_parameters.wall_temperatures[0] = 
        config_.boundary_conditions.wall_temperature;
  }
  
  return temp_config;
}

} // namespace blast::boundary_layer