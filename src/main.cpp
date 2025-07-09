#include "blast/io/config_manager.hpp"
#include "blast/thermophysics/mixture_interface.hpp"
#include "blast/boundary_layer/solver/boundary_layer_solver.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file.yaml>\n";
        return 1;
    }
    
    try {
        std::cout << "=== BLAST Boundary Layer Solver ===" << std::endl;
        std::cout << "Loading configuration from: " << argv[1] << std::endl;
        
        // Load configuration
        blast::io::ConfigurationManager config_manager;
        auto config_result = config_manager.load(argv[1]);
        if (!config_result) {
            std::cerr << "Failed to load config: " << config_result.error().message() << "\n";
            return 1;
        }
        auto config = config_result.value();
        std::cout << "✓ Configuration loaded successfully" << std::endl;
        
        // Create mixture
        std::cout << "Creating mixture: " << config.mixture.name << std::endl;
        auto mixture_result = blast::thermophysics::create_mixture(config.mixture);
        if (!mixture_result) {
            std::cerr << "Failed to create mixture: " << mixture_result.error().message() << "\n";
            return 1;
        }
        auto& mixture = *mixture_result.value();
        std::cout << "✓ Mixture created (" << mixture.n_species() << " species)" << std::endl;
        
        // Print mixture information
        std::cout << "\nMixture species:" << std::endl;
        for (std::size_t i = 0; i < mixture.n_species(); ++i) {
            std::cout << "  [" << std::setw(2) << i << "] " 
                      << std::setw(8) << mixture.species_name(i) 
                      << " (MW: " << std::setw(8) << std::fixed << std::setprecision(3) 
                      << mixture.species_molecular_weight(i) << " kg/kmol)" << std::endl;
        }
        
        // Print simulation configuration
        std::cout << "\nSimulation setup:" << std::endl;
        std::cout << "  Body type: " << (config.simulation.body_type == blast::io::SimulationConfig::BodyType::Axisymmetric ? "Axisymmetric" : "Other") << std::endl;
        std::cout << "  Stagnation only: " << (config.simulation.only_stagnation_point ? "Yes" : "No") << std::endl;
        std::cout << "  Chemical non-eq: " << (config.simulation.chemical_non_equilibrium ? "Yes" : "No") << std::endl;
        std::cout << "  Thermal diffusion: " << (config.simulation.consider_thermal_diffusion ? "Yes" : "No") << std::endl;
        std::cout << "  Grid points (η): " << config.numerical.n_eta << std::endl;
        std::cout << "  η_max: " << config.numerical.eta_max << std::endl;
        std::cout << "  Convergence tol: " << config.numerical.convergence_tolerance << std::endl;
        
        // Print edge conditions
        if (!config.outer_edge.edge_points.empty()) {
            const auto& edge = config.outer_edge.edge_points[0];
            std::cout << "\nEdge conditions:" << std::endl;
            std::cout << "  Pressure: " << edge.pressure << " Pa" << std::endl;
            std::cout << "  Temperature: " << edge.enthalpy / 1000.0 << " kJ/kg (enthalpy)" << std::endl;
            std::cout << "  Density: " << edge.density << " kg/m³" << std::endl;
            std::cout << "  Viscosity: " << edge.viscosity << " Pa·s" << std::endl;
            std::cout << "  Wall temp: " << config.wall_parameters.wall_temperatures[0] << " K" << std::endl;
        }
        
        // Create and run solver
        std::cout << "\n=== STARTING BOUNDARY LAYER SOLUTION ===" << std::endl;
        blast::boundary_layer::solver::BoundaryLayerSolver solver(mixture, config);
        
        auto solution_result = solver.solve();
        if (!solution_result) {
            std::cerr << "Solver failed: " << solution_result.error().message() << "\n";
            return 1;
        }
        
        auto& solution = solution_result.value();
        std::cout << "✓ Solution completed successfully!" << std::endl;
        
        // Print solution summary
        std::cout << "\n=== SOLUTION SUMMARY ===" << std::endl;
        std::cout << "Stations solved: " << solution.stations.size() << std::endl;
        std::cout << "Total iterations: " << solution.total_iterations << std::endl;
        std::cout << "Converged: " << (solution.converged ? "Yes" : "No") << std::endl;
        
        if (!solution.stations.empty()) {
            const auto& final_station = solution.stations.back();
            const auto n_eta = final_station.F.size();
            
            std::cout << "\nFinal station results:" << std::endl;
            std::cout << "  ξ = " << solution.xi_solved.back() << std::endl;
            std::cout << "  F at edge: " << final_station.F.back() << std::endl;
            std::cout << "  g at wall: " << final_station.g[0] << std::endl;
            std::cout << "  g at edge: " << final_station.g.back() << std::endl;
            
            // Print species concentrations at wall and edge
            std::cout << "\nSpecies concentrations:" << std::endl;
            std::cout << "  Wall → Edge" << std::endl;
            for (std::size_t i = 0; i < mixture.n_species(); ++i) {
                std::cout << "  " << std::setw(8) << mixture.species_name(i) 
                          << ": " << std::setw(8) << std::fixed << std::setprecision(4)
                          << final_station.c(i, 0) << " → " 
                          << std::setw(8) << final_station.c(i, n_eta-1) << std::endl;
            }
            
            // Print some profiles
            std::cout << "\nBoundary layer profiles (selected points):" << std::endl;
            std::cout << std::setw(6) << "η/η_max" 
                      << std::setw(10) << "F" 
                      << std::setw(10) << "g" 
                      << std::setw(10) << "V" << std::endl;
            std::cout << std::string(36, '-') << std::endl;
            
            const std::size_t n_print = std::min(static_cast<std::size_t>(11), n_eta);
            for (std::size_t i = 0; i < n_print; ++i) {
                const std::size_t idx = (i * (n_eta - 1)) / (n_print - 1);
                const double eta_norm = static_cast<double>(idx) / (n_eta - 1);
                
                std::cout << std::setw(6) << std::fixed << std::setprecision(2) << eta_norm
                          << std::setw(10) << std::fixed << std::setprecision(4) << final_station.F[idx]
                          << std::setw(10) << std::fixed << std::setprecision(4) << final_station.g[idx]
                          << std::setw(10) << std::fixed << std::setprecision(4) << final_station.V[idx]
                          << std::endl;
            }
        }
        
        std::cout << "\n=== CALCULATION COMPLETED SUCCESSFULLY ===" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
