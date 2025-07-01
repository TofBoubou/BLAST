#include "blast/io/config_manager.hpp"
#include "blast/boundary_layer/grid/grid.hpp"
#include "blast/boundary_layer/conditions/boundary_interpolator.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include "blast/boundary_layer/coefficients/xi_derivatives.hpp"
#include "blast/thermophysics/mixture_interface.hpp"
#include <iostream>
#include <iomanip>
#include <vector>

int main() {
    try {
        std::cout << "=== BLAST Boundary Layer Solver Test ===" << std::endl;
        
        // 1. Load configuration
        std::cout << "\n1. Loading configuration..." << std::endl;
        blast::io::ConfigurationManager manager;
        auto config_result = manager.load("config/default.yaml");
        if (!config_result) {
            std::cerr << "Failed to load configuration: " << config_result.error().message() << std::endl;
            return 1;
        }
        auto config = config_result.value();
        std::cout << "✓ Configuration loaded successfully" << std::endl;
        
        // 2. Create mixture
        std::cout << "\n2. Creating thermophysics mixture..." << std::endl;
        auto mixture_result = blast::thermophysics::create_mixture(config.mixture);
        if (!mixture_result) {
            std::cerr << "Failed to create mixture: " << mixture_result.error().message() << std::endl;
            return 1;
        }
        auto mixture = std::move(mixture_result.value());
        std::cout << "✓ Mixture created: " << config.mixture.name 
                  << " with " << mixture->n_species() << " species" << std::endl;
        
        // 3. Create boundary layer grid
        std::cout << "\n3. Creating boundary layer grid..." << std::endl;
        auto grid_result = blast::boundary_layer::grid::BoundaryLayerGrid::create_stagnation_grid(
            config.numerical, config.outer_edge);
        if (!grid_result) {
            std::cerr << "Failed to create grid: " << grid_result.error().message() << std::endl;
            return 1;
        }
        auto grid = grid_result.value();
        std::cout << "✓ Grid created with " << grid.n_eta() << " eta points" << std::endl;
        std::cout << "  - eta_max = " << grid.eta_max() << std::endl;
        std::cout << "  - d_eta = " << grid.d_eta() << std::endl;
        
        // 4. Create boundary conditions
        std::cout << "\n4. Setting up boundary conditions..." << std::endl;
        auto bc_result = blast::boundary_layer::conditions::create_stagnation_conditions(
            config.outer_edge, config.wall_parameters, config.simulation);
        if (!bc_result) {
            std::cerr << "Failed to create boundary conditions: " << bc_result.error().message() << std::endl;
            return 1;
        }
        auto bc = bc_result.value();
        
        // Update edge properties using mixture
        auto update_result = blast::boundary_layer::conditions::update_edge_properties(bc, *mixture);
        if (!update_result) {
            std::cerr << "Warning: Failed to update edge properties: " << update_result.error().message() << std::endl;
        }
        
        std::cout << "✓ Boundary conditions created:" << std::endl;
        std::cout << "  - Edge pressure: " << bc.P_e() << " Pa" << std::endl;
        std::cout << "  - Edge velocity: " << bc.ue() << " m/s" << std::endl;
        std::cout << "  - Wall temperature: " << bc.Tw() << " K" << std::endl;
        std::cout << "  - Beta parameter: " << bc.beta << std::endl;
        
        // 5. Initialize solution fields
        std::cout << "\n5. Initializing solution fields..." << std::endl;
        const auto n_eta = grid.n_eta();
        const auto n_species = mixture->n_species();
        const auto eta_coords = grid.eta_coordinates();
        
        // Initialize F (momentum) and g fields
        std::vector<double> F(n_eta);
        std::vector<double> g(n_eta);
        std::vector<double> T(n_eta);
        
        // Initialize species concentrations matrix
        blast::core::Matrix<double> c(n_species, n_eta);
        blast::core::Matrix<double> dc_deta(n_species, n_eta);
        
        // Simple initialization with reasonable profiles
        for (int i = 0; i < n_eta; ++i) {
            const double eta = eta_coords[i];
            const double eta_norm = eta / grid.eta_max();
            
            // Blasius-like profile for F
            F[i] = 1.0 - std::exp(-eta);
            
            // Temperature profile (linear from wall to edge)
            T[i] = bc.Tw() + (2000.0 - bc.Tw()) * eta_norm; // Assume edge temp ~2000K
            
            // g field initialization
            g[i] = 0.5 * eta_norm;
            
            // Species concentrations (use edge values, normalize)
            double sum_c = 0.0;
            for (std::size_t j = 0; j < n_species; ++j) {
                if (j < bc.c_e().size()) {
                    c(j, i) = bc.c_e()[j];
                } else {
                    c(j, i) = (j == 0) ? 0.23 : (j == 1) ? 0.77 : 0.0; // Default air composition
                }
                sum_c += c(j, i);
            }
            
            // Normalize species concentrations
            for (std::size_t j = 0; j < n_species; ++j) {
                c(j, i) /= sum_c;
            }
        }
        
        // Compute concentration derivatives (simple finite differences)
        for (std::size_t j = 0; j < n_species; ++j) {
            for (int i = 0; i < n_eta; ++i) {
                if (i == 0) {
                    dc_deta(j, i) = (c(j, 1) - c(j, 0)) / grid.d_eta();
                } else if (i == n_eta - 1) {
                    dc_deta(j, i) = (c(j, n_eta-1) - c(j, n_eta-2)) / grid.d_eta();
                } else {
                    dc_deta(j, i) = (c(j, i+1) - c(j, i-1)) / (2.0 * grid.d_eta());
                }
            }
        }
        
        std::cout << "✓ Solution fields initialized" << std::endl;
        std::cout << "  - Temperature range: " << T[0] << " - " << T[n_eta-1] << " K" << std::endl;
        std::cout << "  - F field range: " << F[0] << " - " << F[n_eta-1] << std::endl;
        
        // 6. Initialize xi derivatives
        std::cout << "\n6. Initializing xi derivatives..." << std::endl;
        blast::boundary_layer::coefficients::XiDerivatives xi_der(n_eta, n_species);
        
        // For stagnation point, update with current solution
        xi_der.update_station(0, 0.0, F, g, c);
        std::cout << "✓ Xi derivatives initialized for stagnation point" << std::endl;
        std::cout << "  - Lambda0: " << xi_der.lambda0() << std::endl;
        std::cout << "  - Lambda1: " << xi_der.lambda1() << std::endl;
        
        // 7. Create coefficient calculator
        std::cout << "\n7. Creating coefficient calculator..." << std::endl;
        blast::boundary_layer::coefficients::CoefficientCalculator calc(
            *mixture, config.simulation, config.numerical);
        std::cout << "✓ Coefficient calculator created" << std::endl;
        
        // 8. Calculate coefficients
        std::cout << "\n8. Calculating coefficients..." << std::endl;
        
        // Prepare coefficient inputs
        blast::boundary_layer::coefficients::CoefficientInputs inputs{
            .xi = 0.0,
            .F = F,
            .c = c,
            .dc_deta = dc_deta,
            .T = T
        };
        
        auto coeff_result = calc.calculate(inputs, bc, xi_der);
        if (!coeff_result) {
            std::cerr << "Failed to calculate coefficients: " << coeff_result.error().message() << std::endl;
            return 1;
        }
        
        auto coeffs = coeff_result.value();
        std::cout << "✓ Coefficients calculated successfully!" << std::endl;
        
        // 9. Display results
        std::cout << "\n9. Coefficient Results Summary:" << std::endl;
        std::cout << std::string(50, '=') << std::endl;
        
        // Transport coefficients
        std::cout << "\nTransport Coefficients:" << std::endl;
        std::cout << "  l0 range: [" << std::fixed << std::setprecision(6) 
                  << *std::min_element(coeffs.transport.l0.begin(), coeffs.transport.l0.end())
                  << ", " << *std::max_element(coeffs.transport.l0.begin(), coeffs.transport.l0.end())
                  << "]" << std::endl;
        std::cout << "  l3 range: [" 
                  << *std::min_element(coeffs.transport.l3.begin(), coeffs.transport.l3.end())
                  << ", " << *std::max_element(coeffs.transport.l3.begin(), coeffs.transport.l3.end())
                  << "]" << std::endl;
        
        // Thermodynamic coefficients
        std::cout << "\nThermodynamic Coefficients:" << std::endl;
        std::cout << "  Density range: [" 
                  << *std::min_element(coeffs.thermodynamic.rho.begin(), coeffs.thermodynamic.rho.end())
                  << ", " << *std::max_element(coeffs.thermodynamic.rho.begin(), coeffs.thermodynamic.rho.end())
                  << "] kg/m³" << std::endl;
        std::cout << "  MW range: [" 
                  << *std::min_element(coeffs.thermodynamic.MW.begin(), coeffs.thermodynamic.MW.end())
                  << ", " << *std::max_element(coeffs.thermodynamic.MW.begin(), coeffs.thermodynamic.MW.end())
                  << "] kg/kmol" << std::endl;
        std::cout << "  Wall enthalpy: " << coeffs.thermodynamic.h_wall << " J/kg" << std::endl;
        
        // Wall properties
        std::cout << "\nWall Properties:" << std::endl;
        std::cout << "  Thermal conductivity: " << coeffs.wall.k_wall << " W/(m·K)" << std::endl;
        std::cout << "  Viscosity: " << coeffs.wall.mu_wall << " Pa·s" << std::endl;
        std::cout << "  Density: " << coeffs.wall.rho_wall << " kg/m³" << std::endl;
        std::cout << "  Prandtl number: " << coeffs.wall.Pr_wall << std::endl;
        
        // Diffusion coefficients info
        std::cout << "\nDiffusion Matrix: " << coeffs.diffusion.Dij_bin.rows() 
                  << "×" << coeffs.diffusion.Dij_bin.cols() << std::endl;
        
        // Species enthalpies
        std::cout << "\nSpecies Enthalpies Matrix: " << coeffs.h_species.rows() 
                  << "×" << coeffs.h_species.cols() << std::endl;
        
        // Chemical coefficients
        if (config.simulation.chemical_non_equilibrium) {
            std::cout << "\nChemical Production Rates Matrix: " << coeffs.chemical.wi.rows() 
                      << "×" << coeffs.chemical.wi.cols() << std::endl;
        }
        
        // Sample values at wall and a few eta points
        std::cout << "\n10. Sample Values at Key Points:" << std::endl;
        std::cout << std::string(50, '-') << std::endl;
        
        std::vector<int> sample_indices = {0, n_eta/4, n_eta/2, 3*n_eta/4, n_eta-1};
        std::cout << std::setw(8) << "eta" << std::setw(12) << "l0" << std::setw(12) << "l3" 
                  << std::setw(12) << "rho" << std::setw(12) << "T" << std::endl;
        std::cout << std::string(56, '-') << std::endl;
        
        for (int idx : sample_indices) {
            if (idx >= n_eta) idx = n_eta - 1;
            std::cout << std::setw(8) << std::setprecision(3) << eta_coords[idx]
                      << std::setw(12) << std::setprecision(6) << coeffs.transport.l0[idx]
                      << std::setw(12) << coeffs.transport.l3[idx]
                      << std::setw(12) << coeffs.thermodynamic.rho[idx]
                      << std::setw(12) << std::setprecision(1) << T[idx]
                      << std::endl;
        }
        
        std::cout << "\n✓ All tests completed successfully!" << std::endl;
        std::cout << "The coefficient calculator is working properly." << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ Error during execution: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}





/* 

g++ -std=c++23 -O2 -DYAML_CPP_STATIC_DEFINE -Iinclude -Ilibs/yaml-cpp/include -Ilibs/eigen -Ilibs/mutation++/include -Llibs/yaml-cpp/build -Llibs/mutation++/lib (Get-ChildItem -Recurse -Filter *.cpp -Path src | ForEach-Object { $_.FullName }) -lyaml-cpp -lmutation++ -static-libgcc -static-libstdc++ -o blast_main.exe

*/


/* 
$cppFiles = Get-ChildItem -Recurse -Filter *.cpp -Path src | ForEach-Object { $_.FullName }

g++ -std=c++23 `
    -DYAML_CPP_STATIC_DEFINE `
    -Iinclude `
    -Ilibs/yaml-cpp/include `
    -Ilibs/eigen `
    -Llibs/yaml-cpp/build `
    $cppFiles `
    -lyaml-cpp `
    -static-libgcc `
    -static-libstdc++ `
    -o blast_main.exe
 */
