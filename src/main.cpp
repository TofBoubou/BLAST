#include "blast/io/config_manager.hpp"
#include "blast/thermophysics/mixture_interface.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include "blast/boundary_layer/coefficients/xi_derivatives.hpp"
#include "blast/boundary_layer/conditions/boundary_interpolator.hpp"
#include <iostream>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file.yaml>\n";
        return 1;
    }
    
    try {
        // Load configuration
        blast::io::ConfigurationManager config_manager;
        auto config_result = config_manager.load(argv[1]);
        if (!config_result) {
            std::cerr << "Failed to load config: " << config_result.error().message() << "\n";
            return 1;
        }
        auto config = config_result.value();
        
        // Create mixture
        auto mixture_result = blast::thermophysics::create_mixture(config.mixture);
        if (!mixture_result) {
            std::cerr << "Failed to create mixture: " << mixture_result.error().message() << "\n";
            return 1;
        }
        auto& mixture = *mixture_result.value();
        
        // Create coefficient calculator
        blast::boundary_layer::coefficients::CoefficientCalculator calculator(
            mixture, config.simulation, config.numerical
        );
        
        // Initialize test data
        const std::size_t n_eta = config.numerical.n_eta;
        const std::size_t n_species = mixture.n_species();
        
        // Create boundary conditions for stagnation point
        auto bc_result = blast::boundary_layer::conditions::create_stagnation_conditions(
            config.outer_edge, config.wall_parameters, config.simulation
        );
        if (!bc_result) {
            std::cerr << "Failed to create boundary conditions: " << bc_result.error().message() << "\n";
            return 1;
        }
        auto bc = bc_result.value();
        
        // Initialize solution fields with test data
        std::vector<double> F(n_eta, 0.0);
        std::vector<double> T(n_eta);
        blast::core::Matrix<double> c(n_species, n_eta);
        blast::core::Matrix<double> dc_deta(n_species, n_eta);
        
        // Fill with linear temperature profile from wall to edge
        double T_wall = bc.Tw();
        double T_edge = 300.0; // Assumed edge temperature
        for (std::size_t i = 0; i < n_eta; ++i) {
            double eta_norm = static_cast<double>(i) / (n_eta - 1);
            T[i] = T_wall + (T_edge - T_wall) * eta_norm;
            F[i] = eta_norm; // Simple linear profile
        }
        
        // Initialize species with edge composition
        for (std::size_t j = 0; j < n_species; ++j) {
            for (std::size_t i = 0; i < n_eta; ++i) {
                c(j, i) = (j < bc.c_e().size()) ? bc.c_e()[j] : 0.0;
                dc_deta(j, i) = 0.0; // Zero gradient for test
            }
        }
        
        // Create xi derivatives
        blast::boundary_layer::coefficients::XiDerivatives xi_der(n_eta, n_species);
        
        // Create coefficient inputs
        blast::boundary_layer::coefficients::CoefficientInputs inputs{
            .xi = 0.0,
            .F = F,
            .c = c,
            .dc_deta = dc_deta,
            .T = T
        };
        
        // Calculate coefficients
        auto coeff_result = calculator.calculate(inputs, bc, xi_der);
        if (!coeff_result) {
            std::cerr << "Failed to calculate coefficients: " << coeff_result.error().message() << "\n";
            return 1;
        }
        
        auto& coeffs = coeff_result.value();
        
        // Print some results to verify
        std::cout << "Coefficient calculation successful!\n";
        std::cout << "Number of eta points: " << n_eta << "\n";
        std::cout << "Number of species: " << n_species << "\n";
        std::cout << "Wall density: " << coeffs.wall.rho_wall << " kg/m³\n";
        std::cout << "Wall viscosity: " << coeffs.wall.mu_wall << " Pa·s\n";
        std::cout << "Wall thermal conductivity: " << coeffs.wall.k_wall << " W/(m·K)\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}



/* 
g++ -std=c++23 -Wall -Wextra -pedantic -fmax-errors=3 -I./include -I./libs/eigen -I./libs/yaml-cpp/include -I./libs/mutationpp/include src/main.cpp src/io/config_manager.cpp src/io/yaml_parser.cpp src/thermophysics/mutation_mixture.cpp src/boundary_layer/coefficients/coefficient_calculator.cpp src/boundary_layer/coefficients/xi_derivatives.cpp src/boundary_layer/conditions/boundary_interpolator.cpp src/boundary_layer/grid/coordinate_transform.cpp src/boundary_layer/grid/grid.cpp -L./libs/yaml-cpp/lib -L./libs/mutationpp/lib -lyaml-cpp -lmutation++ -o blast_test.exe


 */
