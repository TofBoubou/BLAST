#include "blast/io/config_manager.hpp"
#include "blast/thermophysics/mixture_interface.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include "blast/boundary_layer/coefficients/xi_derivatives.hpp"
#include "blast/boundary_layer/conditions/boundary_interpolator.hpp"
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

// Helper function to create realistic air composition based on species names
std::vector<double> create_air_composition(const blast::thermophysics::MixtureInterface& mixture) {
    const auto n_species = mixture.n_species();
    std::vector<double> composition(n_species, 0.0);
    
    std::cout << "DEBUG: Auto-generating air composition for species:" << std::endl;
    
    // Find N2 and O2 indices
    int n2_index = -1, o2_index = -1;
    
    for (std::size_t i = 0; i < n_species; ++i) {
        std::string name(mixture.species_name(i));
        std::cout << "  [" << i << "] " << name << std::endl;
        
        if (name == "N2") {
            n2_index = static_cast<int>(i);
        } else if (name == "O2") {
            o2_index = static_cast<int>(i);
        }
    }
    
    if (n2_index >= 0 && o2_index >= 0) {
        // Standard air composition by mass
        composition[n2_index] = 0.767;  // ~76.7% N2
        composition[o2_index] = 0.233;  // ~23.3% O2
        std::cout << "DEBUG: Set N2[" << n2_index << "] = 0.767, O2[" << o2_index << "] = 0.233" << std::endl;
    } else {
        std::cout << "WARNING: Could not find N2 and O2 species, using uniform distribution" << std::endl;
        // Fallback: distribute equally among all species
        std::fill(composition.begin(), composition.end(), 1.0 / n_species);
    }
    
    return composition;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file.yaml>\n";
        return 1;
    }
    
    try {
        std::cout << "DEBUG: Starting program..." << std::endl;
        
        // Load configuration
        std::cout << "DEBUG: Loading configuration..." << std::endl;
        blast::io::ConfigurationManager config_manager;
        auto config_result = config_manager.load(argv[1]);
        if (!config_result) {
            std::cerr << "Failed to load config: " << config_result.error().message() << "\n";
            return 1;
        }
        auto config = config_result.value();
        std::cout << "DEBUG: Configuration loaded successfully" << std::endl;
        
        // Create mixture
        std::cout << "DEBUG: Creating mixture with name: " << config.mixture.name << std::endl;
        auto mixture_result = blast::thermophysics::create_mixture(config.mixture);
        if (!mixture_result) {
            std::cerr << "Failed to create mixture: " << mixture_result.error().message() << "\n";
            return 1;
        }
        auto& mixture = *mixture_result.value();
        std::cout << "DEBUG: Mixture created successfully, n_species = " << mixture.n_species() << std::endl;
        
        // Print species information
        std::cout << "DEBUG: Species in mixture:" << std::endl;
        for (std::size_t i = 0; i < mixture.n_species(); ++i) {
            std::cout << "  " << i << ": " << mixture.species_name(i) 
                      << " (MW: " << mixture.species_molecular_weight(i) << ")" << std::endl;
        }
        
        // Create coefficient calculator
        std::cout << "DEBUG: Creating coefficient calculator..." << std::endl;
        blast::boundary_layer::coefficients::CoefficientCalculator calculator(
            mixture, config.simulation, config.numerical
        );
        std::cout << "DEBUG: Coefficient calculator created" << std::endl;
        
        // Initialize test data
        const std::size_t n_eta = config.numerical.n_eta;
        const std::size_t n_species = mixture.n_species();
        std::cout << "DEBUG: n_eta = " << n_eta << ", n_species = " << n_species << std::endl;
        
        // Create boundary conditions for stagnation point
        std::cout << "DEBUG: Creating boundary conditions..." << std::endl;
        auto bc_result = blast::boundary_layer::conditions::create_stagnation_conditions(
            config.outer_edge, config.wall_parameters, config.simulation
        );
        if (!bc_result) {
            std::cerr << "Failed to create boundary conditions: " << bc_result.error().message() << "\n";
            return 1;
        }
        auto bc = bc_result.value();
        std::cout << "DEBUG: Boundary conditions created, species_fractions.size() = " << bc.c_e().size() << std::endl;
        
        // Print original edge conditions
        std::cout << "DEBUG: Original edge conditions:" << std::endl;
        std::cout << "  Pressure: " << bc.P_e() << " Pa" << std::endl;
        std::cout << "  Density: " << bc.rho_e() << " kg/m³" << std::endl;
        std::cout << "  Viscosity: " << bc.mu_e() << " Pa·s" << std::endl;
        std::cout << "  Original species fractions:" << std::endl;
        for (std::size_t i = 0; i < bc.c_e().size(); ++i) {
            std::cout << "    [" << i << "] " << mixture.species_name(i) << " = " << bc.c_e()[i] << std::endl;
        }
        
        // Create realistic air composition
        auto air_composition = create_air_composition(mixture);
        
        // Check if composition contains atomic species at low temperature
        bool has_atomic_species = false;
        double atomic_fraction = 0.0;
        
        for (std::size_t i = 0; i < mixture.n_species(); ++i) {
            std::string name(mixture.species_name(i));
            if (name == "N" || name == "O") {
                atomic_fraction += air_composition[i];
                if (air_composition[i] > 0.01) { // More than 1% atomic species
                    has_atomic_species = true;
                }
            }
        }
        
        if (has_atomic_species && bc.Tw() < 1000.0) {
            std::cout << "WARNING: Detected " << atomic_fraction * 100 << "% atomic species at low temperature!" << std::endl;
            std::cout << "         Using corrected air composition instead." << std::endl;
            // Use the auto-generated composition instead
        } else {
            // Use original composition if it's reasonable
            air_composition = bc.c_e();
            if (air_composition.size() != n_species) {
                air_composition.resize(n_species, 0.0);
            }
        }
        
        // Normalize composition
        double sum = std::accumulate(air_composition.begin(), air_composition.end(), 0.0);
        if (sum > 1e-12) {
            for (auto& frac : air_composition) {
                frac /= sum;
            }
        } else {
            air_composition = create_air_composition(mixture);
        }
        
        std::cout << "DEBUG: Final species composition:" << std::endl;
        for (std::size_t i = 0; i < air_composition.size(); ++i) {
            std::cout << "  [" << i << "] " << mixture.species_name(i) << " = " << air_composition[i] << std::endl;
        }
        
        // Test mixture with final composition
        std::cout << "DEBUG: Testing mixture with final composition..." << std::endl;
        auto mw_test = mixture.mixture_molecular_weight(air_composition);
        if (!mw_test) {
            std::cerr << "ERROR: Mixture molecular weight calculation failed: " 
                      << mw_test.error().message() << std::endl;
            return 1;
        }
        std::cout << "DEBUG: Mixture MW = " << mw_test.value() << " kg/kmol" << std::endl;
        
        // Initialize solution fields
        std::cout << "DEBUG: Initializing solution fields..." << std::endl;
        std::vector<double> F(n_eta, 0.0);
        std::vector<double> T(n_eta);
        blast::core::Matrix<double> c(n_species, n_eta);
        blast::core::Matrix<double> dc_deta(n_species, n_eta);
        
        // More realistic temperature profile
        double T_wall = bc.Tw();
        double T_edge = std::max(1000.0, T_wall + 700.0); // Ensure reasonable temperature difference
        std::cout << "DEBUG: T_wall = " << T_wall << ", T_edge = " << T_edge << std::endl;
        
        for (std::size_t i = 0; i < n_eta; ++i) {
            double eta = static_cast<double>(i) * config.numerical.eta_max / (n_eta - 1);
            double eta_norm = static_cast<double>(i) / (n_eta - 1);
            
            // Boundary layer temperature profile
            T[i] = T_wall + (T_edge - T_wall) * (1.0 - std::exp(-eta));
            F[i] = eta_norm * std::exp(-0.1 * eta);
        }
        std::cout << "DEBUG: Temperature range: " << T[0] << " to " << T[n_eta-1] << " K" << std::endl;
        
        // Initialize species with corrected composition
        std::cout << "DEBUG: Initializing species concentrations..." << std::endl;
        for (std::size_t j = 0; j < n_species; ++j) {
            for (std::size_t i = 0; i < n_eta; ++i) {
                c(j, i) = air_composition[j];
                dc_deta(j, i) = 0.0;
            }
        }
        
        // Create and initialize xi derivatives
        std::cout << "DEBUG: Creating xi derivatives..." << std::endl;
        blast::boundary_layer::coefficients::XiDerivatives xi_der(n_eta, n_species);
        std::vector<double> g_dummy(n_eta, 0.0);
        xi_der.update_station(0, 0.0, F, g_dummy, c);
        std::cout << "DEBUG: Xi derivatives initialized" << std::endl;
        
        // Test individual mixture calls
        std::cout << "DEBUG: Testing individual mixture calls..." << std::endl;
        auto viscosity_test = mixture.viscosity(air_composition, T[0], bc.P_e());
        if (!viscosity_test) {
            std::cerr << "ERROR: Viscosity calculation failed: " << viscosity_test.error().message() << std::endl;
            return 1;
        }
        std::cout << "DEBUG: Test viscosity = " << viscosity_test.value() << " Pa·s" << std::endl;
        
        auto cp_test = mixture.frozen_cp(air_composition, T[0], bc.P_e());
        if (!cp_test) {
            std::cerr << "ERROR: Cp calculation failed: " << cp_test.error().message() << std::endl;
            return 1;
        }
        std::cout << "DEBUG: Test Cp = " << cp_test.value() << " J/(kg·K)" << std::endl;
        
        // Create coefficient inputs
        std::cout << "DEBUG: Creating coefficient inputs..." << std::endl;
        blast::boundary_layer::coefficients::CoefficientInputs inputs{
            .xi = 0.0,
            .F = F,
            .c = c,
            .dc_deta = dc_deta,
            .T = T
        };
        
        // Calculate coefficients
        std::cout << "DEBUG: Calculating coefficients..." << std::endl;
        auto coeff_result = calculator.calculate(inputs, bc, xi_der);
        if (!coeff_result) {
            std::cerr << "Failed to calculate coefficients: " << coeff_result.error().message() << "\n";
            return 1;
        }
        
        auto& coeffs = coeff_result.value();
        std::cout << "DEBUG: Coefficients calculated successfully!" << std::endl;
        
        // Print results
        std::cout << "\n=== CALCULATION SUCCESSFUL ===" << std::endl;
        std::cout << "Number of eta points: " << n_eta << "\n";
        std::cout << "Number of species: " << n_species << "\n";
        std::cout << "Wall density: " << coeffs.wall.rho_wall << " kg/m³\n";
        std::cout << "Wall viscosity: " << coeffs.wall.mu_wall << " Pa·s\n";
        std::cout << "Wall thermal conductivity: " << coeffs.wall.k_wall << " W/(m·K)\n";
        std::cout << "Wall Prandtl number: " << coeffs.wall.Pr_wall << "\n";
        
        // Print some transport coefficient ranges
        auto [l0_min, l0_max] = std::minmax_element(coeffs.transport.l0.begin(), coeffs.transport.l0.end());
        auto [l3_min, l3_max] = std::minmax_element(coeffs.transport.l3.begin(), coeffs.transport.l3.end());
        
        std::cout << "l0 coefficient range: [" << *l0_min << ", " << *l0_max << "]\n";
        std::cout << "l3 coefficient range: [" << *l3_min << ", " << *l3_max << "]\n";
        
        // Test case: Display diffusion coefficient values
        std::cout << "\n=== DIFFUSION TEST CASE ===" << std::endl;
        std::cout << "Binary diffusion coefficients (first eta point):\n";

        for (std::size_t i = 0; i < n_species; ++i) {
            for (std::size_t j = 0; j < n_species; ++j) {
                double dij = coeffs.diffusion.Dij_bin(i, j);
                std::cout << "  D[" << mixture.species_name(i) << "][" << mixture.species_name(j) << "] = " << dij << " m²/s\n";
            }
        }

        std::cout << "\nY coordinate values:\n";
        for (std::size_t i = 0; i < std::min(n_eta, static_cast<std::size_t>(5)); ++i) {
            std::cout << "  y[" << i << "] = " << coeffs.diffusion.y[i] << "\n";
        }
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
