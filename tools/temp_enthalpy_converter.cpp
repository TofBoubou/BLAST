#include <cassert>
#include <mutation++/mutation++.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <cmath>
#include <filesystem>
#include <cstdlib>

void print_usage() {
    std::cout << "Usage:\n";
    std::cout << "  temp_enthalpy_converter <mixture_name> <pressure_Pa> <mode> <value>\n";
    std::cout << "\nModes:\n";
    std::cout << "  T2H  - Convert temperature to enthalpy\n";
    std::cout << "  H2T  - Convert enthalpy to temperature\n";
    std::cout << "\nExample:\n";
    std::cout << "  temp_enthalpy_converter air_5 101325 T2H 300\n";
    std::cout << "  temp_enthalpy_converter air_5 101325 H2T 300000\n";
}

double enthalpy_to_temperature(Mutation::Mixture& mix, std::vector<double>& mass_fractions, 
                              double target_h, double pressure) {
    const double tol = 1e-6;
    const int max_iter = 100;
    
    double T_min = 200.0;
    double T_max = 15000.0;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        double T_mid = (T_min + T_max) / 2.0;
        
        double vars[2] = {pressure, T_mid};
        mix.setState(mass_fractions.data(), vars, 2);
        double h_mid = mix.mixtureHMass();
        
        if (std::abs(h_mid - target_h) < tol) {
            return T_mid;
        }
        
        if (h_mid < target_h) {
            T_min = T_mid;
        } else {
            T_max = T_mid;
        }
    }
    
    return (T_min + T_max) / 2.0;
}

void setup_mutation_data_path() {
    // Check if MPP_DATA_DIRECTORY is already set
    if (std::getenv("MPP_DATA_DIRECTORY") != nullptr) {
        return;
    }
    
    // Try to find the data directory relative to the executable
    namespace fs = std::filesystem;
    
    // Get the path to this executable
    fs::path exe_path = fs::current_path();
    
    // Try common relative paths
    std::vector<fs::path> possible_paths = {
        exe_path / ".." / "libs" / "mutationpp" / "data",
        exe_path / "libs" / "mutationpp" / "data",
        fs::path("../libs/mutationpp/data"),
        fs::path("libs/mutationpp/data"),
        fs::path("../../libs/mutationpp/data")
    };
    
    for (const auto& path : possible_paths) {
        if (fs::exists(path / "mixtures")) {
            setenv("MPP_DATA_DIRECTORY", path.string().c_str(), 1);
            std::cout << "Found Mutation++ data at: " << path << "\n";
            return;
        }
    }
    
    std::cerr << "Warning: Could not find Mutation++ data directory.\n";
    std::cerr << "Please set MPP_DATA_DIRECTORY environment variable.\n";
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        print_usage();
        return 1;
    }
    
    // Auto-setup Mutation++ data path
    setup_mutation_data_path();
    
    std::string mixture_name = argv[1];
    double pressure = std::stod(argv[2]);
    std::string mode = argv[3];
    double value = std::stod(argv[4]);
    
    bool use_custom_composition = false;
    std::vector<double> custom_mass_fractions;
    if (argc > 5) {
        use_custom_composition = true;
        for (int i = 5; i < argc; ++i) {
            custom_mass_fractions.push_back(std::stod(argv[i]));
        }
    }
    
    try {
        Mutation::MixtureOptions opts(mixture_name);
        opts.setStateModel("ChemNonEq1T");
        Mutation::Mixture mix(opts);
        
        std::vector<double> mass_fractions(mix.nSpecies());
        if (use_custom_composition && custom_mass_fractions.size() == mix.nSpecies()) {
            mass_fractions = custom_mass_fractions;
            double sum = 0.0;
            for (double& f : mass_fractions) sum += f;
            for (double& f : mass_fractions) f /= sum;
        } else {
            const double* default_comp = mix.getDefaultComposition();
            for (int i = 0; i < mix.nSpecies(); ++i) {
                mass_fractions[i] = default_comp[i];
            }
        }
        
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "\n=== Temperature-Enthalpy Converter ===\n";
        std::cout << "Mixture: " << mixture_name << "\n";
        std::cout << "Pressure: " << pressure << " Pa\n";
        std::cout << "Species (" << mix.nSpecies() << "): ";
        for (int i = 0; i < mix.nSpecies(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << mix.speciesName(i);
        }
        std::cout << "\n";
        
        std::cout << "Composition (mass fractions):\n";
        for (int i = 0; i < mix.nSpecies(); ++i) {
            std::cout << "  " << std::setw(10) << mix.speciesName(i) 
                     << ": " << std::setw(12) << mass_fractions[i] << "\n";
        }
        std::cout << "\n";
        
        if (mode == "T2H") {
            double temperature = value;
            double vars[2] = {pressure, temperature};
            mix.setState(mass_fractions.data(), vars, 2);
            
            double enthalpy = mix.mixtureHMass();
            
            std::cout << "Input Temperature: " << temperature << " K\n";
            std::cout << "Computed Enthalpy: " << enthalpy << " J/kg\n";
            
            std::vector<double> species_h(mix.nSpecies());
            mix.speciesHOverRT(temperature, species_h.data());
            
            const double R_universal = 8314.46261815324;
            std::cout << "\nSpecies Enthalpies (J/kg):\n";
            for (int i = 0; i < mix.nSpecies(); ++i) {
                double h_i = species_h[i] * R_universal * temperature / mix.speciesMw(i);
                std::cout << "  " << std::setw(10) << mix.speciesName(i) 
                         << ": " << std::setw(15) << h_i << "\n";
            }
            
        } else if (mode == "H2T") {
            double target_enthalpy = value;
            
            double temperature = enthalpy_to_temperature(mix, mass_fractions, target_enthalpy, pressure);
            
            double vars[2] = {pressure, temperature};
            mix.setState(mass_fractions.data(), vars, 2);
            double computed_h = mix.mixtureHMass();
            
            std::cout << "Input Enthalpy: " << target_enthalpy << " J/kg\n";
            std::cout << "Computed Temperature: " << temperature << " K\n";
            std::cout << "Verification - Computed Enthalpy: " << computed_h << " J/kg\n";
            std::cout << "Error: " << std::abs(computed_h - target_enthalpy) << " J/kg\n";
            
        } else {
            std::cerr << "Error: Unknown mode '" << mode << "'. Use T2H or H2T.\n";
            return 1;
        }
        
        double density = mix.density();
        double cp = mix.mixtureFrozenCpMass();
        double cv = mix.mixtureFrozenCvMass();
        double gamma = cp / cv;
        
        std::cout << "\nAdditional Properties:\n";
        std::cout << "  Density: " << density << " kg/m³\n";
        std::cout << "  Cp (frozen): " << cp << " J/(kg·K)\n";
        std::cout << "  Cv (frozen): " << cv << " J/(kg·K)\n";
        std::cout << "  Gamma: " << gamma << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

/*
USAGE INSTRUCTIONS:

1. COMPILATION:
   make

2. BASIC USAGE:
   ./temp_enthalpy_converter <mixture> <pressure_Pa> <mode> <value> [mass_fractions]

   Modes:
   - T2H: Convert temperature (K) to enthalpy (J/kg)
   - H2T: Convert enthalpy (J/kg) to temperature (K)

3. EXAMPLES:

   # Temperature to enthalpy conversion at 300K, 1 atm
   ./temp_enthalpy_converter air_5 101325 T2H 300

   # Enthalpy to temperature conversion
   ./temp_enthalpy_converter air_5 101325 H2T 300000

   # With custom composition (mass fractions for each species)
   # For air_5: N, O, NO, N2, O2
   # Standard air: 76.7% N2, 23.3% O2
   ./temp_enthalpy_converter air_5 101325 T2H 300 0 0 0 0.767 0.233

   # High temperature example
   ./temp_enthalpy_converter air_5 101325 T2H 2000

4. NOTES:
   - Pressure is in Pascals (101325 Pa = 1 atm)
   - Temperature is in Kelvin
   - Enthalpy is in J/kg
   - Mass fractions must sum to 1 (will be normalized if not)
   - The program automatically finds the Mutation++ data directory
*/