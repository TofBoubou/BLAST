#pragma once
#include <string>
#include <vector>
#include <optional>
#include <variant>
#include <unordered_map>
#include <concepts>
#include <expected>
#include <source_location>
#include "../core/exceptions.hpp"

namespace blast::io {


using ConfigValue = std::variant<bool, int, double, std::string, std::vector<double>, std::vector<int>, std::vector<std::string>>;


struct SimulationConfig {
    enum class BodyType { Axisymmetric, Cone, TwoD, FlatPlate };
    enum class DiffusionType { Ramshaw, StefanMaxwell, MPP };
    
    BodyType body_type = BodyType::Axisymmetric; // create a type that is BodyType that is an enum_class, we say to him that it's a constant and the constant is BodyType::Axisymmetric
    bool only_stagnation_point = true;
    bool finite_thickness = false;
    DiffusionType diffusion_type = DiffusionType::StefanMaxwell;
    bool consider_thermal_diffusion = false;
    bool chemical_non_equilibrium = true;
    bool catalysis = false;

};

struct NumericalConfig {
    int n_eta = 101;
    double eta_max = 8.0;
    double convergence_tolerance = 1e-6;
    int max_iterations = 1000;
    double under_relaxation = 0.5;
    
    struct StepControl {
        int lower_bound = 10;
        int upper_bound = 50;
    } step_control;
    
    struct Solvers {
        double h2t_tolerance = 1e-8;
        int h2t_max_iterations = 100;
        double stefan_tolerance = 1e-6;
        int stefan_max_iterations = 50;
    } solvers;
    
};

struct MixtureConfig {
    std::string name = "air7";
    enum class Database { RRHO, NASA7, NASA9 };
    Database thermodynamic_database = Database::NASA9;
    
    enum class ViscosityAlgorithm { ChapmanEnskog_CG, GuptaYos, ChapmanEnskog_LDLT, Wilke };
    ViscosityAlgorithm viscosity_algorithm = ViscosityAlgorithm::ChapmanEnskog_CG;
    
    double reference_temperature = 0.0;
    
};

struct OutputConfig {
    std::vector<double> x_stations;
    std::string output_directory = "BLAST_outputs";
    bool generate_lookup_table = false;
    
    struct LookupTable {
        double temperature_min = 300.0;
        double temperature_max = 3000.0;
        double temperature_step = 100.0;
        double gamma_min = 0.0;
        double gamma_max = 1.0;
        double gamma_step = 0.1;
    } lookup_table;
    
};


struct InitialGuessConfig {
    bool use_initial_guess = false;
    std::optional<std::vector<double>> temperature_profile;
    std::optional<std::vector<double>> enthalpy_profile;
    std::optional<std::vector<std::vector<double>>> species_profiles;

};


struct OuterEdgeConfig {
    struct EdgePoint {
        double x;
        double radius;
        double velocity;
        double enthalpy;
        double pressure;
        double density;
        double viscosity;
        std::vector<double> species_fractions;
    };
    
    std::vector<EdgePoint> edge_points;
    double velocity_gradient_stagnation;
    double freestream_density;
    double freestream_velocity;
    
    // Pour finite thickness
    std::optional<double> edge_normal_velocity;
    std::optional<double> d2_ue_dxdy;
    std::optional<double> boundary_layer_thickness;

};


struct WallParametersConfig {
    std::vector<double> wall_temperatures;
};


struct Configuration {
    SimulationConfig simulation;
    NumericalConfig numerical;
    MixtureConfig mixture;
    OutputConfig output;
    InitialGuessConfig initial_guess;
    OuterEdgeConfig outer_edge;
    WallParametersConfig wall_parameters;

};


template<typename T>
concept DefaultApplicable = requires(T& t) {
    { t.apply_defaults() } -> std::same_as<void>;
};

}