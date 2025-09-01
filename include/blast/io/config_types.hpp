#pragma once
#include "../core/exceptions.hpp"
#include "../core/constants.hpp"
#include <concepts>
#include <expected>
#include <optional>
#include <source_location>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace blast::io {

using ConfigValue =
    std::variant<bool, int, double, std::string, std::vector<double>, std::vector<int>, std::vector<std::string>>;

struct SimulationConfig {
  enum class BodyType { Axisymmetric, Cone, TwoD, FlatPlate };
  enum class DiffusionType { Ramshaw, StefanMaxwell, MPP };
  enum class ChemicalMode { Equilibrium, Frozen, NonEquilibrium };
  enum class WallMode { Adiabatic, ImposedTemperature, Radiative };
  BodyType body_type = BodyType::Axisymmetric;
  bool only_stagnation_point = true;
  bool finite_thickness = false;
  DiffusionType diffusion_type = DiffusionType::StefanMaxwell;
  bool consider_thermal_diffusion = false;
  bool consider_dufour_effect = false;
  ChemicalMode chemical_mode = ChemicalMode::NonEquilibrium;
  bool catalytic_wall = false;
  WallMode wall_mode = WallMode::ImposedTemperature;  // Default to imposed temperature
};

struct NumericalConfig {
  int n_eta = constants::grid::default_n_eta;
  double eta_max = constants::grid::default_eta_max;
  double convergence_tolerance = constants::tolerance::standard;
  int max_iterations = constants::iteration_limits::boundary_layer_max;

  struct Solvers {
    double h2t_tolerance = constants::tolerance::high_precision;
    int h2t_max_iterations = constants::iteration_limits::enthalpy_temperature_max;
    double stefan_tolerance = constants::tolerance::high_precision;
    int stefan_max_iterations = constants::iteration_limits::stefan_max;
  } solvers{};
};

struct MixtureConfig {
  std::string name = "air7";
  enum class Database { RRHO, NASA7, NASA9 };
  Database thermodynamic_database = Database::NASA9;

  enum class ViscosityAlgorithm { chapmanEnskog_CG, GuptaYos, chapmanEnskog_LDLT, Wilke };
  ViscosityAlgorithm viscosity_algorithm = ViscosityAlgorithm::chapmanEnskog_CG;

  enum class ThermalConductivityAlgorithm { chapmanEnskog_CG, chapmanEnskog_LDLT, Wilke };
  ThermalConductivityAlgorithm thermal_conductivity_algorithm = ThermalConductivityAlgorithm::chapmanEnskog_CG;

  // State model is now fixed to ChemNonEq1T only
};

struct OutputConfig {
  // x_stations are now automatically derived from edge_points.x
  // No need to specify them separately
  std::string output_directory = "BLAST_outputs";
};

struct InitialGuessConfig {
  bool use_initial_guess = false;
  std::optional<std::vector<double>> temperature_profile;
  std::optional<std::vector<double>> enthalpy_profile;
  std::optional<std::vector<std::vector<double>>> species_profiles;
};

struct OuterEdgeConfig {
  struct BoundaryOverride {
    bool enabled = false;
    std::optional<std::vector<double>> mass_fraction_condition;
  };

  struct EdgePoint {
    double x;
    double radius;
    double velocity;
    double enthalpy;
    double temperature;
    double pressure;

    BoundaryOverride boundary_override;

    // Legacy compatibility properties
    [[nodiscard]] bool boundary_override_enabled() const { return boundary_override.enabled; }
    [[nodiscard]] const std::optional<std::vector<double>>& mass_fraction_condition() const {
      return boundary_override.mass_fraction_condition;
    }
  };

  struct FiniteThicknessParams {
    double v_edge = constants::defaults::default_emissivity;           // Edge velocity normal to surface
    double d2_ue_dxdy = constants::defaults::default_emissivity;       // Second derivative of edge velocity
    double delta_bl = constants::defaults::default_emissivity;         // Boundary layer thickness
  };

  std::vector<EdgePoint> edge_points;
  double velocity_gradient_stagnation;
  double freestream_density;
  double freestream_velocity;
  FiniteThicknessParams finite_thickness_params;
};

struct WallParametersConfig {
  std::vector<double> wall_temperatures;
  bool catalytic_wall = false;

  double emissivity = constants::defaults::default_emissivity;                   
  double environment_temperature = constants::defaults::ambient_temperature;     
};

struct AbacusConfig {
  bool enabled = false;
  std::vector<double> catalyticity_values;
  double temperature_min = constants::defaults::abacus::temperature_min;
  double temperature_max = constants::defaults::abacus::temperature_max;
  int temperature_points = constants::grid::default_temperature_points;
  double radius = 1.0;     // [m] Radius for abacus
  double velocity = 0.0;   // [m/s] Velocity for abacus
};

struct ContinuationConfig {
  double wall_temperature_stable = constants::defaults::stable_wall_temperature;
  double edge_temperature_stable = constants::defaults::stable_edge_temperature;
  double pressure_stable = constants::defaults::stable_pressure;
  bool use_linear_predictor = true;  // Enable linear predictor for better convergence
};

struct EdgeReconstructionConfig {
  bool enabled = false;
  
  // Target values
  double target_heat_flux = 0.0;  // [W/m²] Target total heat flux at wall
  
  // Boundary conditions
  struct BoundaryConditions {
    double pressure = 0.0;           // [Pa] Pressure at edge (fixed)
    double wall_temperature = 0.0;   // [K] Wall temperature
    double catalyticity = 0.0;       // Catalytic efficiency (0-1)
    double radius = 1.0;              // [m] Radius at stagnation point
  } boundary_conditions;
  
  // Flow parameters
  struct FlowParameters {
    double velocity_gradient_stagnation = 0.0;  // [1/s]
    double freestream_density = 0.0;            // [kg/m³]
    double freestream_velocity = 0.0;           // [m/s]
  } flow_parameters;
  
  // Finite thickness parameters (required for reconstruction)
  struct FiniteThicknessParams {
    double v_edge = 0.0;              // [m/s] Edge normal velocity
    double d2_ue_dxdy = 0.0;          // [1/s²] Second derivative
    double delta_bl = 0.0;            // [m] Boundary layer thickness
  } finite_thickness_params;
  
  // Solver configuration
  struct SolverSettings {
    double initial_temperature_guess = 4000.0;  // [K] Initial guess
    double temperature_min = 2000.0;            // [K] Lower bound
    double temperature_max = 8000.0;            // [K] Upper bound
    double tolerance = 1.0e-6;                  // Convergence tolerance
    int max_iterations = 50;                    // Maximum iterations
  } solver;
};

struct Configuration {
  SimulationConfig simulation;
  NumericalConfig numerical;
  MixtureConfig mixture;
  OutputConfig output;
  InitialGuessConfig initial_guess;
  OuterEdgeConfig outer_edge;
  WallParametersConfig wall_parameters;
  AbacusConfig abacus;
  ContinuationConfig continuation;
  EdgeReconstructionConfig edge_reconstruction;
};

template <typename T>
concept DefaultApplicable = requires(T& t) {
  { t.apply_defaults() } -> std::same_as<void>;
};

} // namespace blast::io