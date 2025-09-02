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
  
  // Catalysis provider configuration
  enum class CatalysisProvider { MutationPP, GASP2 };
  CatalysisProvider catalysis_provider = CatalysisProvider::MutationPP;
};

struct NumericalConfig {
  int n_eta = constants::grid::default_n_eta;
  double eta_max = constants::grid::default_eta_max;
  double convergence_tolerance = constants::tolerance::standard;
  int max_iterations = constants::iteration_limits::boundary_layer_max;

  // Global numerical guards and policies
  double divergence_threshold = 1.0e6;   // Max residual before declaring divergence
  double residual_guard = 1.0e10;        // Residual cap to allow continuation attempt

  enum class NanPolicy { ReduceStep, Fail };
  NanPolicy nan_policy = NanPolicy::ReduceStep;

  enum class ContinuationAttemptPolicy { Always, OnlyIfResidualBelowGuard, Never };
  ContinuationAttemptPolicy continuation_attempt_policy = ContinuationAttemptPolicy::OnlyIfResidualBelowGuard;

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

  double emissivity = constants::defaults::default_emissivity;                   
  double environment_temperature = constants::defaults::ambient_temperature;     
};

struct AbacusConfig {
  enum class Mode { TemperatureSweep, GammaSweepAtTw };
  bool enabled = false;
  Mode mode = Mode::TemperatureSweep;
  std::vector<double> catalyticity_values;
  double temperature_min = constants::defaults::abacus::temperature_min;
  double temperature_max = constants::defaults::abacus::temperature_max;
  int temperature_points = constants::grid::default_temperature_points;
  double radius = 1.0;     // [m] Radius for abacus
  double velocity = 0.0;   // [m/s] Velocity for abacus
  // Optional explicit temperatures override (used for gamma_sweep with multiple Tw)
  std::vector<double> temperatures_override;
};

struct ContinuationConfig {
  double wall_temperature_stable = constants::defaults::stable_wall_temperature;
  double edge_temperature_stable = constants::defaults::stable_edge_temperature;
  double pressure_stable = constants::defaults::stable_pressure;
  bool use_linear_predictor = true;  // Back-compat alias for predictor_enabled

  // Step and adaptation parameters
  double step_initial = 0.01;
  double step_min = 1.0e-4;
  double step_max = 0.5;
  double step_increase_factor = 1.3;
  double step_decrease_factor = 0.5;
  int max_steps = 10000;

  // Mode switching thresholds
  int failure_threshold = 4;   // consecutive failures before switching to equilibrium
  int success_threshold = 3;   // consecutive successes before switching back

  // Predictor parameters
  bool predictor_enabled = true;
  int predictor_max_step_reductions = 3;
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

// GASP2 specific configuration (only used when catalysis_provider = GASP2)
struct Gasp2Config {
  std::string xml_file = "config/surface_chemistry/CO2_5_finite_rate.xml";
};

// Mutation++ specific configuration (only used when catalysis_provider = MutationPP)
struct MutationConfig {
  std::string gsi_file = "";  // If empty, uses default from mixture
  bool reload_gsi_each_step = false;
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
  Gasp2Config gasp2;
  MutationConfig mutation;
  // Global verbosity flag for logging
  bool verbose = false;
};

template <typename T>
concept DefaultApplicable = requires(T& t) {
  { t.apply_defaults() } -> std::same_as<void>;
};

} // namespace blast::io
