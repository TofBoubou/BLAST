#pragma once
#include "../core/exceptions.hpp"
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

  BodyType body_type = BodyType::Axisymmetric;
  bool only_stagnation_point = true;
  DiffusionType diffusion_type = DiffusionType::StefanMaxwell;
  bool consider_thermal_diffusion = false;
  bool consider_dufour_effect = false;
  ChemicalMode chemical_mode = ChemicalMode::NonEquilibrium;
  bool catalytic_wall = false; // Enable surface catalysis
};

struct NumericalConfig {
  int n_eta = 101;
  double eta_max = 8.0;
  double convergence_tolerance = 1e-6;
  int max_iterations = 1000;

  struct Solvers {
    double h2t_tolerance = 1.0e-12;
    int h2t_max_iterations = 50000;
    double stefan_tolerance = 1.0e-12;
    int stefan_max_iterations = 50000;
  } solvers{};
};

struct MixtureConfig {
  std::string name = "air7";
  enum class Database { RRHO, NASA7, NASA9 };
  Database thermodynamic_database = Database::NASA9;

  enum class ViscosityAlgorithm { chapmanEnskog_CG, GuptaYos, chapmanEnskog_LDLT, Wilke };
  ViscosityAlgorithm viscosity_algorithm = ViscosityAlgorithm::chapmanEnskog_CG;

  // State model controlling energy mode treatment
  enum class StateModel { ChemNonEq1T, ChemNonEqTTv };
  StateModel state_model = StateModel::ChemNonEq1T;
};

struct OutputConfig {
  std::vector<double> x_stations;
  std::string output_directory = "BLAST_outputs";
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
    double temperature;
    double pressure;

    bool boundary_override = false;
    std::optional<std::vector<double>> mass_fraction_condition;
  };

  std::vector<EdgePoint> edge_points;
  double velocity_gradient_stagnation;
  double freestream_density;
  double freestream_velocity;
};

struct WallParametersConfig {
  std::vector<double> wall_temperatures;
  bool catalytic_wall = false;
};

struct AbaqueConfig {
  bool enabled = false;
  std::vector<double> catalyticity_values;
  double temperature_min = 300.0;
  double temperature_max = 2000.0;
  int temperature_points = 100;
};

struct ContinuationConfig {
  double wall_temperature_stable = 3100.0;
  double edge_temperature_stable = 3100.0;
  double pressure_stable = 7000.0;
};

struct Configuration {
  SimulationConfig simulation;
  NumericalConfig numerical;
  MixtureConfig mixture;
  OutputConfig output;
  InitialGuessConfig initial_guess;
  OuterEdgeConfig outer_edge;
  WallParametersConfig wall_parameters;
  AbaqueConfig abaque;
  ContinuationConfig continuation;
};

template <typename T>
concept DefaultApplicable = requires(T& t) {
  { t.apply_defaults() } -> std::same_as<void>;
};

} // namespace blast::io