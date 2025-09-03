#pragma once
#include "../core/exceptions.hpp"
#include "config_types.hpp"
#include <algorithm>
#include <concepts>
#include <expected>
#include <format>
#include <source_location>
#include <yaml-cpp/yaml.h>

namespace blast::io {

class YamlParser {
private:
  YAML::Node root_;
  std::string file_path_;

  template <typename T>
  [[nodiscard]] auto extract_value(const YAML::Node& node,
                                   std::string_view key) const -> std::expected<T, core::ConfigurationError>;

  template <typename EnumType>
  [[nodiscard]] auto extract_enum(const YAML::Node& node, std::string_view key,
                                  const std::unordered_map<std::string, EnumType>& mapping) const
      -> std::expected<EnumType, core::ConfigurationError>;

  [[nodiscard]] auto
  parse_simulation_config(const YAML::Node& node, bool abacus_mode = false, bool edge_reconstruction_mode = false, bool verbose = false) const -> std::expected<SimulationConfig, core::ConfigurationError>;

  [[nodiscard]] auto
  parse_numerical_config(const YAML::Node& node) const -> std::expected<NumericalConfig, core::ConfigurationError>;

  [[nodiscard]] auto
  parse_mixture_config(const YAML::Node& node) const -> std::expected<MixtureConfig, core::ConfigurationError>;

  [[nodiscard]] auto
  parse_output_config(const YAML::Node& node) const -> std::expected<OutputConfig, core::ConfigurationError>;

  [[nodiscard]] auto parse_initial_guess_config(const YAML::Node& node) const
      -> std::expected<InitialGuessConfig, core::ConfigurationError>;

  [[nodiscard]] auto
  parse_outer_edge_config(const YAML::Node& node, bool edge_reconstruction_mode = false, bool abacus_mode = false, bool verbose = false) const -> std::expected<OuterEdgeConfig, core::ConfigurationError>;

  [[nodiscard]] auto parse_wall_parameters_config(const YAML::Node& node) const
      -> std::expected<WallParametersConfig, core::ConfigurationError>;

  [[nodiscard]] auto
  parse_abacus_config(const YAML::Node& node) const -> std::expected<AbacusConfig, core::ConfigurationError>;

  [[nodiscard]] auto parse_continuation_config(const YAML::Node& node) const
      -> std::expected<ContinuationConfig, core::ConfigurationError>;

  [[nodiscard]] auto parse_edge_reconstruction_config(const YAML::Node& node) const
      -> std::expected<EdgeReconstructionConfig, core::ConfigurationError>;

  [[nodiscard]] auto parse_catalysis_reconstruction_config(const YAML::Node& node) const
      -> std::expected<CatalysisReconstructionConfig, core::ConfigurationError>;

  [[nodiscard]] auto parse_gasp2_config(const YAML::Node& node) const
      -> std::expected<Gasp2Config, core::ConfigurationError>;

  [[nodiscard]] auto parse_mutation_config(const YAML::Node& node) const
      -> std::expected<MutationConfig, core::ConfigurationError>;

public:
  explicit YamlParser(std::string file_path) : file_path_(std::move(file_path)) {} // constructor so no ->

  [[nodiscard]] auto load() -> std::expected<void, core::FileError>;

  [[nodiscard]] auto parse() const -> std::expected<Configuration, core::ConfigurationError>;
};

// Implementation of template methods
template <typename T>
auto YamlParser::extract_value(const YAML::Node& node,
                               std::string_view key) const -> std::expected<T, core::ConfigurationError> {
  try {
    if (!node[std::string(key)]) {
      return std::unexpected(core::ConfigurationError(std::format("Required field '{}' is missing", key)));
    }

    if constexpr (std::same_as<T, std::vector<double>>) {
      auto sequence = node[std::string(key)];
      std::vector<double> result;
      result.reserve(sequence.size());

      for (const auto& item : sequence) {
        result.push_back(item.as<double>());
      }
      return result;
    } else {
      return node[std::string(key)].as<T>();
    }
  } catch (const YAML::Exception& e) {
    return std::unexpected(core::ConfigurationError(std::format("Failed to parse field '{}': {}", key, e.what())));
  }
}

template <typename EnumType>
auto YamlParser::extract_enum(const YAML::Node& node, std::string_view key,
                              const std::unordered_map<std::string, EnumType>& mapping) const
    -> std::expected<EnumType, core::ConfigurationError> {
  auto str_result = extract_value<std::string>(node, key);
  if (!str_result) {
    return std::unexpected(str_result.error());
  }

  auto str_value = str_result.value();
  std::ranges::transform(str_value, str_value.begin(), ::tolower);

  auto it = mapping.find(str_value);
  if (it == mapping.end()) {
    std::string valid_options;
    for (const auto& [option, _] : mapping) {
      valid_options += option + ", ";
    }
    valid_options = valid_options.substr(0, valid_options.length() - 2);

    return std::unexpected(core::ConfigurationError(
        std::format("Invalid value '{}' for field '{}'. Valid options: {}", str_value, key, valid_options)));
  }

  return it->second;
}

// Enum mappings
namespace enum_mappings {

inline const std::unordered_map<std::string, SimulationConfig::BodyType> body_types = {
    {"axisymmetric", SimulationConfig::BodyType::Axisymmetric},
    {"cone", SimulationConfig::BodyType::Cone},
    {"twod", SimulationConfig::BodyType::TwoD},
    {"2d", SimulationConfig::BodyType::TwoD},
    {"flatplate", SimulationConfig::BodyType::FlatPlate},
    {"flat_plate", SimulationConfig::BodyType::FlatPlate}};

inline const std::unordered_map<std::string, SimulationConfig::DiffusionType> diffusion_types = {
    {"ramshaw", SimulationConfig::DiffusionType::Ramshaw},
    {"stefanmaxwell", SimulationConfig::DiffusionType::StefanMaxwell},
    {"stefan_maxwell", SimulationConfig::DiffusionType::StefanMaxwell},
    {"mpp", SimulationConfig::DiffusionType::MPP}};

inline const std::unordered_map<std::string, MixtureConfig::ViscosityAlgorithm> viscosity_algorithms = {
    {"chapman_enskog_cg", MixtureConfig::ViscosityAlgorithm::chapmanEnskog_CG},
    {"gupta_yos", MixtureConfig::ViscosityAlgorithm::GuptaYos},
    {"chapman_enskog_ldlt", MixtureConfig::ViscosityAlgorithm::chapmanEnskog_LDLT},
    {"wilke", MixtureConfig::ViscosityAlgorithm::Wilke}};

inline const std::unordered_map<std::string, MixtureConfig::ThermalConductivityAlgorithm> thermal_conductivity_algorithms = {
    {"chapman_enskog_cg", MixtureConfig::ThermalConductivityAlgorithm::chapmanEnskog_CG},
    {"chapman_enskog_ldlt", MixtureConfig::ThermalConductivityAlgorithm::chapmanEnskog_LDLT},
    {"wilke", MixtureConfig::ThermalConductivityAlgorithm::Wilke}};

inline const std::unordered_map<std::string, SimulationConfig::CatalysisProvider> catalysis_providers = {
    {"mutation++", SimulationConfig::CatalysisProvider::MutationPP},
    {"mutationpp", SimulationConfig::CatalysisProvider::MutationPP},
    {"gasp2", SimulationConfig::CatalysisProvider::GASP2}};

inline const std::unordered_map<std::string, MixtureConfig::Database> databases = {
    {"rrho", MixtureConfig::Database::RRHO},
    {"nasa7", MixtureConfig::Database::NASA7},
    {"nasa-7", MixtureConfig::Database::NASA7},
    {"nasa9", MixtureConfig::Database::NASA9},
    {"nasa-9", MixtureConfig::Database::NASA9}};


inline const std::unordered_map<std::string, SimulationConfig::ChemicalMode> chemical_modes = {
    {"equilibrium", SimulationConfig::ChemicalMode::Equilibrium},
    {"frozen", SimulationConfig::ChemicalMode::Frozen},
    {"non_equilibrium", SimulationConfig::ChemicalMode::NonEquilibrium},
    {"nonequilibrium", SimulationConfig::ChemicalMode::NonEquilibrium}};

inline const std::unordered_map<std::string, SimulationConfig::WallMode> wall_modes = {
    {"adiabatic", SimulationConfig::WallMode::Adiabatic},
    {"imposed_temperature", SimulationConfig::WallMode::ImposedTemperature},
    {"imposed", SimulationConfig::WallMode::ImposedTemperature},
    {"radiative", SimulationConfig::WallMode::Radiative},
    {"radiative_wall", SimulationConfig::WallMode::Radiative}};

inline const std::unordered_map<std::string, AbacusConfig::Mode> abacus_modes = {
    {"temperature_sweep", AbacusConfig::Mode::TemperatureSweep},
    {"gamma_sweep_at_tw", AbacusConfig::Mode::GammaSweepAtTw},
    {"gamma_sweep", AbacusConfig::Mode::GammaSweepAtTw},
    {"tw_fixed", AbacusConfig::Mode::GammaSweepAtTw}};

} // namespace enum_mappings

// Additional enum mappings for numerical policies
inline const std::unordered_map<std::string, NumericalConfig::NanPolicy> nan_policies = {
    {"reduce_step", NumericalConfig::NanPolicy::ReduceStep},
    {"fail", NumericalConfig::NanPolicy::Fail}};

inline const std::unordered_map<std::string, NumericalConfig::ContinuationAttemptPolicy> continuation_policies = {
    {"always", NumericalConfig::ContinuationAttemptPolicy::Always},
    {"only_if_residual_below_threshold", NumericalConfig::ContinuationAttemptPolicy::OnlyIfResidualBelowGuard},
    {"only_if_residual_below_guard", NumericalConfig::ContinuationAttemptPolicy::OnlyIfResidualBelowGuard},
    {"never", NumericalConfig::ContinuationAttemptPolicy::Never}};

} // namespace blast::io
