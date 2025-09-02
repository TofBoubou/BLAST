#include "blast/io/yaml_parser.hpp"
#include "blast/core/exceptions.hpp"
#include "blast/core/constants.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>

namespace blast::io {

auto YamlParser::load() -> std::expected<void, core::FileError> {
  try {
    root_ = YAML::LoadFile(file_path_);
    return {};
  } catch (const YAML::BadFile& e) {
    return std::unexpected(core::FileError{"Failed to open YAML file", file_path_});
  } catch (const YAML::ParserException& e) {
    return std::unexpected(core::FileError{std::format("YAML parsing error: {}", e.what()), file_path_});
  } catch (const std::exception& e) { // YAML::BadFile inherit from std::exception so no problem
    return std::unexpected(core::FileError{std::format("Unexpected error during YAML load: {}", e.what()), file_path_});
  }
}

auto YamlParser::parse() const -> std::expected<Configuration, core::ConfigurationError> {
  try {
    if (!root_ || root_.IsNull()) { // !root_ is override in YAML so we check if it's empty
      return std::unexpected(core::ConfigurationError("No YAML content loaded. Call load() first."));
    }

    Configuration config;

    // Check which modes are enabled to determine which config sections to use
    bool base_enabled = false;
    bool edge_reconstruction_enabled = false;
    bool abacus_enabled = false;
    
    // Debug markers to trace parsing (disabled)

    if (root_["base"] && root_["base"]["enabled"]) {
      base_enabled = root_["base"]["enabled"].as<bool>();
    }
    
    if (root_["edge_reconstruction"] && root_["edge_reconstruction"]["enabled"]) {
      edge_reconstruction_enabled = root_["edge_reconstruction"]["enabled"].as<bool>();
    }
    
    if (root_["abacus"] && root_["abacus"]["enabled"]) {
      abacus_enabled = root_["abacus"]["enabled"].as<bool>();
    }
    
    // Check for conflicting modes - only one can be enabled at a time
    int enabled_count = (base_enabled ? 1 : 0) + (edge_reconstruction_enabled ? 1 : 0) + (abacus_enabled ? 1 : 0);
    
    if (enabled_count == 0) {
      return std::unexpected(core::ConfigurationError(
          "No simulation mode enabled. Please enable exactly one mode: base.enabled, edge_reconstruction.enabled, or abacus.enabled."));
    }
    
    if (enabled_count > 1) {
      return std::unexpected(core::ConfigurationError(
          "Multiple simulation modes enabled. Please enable only one mode at a time: base, edge_reconstruction, or abacus."));
    }

    // Determine which configuration sections to use based on active mode (strict mode - no fallbacks)
    YAML::Node sim_node, num_node, mix_node, out_node, edge_node, wall_node;
    
    if (base_enabled) {
      // Base mode - use base sections
      if (!root_["base"]["simulation"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'simulation' section in base mode."));
      }
      if (!root_["base"]["numerical"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'numerical' section in base mode."));
      }
      if (!root_["base"]["mixture"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'mixture' section in base mode."));
      }
      if (!root_["base"]["output"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'output' section in base mode."));
      }
      if (!root_["base"]["outer_edge"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'outer_edge' section in base mode."));
      }
      if (!root_["base"]["wall_parameters"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'wall_parameters' section in base mode."));
      }
      
      sim_node = root_["base"]["simulation"];
      num_node = root_["base"]["numerical"];
      mix_node = root_["base"]["mixture"];
      out_node = root_["base"]["output"];
      edge_node = root_["base"]["outer_edge"];
      wall_node = root_["base"]["wall_parameters"];
      
    } else if (edge_reconstruction_enabled) {
      // Edge reconstruction mode - only basic sections, no outer_edge (uses specialized config)
      if (!root_["edge_reconstruction"]["simulation"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'simulation' section in edge_reconstruction mode."));
      }
      if (!root_["edge_reconstruction"]["numerical"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'numerical' section in edge_reconstruction mode."));
      }
      if (!root_["edge_reconstruction"]["mixture"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'mixture' section in edge_reconstruction mode."));
      }
      if (!root_["edge_reconstruction"]["output"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'output' section in edge_reconstruction mode."));
      }
      // edge_reconstruction uses boundary_conditions instead of wall_parameters
      
      sim_node = root_["edge_reconstruction"]["simulation"];
      num_node = root_["edge_reconstruction"]["numerical"];
      mix_node = root_["edge_reconstruction"]["mixture"];
      out_node = root_["edge_reconstruction"]["output"];
      wall_node = YAML::Node(); // Empty - edge_reconstruction uses boundary_conditions
      edge_node = YAML::Node(); // Empty - edge reconstruction uses specialized config
      
    } else if (abacus_enabled) {
      // Abacus mode - only basic sections, no outer_edge (uses specialized config)
      if (!root_["abacus"]["simulation"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'simulation' section in abacus mode."));
      }
      if (!root_["abacus"]["numerical"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'numerical' section in abacus mode."));
      }
      if (!root_["abacus"]["mixture"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'mixture' section in abacus mode."));
      }
      if (!root_["abacus"]["output"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'output' section in abacus mode."));
      }
      if (!root_["abacus"]["boundary_conditions"]) {
        return std::unexpected(core::ConfigurationError(
            "Missing required 'boundary_conditions' section in abacus mode."));
      }
      
      sim_node = root_["abacus"]["simulation"];
      num_node = root_["abacus"]["numerical"];
      mix_node = root_["abacus"]["mixture"];
      out_node = root_["abacus"]["output"];
      wall_node = YAML::Node(); // Empty - abacus uses boundary_conditions instead
      edge_node = YAML::Node(); // Empty - abacus uses specialized config
    }

    // Parse configuration from the determined nodes
    
    auto sim_result = parse_simulation_config(sim_node, abacus_enabled, edge_reconstruction_enabled);
    if (!sim_result) {
      return std::unexpected(sim_result.error());
    }
    config.simulation = std::move(sim_result.value());

    
    auto num_result = parse_numerical_config(num_node);
    if (!num_result) {
      return std::unexpected(num_result.error());
    }
    config.numerical = std::move(num_result.value());

    
    auto mix_result = parse_mixture_config(mix_node);
    if (!mix_result) {
      return std::unexpected(mix_result.error());
    }
    config.mixture = std::move(mix_result.value());

    
    auto out_result = parse_output_config(out_node);
    if (!out_result) {
      return std::unexpected(out_result.error());
    }
    config.output = std::move(out_result.value());

    // Parse outer_edge only for base mode (edge_reconstruction and abacus use specialized configs)
    if (base_enabled) {
      
      auto edge_result = parse_outer_edge_config(edge_node, false, false);
      if (!edge_result) {
        return std::unexpected(edge_result.error());
      }
      config.outer_edge = std::move(edge_result.value());
    } else {
      // For specialized modes, outer_edge will be populated from their specific configs
      config.outer_edge = OuterEdgeConfig{};
    }

    // Parse wall_parameters for modes that need it (only base mode)
    if (base_enabled) {
      
      auto wall_result = parse_wall_parameters_config(wall_node);
      if (!wall_result) {
        return std::unexpected(wall_result.error());
      }
      config.wall_parameters = std::move(wall_result.value());
    } else if (edge_reconstruction_enabled) {
      // Default wall parameters for edge_reconstruction (uses boundary_conditions)
      config.wall_parameters = WallParametersConfig{};
      // Set temperature from boundary_conditions for edge_reconstruction
      if (edge_reconstruction_enabled && root_["edge_reconstruction"]["boundary_conditions"]["wall_temperature"]) {
        config.wall_parameters.wall_temperatures = {root_["edge_reconstruction"]["boundary_conditions"]["wall_temperature"].as<double>()};
      }
    } else if (abacus_enabled) {
      // Default wall parameters for abacus (uses boundary_conditions)
      config.wall_parameters = WallParametersConfig{};
    }
    
    // Force emissivity based on mode and wall_mode
    if (edge_reconstruction_enabled || abacus_enabled) {
      // Edge reconstruction and abacus modes always use imposed temperature (no radiative)
      config.wall_parameters.emissivity = 0.0;
    } else if (base_enabled) {
      // Base mode: emissivity depends on wall_mode
      if (config.simulation.wall_mode == SimulationConfig::WallMode::Adiabatic ||
          config.simulation.wall_mode == SimulationConfig::WallMode::ImposedTemperature) {
        if (config.wall_parameters.emissivity > 0.0) {
          std::cerr << constants::string_processing::colors::red 
                    << "Warning: emissivity=" << config.wall_parameters.emissivity 
                    << " is ignored for wall_mode='" 
                    << (config.simulation.wall_mode == SimulationConfig::WallMode::Adiabatic ? "adiabatic" : "imposed_temperature")
                    << "'. Setting emissivity=0." 
                    << constants::string_processing::colors::reset << std::endl;
        }
        config.wall_parameters.emissivity = 0.0;
      } else if (config.simulation.wall_mode == SimulationConfig::WallMode::Radiative) {
        // Check that emissivity is valid for radiative mode
        if (config.wall_parameters.emissivity <= 0.0) {
          return std::unexpected(core::ConfigurationError("Radiative wall mode requires emissivity > 0"));
        }
      }
    }

    // Parse abacus config only if in abacus mode
    if (abacus_enabled) {
      
      auto abacus_result = parse_abacus_config(root_["abacus"]);
      if (!abacus_result) {
        return std::unexpected(abacus_result.error());
      }
      config.abacus = std::move(abacus_result.value());
    } else {
      // Default abacus config for non-abacus modes
      config.abacus = AbacusConfig{};
    }

    // Force parameters for abacus mode
    if (config.abacus.enabled) {
      config.simulation.catalytic_wall = true;
      config.simulation.only_stagnation_point = true;
      config.simulation.wall_mode = SimulationConfig::WallMode::ImposedTemperature;
      config.simulation.chemical_mode = SimulationConfig::ChemicalMode::NonEquilibrium;
      // x_stations are automatically derived from edge_points (single point at x=0 for stagnation)
      
      // Force edge_points with configured radius and velocity
      if (config.outer_edge.edge_points.empty()) {
        OuterEdgeConfig::EdgePoint point;
        point.x = 0.0;
        point.radius = config.abacus.radius;
        point.velocity = config.abacus.velocity;
        // Temperature and pressure will be set during abacus generation
        point.temperature = root_["abacus"]["boundary_conditions"]["temperature"] ? 
                           root_["abacus"]["boundary_conditions"]["temperature"].as<double>() : 5000.0;
        point.pressure = root_["abacus"]["boundary_conditions"]["pressure"] ? 
                        root_["abacus"]["boundary_conditions"]["pressure"].as<double>() : 7000.0;
        config.outer_edge.edge_points.push_back(point);
      } else {
        // Update existing edge_point with abacus parameters
        config.outer_edge.edge_points[0].radius = config.abacus.radius;
        config.outer_edge.edge_points[0].velocity = config.abacus.velocity;
        if (root_["abacus"]["boundary_conditions"]["temperature"]) {
          config.outer_edge.edge_points[0].temperature = root_["abacus"]["boundary_conditions"]["temperature"].as<double>();
        }
        if (root_["abacus"]["boundary_conditions"]["pressure"]) {
          config.outer_edge.edge_points[0].pressure = root_["abacus"]["boundary_conditions"]["pressure"].as<double>();
        }
      }
      
      // Extract catalyticity values for abacus (new schema preferred, fallback to legacy)
      if (root_["abacus"]["catalyticity_values"]) {
        config.abacus.catalyticity_values = root_["abacus"]["catalyticity_values"].as<std::vector<double>>();
        if (config.abacus.catalyticity_values.empty()) {
          return std::unexpected(core::ConfigurationError("catalyticity_values cannot be empty in abacus mode"));
        }
      } else if (root_["abacus"]["boundary_conditions"]["catalyticity_values"]) {
        config.abacus.catalyticity_values = root_["abacus"]["boundary_conditions"]["catalyticity_values"].as<std::vector<double>>();
        if (config.abacus.catalyticity_values.empty()) {
          return std::unexpected(core::ConfigurationError("catalyticity_values cannot be empty in abacus mode"));
        }
      }
      
      // Extract flow parameters for abacus
      if (root_["abacus"]["flow_parameters"]) {
        const auto flow_node = root_["abacus"]["flow_parameters"];
        if (flow_node["velocity_gradient_stagnation"]) {
          config.outer_edge.velocity_gradient_stagnation = flow_node["velocity_gradient_stagnation"].as<double>();
        }
        if (flow_node["freestream_density"]) {
          config.outer_edge.freestream_density = flow_node["freestream_density"].as<double>();
        }
        if (flow_node["freestream_velocity"]) {
          config.outer_edge.freestream_velocity = flow_node["freestream_velocity"].as<double>();
        }
      }
    }

    // Force non-equilibrium chemistry when boundary override is used
    bool has_boundary_override = false;
    for (const auto& point : config.outer_edge.edge_points) {
      if (point.boundary_override_enabled()) {
        has_boundary_override = true;
        break;
      }
    }
    if (has_boundary_override) {
      if (config.simulation.chemical_mode != SimulationConfig::ChemicalMode::NonEquilibrium) {
        std::cout << "INFO: Forcing chemical_mode to 'non_equilibrium' because boundary_override is used" << std::endl;
        config.simulation.chemical_mode = SimulationConfig::ChemicalMode::NonEquilibrium;
      }
    }

    // Parse continuation from mode-specific section if available, otherwise from root
    YAML::Node cont_node;
    if (base_enabled && root_["base"]["continuation"]) {
      cont_node = root_["base"]["continuation"];
    } else if (edge_reconstruction_enabled && root_["edge_reconstruction"]["continuation"]) {
      cont_node = root_["edge_reconstruction"]["continuation"];
    } else if (abacus_enabled && root_["abacus"]["continuation"]) {
      cont_node = root_["abacus"]["continuation"];
    } else {
      cont_node = root_["continuation"];
    }
    
    
    auto cont_result = parse_continuation_config(cont_node);
    if (!cont_result) {
      return std::unexpected(cont_result.error());
    }
    config.continuation = std::move(cont_result.value());

    // Parse edge reconstruction config only if in edge_reconstruction mode
    if (edge_reconstruction_enabled) {
      
      auto edge_recon_result = parse_edge_reconstruction_config(root_["edge_reconstruction"]);
      if (!edge_recon_result) {
        return std::unexpected(edge_recon_result.error());
      }
      config.edge_reconstruction = std::move(edge_recon_result.value());
      
      // Force parameters for edge reconstruction mode
      if (config.edge_reconstruction.enabled) {
        config.simulation.finite_thickness = true;
        config.simulation.only_stagnation_point = true;
        config.simulation.catalytic_wall = true;
        config.simulation.wall_mode = SimulationConfig::WallMode::ImposedTemperature;
        config.simulation.chemical_mode = SimulationConfig::ChemicalMode::NonEquilibrium;
        // x_stations are automatically derived from edge_points (single point at x=0 for stagnation)
        
        // Edge reconstruction has its own specialized config, no outer_edge needed
      }
    } else {
      // Default edge reconstruction config for non-edge_reconstruction modes
      config.edge_reconstruction = EdgeReconstructionConfig{};
    }
    
    // Parse GASP2 configuration
    YAML::Node gasp2_node;
    if (base_enabled && root_["base"]["gasp2"]) {
      gasp2_node = root_["base"]["gasp2"];
    } else if (edge_reconstruction_enabled && root_["edge_reconstruction"]["gasp2"]) {
      gasp2_node = root_["edge_reconstruction"]["gasp2"];
    } else if (abacus_enabled && root_["abacus"]["gasp2"]) {
      gasp2_node = root_["abacus"]["gasp2"];
    } else {
      gasp2_node = root_["gasp2"];
    }
    
    
    auto gasp2_result = parse_gasp2_config(gasp2_node);
    if (!gasp2_result) {
      return std::unexpected(gasp2_result.error());
    }
    config.gasp2 = std::move(gasp2_result.value());

    // Parse Mutation++ configuration
    if (abacus_enabled) {
      // In abacus mode, mutation config is not required; default it
      config.mutation = MutationConfig{};
    } else {
      YAML::Node mutation_node;
      if (base_enabled && root_["base"]["mutation"]) {
        mutation_node = root_["base"]["mutation"];
      } else if (edge_reconstruction_enabled && root_["edge_reconstruction"]["mutation"]) {
        mutation_node = root_["edge_reconstruction"]["mutation"];
      } else {
        mutation_node = root_["mutation"];
      }

      
      auto mutation_result = parse_mutation_config(mutation_node);
      if (!mutation_result) {
        return std::unexpected(mutation_result.error());
      }
      config.mutation = std::move(mutation_result.value());
    }

    
    return config;

  } catch (const YAML::Exception& e) {
    return std::unexpected(
        core::ConfigurationError(std::format("YAML parsing error: {} at line {}", e.what(), e.mark.line)));
  } catch (const std::exception& e) {
    return std::unexpected(core::ConfigurationError(std::format("Unexpected error during parsing: {}", e.what())));
  }
}

auto YamlParser::parse_simulation_config(const YAML::Node& node, bool abacus_mode, bool edge_reconstruction_mode) const
    -> std::expected<SimulationConfig, core::ConfigurationError> {

  SimulationConfig config;

  try {
    auto body_type_result = extract_enum(node, "body_type", enum_mappings::body_types);
    if (!body_type_result) {
      return std::unexpected(body_type_result.error());
    }
    config.body_type = body_type_result.value();

    // only_stagnation_point is required for base mode, but forced to true for abacus and edge_reconstruction modes
    if (abacus_mode || edge_reconstruction_mode) {
      // Force to true for specialized modes
      config.only_stagnation_point = true;
      // Warn if user specified a different value
      if (node["only_stagnation_point"]) {
        bool user_value = node["only_stagnation_point"].as<bool>();
        if (!user_value) {
          std::cout << "INFO: Ignoring 'only_stagnation_point: false' in " 
                    << (abacus_mode ? "abacus" : "edge_reconstruction") 
                    << " mode. Forcing to true." << std::endl;
        }
      }
    } else {
      // Base mode - only_stagnation_point is required
      auto only_stag_result = extract_value<bool>(node, "only_stagnation_point");
      if (!only_stag_result) {
        return std::unexpected(only_stag_result.error());
      }
      config.only_stagnation_point = only_stag_result.value();
    }

    // Parse finite thickness - optional parameter, defaults to false
    if (node["finite_thickness"]) {
      auto finite_thickness_result = extract_value<bool>(node, "finite_thickness");
      if (!finite_thickness_result) {
        return std::unexpected(finite_thickness_result.error());
      }
      config.finite_thickness = finite_thickness_result.value();
    }

    auto diffusion_result = extract_enum(node, "diffusion_type", enum_mappings::diffusion_types);
    if (!diffusion_result) {
      return std::unexpected(diffusion_result.error());
    }
    config.diffusion_type = diffusion_result.value();

    auto thermal_diff_result = extract_value<bool>(node, "consider_thermal_diffusion");
    if (!thermal_diff_result) {
      return std::unexpected(thermal_diff_result.error());
    }
    config.consider_thermal_diffusion = thermal_diff_result.value();

    auto dufour_result = extract_value<bool>(node, "consider_dufour_effect");
    if (!dufour_result) {
      return std::unexpected(dufour_result.error());
    }
    config.consider_dufour_effect = dufour_result.value();

    // chemical_mode is required for base mode, but forced to non_equilibrium for abacus and edge_reconstruction modes
    if (abacus_mode || edge_reconstruction_mode) {
      // Force to NonEquilibrium for specialized modes
      config.chemical_mode = SimulationConfig::ChemicalMode::NonEquilibrium;
      // Warn if user specified a different value
      if (node["chemical_mode"]) {
        std::string user_value = node["chemical_mode"].as<std::string>();
        if (user_value != "non_equilibrium" && user_value != "nonequilibrium") {
          std::cout << "INFO: Ignoring 'chemical_mode: " << user_value << "' in " 
                    << (abacus_mode ? "abacus" : "edge_reconstruction") 
                    << " mode. Forcing to non_equilibrium." << std::endl;
        }
      }
    } else {
      // Base mode - chemical_mode is required
      auto chemical_result = extract_enum(node, "chemical_mode", enum_mappings::chemical_modes);
      if (!chemical_result) {
        return std::unexpected(chemical_result.error());
      }
      config.chemical_mode = chemical_result.value();
    }

    // catalytic_wall is optional for base mode (defaults to true), but forced to true for abacus and edge_reconstruction modes
    if (abacus_mode || edge_reconstruction_mode) {
      // Force to true for specialized modes
      config.catalytic_wall = true;
      // Warn if user specified false
      if (node["catalytic_wall"]) {
        bool user_value = node["catalytic_wall"].as<bool>();
        if (!user_value) {
          std::cout << "INFO: Ignoring 'catalytic_wall: false' in " 
                    << (abacus_mode ? "abacus" : "edge_reconstruction") 
                    << " mode. Forcing to true." << std::endl;
        }
      }
      
      // Force catalysis_provider to MutationPP for specialized modes
      config.catalysis_provider = SimulationConfig::CatalysisProvider::MutationPP;
    } else {
      // Base mode - catalytic_wall is optional, defaults to true
      if (node["catalytic_wall"]) {
        auto catalytic_result = extract_value<bool>(node, "catalytic_wall");
        if (!catalytic_result) {
          return std::unexpected(catalytic_result.error());
        }
        config.catalytic_wall = catalytic_result.value();
      } else {
        config.catalytic_wall = true;  // Default to catalytic
      }
    }

    // Parse wall_mode - optional parameter, defaults to ImposedTemperature
    if (node["wall_mode"]) {
      auto wall_mode_result = extract_enum(node, "wall_mode", enum_mappings::wall_modes);
      if (!wall_mode_result) {
        return std::unexpected(wall_mode_result.error());
      }
      config.wall_mode = wall_mode_result.value();
    }
    
    // Parse catalysis_provider - optional parameter, defaults to MutationPP
    if (node["catalysis_provider"]) {
      auto provider_result = extract_enum(node, "catalysis_provider", enum_mappings::catalysis_providers);
      if (!provider_result) {
        return std::unexpected(provider_result.error());
      }
      config.catalysis_provider = provider_result.value();
    }

    // Legacy support: if old "adiabatic" parameter is present, convert to wall_mode
    if (node["adiabatic"]) {
      auto adiabatic_result = extract_value<bool>(node, "adiabatic");
      if (!adiabatic_result) {
        return std::unexpected(adiabatic_result.error());
      }
      // Convert legacy adiabatic flag to wall_mode
      if (adiabatic_result.value()) {
        config.wall_mode = SimulationConfig::WallMode::Adiabatic;
      } else {
        config.wall_mode = SimulationConfig::WallMode::ImposedTemperature;
      }
    }

    return config;

  } catch (const core::ConfigurationError& e) {
    return std::unexpected(core::ConfigurationError(std::format("In 'simulation' section: {}", e.message())));
  }
}

auto YamlParser::parse_numerical_config(const YAML::Node& node) const
    -> std::expected<NumericalConfig, core::ConfigurationError> {

  NumericalConfig config;

  try {
    auto n_eta_result = extract_value<int>(node, "n_eta");
    if (!n_eta_result)
      return std::unexpected(n_eta_result.error());
    config.n_eta = n_eta_result.value();

    auto eta_max_result = extract_value<double>(node, "eta_max");
    if (!eta_max_result)
      return std::unexpected(eta_max_result.error());
    config.eta_max = eta_max_result.value();

    auto conv_tol_result = extract_value<double>(node, "convergence_tolerance");
    if (!conv_tol_result)
      return std::unexpected(conv_tol_result.error());
    config.convergence_tolerance = conv_tol_result.value();

    auto max_iter_result = extract_value<int>(node, "max_iterations");
    if (!max_iter_result)
      return std::unexpected(max_iter_result.error());
    config.max_iterations = max_iter_result.value();

    return config;

  } catch (const core::ConfigurationError& e) {
    return std::unexpected(core::ConfigurationError(std::format("In 'numerical' section: {}", e.message())));
  }
}

auto YamlParser::parse_mixture_config(const YAML::Node& node) const
    -> std::expected<MixtureConfig, core::ConfigurationError> {

  MixtureConfig config;

  try {
    auto name_result = extract_value<std::string>(node, "name");
    if (!name_result)
      return std::unexpected(name_result.error());
    config.name = name_result.value();

    auto db_result = extract_enum(node, "thermodynamic_database", enum_mappings::databases);
    if (!db_result)
      return std::unexpected(db_result.error());
    config.thermodynamic_database = db_result.value();

    auto visc_result = extract_enum(node, "viscosity_algorithm", enum_mappings::viscosity_algorithms);
    if (!visc_result)
      return std::unexpected(visc_result.error());
    config.viscosity_algorithm = visc_result.value();

    // Optional thermal conductivity algorithm (defaults to chapmanEnskog_CG)
    if (node["thermal_conductivity_algorithm"]) {
      auto thermal_result = extract_enum(node, "thermal_conductivity_algorithm", enum_mappings::thermal_conductivity_algorithms);
      if (!thermal_result)
        return std::unexpected(thermal_result.error());
      config.thermal_conductivity_algorithm = thermal_result.value();
    }


    return config;

  } catch (const core::ConfigurationError& e) {
    return std::unexpected(core::ConfigurationError(std::format("In 'mixture' section: {}", e.message())));
  }
}

auto YamlParser::parse_output_config(const YAML::Node& node) const
    -> std::expected<OutputConfig, core::ConfigurationError> {

  OutputConfig config;

  try {
    // x_stations are now automatically derived from edge_points, no need to parse them

    auto output_dir_result = extract_value<std::string>(node, "output_directory");
    if (!output_dir_result)
      return std::unexpected(output_dir_result.error());
    config.output_directory = output_dir_result.value();

    return config;

  } catch (const core::ConfigurationError& e) {
    return std::unexpected(core::ConfigurationError(std::format("In 'output' section: {}", e.message())));
  }
}

auto YamlParser::parse_outer_edge_config(const YAML::Node& node, bool edge_reconstruction_mode, bool abacus_mode) const
    -> std::expected<OuterEdgeConfig, core::ConfigurationError> {

  OuterEdgeConfig config;

  try {
    if (!node["edge_points"]) {
      return std::unexpected(core::ConfigurationError("Missing required 'edge_points' in outer_edge"));
    }

    auto points_node = node["edge_points"];
    for (const auto& point_node : points_node) {
      OuterEdgeConfig::EdgePoint point;

      // Handle x coordinate with mode-specific implicit values
      if (edge_reconstruction_mode || abacus_mode) {
        // For edge reconstruction and abacus modes, always use stagnation point (x=0)
        point.x = 0.0;
        // Warn if x is specified in config but will be ignored
        if (point_node["x"]) {
          std::cout << "INFO: Ignoring 'x' value in edge_points for " << (edge_reconstruction_mode ? "edge_reconstruction" : "abacus") << " mode. Using x=0.0 (stagnation point)." << std::endl;
        }
      } else {
        // Base mode - x is required
        auto x_result = extract_value<double>(point_node, "x");
        if (!x_result)
          return std::unexpected(x_result.error());
        point.x = x_result.value();
      }

      auto radius_result = extract_value<double>(point_node, "radius");
      if (!radius_result)
        return std::unexpected(radius_result.error());
      point.radius = radius_result.value();

      auto velocity_result = extract_value<double>(point_node, "velocity");
      if (!velocity_result)
        return std::unexpected(velocity_result.error());
      point.velocity = velocity_result.value();

      auto temperature_result = extract_value<double>(point_node, "temperature");
      if (!temperature_result)
        return std::unexpected(temperature_result.error());
      point.temperature = temperature_result.value();

      auto pressure_result = extract_value<double>(point_node, "pressure");
      if (!pressure_result)
        return std::unexpected(pressure_result.error());
      point.pressure = pressure_result.value();

      if (point_node["boundary_override"]) {
        auto boundary_override_node = point_node["boundary_override"];

        // Handle both old format (boolean) and new format (object)
        if (boundary_override_node.IsScalar()) {
          // Legacy format: boundary_override: true
          point.boundary_override.enabled = boundary_override_node.as<bool>();

          // Look for mass_fraction_condition at the same level (legacy)
          if (point_node["mass_fraction_condition"]) {
            if (point.boundary_override.enabled) {
              auto mass_fractions = point_node["mass_fraction_condition"].as<std::vector<double>>();

              double sum = std::accumulate(mass_fractions.begin(), mass_fractions.end(), constants::defaults::default_emissivity);
              if (std::abs(sum - constants::hermite::basis_functions::h00_constant) > constants::tolerance::mass_fraction_sum) {
                return std::unexpected(core::ConfigurationError(
                    std::format("mass_fraction_condition must sum to {} (current sum: {})", constants::hermite::basis_functions::h00_constant, sum)));
              }

              if (std::any_of(mass_fractions.begin(), mass_fractions.end(), [](double val) { return val < constants::defaults::default_emissivity; })) {
                return std::unexpected(core::ConfigurationError("All mass fractions must be non-negative"));
              }

              point.boundary_override.mass_fraction_condition = mass_fractions;
            }
          }
        } else {
          // New format: boundary_override: { enabled: true, mass_fraction_condition: [...] }
          if (boundary_override_node["enabled"]) {
            point.boundary_override.enabled = boundary_override_node["enabled"].as<bool>();
          }

          if (point.boundary_override.enabled) {
            if (!boundary_override_node["mass_fraction_condition"]) {
              return std::unexpected(core::ConfigurationError(
                  "mass_fraction_condition is required when boundary_override.enabled is true"));
            }

            auto mass_fractions = boundary_override_node["mass_fraction_condition"].as<std::vector<double>>();

            double sum = std::accumulate(mass_fractions.begin(), mass_fractions.end(), constants::defaults::default_emissivity);
            if (std::abs(sum - constants::hermite::basis_functions::h00_constant) > constants::tolerance::mass_fraction_sum) {
              return std::unexpected(core::ConfigurationError(
                  std::format("mass_fraction_condition must sum to {} (current sum: {})", constants::hermite::basis_functions::h00_constant, sum)));
            }

            if (std::any_of(mass_fractions.begin(), mass_fractions.end(), [](double val) { return val < constants::defaults::default_emissivity; })) {
              return std::unexpected(core::ConfigurationError("All mass fractions must be non-negative"));
            }

            point.boundary_override.mass_fraction_condition = mass_fractions;
          }
        }
      }

      config.edge_points.push_back(point);
    }

    auto vel_grad_result = extract_value<double>(node, "velocity_gradient_stagnation");
    if (!vel_grad_result)
      return std::unexpected(vel_grad_result.error());
    config.velocity_gradient_stagnation = vel_grad_result.value();

    auto freestream_dens_result = extract_value<double>(node, "freestream_density");
    if (!freestream_dens_result)
      return std::unexpected(freestream_dens_result.error());
    config.freestream_density = freestream_dens_result.value();

    auto freestream_vel_result = extract_value<double>(node, "freestream_velocity");
    if (!freestream_vel_result)
      return std::unexpected(freestream_vel_result.error());
    config.freestream_velocity = freestream_vel_result.value();

    // Parse finite thickness parameters if present
    if (node["finite_thickness_params"]) {
      auto ft_node = node["finite_thickness_params"];
      
      auto v_edge_result = extract_value<double>(ft_node, "v_edge");
      if (!v_edge_result)
        return std::unexpected(v_edge_result.error());
      config.finite_thickness_params.v_edge = v_edge_result.value();
      
      auto d2_ue_dxdy_result = extract_value<double>(ft_node, "d2_ue_dxdy");
      if (!d2_ue_dxdy_result)
        return std::unexpected(d2_ue_dxdy_result.error());
      config.finite_thickness_params.d2_ue_dxdy = d2_ue_dxdy_result.value();
      
      auto delta_bl_result = extract_value<double>(ft_node, "delta_bl");
      if (!delta_bl_result)
        return std::unexpected(delta_bl_result.error());
      config.finite_thickness_params.delta_bl = delta_bl_result.value();
    }

    return config;

  } catch (const core::ConfigurationError& e) {
    return std::unexpected(core::ConfigurationError(std::format("In 'outer_edge' section: {}", e.message())));
  } catch (const YAML::Exception& e) {
    return std::unexpected(core::ConfigurationError(std::format("In 'outer_edge' section: YAML error - {}", e.what())));
  }
}

auto YamlParser::parse_wall_parameters_config(const YAML::Node& node) const
    -> std::expected<WallParametersConfig, core::ConfigurationError> {

  WallParametersConfig config;

  try {
    auto temp_result = extract_value<std::vector<double>>(node, "temperatures");
    if (!temp_result)
      return std::unexpected(temp_result.error());
    config.wall_temperatures = temp_result.value();

    if (config.wall_temperatures.empty()) {
      return std::unexpected(core::ConfigurationError("Wall temperatures cannot be empty"));
    }
    
    if (node["emissivity"]) {
      auto emiss_result = extract_value<double>(node, "emissivity");
      if (!emiss_result)
        return std::unexpected(emiss_result.error());
      config.emissivity = emiss_result.value();
      
      if (config.emissivity < constants::defaults::default_emissivity || config.emissivity > constants::defaults::max_emissivity) {
        return std::unexpected(core::ConfigurationError(std::format("emissivity must be between {} and {}", constants::defaults::default_emissivity, constants::defaults::max_emissivity)));
      }
    }

    if (node["environment_temperature"]) {
      auto env_temp_result = extract_value<double>(node, "environment_temperature");
      if (!env_temp_result)
        return std::unexpected(env_temp_result.error());
      config.environment_temperature = env_temp_result.value();
      
      if (config.environment_temperature <= constants::defaults::default_emissivity) {
        return std::unexpected(core::ConfigurationError("environment_temperature must be positive"));
      }
    }
    
    return config;

  } catch (const core::ConfigurationError& e) {
    return std::unexpected(core::ConfigurationError(std::format("In 'wall_parameters' section: {}", e.message())));
  }
}

auto YamlParser::parse_abacus_config(const YAML::Node& node) const
    -> std::expected<AbacusConfig, core::ConfigurationError> {

  AbacusConfig config;

  // If no abacus section, return default (disabled)
  if (!node) {
    return config;
  }

  try {
    // Parse enabled flag
    if (node["enabled"]) {
      config.enabled = node["enabled"].as<bool>();
    }

    // Only parse other fields if enabled
    if (config.enabled) {
      // Default catalyticity values (overridden below if provided)
      config.catalyticity_values = {0.0, 0.1};

      // Parse mode if present (supports both new and legacy naming)
      if (node["mode"]) {
        auto mode_result = extract_enum(node, "mode", enum_mappings::abacus_modes);
        if (!mode_result) {
          return std::unexpected(core::ConfigurationError(std::format("In 'abacus' section: {}", mode_result.error().message())));
        }
        config.mode = mode_result.value();
      }

      // Respect selected mode and ignore unrelated blocks
      if (config.mode == AbacusConfig::Mode::TemperatureSweep) {
        // New schema: temperature_sweep { range: [min, max], points: N }
        if (node["temperature_sweep"]) {
          auto ts_node = node["temperature_sweep"];
          if (ts_node["range"]) {
            auto range = ts_node["range"].as<std::vector<double>>();
            if (range.size() != constants::indexing::min_range_size) {
              return std::unexpected(core::ConfigurationError(std::format("temperature_sweep.range must have exactly {} values", constants::indexing::min_range_size)));
            }
            config.temperature_min = range[constants::indexing::first];
            config.temperature_max = range[constants::indexing::second];
            if (config.temperature_min >= config.temperature_max) {
              return std::unexpected(core::ConfigurationError("temperature_min must be less than temperature_max"));
            }
          }
          if (ts_node["points"]) {
            config.temperature_points = ts_node["points"].as<int>();
            if (config.temperature_points < constants::indexing::min_range_size) {
              return std::unexpected(core::ConfigurationError(std::format("temperature_sweep.points must be at least {}", constants::indexing::min_range_size)));
            }
          }
        }

        // Legacy fallbacks only for temperature_sweep mode
        if (node["temperature_points"]) {
          config.temperature_points = node["temperature_points"].as<int>();
          if (config.temperature_points < constants::indexing::min_range_size) {
            return std::unexpected(core::ConfigurationError(std::format("temperature_points must be at least {}", constants::indexing::min_range_size)));
          }
        }
        if (node["boundary_conditions"]) {
          auto bc_node = node["boundary_conditions"];
          if (bc_node["temperature_range"]) {
            auto range = bc_node["temperature_range"].as<std::vector<double>>();
            if (range.size() != constants::indexing::min_range_size) {
              return std::unexpected(core::ConfigurationError(std::format("temperature_range must have exactly {} values", constants::indexing::min_range_size)));
            }
            config.temperature_min = range[constants::indexing::first];
            config.temperature_max = range[constants::indexing::second];
            if (config.temperature_min >= config.temperature_max) {
              return std::unexpected(core::ConfigurationError("temperature_min must be less than temperature_max"));
            }
          }
        }
      } else if (config.mode == AbacusConfig::Mode::GammaSweepAtTw) {
        // New schema: gamma_sweep { wall_temperature: Tw }
        if (node["gamma_sweep"]) {
          auto gs_node = node["gamma_sweep"];
          if (gs_node["wall_temperature"]) {
            // Accept scalar or sequence of temperatures
            if (gs_node["wall_temperature"].IsSequence()) {
              auto temps = gs_node["wall_temperature"].as<std::vector<double>>();
              if (temps.empty()) {
                return std::unexpected(core::ConfigurationError("gamma_sweep.wall_temperature list cannot be empty"));
              }
              config.temperatures_override = temps;
              config.temperature_points = static_cast<int>(temps.size());
              config.temperature_min = temps.front();
              config.temperature_max = temps.back();
            } else {
              const double Tw = gs_node["wall_temperature"].as<double>();
              config.temperatures_override = {Tw};
              config.temperature_min = Tw;
              config.temperature_max = Tw;
              config.temperature_points = 1; // single temperature
            }
          } else {
            return std::unexpected(core::ConfigurationError("gamma_sweep.wall_temperature is required for gamma_sweep_at_Tw mode"));
          }
        } else {
          return std::unexpected(core::ConfigurationError("gamma_sweep section is required for gamma_sweep_at_Tw mode"));
        }
      }

      // Geometry/flow optional parameters shared across modes
      if (node["boundary_conditions"]) {
        auto bc_node = node["boundary_conditions"];
        // Optional radius (defaults to 1.0)
        if (bc_node["radius"]) {
          config.radius = bc_node["radius"].as<double>();
        }
        // Optional velocity (defaults to 0.0)
        if (bc_node["velocity"]) {
          config.velocity = bc_node["velocity"].as<double>();
        }
      }

      // Catalyticity values (new schema preferred, fallback to legacy)
      if (node["catalyticity_values"]) {
        config.catalyticity_values = node["catalyticity_values"].as<std::vector<double>>();
        if (config.catalyticity_values.empty()) {
          return std::unexpected(core::ConfigurationError("catalyticity_values cannot be empty in abacus mode"));
        }
      } else if (node["boundary_conditions"] && node["boundary_conditions"]["catalyticity_values"]) {
        config.catalyticity_values = node["boundary_conditions"]["catalyticity_values"].as<std::vector<double>>();
        if (config.catalyticity_values.empty()) {
          return std::unexpected(core::ConfigurationError("catalyticity_values cannot be empty in abacus mode"));
        }
      }
    }

    return config;

  } catch (const YAML::Exception& e) {
    return std::unexpected(core::ConfigurationError(std::format("In 'abacus' section: {}", e.what())));
  }
}

auto YamlParser::parse_continuation_config(const YAML::Node& node) const
    -> std::expected<ContinuationConfig, core::ConfigurationError> {

  ContinuationConfig config;

  if (!node) {
    return config;
  }

  try {
    if (node["wall_temperature_stable"]) {
      config.wall_temperature_stable = node["wall_temperature_stable"].as<double>();
    }
    if (node["edge_temperature_stable"]) {
      config.edge_temperature_stable = node["edge_temperature_stable"].as<double>();
    }
    if (node["pressure_stable"]) {
      config.pressure_stable = node["pressure_stable"].as<double>();
    }
    if (node["use_linear_predictor"]) {
      config.use_linear_predictor = node["use_linear_predictor"].as<bool>();
    }
    return config;
  } catch (const YAML::Exception& e) {
    return std::unexpected(
        core::ConfigurationError(std::format("In 'continuation' section: failed to parse - {}", e.what())));
  }
}

auto YamlParser::parse_edge_reconstruction_config(const YAML::Node& node) const
    -> std::expected<EdgeReconstructionConfig, core::ConfigurationError> {

  EdgeReconstructionConfig config;

  if (!node) {
    return config;
  }

  try {
    if (node["enabled"]) {
      config.enabled = node["enabled"].as<bool>();
    }

    if (config.enabled) {
      // Target values
      if (!node["target_heat_flux"]) {
        return std::unexpected(
            core::ConfigurationError("Missing required 'target_heat_flux' in edge_reconstruction"));
      }
      config.target_heat_flux = node["target_heat_flux"].as<double>();

      // Boundary conditions
      if (!node["boundary_conditions"]) {
        return std::unexpected(
            core::ConfigurationError("Missing required 'boundary_conditions' in edge_reconstruction"));
      }
      auto bc_node = node["boundary_conditions"];
      if (!bc_node["pressure"]) {
        return std::unexpected(
            core::ConfigurationError("Missing required 'pressure' in edge_reconstruction.boundary_conditions"));
      }
      config.boundary_conditions.pressure = bc_node["pressure"].as<double>();
      
      if (!bc_node["wall_temperature"]) {
        return std::unexpected(
            core::ConfigurationError("Missing required 'wall_temperature' in edge_reconstruction.boundary_conditions"));
      }
      config.boundary_conditions.wall_temperature = bc_node["wall_temperature"].as<double>();
      
      if (!bc_node["catalyticity"]) {
        return std::unexpected(
            core::ConfigurationError("Missing required 'catalyticity' in edge_reconstruction.boundary_conditions"));
      }
      config.boundary_conditions.catalyticity = bc_node["catalyticity"].as<double>();
      
      // Optional radius (defaults to 1.0)
      if (bc_node["radius"]) {
        config.boundary_conditions.radius = bc_node["radius"].as<double>();
      }

      // Flow parameters
      if (!node["flow_parameters"]) {
        return std::unexpected(
            core::ConfigurationError("Missing required 'flow_parameters' in edge_reconstruction"));
      }
      auto flow_node = node["flow_parameters"];
      if (!flow_node["velocity_gradient_stagnation"]) {
        return std::unexpected(
            core::ConfigurationError("Missing required 'velocity_gradient_stagnation' in edge_reconstruction.flow_parameters"));
      }
      config.flow_parameters.velocity_gradient_stagnation = flow_node["velocity_gradient_stagnation"].as<double>();
      
      if (!flow_node["freestream_density"]) {
        return std::unexpected(
            core::ConfigurationError("Missing required 'freestream_density' in edge_reconstruction.flow_parameters"));
      }
      config.flow_parameters.freestream_density = flow_node["freestream_density"].as<double>();
      
      if (!flow_node["freestream_velocity"]) {
        return std::unexpected(
            core::ConfigurationError("Missing required 'freestream_velocity' in edge_reconstruction.flow_parameters"));
      }
      config.flow_parameters.freestream_velocity = flow_node["freestream_velocity"].as<double>();

      // Finite thickness params
      if (!node["finite_thickness_params"]) {
        return std::unexpected(
            core::ConfigurationError("Missing required 'finite_thickness_params' in edge_reconstruction"));
      }
      auto ft_node = node["finite_thickness_params"];
      if (!ft_node["v_edge"]) {
        return std::unexpected(
            core::ConfigurationError("Missing required 'v_edge' in edge_reconstruction.finite_thickness_params"));
      }
      config.finite_thickness_params.v_edge = ft_node["v_edge"].as<double>();
      
      if (!ft_node["d2_ue_dxdy"]) {
        return std::unexpected(
            core::ConfigurationError("Missing required 'd2_ue_dxdy' in edge_reconstruction.finite_thickness_params"));
      }
      config.finite_thickness_params.d2_ue_dxdy = ft_node["d2_ue_dxdy"].as<double>();
      
      if (!ft_node["delta_bl"]) {
        return std::unexpected(
            core::ConfigurationError("Missing required 'delta_bl' in edge_reconstruction.finite_thickness_params"));
      }
      config.finite_thickness_params.delta_bl = ft_node["delta_bl"].as<double>();

      // Solver settings (optional with defaults)
      if (node["solver"]) {
        auto solver_node = node["solver"];
        if (solver_node["initial_temperature_guess"]) {
          config.solver.initial_temperature_guess = solver_node["initial_temperature_guess"].as<double>();
        }
        if (solver_node["temperature_min"]) {
          config.solver.temperature_min = solver_node["temperature_min"].as<double>();
        }
        if (solver_node["temperature_max"]) {
          config.solver.temperature_max = solver_node["temperature_max"].as<double>();
        }
        if (solver_node["tolerance"]) {
          config.solver.tolerance = solver_node["tolerance"].as<double>();
        }
        if (solver_node["max_iterations"]) {
          config.solver.max_iterations = solver_node["max_iterations"].as<int>();
        }
      }
    }

    return config;
  } catch (const YAML::Exception& e) {
    return std::unexpected(
        core::ConfigurationError(std::format("In 'edge_reconstruction' section: failed to parse - {}", e.what())));
  }
}

auto YamlParser::parse_gasp2_config(const YAML::Node& node) const
    -> std::expected<Gasp2Config, core::ConfigurationError> {
  try {
    Gasp2Config config;
    
    if (!node || node.IsNull()) {
      // No gasp2 section - return defaults
      return config;
    }

    if (node["xml_file"]) {
      config.xml_file = node["xml_file"].as<std::string>();
    }

    return config;
  } catch (const YAML::Exception& e) {
    return std::unexpected(
        core::ConfigurationError(std::format("In 'gasp2' section: failed to parse - {}", e.what())));
  }
}

auto YamlParser::parse_mutation_config(const YAML::Node& node) const
    -> std::expected<MutationConfig, core::ConfigurationError> {
  try {
    MutationConfig config;
    
    if (!node || node.IsNull()) {
      // No mutation section - return defaults
      return config;
    }

    if (node["gsi_file"]) {
      config.gsi_file = node["gsi_file"].as<std::string>();
    }
    
    if (node["reload_gsi_each_step"]) {
      config.reload_gsi_each_step = node["reload_gsi_each_step"].as<bool>();
    }

    return config;
  } catch (const YAML::Exception& e) {
    return std::unexpected(
        core::ConfigurationError(std::format("In 'mutation' section: failed to parse - {}", e.what())));
  }
}

} // namespace blast::io
