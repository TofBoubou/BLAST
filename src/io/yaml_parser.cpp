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

    // Simulation
    if (!root_["simulation"]) { // !root_ is override in YAML so we check if
                                // it's empty
      return std::unexpected(core::ConfigurationError("Missing required 'simulation' section"));
    }
    auto sim_result = parse_simulation_config(root_["simulation"]);
    if (!sim_result) {
      return std::unexpected(sim_result.error());
    }
    config.simulation = std::move(sim_result.value()); // emptying sim_result, avoid a copy

    // Numerical parameters
    if (!root_["numerical"]) {
      return std::unexpected(core::ConfigurationError("Missing required 'numerical' section"));
    }
    auto num_result = parse_numerical_config(root_["numerical"]); // The root changes
    if (!num_result) {
      return std::unexpected(num_result.error());
    }
    config.numerical = std::move(num_result.value());

    // Mixture
    if (!root_["mixture"]) {
      return std::unexpected(core::ConfigurationError("Missing required 'mixture' section"));
    }
    auto mix_result = parse_mixture_config(root_["mixture"]);
    if (!mix_result) {
      return std::unexpected(mix_result.error());
    }
    config.mixture = std::move(mix_result.value());

    // Output
    if (!root_["output"]) {
      return std::unexpected(core::ConfigurationError("Missing required 'output' section"));
    }
    auto out_result = parse_output_config(root_["output"]);
    if (!out_result) {
      return std::unexpected(out_result.error());
    }
    config.output = std::move(out_result.value());

    // outer_edge
    if (!root_["outer_edge"]) {
      return std::unexpected(core::ConfigurationError("Missing required 'outer_edge' section"));
    }
    auto edge_result = parse_outer_edge_config(root_["outer_edge"]);
    if (!edge_result) {
      return std::unexpected(edge_result.error());
    }
    config.outer_edge = std::move(edge_result.value());

    if (!root_["wall_parameters"]) {
      return std::unexpected(core::ConfigurationError("Missing required 'wall_parameters' section"));
    }
    auto wall_result = parse_wall_parameters_config(root_["wall_parameters"]);
    if (!wall_result) {
      return std::unexpected(wall_result.error());
    }
    config.wall_parameters = std::move(wall_result.value());

    auto abacus_result = parse_abacus_config(root_["abacus"]);
    if (!abacus_result) {
      return std::unexpected(abacus_result.error());
    }
    config.abacus = std::move(abacus_result.value());

    // Force catalytic wall when abacus generation is enabled
    if (config.abacus.enabled) {
      config.simulation.catalytic_wall = true;
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

    auto cont_result = parse_continuation_config(root_["continuation"]);
    if (!cont_result) {
      return std::unexpected(cont_result.error());
    }
    config.continuation = std::move(cont_result.value());

    // Edge reconstruction (optional)
    if (root_["edge_reconstruction"]) {
      auto edge_recon_result = parse_edge_reconstruction_config(root_["edge_reconstruction"]);
      if (!edge_recon_result) {
        return std::unexpected(edge_recon_result.error());
      }
      config.edge_reconstruction = std::move(edge_recon_result.value());
      
      // Force finite thickness when edge reconstruction is enabled
      if (config.edge_reconstruction.enabled) {
        config.simulation.finite_thickness = true;
      }
    }

    return config;

  } catch (const YAML::Exception& e) {
    return std::unexpected(
        core::ConfigurationError(std::format("YAML parsing error: {} at line {}", e.what(), e.mark.line)));
  } catch (const std::exception& e) {
    return std::unexpected(core::ConfigurationError(std::format("Unexpected error during parsing: {}", e.what())));
  }
}

auto YamlParser::parse_simulation_config(const YAML::Node& node) const
    -> std::expected<SimulationConfig, core::ConfigurationError> {

  SimulationConfig config;

  try {
    auto body_type_result = extract_enum(node, "body_type", enum_mappings::body_types);
    if (!body_type_result) {
      return std::unexpected(body_type_result.error());
    }
    config.body_type = body_type_result.value();

    auto only_stag_result = extract_value<bool>(node, "only_stagnation_point");
    if (!only_stag_result) {
      return std::unexpected(only_stag_result.error());
    }
    config.only_stagnation_point = only_stag_result.value();

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

    auto chemical_result = extract_enum(node, "chemical_mode", enum_mappings::chemical_modes);
    if (!chemical_result) {
      return std::unexpected(chemical_result.error());
    }
    config.chemical_mode = chemical_result.value();

    auto catalytic_result = extract_value<bool>(node, "catalytic_wall");
    if (!catalytic_result) {
      return std::unexpected(catalytic_result.error());
    }
    config.catalytic_wall = catalytic_result.value();

    if (node["adiabatic"]) {
      auto adiabatic_result = extract_value<bool>(node, "adiabatic");
      if (!adiabatic_result) {
        return std::unexpected(adiabatic_result.error());
      }
      config.adiabatic = adiabatic_result.value();
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

    // Optional state model selection (defaults to ChemNonEq1T)
    if (node["state_model"]) {
      auto state_result = extract_enum(node, "state_model", enum_mappings::state_models);
      if (!state_result)
        return std::unexpected(state_result.error());
      config.state_model = state_result.value();
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
    auto x_stations_result = extract_value<std::vector<double>>(node, "x_stations");
    if (!x_stations_result)
      return std::unexpected(x_stations_result.error());
    config.x_stations = x_stations_result.value();
    if (config.x_stations.empty()) {
      return std::unexpected(core::ConfigurationError("x_stations cannot be empty"));
    }
    if (std::any_of(config.x_stations.begin(), config.x_stations.end(), [](double x) { return x < 0; })) {
      return std::unexpected(core::ConfigurationError("x_stations must be non-negative"));
    }
    for (std::size_t i = 1; i < config.x_stations.size(); ++i) {
      if (config.x_stations[i] <= config.x_stations[i - 1]) {
        return std::unexpected(core::ConfigurationError("x_stations must be in strictly increasing order"));
      }
    }

    auto output_dir_result = extract_value<std::string>(node, "output_directory");
    if (!output_dir_result)
      return std::unexpected(output_dir_result.error());
    config.output_directory = output_dir_result.value();

    return config;

  } catch (const core::ConfigurationError& e) {
    return std::unexpected(core::ConfigurationError(std::format("In 'output' section: {}", e.message())));
  }
}

auto YamlParser::parse_outer_edge_config(const YAML::Node& node) const
    -> std::expected<OuterEdgeConfig, core::ConfigurationError> {

  OuterEdgeConfig config;

  try {
    if (!node["edge_points"]) {
      return std::unexpected(core::ConfigurationError("Missing required 'edge_points' in outer_edge"));
    }

    auto points_node = node["edge_points"];
    for (const auto& point_node : points_node) {
      OuterEdgeConfig::EdgePoint point;

      auto x_result = extract_value<double>(point_node, "x");
      if (!x_result)
        return std::unexpected(x_result.error());
      point.x = x_result.value();

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
      // Parse catalyticity values
      if (!node["catalyticity_values"]) {
        return std::unexpected(
            core::ConfigurationError("Missing required 'catalyticity_values' when abacus is enabled"));
      }
      config.catalyticity_values = node["catalyticity_values"].as<std::vector<double>>();

      if (config.catalyticity_values.empty()) {
        return std::unexpected(core::ConfigurationError("catalyticity_values cannot be empty"));
      }

      // Parse temperature range
      if (node["temperature_range"]) {
        auto range = node["temperature_range"].as<std::vector<double>>();
        if (range.size() != constants::indexing::min_range_size) {
          return std::unexpected(core::ConfigurationError(std::format("temperature_range must have exactly {} values", constants::indexing::min_range_size)));
        }
        config.temperature_min = range[constants::indexing::first];
        config.temperature_max = range[constants::indexing::second];

        if (config.temperature_min >= config.temperature_max) {
          return std::unexpected(core::ConfigurationError("temperature_min must be less than temperature_max"));
        }
      }

      // Parse temperature points
      if (node["temperature_points"]) {
        config.temperature_points = node["temperature_points"].as<int>();
        if (config.temperature_points < constants::indexing::min_range_size) {
          return std::unexpected(core::ConfigurationError(std::format("temperature_points must be at least {}", constants::indexing::min_range_size)));
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

} // namespace blast::io
