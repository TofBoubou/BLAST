#include "blast/io/yaml_parser.hpp"
#include "blast/core/exceptions.hpp"

#include <fstream>
#include <filesystem>

namespace blast::io {

auto YamlParser::load() -> std::expected<void, core::FileError> {
    try {
        root_ = YAML::LoadFile(file_path_);  
        return {}; 
    } catch (const YAML::BadFile& e) {
        return std::unexpected(core::FileError{"Failed to open YAML file", file_path_});
    } catch (const YAML::ParserException& e) {
        return std::unexpected(core::FileError{std::format("YAML parsing error: {}", e.what()), file_path_});
    } catch (const std::exception& e) { //YAML::BadFile inherit from std::exception so no problem
        return std::unexpected(core::FileError{std::format("Unexpected error during YAML load: {}", e.what()), file_path_});
    }
}


auto YamlParser::parse() const -> std::expected<Configuration, core::ConfigurationError> {
    try {
        if (!root_ || root_.IsNull()) { // !root_ is override in YAML so we check if it's empty
            return std::unexpected(
                core::ConfigurationError("No YAML content loaded. Call load() first.")
            );
        }
        
        Configuration config;
        
        // Simulation
        if (!root_["simulation"]) { // !root_ is override in YAML so we check if it's empty
            return std::unexpected(
                core::ConfigurationError("Missing required 'simulation' section")
            );
        }
        auto sim_result = parse_simulation_config(root_["simulation"]);
        if (!sim_result) {
            return std::unexpected(sim_result.error());
        }
        config.simulation = std::move(sim_result.value()); // emptying sim_result, avoid a copy 
        
        // Numerical parameters
        if (!root_["numerical"]) {
            return std::unexpected(
                core::ConfigurationError("Missing required 'numerical' section")
            );
        }
        auto num_result = parse_numerical_config(root_["numerical"]); // The root changes
        if (!num_result) {
            return std::unexpected(num_result.error());
        }
        config.numerical = std::move(num_result.value());
        
        // Mixture
        if (!root_["mixture"]) {
            return std::unexpected(
                core::ConfigurationError("Missing required 'mixture' section")
            );
        }
        auto mix_result = parse_mixture_config(root_["mixture"]);
        if (!mix_result) {
            return std::unexpected(mix_result.error());
        }
        config.mixture = std::move(mix_result.value());
        
        //  Output
        if (root_["output"]) {
            auto out_result = parse_output_config(root_["output"]);
            if (!out_result) {
                return std::unexpected(out_result.error());
            }
            config.output = std::move(out_result.value());
        }
        
        // outer_edge
        if (!root_["outer_edge"]) {
            return std::unexpected(
                core::ConfigurationError("Missing required 'outer_edge' section")
            );
        }
        auto edge_result = parse_outer_edge_config(root_["outer_edge"]);
        if (!edge_result) {
            return std::unexpected(edge_result.error());
        }
        config.outer_edge = std::move(edge_result.value());
        
        // 6. Parser la section wall_parameters (obligatoire)
        if (!root_["wall_parameters"]) {
            return std::unexpected(
                core::ConfigurationError("Missing required 'wall_parameters' section")
            );
        }
        auto wall_result = parse_wall_parameters_config(root_["wall_parameters"]);
        if (!wall_result) {
            return std::unexpected(wall_result.error());
        }
        config.wall_parameters = std::move(wall_result.value());
        
        // 7. Parser la section initial_guess (optionnelle)
        if (root_["initial_guess"]) {
            auto guess_result = parse_initial_guess_config(root_["initial_guess"]);
            if (!guess_result) {
                return std::unexpected(guess_result.error());
            }
            config.initial_guess = std::move(guess_result.value());
        }
        // Sinon, use_initial_guess reste à false par défaut
        
        return config;
        
    } catch (const YAML::Exception& e) {
        return std::unexpected(
            core::ConfigurationError(
                std::format("YAML parsing error: {} at line {}", 
                    e.what(), e.mark.line)
            )
        );
    } catch (const std::exception& e) {
        return std::unexpected(
            core::ConfigurationError(
                std::format("Unexpected error during parsing: {}", e.what())
            )
        );
    }
}


auto YamlParser::parse_simulation_config(const YAML::Node& node) const 
    -> std::expected<SimulationConfig, core::ConfigurationError> {
    
    SimulationConfig config;
    
    try {
        config.body_type = extract_enum(
            node, "body_type", 
            enum_mappings::body_types,
            SimulationConfig::BodyType::Axisymmetric
        );
        
        config.only_stagnation_point = extract_value<bool>(
            node, "only_stagnation_point", true
        );
        
        config.diffusion_type = extract_enum(
            node, "diffusion_type",
            enum_mappings::diffusion_types,
            SimulationConfig::DiffusionType::StefanMaxwell
        );
        
        config.consider_thermal_diffusion = extract_value<bool>(
            node, "consider_thermal_diffusion", false
        );
        
        config.chemical_non_equilibrium = extract_value<bool>(
            node, "chemical_non_equilibrium", true
        );
        
        config.catalysis = extract_value<bool>(
            node, "catalysis", false
        );
        
        return config;
        
    } catch (const core::ConfigurationError& e) {
        return std::unexpected(
            core::ConfigurationError(
                std::format("In 'simulation' section: {}", e.message())
            )
        );
    }
}


auto YamlParser::parse_numerical_config(const YAML::Node& node) const 
    -> std::expected<NumericalConfig, core::ConfigurationError> {
    
    NumericalConfig config;
    
    try {
        config.n_eta = extract_value<int>(node, "n_eta", 101);
        config.eta_max = extract_value<double>(node, "eta_max", 8.0);
        config.convergence_tolerance = extract_value<double>(
            node, "convergence_tolerance", 1e-6
        );
        config.max_iterations = extract_value<int>(node, "max_iterations", 1000);
        config.under_relaxation = extract_value<double>(
            node, "under_relaxation", 0.5
        );
        
        if (node["step_control"]) {
            auto sc_node = node["step_control"];
            config.step_control.lower_bound = extract_value<int>(
                sc_node, "lower_bound", 10
            );
            config.step_control.upper_bound = extract_value<int>(
                sc_node, "upper_bound", 50
            );
        }
        
        if (node["solvers"]) {
            auto solv_node = node["solvers"];
            config.solvers.h2t_tolerance = extract_value<double>(
                solv_node, "h2t_tolerance", 1e-8
            );
            config.solvers.h2t_max_iterations = extract_value<int>(
                solv_node, "h2t_max_iterations", 100
            );
            config.solvers.stefan_tolerance = extract_value<double>(
                solv_node, "stefan_tolerance", 1e-6
            );
            config.solvers.stefan_max_iterations = extract_value<int>(
                solv_node, "stefan_max_iterations", 50
            );
        }
        
        return config;
        
    } catch (const core::ConfigurationError& e) {
        return std::unexpected(
            core::ConfigurationError(
                std::format("In 'numerical' section: {}", e.message())
            )
        );
    }
}


auto YamlParser::parse_mixture_config(const YAML::Node& node) const 
    -> std::expected<MixtureConfig, core::ConfigurationError> {
    
    MixtureConfig config;
    
    try {
        config.name = extract_value<std::string>(node, "name", "air7");
        
        config.thermodynamic_database = extract_enum(
            node, "thermodynamic_database",
            enum_mappings::databases,
            MixtureConfig::Database::NASA9
        );
        
        config.viscosity_algorithm = extract_enum(
            node, "viscosity_algorithm",
            enum_mappings::viscosity_algorithms,
            MixtureConfig::ViscosityAlgorithm::ChapmanEnskog_CG
        );
        
        config.reference_temperature = extract_value<double>(
            node, "reference_temperature", 0.0
        );
        
        return config;
        
    } catch (const core::ConfigurationError& e) {
        return std::unexpected(
            core::ConfigurationError(
                std::format("In 'mixture' section: {}", e.message())
            )
        );
    }
}


auto YamlParser::parse_output_config(const YAML::Node& node) const 
    -> std::expected<OutputConfig, core::ConfigurationError> {
    
    OutputConfig config;
    
    try {
        config.x_stations = extract_value<std::vector<double>>(
            node, "x_stations", std::vector<double>{}
        );
        
        config.output_directory = extract_value<std::string>(
            node, "output_directory", "BLAST_outputs"
        );
        
        config.generate_lookup_table = extract_value<bool>(
            node, "generate_lookup_table", false
        );
        
        // lookup_table
        if (node["lookup_table"]) {
            auto lt_node = node["lookup_table"];
            
            config.lookup_table.temperature_min = extract_value<double>(
                lt_node, "temperature_min", 300.0
            );
            config.lookup_table.temperature_max = extract_value<double>(
                lt_node, "temperature_max", 3000.0
            );
            config.lookup_table.temperature_step = extract_value<double>(
                lt_node, "temperature_step", 100.0
            );
            config.lookup_table.gamma_min = extract_value<double>(
                lt_node, "gamma_min", 0.0
            );
            config.lookup_table.gamma_max = extract_value<double>(
                lt_node, "gamma_max", 1.0
            );
            config.lookup_table.gamma_step = extract_value<double>(
                lt_node, "gamma_step", 0.1
            );
        }
        
        return config;
        
    } catch (const core::ConfigurationError& e) {
        return std::unexpected(
            core::ConfigurationError(
                std::format("In 'output' section: {}", e.message())
            )
        );
    }
}


auto YamlParser::parse_initial_guess_config(const YAML::Node& node) const 
    -> std::expected<InitialGuessConfig, core::ConfigurationError> {
    
    InitialGuessConfig config;
    
    try {
        config.use_initial_guess = extract_value<bool>(
            node, "use_initial_guess", false
        );
        
        if (config.use_initial_guess) {
            // Parser les profils seulement si on utilise l'initial guess
            if (node["temperature_profile"]) {
                config.temperature_profile = extract_value<std::vector<double>>(
                    node, "temperature_profile", std::vector<double>{}
                );
            }
            
            if (node["enthalpy_profile"]) {
                config.enthalpy_profile = extract_value<std::vector<double>>(
                    node, "enthalpy_profile", std::vector<double>{}
                );
            }
            
            if (node["species_profiles"]) {
                auto sp_node = node["species_profiles"];
                std::vector<std::vector<double>> profiles;
                
                for (const auto& profile_node : sp_node) {
                    profiles.push_back(profile_node.as<std::vector<double>>());
                }
                
                config.species_profiles = profiles;
            }
        }
        
        return config;
        
    } catch (const core::ConfigurationError& e) {
        return std::unexpected(
            core::ConfigurationError(
                std::format("In 'initial_guess' section: {}", e.message())
            )
        );
    }
}


auto YamlParser::parse_outer_edge_config(const YAML::Node& node) const 
    -> std::expected<OuterEdgeConfig, core::ConfigurationError> {
    
    OuterEdgeConfig config;
    
    try {
        if (!node["edge_points"]) {
            return std::unexpected(
                core::ConfigurationError("Missing required 'edge_points' in outer_edge")
            );
        }
        
        auto points_node = node["edge_points"];
        for (const auto& point_node : points_node) {
            OuterEdgeConfig::EdgePoint point;
            
            point.x = extract_value<double>(point_node, "x", 0.0);
            point.radius = extract_value<double>(point_node, "radius", 1.0);
            point.velocity = extract_value<double>(point_node, "velocity", 0.0);
            point.enthalpy = extract_value<double>(point_node, "enthalpy", 0.0);
            point.pressure = extract_value<double>(point_node, "pressure", 101325.0);
            point.density = extract_value<double>(point_node, "density", 1.225);
            point.viscosity = extract_value<double>(point_node, "viscosity", 1.8e-5);
            
            if (point_node["species"]) {
                point.species_fractions = extract_value<std::vector<double>>(
                    point_node, "species", std::vector<double>{}
                );
            }
            
            config.edge_points.push_back(point);
        }
        
        // Scalars
        config.velocity_gradient_stagnation = extract_value<double>(
            node, "velocity_gradient_stagnation", 0.0
        );
        config.freestream_density = extract_value<double>(
            node, "freestream_density", 1.225
        );
        config.freestream_velocity = extract_value<double>(
            node, "freestream_velocity", 0.0
        );

        return config;
        
    } catch (const core::ConfigurationError& e) {
        return std::unexpected(
            core::ConfigurationError(
                std::format("In 'outer_edge' section: {}", e.message())
            )
        );
    }
}


auto YamlParser::parse_wall_parameters_config(const YAML::Node& node) const 
    -> std::expected<WallParametersConfig, core::ConfigurationError> {
    
    WallParametersConfig config;
    
    try {
        // Parser les températures de paroi (requis)
        config.wall_temperatures = extract_value<std::vector<double>>(
            node, "temperatures", std::vector<double>{}
        );
        
        if (config.wall_temperatures.empty()) {
            return std::unexpected(
                core::ConfigurationError(
                    "Wall temperatures cannot be empty"
                )
            );
        }
        return config;
        
    } catch (const core::ConfigurationError& e) {
        return std::unexpected(
            core::ConfigurationError(
                std::format("In 'wall_parameters' section: {}", e.message())
            )
        );
    }
}

} // namespace blast::io
