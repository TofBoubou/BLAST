#include "blast/core/configuration_loader.hpp"
#include "blast/core/constants.hpp"
#include "blast/io/config_manager.hpp"
#include "blast/thermophysics/mixture_interface.hpp"
#include <format>
#include <iomanip>
#include <iostream>

namespace blast::core {

auto ConfigurationLoader::load_configuration(const std::string& config_file) 
  -> std::expected<LoadResult, ApplicationError> {
    
  std::cout << "Loading configuration from: " << config_file << std::endl;
  
  auto config_result = load_config_file(config_file);
  if (!config_result) {
    return std::unexpected(config_result.error());
  }
  auto config = std::move(config_result.value());
  
  std::cout << "✓ Configuration loaded successfully" << std::endl;
  
  auto mixture_result = create_mixture(config.mixture);
  if (!mixture_result) {
    return std::unexpected(mixture_result.error());
  }
  auto mixture = std::move(mixture_result.value());
  
  std::cout << "✓ Mixture created (" << mixture->n_species() << " species)" << std::endl;
  
  display_configuration_info(config, *mixture);
  
  return LoadResult{std::move(config), std::move(mixture)};
}

auto ConfigurationLoader::display_configuration_info(const io::Configuration& config, 
                                                     const thermophysics::MixtureInterface& mixture) const -> void {
  display_mixture_info(mixture);
  display_simulation_info(config);
  display_edge_conditions(config);
}

auto ConfigurationLoader::load_config_file(const std::string& config_file) 
  -> std::expected<io::Configuration, ApplicationError> {
    
  io::ConfigurationManager config_manager;
  auto config_result = config_manager.load(config_file);
  
  if (!config_result) {
    return std::unexpected(ApplicationError{
      "Failed to load config: " + config_result.error().message(),
      constants::indexing::second
    });
  }
  
  return std::move(config_result.value());
}

auto ConfigurationLoader::create_mixture(const io::MixtureConfig& mixture_config) 
  -> std::expected<std::unique_ptr<thermophysics::MixtureInterface>, ApplicationError> {
    
  std::cout << "Creating mixture: " << mixture_config.name << std::endl;
  
  auto mixture_result = thermophysics::create_mixture(mixture_config);
  if (!mixture_result) {
    return std::unexpected(ApplicationError{
      "Failed to create mixture: " + mixture_result.error().message(),
      constants::indexing::second
    });
  }
  
  return std::move(mixture_result.value());
}

auto ConfigurationLoader::display_mixture_info(const thermophysics::MixtureInterface& mixture) const -> void {
  std::cout << "\nMixture species:" << std::endl;
  for (std::size_t i = constants::indexing::first; i < mixture.n_species(); ++i) {
    std::cout << std::format("  [{:2}] {:>8} (MW: {:8.3f} kg/kmol)", 
                             i, mixture.species_name(i), mixture.species_molecular_weight(i))
              << std::endl;
  }
}

auto ConfigurationLoader::display_simulation_info(const io::Configuration& config) const -> void {
  std::cout << "\nSimulation setup:" << std::endl;
  std::cout << "  Body type: "
            << (config.simulation.body_type == io::SimulationConfig::BodyType::Axisymmetric ? "Axisymmetric" : "Other")
            << std::endl;
  std::cout << "  Stagnation only: " << (config.simulation.only_stagnation_point ? "Yes" : "No") << std::endl;
  std::cout << "  Chemical mode: ";
  
  switch (config.simulation.chemical_mode) {
  case io::SimulationConfig::ChemicalMode::Equilibrium:
    std::cout << "Equilibrium";
    break;
  case io::SimulationConfig::ChemicalMode::Frozen:
    std::cout << "Frozen";
    break;
  case io::SimulationConfig::ChemicalMode::NonEquilibrium:
    std::cout << "Non-equilibrium";
    break;
  }
  
  std::cout << std::endl;
  std::cout << "  Thermal diffusion: " << (config.simulation.consider_thermal_diffusion ? "Yes" : "No") << std::endl;
  std::cout << "  Grid points (η): " << config.numerical.n_eta << std::endl;
  std::cout << "  η_max: " << config.numerical.eta_max << std::endl;
  std::cout << "  Convergence tol: " << config.numerical.convergence_tolerance << std::endl;
}

auto ConfigurationLoader::display_edge_conditions(const io::Configuration& config) const -> void {
  if (!config.outer_edge.edge_points.empty()) {
    auto edge = config.outer_edge.edge_points[constants::indexing::first];
    std::cout << "\nEdge conditions:" << std::endl;
    std::cout << "  Pressure: " << edge.pressure << " Pa" << std::endl;
    std::cout << "  Temperature: " << edge.temperature << " K" << std::endl;
    
    if (!config.wall_parameters.wall_temperatures.empty()) {
      std::cout << "  Wall temp: " << config.wall_parameters.wall_temperatures[constants::indexing::first] << " K" << std::endl;
    } else {
      std::cout << "  Wall temp: NOT SET" << std::endl;
    }
  }
}

} // namespace blast::core