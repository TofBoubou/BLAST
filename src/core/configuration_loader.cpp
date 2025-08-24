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
  std::cout << "\n" << constants::string_processing::colors::cyan 
            << "┌─ SIMULATION SETUP ────────────────────────┐" 
            << constants::string_processing::colors::reset << std::endl;
  std::cout << "│ Body type       : " << std::setw(20) << std::left
            << (config.simulation.body_type == io::SimulationConfig::BodyType::Axisymmetric ? "Axisymmetric" : "Other")
            << " │" << std::endl;
  std::cout << "│ Stagnation only : " << std::setw(20) << std::left 
            << (config.simulation.only_stagnation_point ? "Yes" : "No") << " │" << std::endl;
  
  std::cout << "│ Chemical mode   : " << std::setw(20) << std::left;
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
  std::cout << " │" << std::endl;
  
  std::cout << "│ Thermal diffusion: " << std::setw(20) << std::left 
            << (config.simulation.consider_thermal_diffusion ? "Yes" : "No") << " │" << std::endl;
  std::cout << "│ Dufour effect   : " << std::setw(20) << std::left 
            << (config.simulation.consider_dufour_effect ? "Yes" : "No") << " │" << std::endl;
  std::cout << "│ Finite thickness: " << std::setw(20) << std::left 
            << (config.simulation.finite_thickness ? "Yes" : "No") << " │" << std::endl;
  std::cout << "│ Catalytic wall  : " << std::setw(20) << std::left 
            << (config.simulation.catalytic_wall ? "Yes" : "No") << " │" << std::endl;
  std::cout << "│ Grid points (η) : " << std::setw(20) << std::left 
            << config.numerical.n_eta << " │" << std::endl;
  std::cout << "│ η_max           : " << std::setw(20) << std::left 
            << config.numerical.eta_max << " │" << std::endl;
  std::cout << "│ Convergence tol : " << std::setw(20) << std::left 
            << config.numerical.convergence_tolerance << " │" << std::endl;
  std::cout << constants::string_processing::colors::cyan 
            << "└───────────────────────────────────────────┘" 
            << constants::string_processing::colors::reset << std::endl;
}

auto ConfigurationLoader::display_edge_conditions(const io::Configuration& config) const -> void {
  if (!config.outer_edge.edge_points.empty()) {
    auto edge = config.outer_edge.edge_points[constants::indexing::first];
    
    std::cout << "\n" << constants::string_processing::colors::blue 
              << "┌─ EDGE CONDITIONS ─────────────────────────┐" 
              << constants::string_processing::colors::reset << std::endl;
    std::cout << "│ Pressure        : " << std::setw(15) << std::right 
              << edge.pressure << " Pa   │" << std::endl;
    std::cout << "│ Temperature     : " << std::setw(15) << std::right 
              << std::fixed << std::setprecision(1) << edge.temperature << " K    │" << std::endl;
    
    if (!config.wall_parameters.wall_temperatures.empty()) {
      std::cout << "│ Wall temp       : " << std::setw(15) << std::right 
                << std::fixed << std::setprecision(0) << config.wall_parameters.wall_temperatures[constants::indexing::first] 
                << " K    │" << std::endl;
    } else {
      std::cout << "│ Wall temp       : " << std::setw(15) << std::right 
                << "NOT SET" << "      │" << std::endl;
    }
    
    std::cout << "│ Velocity        : " << std::setw(15) << std::right 
              << edge.velocity << " m/s  │" << std::endl;
    std::cout << "│ Radius          : " << std::setw(15) << std::right 
              << edge.radius << " m    │" << std::endl;
    std::cout << constants::string_processing::colors::blue 
              << "└───────────────────────────────────────────┘" 
              << constants::string_processing::colors::reset << std::endl;
  }
}

} // namespace blast::core