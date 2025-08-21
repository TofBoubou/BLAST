#pragma once
#include "../io/config_types.hpp"
#include "../thermophysics/mixture_interface.hpp"
#include "application_types.hpp"
#include <expected>
#include <memory>
#include <string>

namespace blast::core {

class ConfigurationLoader {
public:
  struct LoadResult {
    io::Configuration config;
    std::unique_ptr<thermophysics::MixtureInterface> mixture;
  };

  // Load configuration and create mixture
  [[nodiscard]] auto load_configuration(const std::string& config_file) 
    -> std::expected<LoadResult, ApplicationError>;

  // Display configuration information
  auto display_configuration_info(const io::Configuration& config, 
                                  const thermophysics::MixtureInterface& mixture) const -> void;

private:
  // Load YAML configuration
  [[nodiscard]] auto load_config_file(const std::string& config_file) 
    -> std::expected<io::Configuration, ApplicationError>;
    
  // Create thermophysical mixture
  [[nodiscard]] auto create_mixture(const io::MixtureConfig& mixture_config) 
    -> std::expected<std::unique_ptr<thermophysics::MixtureInterface>, ApplicationError>;
    
  // Display mixture information
  auto display_mixture_info(const thermophysics::MixtureInterface& mixture) const -> void;
  
  // Display simulation setup
  auto display_simulation_info(const io::Configuration& config) const -> void;
  
  // Display edge conditions
  auto display_edge_conditions(const io::Configuration& config) const -> void;
};

} // namespace blast::core