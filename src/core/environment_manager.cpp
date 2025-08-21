#include "blast/core/environment_manager.hpp"
#include "blast/core/constants.hpp"
#include <cstdlib>
#include <iostream>

namespace blast::core {

auto EnvironmentManager::configure_environment(const std::filesystem::path& executable_path) 
  -> std::expected<void, ApplicationError> {
    
  if (auto result = configure_mutation_pp(executable_path); !result) {
    return std::unexpected(result.error());
  }
  
  return validate_environment();
}

auto EnvironmentManager::configure_mutation_pp(const std::filesystem::path& exe_path) 
  -> std::expected<void, ApplicationError> {
    
  try {
    std::filesystem::path mpp_data_path = exe_path / "libs" / "mutationpp" / "data";
    
    if (std::filesystem::exists(mpp_data_path)) {
      std::string mpp_data_str = std::filesystem::canonical(mpp_data_path).string();
      setenv("MPP_DATA_DIRECTORY", mpp_data_str.c_str(), constants::indexing::second);
      std::cout << "MPP_DATA_DIRECTORY auto-set to: " << mpp_data_str << std::endl;
    } else {
      std::cerr << "Warning: Mutation++ data directory not found at: " << mpp_data_path << std::endl;
      std::cerr << "Looking for MPP_DATA_DIRECTORY in environment..." << std::endl;
      
      if (const char* env_mpp = std::getenv("MPP_DATA_DIRECTORY")) {
        std::cout << "Using existing MPP_DATA_DIRECTORY: " << env_mpp << std::endl;
      } else {
        std::cerr << "Warning: MPP_DATA_DIRECTORY not set. Mutation++ may fail." << std::endl;
      }
    }
    
    return {};
    
  } catch (const std::exception& e) {
    std::cerr << "Warning: Could not auto-configure MPP_DATA_DIRECTORY: " << e.what() << std::endl;
    return {}; // Not fatal, continue execution
  }
}

auto EnvironmentManager::validate_environment() const 
  -> std::expected<void, ApplicationError> {
  // Could add more environment validation here if needed
  return {};
}

} // namespace blast::core