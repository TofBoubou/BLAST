#pragma once
#include "application_types.hpp"
#include <expected>
#include <filesystem>
#include <string>

namespace blast::core {

class EnvironmentManager {
public:
  // Configure environment for the application
  [[nodiscard]] auto configure_environment(const std::filesystem::path& executable_path) 
    -> std::expected<void, ApplicationError>;

private:
  // Configure Mutation++ data directory
  [[nodiscard]] auto configure_mutation_pp(const std::filesystem::path& exe_path) 
    -> std::expected<void, ApplicationError>;
    
  // Validate environment setup
  [[nodiscard]] auto validate_environment() const 
    -> std::expected<void, ApplicationError>;
};

} // namespace blast::core