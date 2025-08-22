#pragma once

#include "config_types.hpp"
#include <expected>
#include <filesystem>
#include <string>

namespace blast::io {

/**
 * @brief Utility class for managing GSI file modifications
 * 
 * This class handles backing up, modifying, and restoring GSI files,
 * particularly for updating catalyticity values. It's used by both
 * AbacusGenerator and EdgeTemperatureReconstructor.
 */
class GsiManager {
public:
  explicit GsiManager(const Configuration& config);
  ~GsiManager();

  // Non-copyable but movable
  GsiManager(const GsiManager&) = delete;
  GsiManager& operator=(const GsiManager&) = delete;
  GsiManager(GsiManager&&) = default;
  GsiManager& operator=(GsiManager&&) = default;

  /**
   * @brief Backup the original GSI file content
   */
  [[nodiscard]] auto backup_gsi_file() -> std::expected<void, std::string>;

  /**
   * @brief Update GSI file with new catalyticity value
   * @param gamma The catalyticity value to set
   */
  [[nodiscard]] auto update_gsi_catalyticity(double gamma) -> std::expected<void, std::string>;

  /**
   * @brief Restore the original GSI file content
   */
  auto restore_gsi_file() -> void;

  /**
   * @brief Get the path to the GSI file
   */
  [[nodiscard]] auto get_gsi_path() const -> const std::filesystem::path& { return gsi_file_path_; }

private:
  /**
   * @brief Construct the GSI file path from mixture name and environment
   */
  [[nodiscard]] auto construct_gsi_path() const -> std::filesystem::path;

  const Configuration& config_;
  std::filesystem::path gsi_file_path_;
  std::string original_gsi_content_;
  bool gsi_backed_up_ = false;
};

} // namespace blast::io