#pragma once
#include "../core/exceptions.hpp"
#include "config_types.hpp"
#include "yaml_parser.hpp"
#include <expected>
#include <filesystem>
#include <memory>

namespace blast::io {

class ConfigurationManager {
private:
  std::unique_ptr<YamlParser> parser_;
  std::optional<Configuration> current_config_;
  std::filesystem::path config_file_path_;

  [[nodiscard]] auto
  resolve_config_path(std::string_view config_file) const -> std::expected<std::filesystem::path, core::FileError>;

public:
  explicit ConfigurationManager() = default;

  [[nodiscard]] auto load(std::string_view config_file) -> std::expected<Configuration, core::ConfigurationError>;
};
} // namespace blast::io