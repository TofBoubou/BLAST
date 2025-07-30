#include "blast/io/config_manager.hpp"
#include <iostream>

namespace blast::io {

auto ConfigurationManager::resolve_config_path(std::string_view config_file) const
    -> std::expected<std::filesystem::path, core::FileError> {
  std::filesystem::path path_candidate(config_file); // convert, more intelligent variable

  if (!std::filesystem::exists(path_candidate)) {
    return std::unexpected(core::FileError{"could not locate config file", std::string(config_file)});
  }
  return std::filesystem::absolute(path_candidate);
}

auto ConfigurationManager::load(std::string_view config_file)
    -> std::expected<Configuration, core::ConfigurationError> {

  auto path_result = resolve_config_path(config_file);
  if (!path_result) {
    return std::unexpected(
        core::ConfigurationError(std::format("Failed to resolve config path: {}",
                                             path_result.error().message())) // .error() to go to the exception
    );
  }

  config_file_path_ = path_result.value(); // extract the value of std::expected

  parser_ = std::make_unique<YamlParser>(config_file_path_.string()); // free the memory automatically, use the
                                                                      // contructor of YamlParser with the arguments
                                                                      // after, type std::unique_ptr<YamlParser>

  auto load_result = parser_->load(); // normally include nothing because it in
                                      // an expected void, error. Fill the root
  if (!load_result) {
    return std::unexpected(
        core::ConfigurationError(std::format("Failed to load YAML file: {}", load_result.error().message())));
  }

  auto parse_result = parser_->parse();
  if (!parse_result) {
    return std::unexpected(parse_result.error());
  }

  current_config_ = std::move(parse_result.value());

  return *current_config_;
}

} // namespace blast::io
