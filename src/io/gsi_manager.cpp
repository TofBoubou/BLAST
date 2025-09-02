#include "blast/io/gsi_manager.hpp"
#include <cstdlib>
#include <fstream>
#include <format>
#include <iostream>
#include <regex>
#include <sstream>

namespace blast::io {

GsiManager::GsiManager(const Configuration& config) : config_(config) {
  gsi_file_path_ = construct_gsi_path();
}

GsiManager::~GsiManager() {
  // Always restore original GSI on destruction if it was backed up
  if (gsi_backed_up_) {
    restore_gsi_file();
  }
}

auto GsiManager::backup_gsi_file() -> std::expected<void, std::string> {
  if (!std::filesystem::exists(gsi_file_path_)) {
    return std::unexpected(std::format("GSI file not found: {}", gsi_file_path_.string()));
  }

  std::ifstream file(gsi_file_path_);
  if (!file.is_open()) {
    return std::unexpected(std::format("Failed to open GSI file: {}", gsi_file_path_.string()));
  }

  std::stringstream buffer;
  buffer << file.rdbuf();
  original_gsi_content_ = buffer.str();
  file.close();

  gsi_backed_up_ = true;
  return {};
}

auto GsiManager::restore_gsi_file() -> void {
  if (!gsi_backed_up_ || original_gsi_content_.empty()) {
    return;
  }

  std::ofstream file(gsi_file_path_);
  if (file.is_open()) {
    file << original_gsi_content_;
    file.close();
    gsi_backed_up_ = false;  // Mark as restored to avoid double restoration
  }
}

auto GsiManager::update_gsi_catalyticity(double gamma) -> std::expected<void, std::string> {
  if (!gsi_backed_up_) {
    return std::unexpected("GSI file not backed up. Call backup_gsi_file() first.");
  }

  std::string content = original_gsi_content_;

  // Regex to match gamma_const elements
  std::regex gamma_pattern(R"(<gamma_const>\s*([^<]+)\s*</gamma_const>)");
  std::string result_content;

  auto current_pos = content.cbegin();
  auto content_end = content.cend();

  std::smatch match;
  while (std::regex_search(current_pos, content_end, match, gamma_pattern)) {
    // Append everything before the match
    result_content.append(current_pos, current_pos + match.position());

    // Parse the species list and replace gamma values
    std::string species_list = match[1].str();
    std::stringstream new_gamma_content;
    new_gamma_content << "<gamma_const> ";

    // Match individual species:value pairs
    std::regex species_pattern(R"((\w+):[\d.]+(?:\s+|$))");
    std::sregex_iterator species_begin(species_list.begin(), species_list.end(), species_pattern);
    std::sregex_iterator species_end;

    bool first = true;
    for (auto it = species_begin; it != species_end; ++it) {
      if (!first)
        new_gamma_content << " ";

      std::string species_name = (*it)[1].str();
      new_gamma_content << species_name << ":" << gamma;
      first = false;
    }

    new_gamma_content << " </gamma_const>";
    result_content.append(new_gamma_content.str());

    // Move past this match
    current_pos = current_pos + match.position() + match.length();
  }

  // Append any remaining content
  result_content.append(current_pos, content_end);

  // Write updated content to file
  std::ofstream file(gsi_file_path_);
  if (!file.is_open()) {
    return std::unexpected(std::format("Failed to open GSI file for writing: {}", gsi_file_path_.string()));
  }

  file << result_content;
  file.close();

  return {};
}

auto GsiManager::construct_gsi_path() const -> std::filesystem::path {
  // Get MPP_DATA_DIRECTORY from environment
  const char* mpp_data = std::getenv("MPP_DATA_DIRECTORY");

  if (!mpp_data) {
    // Try to find it relative to the executable if available (Linux), non-throwing
    std::error_code ec;
    std::filesystem::path exe_path = std::filesystem::read_symlink("/proc/self/exe", ec);
    if (!ec && !exe_path.empty()) {
      std::filesystem::path base_path = exe_path.parent_path();
      std::filesystem::path mpp_path = base_path / "libs" / "mutationpp" / "data";
      if (std::filesystem::exists(mpp_path)) {
        return mpp_path / "gsi" / (config_.mixture.name + "_cata.xml");
      }
    }

    // Fallback to current working directory (cross-platform)
    return std::filesystem::current_path() / "libs" / "mutationpp" / "data" / "gsi" /
           (config_.mixture.name + "_cata.xml");
  }

  // Use MPP_DATA_DIRECTORY
  std::filesystem::path gsi_path = std::filesystem::path(mpp_data) / "gsi";
  std::string gsi_filename = config_.mixture.name + "_cata.xml";

  return gsi_path / gsi_filename;
}

} // namespace blast::io
