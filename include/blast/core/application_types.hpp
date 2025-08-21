#pragma once
#include <expected>
#include <filesystem>
#include <string>
#include <vector>

namespace blast::core {

// Application-level error handling
struct ApplicationError {
  std::string message;
  int exit_code;
};

// Command line arguments structure
struct CommandLineArgs {
  std::string config_file;
  std::string output_name = "simulation";
  bool help_requested = false;
};

// Application result for clean exit handling
struct ApplicationResult {
  bool success;
  int exit_code;
  std::string message;
};

// Performance metrics for reporting
struct PerformanceMetrics {
  std::chrono::milliseconds total_time{0};
  std::chrono::milliseconds solve_time{0};
  std::chrono::milliseconds output_time{0};
  std::vector<std::filesystem::path> output_files;
};

} // namespace blast::core