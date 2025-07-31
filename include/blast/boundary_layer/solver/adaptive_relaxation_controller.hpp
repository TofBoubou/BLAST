#pragma once
#include "../../core/exceptions.hpp"
#include "../equations/equation_types.hpp"
#include <algorithm>
#include <cmath>
#include <expected>
#include <numeric>
#include <vector>

namespace blast::boundary_layer::solver {

// Forward declarations
class SolverError;

// ConvergenceInfo definition
struct ConvergenceInfo {
  double residual_F = 1e10;
  double residual_g = 1e10;
  double residual_c = 1e10;
  int iterations = 0;
  bool converged = false;

  [[nodiscard]] constexpr auto max_residual() const noexcept -> double {
    return std::max({residual_F, residual_g, residual_c});
  }
};

// Structure to handle adaptive relaxation
class AdaptiveRelaxationController {
public:
  // Configuration parameters
  struct Config {
    double initial_factor = 0.01; // Conservative initial factor
    double min_factor = 0.00001;  // Absolute minimum
    double max_factor = 1;        // Reasonable maximum

    // Thresholds based on your logic
    double excellent_threshold = 1;  // Threshold for excellent convergence (< 0.8)
    double divergence_threshold = 1; // Threshold for divergence (> 1.0)

    // Adaptation factors
    double strong_increase = 1.0;    // Strong increase for ratio < 0.8
    double moderate_increase = 1.0;  // Moderate increase for [0.8, 1.0]
    double decrease_factor = 1;     // Decrease for ratio > X
    double oscillation_penalty = 1; // Penalty factor for oscillations

    // Specialized configurations
    [[nodiscard]] static auto for_stagnation_point() -> Config {
      Config config;
      config.initial_factor = 0.001;
      config.max_factor = 1;           // More conservative
      config.strong_increase = 1.0;   // Slower growth
      config.moderate_increase = 1.0; // Very moderate
      return config;
    }

    [[nodiscard]] static auto for_downstream_station() -> Config {
      Config config;
      return config;
    }
  };

private:
  Config config_;

  // Internal state
  double current_factor_;
  double previous_residual_ = 1e10;
  int consecutive_improvements_ = 0;
  int consecutive_deteriorations_ = 0;
  bool first_iteration_ = true;

  // History for oscillation detection
  std::vector<double> residual_history_;
  static constexpr std::size_t history_size_ = 5;

public:
  AdaptiveRelaxationController() : config_(Config{}), current_factor_(config_.initial_factor) {
    residual_history_.reserve(history_size_);
  }

  explicit AdaptiveRelaxationController(const Config& config)
      : config_(config), current_factor_(config.initial_factor) {
    residual_history_.reserve(history_size_);
  }

  // Main method to adapt the factor
  [[nodiscard]] auto adapt_relaxation_factor(const ConvergenceInfo& conv_info, int iteration) -> double;

  // Accessors
  [[nodiscard]] auto current_factor() const noexcept -> double { return current_factor_; }

  [[nodiscard]] auto config() const noexcept -> const Config& { return config_; }

  // Reset for a new station
  auto reset_for_new_station() -> void {
    current_factor_ = config_.initial_factor;
    previous_residual_ = 1e10;
    consecutive_improvements_ = 0;
    consecutive_deteriorations_ = 0;
    first_iteration_ = true;
    residual_history_.clear();
  }

private:
  [[nodiscard]] auto detect_oscillations() const -> std::expected<bool, SolverError>;
  [[nodiscard]] auto compute_residual_trend() const -> std::expected<double, SolverError>;
};

// Factory functions for different problem types
namespace relaxation_factory {
[[nodiscard]] auto create_for_station(int station_number) -> AdaptiveRelaxationController;
[[nodiscard]] auto create_conservative() -> AdaptiveRelaxationController;
[[nodiscard]] auto create_aggressive() -> AdaptiveRelaxationController;
} // namespace relaxation_factory

} // namespace blast::boundary_layer::solver
