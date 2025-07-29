#pragma once
#include "../../core/containers.hpp"
#include "../../core/exceptions.hpp"
#include "../../io/config_types.hpp"
#include "coordinate_transform.hpp"
#include <concepts>
#include <expected>
#include <ranges>
#include <span>
#include <vector>

namespace blast::boundary_layer::grid {

class GridError : public core::BlastException {
public:
  explicit GridError(std::string_view message, std::source_location location = std::source_location::current()) noexcept
      : BlastException(std::format("Grid Error: {}", message), location) {}
};

template <typename T>
concept NumericRange = std::ranges::range<T> && std::floating_point<std::ranges::range_value_t<T>>;

template <typename T>
concept GridConfigType = requires(const T& config) {
  config.n_eta;
  config.eta_max;
};

// Compile-time constants
namespace constants {
constexpr double default_tolerance = 1e-6;
constexpr size_t default_grid_reserve = 100;
constexpr double step_reduction_factor = 0.5;
} // namespace constants

class BoundaryLayerGrid {
private:
  std::vector<double> xi_;
  std::vector<double> xi_output_;
  std::vector<double> eta_;
  double d_eta_;
  int n_eta_;
  double eta_max_;

public:
  template <GridConfigType NumericalConfig>
  [[nodiscard]] static constexpr auto
  create_stagnation_grid(NumericalConfig&& numerical_config,
                         const io::OuterEdgeConfig& edge_config) -> std::expected<BoundaryLayerGrid, GridError>;

  template <GridConfigType NumericalConfig>
  [[nodiscard]] static constexpr auto
  create_downstream_grid(NumericalConfig&& numerical_config, const io::OuterEdgeConfig& edge_config,
                         const io::OutputConfig& output_config) -> std::expected<BoundaryLayerGrid, GridError>;

  // Modern accessors with noexcept
  [[nodiscard]] constexpr auto xi_coordinates() const noexcept -> std::span<const double> { return xi_; }
  [[nodiscard]] constexpr auto xi_output_coordinates() const noexcept -> std::span<const double> { return xi_output_; }
  [[nodiscard]] constexpr auto eta_coordinates() const noexcept -> std::span<const double> { return eta_; }

  [[nodiscard]] constexpr auto n_eta() const noexcept -> int { return n_eta_; }
  [[nodiscard]] constexpr auto eta_max() const noexcept -> double { return eta_max_; }
  [[nodiscard]] constexpr auto d_eta() const noexcept -> double { return d_eta_; }

  // Grid operations
  [[nodiscard]] auto find_xi_interval(double xi_target) const noexcept -> std::expected<std::pair<int, int>, GridError>;

  template <NumericRange XGrid>
  [[nodiscard]] auto interpolate_x_from_xi(double xi_target, XGrid&& x_grid) const -> std::expected<double, GridError>;

private:
  constexpr BoundaryLayerGrid(int n_eta, double eta_max) noexcept
      : n_eta_(n_eta), eta_max_(eta_max), d_eta_(eta_max / static_cast<double>(n_eta - 1)) {

    xi_.reserve(constants::default_grid_reserve);
    xi_output_.reserve(constants::default_grid_reserve);
    eta_.reserve(n_eta_);
  }

  inline constexpr auto generate_eta_distribution() noexcept -> void;
  auto generate_xi_distribution(const io::OuterEdgeConfig& edge_config) -> std::expected<void, GridError>;
  auto generate_xi_output_distribution(const io::OutputConfig& output_config,
                                       std::span<const double> x_edge) -> std::expected<void, GridError>;
};

// Utility functions with constexpr where possible
[[nodiscard]] constexpr auto
compute_xi_step_size(double current_xi, double d_xi, int iterations,
                     const io::NumericalConfig::StepControl& step_control) noexcept -> double;

[[nodiscard]] constexpr auto reduce_step_size(double current_xi, double& d_xi) noexcept -> double;

template <GridConfigType NumericalConfig>
constexpr auto BoundaryLayerGrid::create_stagnation_grid(NumericalConfig&& numerical_config,
                                                         const io::OuterEdgeConfig& edge_config)
    -> std::expected<BoundaryLayerGrid, GridError> {

  auto grid = BoundaryLayerGrid(numerical_config.n_eta, numerical_config.eta_max);

  grid.generate_eta_distribution();
  grid.xi_.emplace_back(0.0);
  grid.xi_output_.emplace_back(0.0);

  return grid;
}

inline constexpr auto BoundaryLayerGrid::generate_eta_distribution() noexcept -> void {
  eta_.clear();

  // Modern range-based generation
  auto eta_indices = std::views::iota(0, n_eta_);
  std::ranges::transform(eta_indices, std::back_inserter(eta_),
                         [this](int i) constexpr { return static_cast<double>(i) * d_eta_; });
}

} // namespace blast::boundary_layer::grid