#pragma once
#include "../../core/exceptions.hpp"
#include "../../io/config_types.hpp"
#include <algorithm>
#include <cmath>
#include <concepts>
#include <expected>
#include <iostream>
#include <ranges>
#include <span>
#include <vector>

namespace blast::boundary_layer::grid::coordinate_transform {
using blast::core::TransformError;

// Simpson integration constants
namespace simpson_constants {
constexpr double coeff_5pt_1 = 17.0;
constexpr double coeff_5pt_2 = 42.0;
constexpr double coeff_5pt_3 = -16.0;
constexpr double coeff_5pt_4 = 6.0;
constexpr double coeff_5pt_5 = -1.0;
constexpr double divisor_5pt = 48.0;

constexpr double coeff_3pt_1 = 1.0;
constexpr double coeff_3pt_2 = 4.0;
constexpr double coeff_3pt_3 = 1.0;
constexpr double divisor_3pt = 3.0;
} // namespace simpson_constants

// Concepts for type safety
template <typename T>
concept PhysicalQuantity = std::floating_point<T> && requires(T t) {
  { t >= T{0} } -> std::convertible_to<bool>;
};

template <typename Range>
concept PhysicalRange = std::ranges::sized_range<Range> && PhysicalQuantity<std::ranges::range_value_t<Range>>;

template <PhysicalRange Values>
auto simpson_integrate(Values&& f_values, PhysicalQuantity auto dx, PhysicalQuantity auto f0) -> std::vector<double> {

  const auto n = std::ranges::size(f_values);
  std::vector<double> result(n);

  if (n == 0)
    return result;

  result[0] = f0;
  if (n == 1)
    return result;

  using namespace simpson_constants;

  // Special first point using 5-point formula if possible
  if (n >= 5) {
    result[1] = f0 + (coeff_5pt_1 * f_values[0] + coeff_5pt_2 * f_values[1] + coeff_5pt_3 * f_values[2] +
                      coeff_5pt_4 * f_values[3] + coeff_5pt_5 * f_values[4]) *
                         dx / divisor_5pt;
  } else {
    result[1] = f0 + (f_values[0] + f_values[1]) * dx * 0.5; // Trapezoidal fallback
  }

  // Simpson's rule for remaining points using modern range
  for (auto i : std::views::iota(size_t{2}, n)) {
    result[i] = result[i - 2] + (f_values[i - 2] + coeff_3pt_2 * f_values[i - 1] + f_values[i]) * dx / divisor_3pt;
  }

  return result;
}

// Core transformation functions with concepts
template <PhysicalRange EdgeData>
auto compute_xi_from_integration(EdgeData&& x_edge, EdgeData&& rho_edge, EdgeData&& mu_edge, EdgeData&& u_edge,
                                 EdgeData&& r_body) -> std::expected<std::vector<double>, TransformError> {

  constexpr auto min_points = 5;
  const auto n_points = std::ranges::size(x_edge);

  if (n_points < min_points) {
    return std::unexpected(TransformError("Need at least 5 points for xi integration"));
  }

  // Check uniform grid with constexpr tolerance
  constexpr auto tolerance = 1e-6;
  const auto dx = x_edge[1] - x_edge[0];

  auto x_diffs = std::views::zip(x_edge, x_edge | std::views::drop(1)) | std::views::transform([](const auto& pair) {
                   return std::abs(std::get<1>(pair) - std::get<0>(pair));
                 });

  if (!std::ranges::all_of(x_diffs, [dx, tolerance](auto diff) { return std::abs(diff - dx) <= tolerance; })) {
    return std::unexpected(TransformError("Non-uniform x grid not supported"));
  }

  std::vector<double> xi(n_points);
  xi[0] = 0.0;

  auto integrand_range =
      std::views::zip(rho_edge, mu_edge, u_edge, r_body) | std::views::transform([](const auto& tuple) {
        const auto [rho, mu, u, r] = tuple;
        return rho * mu * u * r * r;
      });

  std::vector<double> integrand_values;
  std::ranges::copy(integrand_range, std::back_inserter(integrand_values));

  // Integration using Simpson's rule with named constants
  xi = simpson_integrate(integrand_values, dx, 0.0);

  return xi;
}

template <PhysicalRange RhoData>
[[nodiscard]] auto
compute_physical_y_from_eta(PhysicalQuantity auto eta, PhysicalQuantity auto xi, PhysicalQuantity auto rho_wall,
                            RhoData&& rho_eta, PhysicalQuantity auto d_eta, const io::SimulationConfig& sim_config,
                            PhysicalQuantity auto d_ue_dx = 0.0, PhysicalQuantity auto rho_e = 0.0,
                            PhysicalQuantity auto mu_e = 0.0, PhysicalQuantity auto u_e = 0.0,
                            PhysicalQuantity auto r_body = 1.0) -> std::expected<double, TransformError>;

[[nodiscard]] constexpr auto
compute_derivative_factor(int station, PhysicalQuantity auto xi, const io::SimulationConfig& sim_config,
                          PhysicalQuantity auto d_ue_dx, PhysicalQuantity auto rho_e, PhysicalQuantity auto mu_e,
                          PhysicalQuantity auto u_e = 0.0,
                          PhysicalQuantity auto r_body = 1.0) -> std::expected<double, TransformError>;

[[nodiscard]] constexpr auto linear_interpolate(std::floating_point auto x, std::floating_point auto x1,
                                                std::floating_point auto x2, std::floating_point auto y1,
                                                std::floating_point auto y2) noexcept -> double {
  constexpr auto epsilon = 1e-15;
  return (std::abs(x2 - x1) < epsilon) ? y1 : y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

// Hermite interpolation with derivatives for higher accuracy O(h^4)
[[nodiscard]] constexpr auto hermite_interpolate(std::floating_point auto x, std::floating_point auto x1,
                                                 std::floating_point auto x2, std::floating_point auto y1,
                                                 std::floating_point auto y2,
                                                 std::floating_point auto dy1, // derivative at x1
                                                 std::floating_point auto dy2  // derivative at x2
                                                 ) noexcept -> double {
  constexpr auto epsilon = 1e-15;

  // Handle degenerate case
  if (std::abs(x2 - x1) < epsilon) {
    return y1;
  }

  const double h = x2 - x1;
  const double t = (x - x1) / h;
  const double t2 = t * t;
  const double t3 = t2 * t;

  // Hermite basis functions
  const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0; // (1,0) at x1
  const double h10 = t3 - 2.0 * t2 + t;         // (0,1) at x1
  const double h01 = -2.0 * t3 + 3.0 * t2;      // (1,0) at x2
  const double h11 = t3 - t2;                   // (0,1) at x2

  return h00 * y1 + h10 * h * dy1 + h01 * y2 + h11 * h * dy2;
}

template <std::ranges::random_access_range Grid>
[[nodiscard]] auto search_interval(Grid&& grid, std::floating_point auto target) noexcept
    -> std::expected<std::pair<int, int>, TransformError> {
  if (std::ranges::empty(grid)) {
    return std::unexpected(TransformError("Empty grid"));
  }

  constexpr auto epsilon = 1e-15;
  if (auto exact_it = std::ranges::find_if(grid, [target](auto val) { return std::abs(val - target) < epsilon; });
      exact_it != std::ranges::end(grid)) {
    const auto index = static_cast<int>(std::ranges::distance(std::ranges::begin(grid), exact_it));
    return std::make_pair(index, index);
  }

  auto upper_it = std::ranges::upper_bound(grid, target);

  if (upper_it == std::ranges::begin(grid)) {
    return std::unexpected(TransformError("Target below grid minimum"));
  }

  if (upper_it == std::ranges::end(grid)) {
    return std::unexpected(TransformError("Target above grid maximum"));
  }

  const auto i2 = static_cast<int>(std::ranges::distance(std::ranges::begin(grid), upper_it));
  const auto i1 = i2 - 1;

  return std::make_pair(i1, i2);
}

constexpr auto compute_stagnation_y_factor(const io::SimulationConfig& sim_config, PhysicalQuantity auto d_ue_dx,
                                           PhysicalQuantity auto rho_e,
                                           PhysicalQuantity auto mu_e) -> std::expected<double, TransformError>;

} // namespace blast::boundary_layer::grid::coordinate_transform