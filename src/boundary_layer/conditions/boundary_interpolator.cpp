#include "blast/boundary_layer/conditions/boundary_interpolator.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include <algorithm>
#include <expected>
#include <format>
#include <numeric>

namespace blast::boundary_layer::conditions {

using blast::boundary_layer::grid::coordinate_transform::hermite_interpolate;
using blast::boundary_layer::grid::coordinate_transform::linear_interpolate;
using blast::boundary_layer::grid::coordinate_transform::search_interval;

namespace {

template <typename Container>
[[nodiscard]] auto interpolate_property(const Container& values, std::size_t i1, std::size_t i2,
                                        const std::vector<double>& x_grid, double x_target) -> double {

  if (i1 == i2)
    return values[i1];

  return linear_interpolate(x_target, x_grid[i1], x_grid[i2], values[i1], values[i2]);
}

// Overload with optional derivatives for Hermite interpolation
template <typename Container>
[[nodiscard]] auto interpolate_property(const Container& values, std::size_t i1, std::size_t i2,
                                        const std::vector<double>& x_grid, double x_target,
                                        const Container& derivatives) -> double {

  if (i1 == i2)
    return values[i1];

  // Use Hermite interpolation when derivatives are available
  return hermite_interpolate(x_target, x_grid[i1], x_grid[i2], values[i1], values[i2], derivatives[i1],
                             derivatives[i2]);
}

template <typename MemberPtr>
[[nodiscard]] auto extract_edge_property(const std::vector<io::OuterEdgeConfig::EdgePoint>& points,
                                         MemberPtr member_ptr) -> std::vector<double> {

  std::vector<double> result;
  result.reserve(points.size());

  std::ranges::transform(points, std::back_inserter(result),
                         [member_ptr](const auto& point) { return point.*member_ptr; });

  return result;
}

// Compute derivatives reusing existing high-order function
[[nodiscard]] auto compute_property_derivatives(const std::vector<double>& values, const std::vector<double>& x_grid)
    -> std::expected<std::vector<double>, BoundaryConditionError> {

  if (values.size() != x_grid.size()) {
    return std::unexpected(BoundaryConditionError("Values and grid sizes must match for derivative computation"));
  }
  if (values.size() < 2) {
    return std::unexpected(BoundaryConditionError("Need at least 2 points for derivative computation"));
  }

  // Check if grid is uniform (within tolerance)
  constexpr double tolerance = 1e-12;
  bool is_uniform = true;
  const double dx_first = x_grid[1] - x_grid[0];

  for (std::size_t i = 2; i < x_grid.size(); ++i) {
    const double dx_current = x_grid[i] - x_grid[i - 1];
    if (std::abs(dx_current - dx_first) > tolerance) {
      is_uniform = false;
      break;
    }
  }

  if (is_uniform) {
    // Reuse existing high-order O(h^4) function for uniform grids
    auto result = coefficients::derivatives::compute_eta_derivative(values, dx_first);
    if (!result) {
      return std::unexpected(BoundaryConditionError("Failed to compute property derivatives"));
    }
    return result.value();
  } else {
    // Handle non-uniform grids with weighted central differences
    const auto n = values.size();
    std::vector<double> derivatives(n);

    // Forward difference at start
    derivatives[0] = (values[1] - values[0]) / (x_grid[1] - x_grid[0]);

    // Central differences in interior
    for (std::size_t i = 1; i < n - 1; ++i) {
      const double dx_left = x_grid[i] - x_grid[i - 1];
      const double dx_right = x_grid[i + 1] - x_grid[i];
      const double dx_total = x_grid[i + 1] - x_grid[i - 1];

      // Weighted central difference for non-uniform grids
      derivatives[i] =
          ((values[i + 1] - values[i]) / dx_right * dx_left + (values[i] - values[i - 1]) / dx_left * dx_right) /
          dx_total;
    }

    // Backward difference at end
    derivatives[n - 1] = (values[n - 1] - values[n - 2]) / (x_grid[n - 1] - x_grid[n - 2]);

    return derivatives;
  }
}

} // namespace

constexpr auto compute_beta(int station, double xi, const io::SimulationConfig& sim_config, double ue,
                            double d_ue_dxi) noexcept -> double {

  if (station == 0) {
    switch (sim_config.body_type) {
    case io::SimulationConfig::BodyType::Axisymmetric:
      return 0.5;
    case io::SimulationConfig::BodyType::Cone:
      return 0.0;
    case io::SimulationConfig::BodyType::FlatPlate:
      return 0.0;
    case io::SimulationConfig::BodyType::TwoD:
      return 1.0;
    }
  }

  return 2.0 * xi * d_ue_dxi / ue;
}

auto create_stagnation_conditions(const io::OuterEdgeConfig& edge_config, 
                                 const io::WallParametersConfig& wall_config,
                                 const io::SimulationConfig& sim_config,
                                 const thermophysics::MixtureInterface& mixture)
   -> std::expected<BoundaryConditions, BoundaryConditionError> {

 if (edge_config.edge_points.empty()) {
   return std::unexpected(BoundaryConditionError("No edge points defined"));
 }

 if (wall_config.wall_temperatures.empty()) {
   return std::unexpected(BoundaryConditionError("No wall temperatures defined"));
 }

 const auto& edge_point = edge_config.edge_points[0];

 auto eq_result = mixture.equilibrium_composition(edge_point.temperature, edge_point.pressure);
 if (!eq_result) {
   return std::unexpected(BoundaryConditionError(
     std::format("Failed to compute equilibrium composition: {}", eq_result.error().message())));
 }
 auto species_fractions = eq_result.value();

 auto h_result = mixture.mixture_enthalpy(species_fractions, edge_point.temperature, edge_point.pressure);
 if (!h_result) {
   return std::unexpected(BoundaryConditionError("Failed to compute mixture enthalpy"));
 }
 auto enthalpy = h_result.value();

 auto mw_result = mixture.mixture_molecular_weight(species_fractions);
 if (!mw_result) {
   return std::unexpected(BoundaryConditionError("Failed to compute mixture MW"));
 }
 auto MW = mw_result.value();
 auto density = edge_point.pressure * MW / (edge_point.temperature * thermophysics::constants::R_universal);

 auto visc_result = mixture.viscosity(species_fractions, edge_point.temperature, edge_point.pressure);
 if (!visc_result) {
   return std::unexpected(BoundaryConditionError("Failed to compute viscosity"));
 }
 auto viscosity = visc_result.value();

 EdgeConditions edge{.pressure = edge_point.pressure,
                     .viscosity = viscosity,
                     .velocity = edge_point.velocity,
                     .enthalpy = enthalpy,
                     .density = density,
                     .species_fractions = species_fractions,
                     .d_xi_dx = edge_config.velocity_gradient_stagnation,
                     .d_ue_dx = edge_config.velocity_gradient_stagnation,
                     .d_he_dx = 0.0,
                     .d_he_dxi = 0.0,
                     .body_radius = edge_point.radius};

 WallConditions wall{.temperature = wall_config.wall_temperatures[0]};

 return BoundaryConditions{.edge = std::move(edge),
                           .wall = std::move(wall),
                           .beta = compute_beta(0, 0.0, sim_config),
                           .station = 0,
                           .xi = 0.0};
}

auto interpolate_boundary_conditions(
    int station, double xi, std::span<const double> xi_grid, const io::OuterEdgeConfig& edge_config,
    const io::WallParametersConfig& wall_config,
    const io::SimulationConfig& sim_config, const thermophysics::MixtureInterface& mixture) -> std::expected<BoundaryConditions, BoundaryConditionError> {

  // Stagnation point special case
  if (station == 0) {
    return create_stagnation_conditions(edge_config, wall_config, sim_config, mixture);
  }

  // Find interpolation interval in xi space
  auto interval_result = search_interval(xi_grid, xi);
  if (!interval_result) {
    return std::unexpected(
        BoundaryConditionError(std::format("Failed to find xi={} in grid: {}", xi, interval_result.error().message())));
  }

  const auto [i1, i2] = interval_result.value();

  // Extract x coordinates for interpolation
  auto x_grid = extract_edge_property(edge_config.edge_points, &io::OuterEdgeConfig::EdgePoint::x);

  // Compute interpolated x position
  const double x_interp =
      (i1 == i2) ? x_grid[i1] : linear_interpolate(xi, xi_grid[i1], xi_grid[i2], x_grid[i1], x_grid[i2]);

  // Find interval in physical x space
  auto x_interval_result = search_interval(std::span(x_grid), x_interp);
  if (!x_interval_result) {
    return std::unexpected(BoundaryConditionError(std::format("Failed to find x={} in edge grid", x_interp)));
  }

  const auto [ix1, ix2] = x_interval_result.value();

  auto u_edge = extract_edge_property(edge_config.edge_points, &io::OuterEdgeConfig::EdgePoint::velocity);
  auto t_edge = extract_edge_property(edge_config.edge_points, &io::OuterEdgeConfig::EdgePoint::temperature);

  auto du_dx_result = compute_property_derivatives(u_edge, x_grid);
  if (!du_dx_result) {
    return std::unexpected(
        BoundaryConditionError(std::format("Failed to compute d_ue_dx: {}", du_dx_result.error().message())));
  }
  auto du_dx = du_dx_result.value();

  auto dt_dx_result = compute_property_derivatives(t_edge, x_grid);
  if (!dt_dx_result) {
    return std::unexpected(
        BoundaryConditionError(std::format("Failed to compute d_te_dx: {}", dt_dx_result.error().message())));
  }
  auto dt_dx = dt_dx_result.value();

  const double d_ue_dx_interp = interpolate_property(du_dx, ix1, ix2, x_grid, x_interp);
  const double d_te_dx_interp = interpolate_property(dt_dx, ix1, ix2, x_grid, x_interp);

  auto interp_linear = [&](auto member_ptr) {
    auto values = extract_edge_property(edge_config.edge_points, member_ptr);
    return interpolate_property(values, ix1, ix2, x_grid, x_interp);
  };

  auto interp_hermite = [&](auto member_ptr) -> std::expected<double, BoundaryConditionError> {
    auto values = extract_edge_property(edge_config.edge_points, member_ptr);
    auto derivatives_result = compute_property_derivatives(values, x_grid);
    if (!derivatives_result) {
      return std::unexpected(derivatives_result.error());
    }
    return interpolate_property(values, ix1, ix2, x_grid, x_interp, derivatives_result.value());
  };

  auto pressure_result = interp_hermite(&io::OuterEdgeConfig::EdgePoint::pressure);
  if (!pressure_result) {
    return std::unexpected(pressure_result.error());
  }

  auto velocity_result = interp_hermite(&io::OuterEdgeConfig::EdgePoint::velocity);
  if (!velocity_result) {
    return std::unexpected(velocity_result.error());
  }

  auto temperature_result = interp_hermite(&io::OuterEdgeConfig::EdgePoint::temperature);
  if (!temperature_result) {
    return std::unexpected(temperature_result.error());
  }

  const double radius_interp = interp_linear(&io::OuterEdgeConfig::EdgePoint::radius);

  // Calculate derived properties from interpolated temperature and pressure
  const double temperature_interp = temperature_result.value();
  const double pressure_interp = pressure_result.value();

  auto eq_result = mixture.equilibrium_composition(temperature_interp, pressure_interp);
  if (!eq_result) {
    return std::unexpected(BoundaryConditionError("Failed to compute equilibrium composition"));
  }
  auto species_fractions = eq_result.value();

  auto h_result = mixture.mixture_enthalpy(species_fractions, temperature_interp, pressure_interp);
  if (!h_result) {
    return std::unexpected(BoundaryConditionError("Failed to compute enthalpy"));
  }
  auto enthalpy_interp = h_result.value();

  auto mw_result = mixture.mixture_molecular_weight(species_fractions);
  if (!mw_result) {
    return std::unexpected(BoundaryConditionError("Failed to compute MW"));
  }
  auto density_interp = pressure_interp * mw_result.value() / (temperature_interp * thermophysics::constants::R_universal);

  auto visc_result = mixture.viscosity(species_fractions, temperature_interp, pressure_interp);
  if (!visc_result) {
    return std::unexpected(BoundaryConditionError("Failed to compute viscosity"));
  }
  auto viscosity_interp = visc_result.value();

  // Calculate dh/dx from dt/dx
  auto cp_result = mixture.frozen_cp(species_fractions, temperature_interp, pressure_interp);
  if (!cp_result) {
    return std::unexpected(BoundaryConditionError("Failed to compute Cp"));
  }
  const double d_he_dx_interp = cp_result.value() * d_te_dx_interp;

  // By definition: dξ/dx = ρ_e(x) * μ_e(x) * u_e(x) * r_body(x)²
  const double d_xi_dx_calculated =
      density_interp * viscosity_interp * velocity_result.value() * radius_interp * radius_interp;

  if (d_xi_dx_calculated <= 0.0 || !std::isfinite(d_xi_dx_calculated)) {
    return std::unexpected(
        BoundaryConditionError(std::format("Invalid calculated d_xi_dx={} at x={}", d_xi_dx_calculated, x_interp)));
  }

  EdgeConditions edge{.pressure = pressure_interp,
                      .viscosity = viscosity_interp,
                      .velocity = velocity_result.value(),
                      .enthalpy = enthalpy_interp,
                      .density = density_interp,
                      .species_fractions = species_fractions,
                      .d_xi_dx = d_xi_dx_calculated,
                      .d_ue_dx = d_ue_dx_interp,
                      .d_he_dx = d_he_dx_interp,
                      .d_he_dxi = 0.0,
                      .body_radius = radius_interp};

  if (wall_config.wall_temperatures.size() != edge_config.edge_points.size()) {
    return std::unexpected(
        BoundaryConditionError(std::format("Wall temperatures array size ({}) doesn't match edge points size ({})",
                                           wall_config.wall_temperatures.size(), edge_config.edge_points.size())));
  }

  auto wall_temp_derivatives_result = compute_property_derivatives(wall_config.wall_temperatures, x_grid);
  if (!wall_temp_derivatives_result) {
    return std::unexpected(BoundaryConditionError(std::format("Failed to compute wall temperature derivatives: {}",
                                                              wall_temp_derivatives_result.error().message())));
  }

  const double wall_temp = interpolate_property(wall_config.wall_temperatures, ix1, ix2, x_grid, x_interp,
                                                wall_temp_derivatives_result.value());

  if (wall_temp <= 0.0 || !std::isfinite(wall_temp)) {
    return std::unexpected(BoundaryConditionError(std::format("Invalid wall temperature: {} K", wall_temp)));
  }

  WallConditions wall{.temperature = wall_temp};

  const double d_ue_dxi = edge.d_ue_dx / edge.d_xi_dx;
  edge.d_he_dxi = edge.d_he_dx / edge.d_xi_dx;

  if (!std::isfinite(d_ue_dxi) || !std::isfinite(edge.d_he_dxi)) {
    return std::unexpected(BoundaryConditionError(
        std::format("Invalid computed derivatives: d_ue_dxi={}, d_he_dxi={}", d_ue_dxi, edge.d_he_dxi)));
  }

  return BoundaryConditions{.edge = std::move(edge),
                            .wall = std::move(wall),
                            .beta = compute_beta(station, xi, sim_config, edge.velocity, d_ue_dxi),
                            .station = station,
                            .xi = xi};
}

} // namespace blast::boundary_layer::conditions