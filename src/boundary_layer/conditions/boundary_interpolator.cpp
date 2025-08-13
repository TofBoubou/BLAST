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

[[nodiscard]] static auto get_species_fractions(
    const io::OuterEdgeConfig::EdgePoint& edge_point,
    const thermophysics::MixtureInterface& mixture)
    -> std::expected<std::vector<double>, BoundaryConditionError> {
  
  if (edge_point.boundary_override_enabled()) {
    // Utiliser les fractions spécifiées par l'utilisateur
    if (!edge_point.mass_fraction_condition()) {
      return std::unexpected(BoundaryConditionError(
        "mass_fraction_condition missing despite boundary_override being enabled"));
    }
    
    // Vérifier que le nombre d'espèces correspond
    if (edge_point.mass_fraction_condition()->size() != mixture.n_species()) {
      return std::unexpected(BoundaryConditionError(
        std::format("mass_fraction_condition size ({}) doesn't match number of species ({})",
                   edge_point.mass_fraction_condition()->size(), mixture.n_species())));
    }
    
    return edge_point.mass_fraction_condition().value();
  } else {
    // Comportement par défaut : calculer l'équilibre
    auto eq_result = mixture.equilibrium_composition(edge_point.temperature, edge_point.pressure);
    if (!eq_result) {
      return std::unexpected(BoundaryConditionError(
        std::format("Failed to compute equilibrium composition: {}", eq_result.error().message())));
    }
    return eq_result.value();
  }
}

[[nodiscard]] auto create_stagnation_conditions(const io::OuterEdgeConfig& edge_config,
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

  // Utilise la nouvelle fonction helper pour obtenir les fractions d'espèces
  auto species_fractions_result = get_species_fractions(edge_point, mixture);
  if (!species_fractions_result) {
    return std::unexpected(species_fractions_result.error());
  }
  auto species_fractions = species_fractions_result.value();

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

  EdgeConditions edge{
      .pressure = edge_point.pressure,
      .viscosity = viscosity,
      .velocity = edge_point.velocity,
      .enthalpy = enthalpy,
      .density = density,
      .species_fractions = species_fractions,
      .d_xi_dx = edge_config.velocity_gradient_stagnation,
      .d_ue_dx = edge_config.velocity_gradient_stagnation,
      .d_he_dx = 0.0,
      .d_he_dxi = 0.0,
      .body_radius = edge_point.radius,
      .boundary_override = edge_point.boundary_override_enabled()
  };

  WallConditions wall{.temperature = wall_config.wall_temperatures[0]};

  return BoundaryConditions{
      .edge = edge,
      .wall = wall,
      .beta = 0.0,
      .station = 0,
      .xi = 0.0
  };
}

[[nodiscard]] auto interpolate_boundary_conditions(
    int station, double xi, std::span<const double> xi_grid,
    const io::OuterEdgeConfig& edge_config,
    const io::WallParametersConfig& wall_config,
    const io::SimulationConfig& sim_config,
    const thermophysics::MixtureInterface& mixture)
    -> std::expected<BoundaryConditions, BoundaryConditionError>
{
  // Stagnation point special case
  if (station == 0) {
    return create_stagnation_conditions(edge_config, wall_config, sim_config, mixture);
  }

  // Locate xi interval
  auto interval_result = search_interval(xi_grid, xi);
  if (!interval_result) {
    return std::unexpected(BoundaryConditionError(
        std::format("Failed to find xi = {} in grid: {}", xi, interval_result.error().message())));
  }
  const auto [i1, i2] = interval_result.value();

  // Extract x-grid from edge config and interpolate x at current xi
  const auto x_grid = extract_edge_property(edge_config.edge_points, &io::OuterEdgeConfig::EdgePoint::x);
  const double x_interp = (i1 == i2)
      ? x_grid[i1]
      : linear_interpolate(xi, xi_grid[i1], xi_grid[i2], x_grid[i1], x_grid[i2]);

  // Find interval in physical x
  auto x_interval_result = search_interval(std::span(x_grid), x_interp);
  if (!x_interval_result) {
    return std::unexpected(BoundaryConditionError(
        std::format("Failed to find x = {} in edge grid", x_interp)));
  }
  const auto [ix1, ix2] = x_interval_result.value();

  // Prepare edge properties and their derivatives
  const auto u_edge = extract_edge_property(edge_config.edge_points, &io::OuterEdgeConfig::EdgePoint::velocity);
  const auto t_edge = extract_edge_property(edge_config.edge_points, &io::OuterEdgeConfig::EdgePoint::temperature);

  auto du_dx_result = compute_property_derivatives(u_edge, x_grid);
  if (!du_dx_result) {
    return std::unexpected(BoundaryConditionError(
        std::format("Failed to compute d_ue_dx: {}", du_dx_result.error().message())));
  }
  const auto du_dx = du_dx_result.value();

  auto dt_dx_result = compute_property_derivatives(t_edge, x_grid);
  if (!dt_dx_result) {
    return std::unexpected(BoundaryConditionError(
        std::format("Failed to compute d_te_dx: {}", dt_dx_result.error().message())));
  }
  const auto dt_dx = dt_dx_result.value();

  const double d_ue_dx_interp = interpolate_property(du_dx, ix1, ix2, x_grid, x_interp);
  const double d_te_dx_interp = interpolate_property(dt_dx, ix1, ix2, x_grid, x_interp);

  // Helpers for (Hermite) interpolation of scalar fields defined on edge points
  const auto interp_linear = [&](auto member_ptr) {
    const auto values = extract_edge_property(edge_config.edge_points, member_ptr);
    return interpolate_property(values, ix1, ix2, x_grid, x_interp);
  };

  const auto interp_hermite = [&](auto member_ptr) -> std::expected<double, BoundaryConditionError> {
    const auto values = extract_edge_property(edge_config.edge_points, member_ptr);
    auto dvals_dx = compute_property_derivatives(values, x_grid);
    if (!dvals_dx) return std::unexpected(dvals_dx.error());
    return interpolate_property(values, ix1, ix2, x_grid, x_interp, dvals_dx.value());
  };

  // Interpolate primitive edge conditions (Hermite for smoothness)
  auto pressure_result    = interp_hermite(&io::OuterEdgeConfig::EdgePoint::pressure);
  if (!pressure_result)   return std::unexpected(pressure_result.error());
  auto velocity_result    = interp_hermite(&io::OuterEdgeConfig::EdgePoint::velocity);
  if (!velocity_result)   return std::unexpected(velocity_result.error());
  auto temperature_result = interp_hermite(&io::OuterEdgeConfig::EdgePoint::temperature);
  if (!temperature_result) return std::unexpected(temperature_result.error());

  const double pressure_interp    = pressure_result.value();
  const double temperature_interp = temperature_result.value();
  const double radius_interp      = interp_linear(&io::OuterEdgeConfig::EdgePoint::radius);

  // === Species fractions: override vs equilibrium ===
  const bool all_override = std::all_of(edge_config.edge_points.begin(), edge_config.edge_points.end(),
                                        [](const auto& p){ return p.boundary_override_enabled(); });
  const bool none_override = std::none_of(edge_config.edge_points.begin(), edge_config.edge_points.end(),
                                          [](const auto& p){ return p.boundary_override_enabled(); });
  
  std::vector<double> species_fractions;
  {

    if (!all_override && !none_override) {
      return std::unexpected(BoundaryConditionError(
        "All edge points must have the same boundary_override setting for interpolation"));
    }

    if (all_override) {
      // Interpolate each species mass fraction provided at edge points
      const std::size_t n_species = mixture.n_species();
      if (n_species == 0) {
        return std::unexpected(BoundaryConditionError("Mixture reports zero species"));
      }
      species_fractions.assign(n_species, 0.0);

      for (std::size_t j = 0; j < n_species; ++j) {
        std::vector<double> yj_values;
        yj_values.reserve(edge_config.edge_points.size());
        for (const auto& pt : edge_config.edge_points) {
          if (!pt.mass_fraction_condition()) {
            return std::unexpected(BoundaryConditionError(
              "Missing mass_fraction_condition in edge point with boundary_override enabled"));
          }
          const auto& y = *pt.mass_fraction_condition();
          if (y.size() != n_species) {
            return std::unexpected(BoundaryConditionError(
              std::format("mass_fraction_condition size mismatch (got {}, expected {})",
                          y.size(), n_species)));
          }
          yj_values.push_back(y[j]);
        }
        // Interpolate species j at x_interp
        species_fractions[j] = interpolate_property(yj_values, ix1, ix2, x_grid, x_interp);
        if (!std::isfinite(species_fractions[j])) {
          return std::unexpected(BoundaryConditionError(
            std::format("Non-finite interpolated mass fraction for species {} at x = {}", j, x_interp)));
        }
      }

      // Clean small negatives (interp artifacts), enforce non-negativity and renormalize
      for (auto& v : species_fractions) {
        if (v < 0.0 && v > -1e-12) v = 0.0;
        if (v < 0.0) {
          return std::unexpected(BoundaryConditionError(
            "Interpolated mass fractions contain negative values beyond tolerance"));
        }
      }
      double sum = std::accumulate(species_fractions.begin(), species_fractions.end(), 0.0);
      if (!(sum > 0.0) || !std::isfinite(sum)) {
        return std::unexpected(BoundaryConditionError("Sum of interpolated mass fractions is not positive/finite"));
      }
      for (auto& v : species_fractions) v /= sum;
    } else {
      // Default: equilibrium composition at interpolated (T, p)
      auto eq_result = mixture.equilibrium_composition(temperature_interp, pressure_interp);
      if (!eq_result) {
        return std::unexpected(BoundaryConditionError("Failed to compute equilibrium composition"));
      }
      species_fractions = std::move(eq_result.value());
      if (species_fractions.empty()) {
        return std::unexpected(BoundaryConditionError("Equilibrium composition returned empty vector"));
      }
    }
  }

  // Thermophysical properties based on chosen composition
  auto h_result = mixture.mixture_enthalpy(species_fractions, temperature_interp, pressure_interp);
  if (!h_result) {
    return std::unexpected(BoundaryConditionError("Failed to compute enthalpy"));
  }
  const double enthalpy_interp = h_result.value();

  auto mw_result = mixture.mixture_molecular_weight(species_fractions);
  if (!mw_result) {
    return std::unexpected(BoundaryConditionError("Failed to compute MW"));
  }
  const double density_interp =
      pressure_interp * mw_result.value() / (temperature_interp * thermophysics::constants::R_universal);

  auto visc_result = mixture.viscosity(species_fractions, temperature_interp, pressure_interp);
  if (!visc_result) {
    return std::unexpected(BoundaryConditionError("Failed to compute viscosity"));
  }
  const double viscosity_interp = visc_result.value();

  // dh/dx from Cp(T,p,Y) * dT/dx
  auto cp_result = mixture.frozen_cp(species_fractions, temperature_interp, pressure_interp);
  if (!cp_result) {
    return std::unexpected(BoundaryConditionError("Failed to compute Cp"));
  }
  const double d_he_dx_interp = cp_result.value() * d_te_dx_interp;

  // By definition: dξ/dx = ρ_e * μ_e * u_e * r(x)^2
  const double d_xi_dx_calculated =
      density_interp * viscosity_interp * velocity_result.value() * radius_interp * radius_interp;

  if (!(d_xi_dx_calculated > 0.0) || !std::isfinite(d_xi_dx_calculated)) {
    return std::unexpected(BoundaryConditionError(
      std::format("Invalid calculated d_xi_dx = {} at x = {}", d_xi_dx_calculated, x_interp)));
  }

  // Wall temperature (Hermite)
  if (wall_config.wall_temperatures.size() != edge_config.edge_points.size()) {
    return std::unexpected(BoundaryConditionError(
        std::format("Wall temperatures size ({}) doesn't match edge points size ({})",
                    wall_config.wall_temperatures.size(), edge_config.edge_points.size())));
  }
  auto dTw_dx = compute_property_derivatives(wall_config.wall_temperatures, x_grid);
  if (!dTw_dx) {
    return std::unexpected(BoundaryConditionError(
        std::format("Failed to compute wall temperature derivatives: {}", dTw_dx.error().message())));
  }
  const double wall_temp = interpolate_property(wall_config.wall_temperatures, ix1, ix2, x_grid, x_interp, dTw_dx.value());
  if (!(wall_temp > 0.0) || !std::isfinite(wall_temp)) {
    return std::unexpected(BoundaryConditionError(std::format("Invalid wall temperature: {} K", wall_temp)));
  }

  EdgeConditions edge{
    .pressure          = pressure_interp,
    .viscosity         = viscosity_interp,
    .velocity          = velocity_result.value(),
    .enthalpy          = enthalpy_interp,
    .density           = density_interp,
    .species_fractions = species_fractions,
    .d_xi_dx           = d_xi_dx_calculated,
    .d_ue_dx           = d_ue_dx_interp,
    .d_he_dx           = d_he_dx_interp,
    .d_he_dxi          = 0.0, // set below
    .body_radius       = radius_interp,
    .boundary_override = all_override
  };

  const double d_ue_dxi = edge.d_ue_dx / edge.d_xi_dx;
  edge.d_he_dxi = edge.d_he_dx / edge.d_xi_dx;

  if (!std::isfinite(d_ue_dxi) || !std::isfinite(edge.d_he_dxi)) {
    return std::unexpected(BoundaryConditionError(
        std::format("Invalid computed derivatives: d_ue_dxi = {}, d_he_dxi = {}",
                    d_ue_dxi, edge.d_he_dxi)));
  }

  WallConditions wall{ .temperature = wall_temp };

  return BoundaryConditions{
    .edge    = std::move(edge),
    .wall    = std::move(wall),
    .beta    = compute_beta(station, xi, sim_config, /*u_e=*/velocity_result.value(), d_ue_dxi),
    .station = station,
    .xi      = xi
  };
}



} // namespace blast::boundary_layer::conditions