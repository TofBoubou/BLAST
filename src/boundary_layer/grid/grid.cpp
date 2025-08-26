
#include "blast/boundary_layer/grid/grid.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include "blast/core/constants.hpp"
#include <algorithm>
#include <numeric>
#include <ranges>

namespace blast::boundary_layer::grid {

namespace {
// Helper function to compute coordinate derivatives for Hermite interpolation
template <typename CoordGrid>
[[nodiscard]] auto compute_coordinate_derivatives(const CoordGrid& coord_values,
                                                  const std::vector<double>& grid_spacing) -> std::vector<double> {

  if (coord_values.size() != grid_spacing.size() || coord_values.size() < 2) {
    return std::vector<double>(coord_values.size(), 0.0);
  }

  // Check if grid spacing is uniform
  constexpr double tolerance = 1e-12;
  bool is_uniform = true;
  const double dx_first = grid_spacing[1] - grid_spacing[0];

  for (std::size_t i = 2; i < grid_spacing.size(); ++i) {
    const double dx_current = grid_spacing[i] - grid_spacing[i - 1];
    if (std::abs(dx_current - dx_first) > tolerance) {
      is_uniform = false;
      break;
    }
  }

  if (is_uniform) {
    // Reuse high-order function for uniform grids
    auto result = coefficients::derivatives::compute_eta_derivative(coord_values, dx_first);
    if (!result) {
      return std::vector<double>(coord_values.size(), 0.0);
    }
    return result.value();
  } else {
    // Handle non-uniform grids
    const auto n = coord_values.size();
    std::vector<double> derivatives(n);

    // Forward difference at start
    derivatives[0] = (coord_values[1] - coord_values[0]) / (grid_spacing[1] - grid_spacing[0]);

    // Central differences in interior
    for (std::size_t i = 1; i < n - 1; ++i) {
      const double dx_left = grid_spacing[i] - grid_spacing[i - 1];
      const double dx_right = grid_spacing[i + 1] - grid_spacing[i];
      const double dx_total = grid_spacing[i + 1] - grid_spacing[i - 1];

      derivatives[i] = ((coord_values[i + 1] - coord_values[i]) / dx_right * dx_left +
                        (coord_values[i] - coord_values[i - 1]) / dx_left * dx_right) /
                       dx_total;
    }

    // Backward difference at end
    derivatives[n - 1] = (coord_values[n - 1] - coord_values[n - 2]) / (grid_spacing[n - 1] - grid_spacing[n - 2]);

    return derivatives;
  }
}
} // namespace

template <GridConfigType NumericalConfig>
constexpr auto BoundaryLayerGrid::create_downstream_grid(
    NumericalConfig&& numerical_config, const io::OuterEdgeConfig& edge_config,
    const thermophysics::MixtureInterface& mixture) -> std::expected<BoundaryLayerGrid, GridError> {

  auto grid = BoundaryLayerGrid(numerical_config.n_eta, numerical_config.eta_max);

  grid.generate_eta_distribution();

  if (auto xi_result = grid.generate_xi_distribution(edge_config, mixture); !xi_result) {
    return std::unexpected(xi_result.error());
  }

  // xi_output is now the same as xi since output is at same locations as edge_points

  return grid;
}

auto BoundaryLayerGrid::generate_xi_distribution(const io::OuterEdgeConfig& edge_config,
                                                 const thermophysics::MixtureInterface& mixture)
    -> std::expected<void, GridError> {

  constexpr auto min_points = 2;
  if (edge_config.edge_points.size() < min_points) {
    return std::unexpected(GridError("Need at least 2 edge points for downstream grid"));
  }

  auto extract_property = [&edge_config](auto member_ptr) {
    auto range =
        edge_config.edge_points | std::views::transform([member_ptr](const auto& point) { return point.*member_ptr; });

    std::vector<double> result;
    std::ranges::copy(range, std::back_inserter(result));
    return result;
  };

  auto x_edge = extract_property(&io::OuterEdgeConfig::EdgePoint::x);
  auto u_edge = extract_property(&io::OuterEdgeConfig::EdgePoint::velocity);
  auto r_body = extract_property(&io::OuterEdgeConfig::EdgePoint::radius);

  std::vector<double> rho_edge, mu_edge;
  rho_edge.reserve(edge_config.edge_points.size());
  mu_edge.reserve(edge_config.edge_points.size());

  for (const auto& point : edge_config.edge_points) {
    auto eq_result = mixture.equilibrium_composition(point.temperature, point.pressure);
    if (!eq_result) {
      return std::unexpected(GridError("Failed to compute equilibrium composition for grid generation"));
    }
    auto species_fractions = eq_result.value();

    auto mw_result = mixture.mixture_molecular_weight(species_fractions);
    if (!mw_result) {
      return std::unexpected(GridError("Failed to compute MW for grid generation"));
    }
    auto density = point.pressure * mw_result.value() / (point.temperature * blast::constants::physical::universal_gas_constant);
    rho_edge.push_back(density);

    auto visc_result = mixture.viscosity(species_fractions, point.temperature, point.pressure);
    if (!visc_result) {
      return std::unexpected(GridError("Failed to compute viscosity for grid generation"));
    }
    mu_edge.push_back(visc_result.value());
  }

  auto xi_result = coordinate_transform::compute_xi_from_integration(x_edge, rho_edge, mu_edge, u_edge, r_body);

  if (!xi_result) {
    return std::unexpected(
        GridError(std::format("Failed to compute xi distribution: {}", xi_result.error().message())));
  }

  xi_ = std::move(xi_result.value());
  return {};
}


auto BoundaryLayerGrid::find_xi_interval(double xi_target) const noexcept
    -> std::expected<std::pair<int, int>, GridError> {
  auto result = coordinate_transform::search_interval(xi_, xi_target);
  if (!result) {
    return std::unexpected(GridError(result.error().message()));
  }
  return result.value();
}


// Explicit instantiations for common use cases
template auto BoundaryLayerGrid::create_downstream_grid<const io::NumericalConfig&>(
    const io::NumericalConfig& numerical_config, const io::OuterEdgeConfig& edge_config,
    const thermophysics::MixtureInterface& mixture) -> std::expected<BoundaryLayerGrid, GridError>;

} // namespace blast::boundary_layer::grid