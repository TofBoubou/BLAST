
#include "blast/boundary_layer/grid/grid.hpp"
#include <algorithm>
#include <numeric>
#include <ranges>

namespace blast::boundary_layer::grid {

constexpr BoundaryLayerGrid::BoundaryLayerGrid(int n_eta, double eta_max) noexcept
    : n_eta_(n_eta), eta_max_(eta_max), d_eta_(eta_max / static_cast<double>(n_eta - 1)) {
    
    xi_.reserve(constants::default_grid_reserve);
    xi_output_.reserve(constants::default_grid_reserve);
    eta_.reserve(n_eta_);
}

template<GridConfigType NumericalConfig>
constexpr auto BoundaryLayerGrid::create_stagnation_grid(
    NumericalConfig&& numerical_config,
    const io::OuterEdgeConfig& edge_config
) -> std::expected<BoundaryLayerGrid, GridError> {
    
    auto grid = BoundaryLayerGrid(numerical_config.n_eta, numerical_config.eta_max);
    
    grid.generate_eta_distribution();
    grid.xi_.emplace_back(0.0);  // Stagnation point
    grid.xi_output_.emplace_back(0.0);
    
    return grid;
}

template<GridConfigType NumericalConfig>
constexpr auto BoundaryLayerGrid::create_downstream_grid(
    NumericalConfig&& numerical_config,
    const io::OuterEdgeConfig& edge_config,
    const io::OutputConfig& output_config
) -> std::expected<BoundaryLayerGrid, GridError> {
    
    auto grid = BoundaryLayerGrid(numerical_config.n_eta, numerical_config.eta_max);
    
    grid.generate_eta_distribution();
    
    if (auto xi_result = grid.generate_xi_distribution(edge_config); !xi_result) {
        return std::unexpected(xi_result.error());
    }
    
    // Extract x coordinates using modern ranges
    auto x_edge = edge_config.edge_points 
                | std::views::transform([](const auto& point) { return point.x; })
                | std::ranges::to<std::vector>();
    
    if (auto output_result = grid.generate_xi_output_distribution(output_config, x_edge); !output_result) {
        return std::unexpected(output_result.error());
    }
    
    return grid;
}

constexpr auto BoundaryLayerGrid::generate_eta_distribution() noexcept -> void {
    eta_.clear();
    
    // Modern range-based generation
    auto eta_indices = std::views::iota(0, n_eta_);
    std::ranges::transform(eta_indices, std::back_inserter(eta_),
        [this](int i) constexpr { return static_cast<double>(i) * d_eta_; });
}

auto BoundaryLayerGrid::generate_xi_distribution(const io::OuterEdgeConfig& edge_config) 
    -> std::expected<void, GridError> {
    
    constexpr auto min_points = 2;
    if (edge_config.edge_points.size() < min_points) {
        return std::unexpected(GridError("Need at least 2 edge points for downstream grid"));
    }
    
    // Extract data using modern transformations
    auto extract_property = [&edge_config](auto member_ptr) {
        return edge_config.edge_points 
             | std::views::transform([member_ptr](const auto& point) { return point.*member_ptr; })
             | std::ranges::to<std::vector>();
    };
    
    auto x_edge = extract_property(&io::OuterEdgeConfig::EdgePoint::x);
    auto rho_edge = extract_property(&io::OuterEdgeConfig::EdgePoint::density);
    auto mu_edge = extract_property(&io::OuterEdgeConfig::EdgePoint::viscosity);
    auto u_edge = extract_property(&io::OuterEdgeConfig::EdgePoint::velocity);
    auto r_body = extract_property(&io::OuterEdgeConfig::EdgePoint::radius);
    
    auto xi_result = coordinate_transform::compute_xi_from_integration(
        x_edge, rho_edge, mu_edge, u_edge, r_body
    );
    
    if (!xi_result) {
        return std::unexpected(GridError(
            std::format("Failed to compute xi distribution: {}", xi_result.error().message())
        ));
    }
    
    xi_ = std::move(xi_result.value());
    return {};
}

auto BoundaryLayerGrid::generate_xi_output_distribution(
    const io::OutputConfig& output_config,
    std::span<const double> x_edge
) -> std::expected<void, GridError> {
    
    xi_output_.clear();
    xi_output_.reserve(output_config.x_stations.size());
    
    for (const auto x_out : output_config.x_stations) {
        auto interval_result = coordinate_transform::search_interval(x_edge, x_out);
        if (!interval_result) {
            return std::unexpected(GridError(
                std::format("Output station x={} is outside edge domain", x_out)
            ));
        }
        
        const auto [i1, i2] = interval_result.value();
        const auto xi_out = (i1 == i2) ? xi_[i1] : 
            coordinate_transform::linear_interpolate(x_out, x_edge[i1], x_edge[i2], xi_[i1], xi_[i2]);
        
        xi_output_.emplace_back(xi_out);
    }
    
    return {};
}

auto BoundaryLayerGrid::find_xi_interval(double xi_target) const noexcept
    -> std::expected<std::pair<int, int>, GridError> {
    return coordinate_transform::search_interval(xi_, xi_target);
}

template<NumericRange XGrid>
auto BoundaryLayerGrid::interpolate_x_from_xi(double xi_target, XGrid&& x_grid) const 
    -> std::expected<double, GridError> {
    
    auto interval_result = find_xi_interval(xi_target);
    if (!interval_result) {
        return std::unexpected(interval_result.error());
    }
    
    const auto [i1, i2] = interval_result.value();
    return (i1 == i2) ? x_grid[i1] : 
        coordinate_transform::linear_interpolate(xi_target, xi_[i1], xi_[i2], x_grid[i1], x_grid[i2]);
}

constexpr auto compute_xi_step_size(
    double current_xi, double d_xi, int iterations,
    const io::NumericalConfig::StepControl& step_control
) noexcept -> double {
    
    if (iterations <= step_control.lower_bound) {
        return d_xi + d_xi / static_cast<double>(iterations);
    }
    
    if (iterations >= step_control.upper_bound) {
        return d_xi * static_cast<double>(step_control.upper_bound) / static_cast<double>(iterations);
    }
    
    return d_xi;
}

constexpr auto reduce_step_size(double current_xi, double& d_xi) noexcept -> double {
    const auto new_xi = current_xi - d_xi;
    d_xi *= constants::step_reduction_factor;
    return new_xi + d_xi;
}

// Explicit template instantiations for common types
template auto BoundaryLayerGrid::create_stagnation_grid(const io::NumericalConfig&, const io::OuterEdgeConfig&);
template auto BoundaryLayerGrid::create_downstream_grid(const io::NumericalConfig&, const io::OuterEdgeConfig&, const io::OutputConfig&);

} // namespace blast::boundary_layer::grid