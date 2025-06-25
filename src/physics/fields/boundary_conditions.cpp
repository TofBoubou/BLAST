#include "blast/physics/fields/boundary_conditions.hpp"
#include <algorithm>
#include <ranges>
#include <cmath>

namespace blast::physics::fields {

BoundaryConditionInterpolator::BoundaryConditionInterpolator(
    const io::OuterEdgeConfig& edge_config,
    const io::WallParametersConfig& wall_config
) : edge_config_(edge_config), wall_config_(wall_config) {
    
    x_coordinates_.reserve(edge_config_.edge_points.size());
    
    std::ranges::transform(edge_config_.edge_points, 
                          std::back_inserter(x_coordinates_),
                          [](const auto& point) { return point.x; });
}

auto BoundaryConditionInterpolator::find_interval(double x_target) const noexcept
    -> std::expected<std::pair<std::size_t, std::size_t>, core::BlastException> {
    
    if (x_coordinates_.empty()) {
        return std::unexpected(core::BlastException("Empty edge coordinates"));
    }
    
    // Check exact match
    constexpr auto epsilon = 1e-15;
    if (auto exact_it = std::ranges::find_if(x_coordinates_, 
            [x_target](auto x) { return std::abs(x - x_target) < epsilon; });
        exact_it != std::ranges::end(x_coordinates_)) {
        
        const auto index = static_cast<std::size_t>(
            std::ranges::distance(std::ranges::begin(x_coordinates_), exact_it)
        );
        return std::make_pair(index, index);
    }
    
    // Find interval
    auto upper_it = std::ranges::upper_bound(x_coordinates_, x_target);
    
    if (upper_it == std::ranges::begin(x_coordinates_)) {
        return std::unexpected(core::BlastException("x_target below domain minimum"));
    }
    
    if (upper_it == std::ranges::end(x_coordinates_)) {
        return std::unexpected(core::BlastException("x_target above domain maximum"));
    }
    
    const auto i2 = static_cast<std::size_t>(
        std::ranges::distance(std::ranges::begin(x_coordinates_), upper_it)
    );
    const auto i1 = i2 - 1;
    
    return std::make_pair(i1, i2);
}

constexpr auto BoundaryConditionInterpolator::linear_interpolate(
    double x, double x1, double x2, double y1, double y2
) const noexcept -> double {
    
    constexpr auto epsilon = 1e-15;
    return (std::abs(x2 - x1) < epsilon) ? y1 : y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

auto BoundaryConditionInterpolator::compute_edge_properties_at_x(double x_station) const
    -> std::expected<EdgeProperties, core::BlastException> {
    
    auto interval_result = find_interval(x_station);
    if (!interval_result) {
        return std::unexpected(interval_result.error());
    }
    
    const auto [i1, i2] = interval_result.value();
    const auto& point1 = edge_config_.edge_points[i1];
    
    if (i1 == i2) {
        // Exact match
        EdgeProperties props(point1.species_fractions.size());
        props.pressure = point1.pressure;
        props.temperature = point1.enthalpy; // Will need thermodynamic conversion later
        props.velocity = point1.velocity;
        props.enthalpy = point1.enthalpy;
        props.density = point1.density;
        props.viscosity = point1.viscosity;
        
        std::ranges::copy(point1.species_fractions, props.species_fractions.begin());
        return props;
    }
    
    // Interpolation
    const auto& point2 = edge_config_.edge_points[i2];
    const auto x1 = x_coordinates_[i1];
    const auto x2 = x_coordinates_[i2];
    
    EdgeProperties props(point1.species_fractions.size());
    props.pressure = linear_interpolate(x_station, x1, x2, point1.pressure, point2.pressure);
    props.velocity = linear_interpolate(x_station, x1, x2, point1.velocity, point2.velocity);
    props.enthalpy = linear_interpolate(x_station, x1, x2, point1.enthalpy, point2.enthalpy);
    props.density = linear_interpolate(x_station, x1, x2, point1.density, point2.density);
    props.viscosity = linear_interpolate(x_station, x1, x2, point1.viscosity, point2.viscosity);
    
    for (std::size_t i = 0; i < point1.species_fractions.size(); ++i) {
        props.species_fractions[i] = linear_interpolate(
            x_station, x1, x2, point1.species_fractions[i], point2.species_fractions[i]
        );
    }
    
    props.temperature = props.enthalpy; // Placeholder - will need proper conversion
    
    return props;
}

auto BoundaryConditionInterpolator::compute_wall_properties_at_x(double x_station) const
    -> std::expected<WallProperties, core::BlastException> {
    
    if (wall_config_.wall_temperatures.empty()) {
        return std::unexpected(core::BlastException("No wall temperatures provided"));
    }
    
    // For now, simple: use first temperature if only one, or interpolate if multiple
    if (wall_config_.wall_temperatures.size() == 1) {
        return WallProperties(wall_config_.wall_temperatures[0]);
    }
    
    // If multiple temperatures, need corresponding x coordinates - placeholder for now
    return WallProperties(wall_config_.wall_temperatures[0]);
}

constexpr auto BoundaryConditionInterpolator::x_min() const noexcept -> double {
    return x_coordinates_.empty() ? 0.0 : x_coordinates_.front();
}

constexpr auto BoundaryConditionInterpolator::x_max() const noexcept -> double {
    return x_coordinates_.empty() ? 0.0 : x_coordinates_.back();
}

} // namespace blast::physics::fields