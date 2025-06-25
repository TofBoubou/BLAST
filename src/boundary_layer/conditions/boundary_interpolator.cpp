#include "blast/boundary_layer/conditions/boundary_interpolator.hpp"
#include <algorithm>
#include <format>
#include <expected>

namespace blast::boundary_layer::conditions {

using blast::boundary_layer::grid::coordinate_transform::search_interval;
using blast::boundary_layer::grid::coordinate_transform::linear_interpolate;

namespace {
    
template<typename Container>
[[nodiscard]] auto interpolate_property(
    const Container& values,
    std::size_t i1,
    std::size_t i2,
    const std::vector<double>& x_grid,
    double x_target
) -> double {
    
    if (i1 == i2) return values[i1];
    
    return linear_interpolate(x_target, 
                             x_grid[i1], x_grid[i2],
                             values[i1], values[i2]);
}

template<typename MemberPtr>
[[nodiscard]] auto extract_edge_property(
    const std::vector<io::OuterEdgeConfig::EdgePoint>& points,
    MemberPtr member_ptr
) -> std::vector<double> {
    
    std::vector<double> result;
    result.reserve(points.size());
    
    std::ranges::transform(points, std::back_inserter(result),
        [member_ptr](const auto& point) { return point.*member_ptr; });
    
    return result;
}

} 

constexpr auto compute_beta(
    int station,
    double xi,
    const io::SimulationConfig& sim_config,
    double ue,
    double d_ue_dxi
) noexcept -> double {
    
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

auto create_stagnation_conditions(
    const io::OuterEdgeConfig& edge_config,
    const io::WallParametersConfig& wall_config,
    const io::SimulationConfig& sim_config
) -> std::expected<BoundaryConditions, BoundaryConditionError> {
    
    if (edge_config.edge_points.empty()) {
        return std::unexpected(BoundaryConditionError("No edge points defined"));
    }
    
    if (wall_config.wall_temperatures.empty()) {
        return std::unexpected(BoundaryConditionError("No wall temperatures defined"));
    }
    
    const auto& edge_point = edge_config.edge_points[0];
    
    EdgeConditions edge{
        .pressure = edge_point.pressure,
        .viscosity = edge_point.viscosity,
        .velocity = edge_point.velocity,
        .enthalpy = edge_point.enthalpy,
        .density = edge_point.density,
        .species_fractions = edge_point.species_fractions,
        .d_xi_dx = edge_config.velocity_gradient_stagnation,  
        .d_ue_dx = edge_config.velocity_gradient_stagnation,
        .d_he_dx = 0.0,  
        .d_he_dxi = 0.0,
        .body_radius = edge_point.radius
    };
    
    WallConditions wall{
        .temperature = wall_config.wall_temperatures[0],
    };
    
    return BoundaryConditions{
        .edge = std::move(edge),
        .wall = std::move(wall),
        .beta = compute_beta(0, 0.0, sim_config),
        .station = 0,
        .xi = 0.0
    };
}

auto interpolate_boundary_conditions(
    int station,
    double xi,
    std::span<const double> xi_grid,
    const io::OuterEdgeConfig& edge_config,
    const io::WallParametersConfig& wall_config,
    const io::SimulationConfig& sim_config
) -> std::expected<BoundaryConditions, BoundaryConditionError> {
    
    // Stagnation point special case
    if (station == 0) {
        return create_stagnation_conditions(edge_config, wall_config, sim_config);
    }
    
    // Find interpolation interval
    auto interval_result = search_interval(xi_grid, xi);
    if (!interval_result) {
        return std::unexpected(BoundaryConditionError(
            std::format("Failed to find xi={} in grid: {}", xi, interval_result.error().message())
        ));
    }
    
    const auto [i1, i2] = interval_result.value();
    
    // Extract x coordinates for interpolation  
    auto x_grid = extract_edge_property(edge_config.edge_points, 
                                       &io::OuterEdgeConfig::EdgePoint::x);
    
    // Compute interpolated x position
    const double x_interp = (i1 == i2) ? x_grid[i1] :
        linear_interpolate(xi, xi_grid[i1], xi_grid[i2], x_grid[i1], x_grid[i2]);
    
    // Find interval in physical x space
    auto x_interval_result = search_interval(std::span(x_grid), x_interp);
    if (!x_interval_result) {
        return std::unexpected(BoundaryConditionError(
            std::format("Failed to find x={} in edge grid", x_interp)
        ));
    }
    
    const auto [ix1, ix2] = x_interval_result.value();
    
    auto interp = [&](auto member_ptr) {
        auto values = extract_edge_property(edge_config.edge_points, member_ptr);
        return interpolate_property(values, ix1, ix2, x_grid, x_interp);
    };
    
    // Interpolate all edge properties
    EdgeConditions edge{
        .pressure = interp(&io::OuterEdgeConfig::EdgePoint::pressure),
        .viscosity = interp(&io::OuterEdgeConfig::EdgePoint::viscosity),
        .velocity = interp(&io::OuterEdgeConfig::EdgePoint::velocity),
        .enthalpy = interp(&io::OuterEdgeConfig::EdgePoint::enthalpy),
        .density = interp(&io::OuterEdgeConfig::EdgePoint::density),
        .species_fractions = {}, // Handle separately
        .d_xi_dx = 0.0,  // TODO: Compute from grid
        .d_ue_dx = 0.0,  // TODO: Compute derivative
        .d_he_dx = 0.0,  // TODO: Compute derivative
        .d_he_dxi = 0.0,
        .body_radius = interp(&io::OuterEdgeConfig::EdgePoint::radius)
    };
    
    // Handle species interpolation
    if (!edge_config.edge_points[0].species_fractions.empty()) {
        const auto n_species = edge_config.edge_points[0].species_fractions.size();
        edge.species_fractions.resize(n_species);
        
        for (std::size_t s = 0; s < n_species; ++s) {
            std::vector<double> species_values;
            for (const auto& point : edge_config.edge_points) {
                species_values.push_back(point.species_fractions[s]);
            }
            edge.species_fractions[s] = interpolate_property(
                species_values, ix1, ix2, x_grid, x_interp
            );
        }
    }
    
    // Interpolate wall temperature
    const double wall_temp = interpolate_property(
        wall_config.wall_temperatures, ix1, ix2, x_grid, x_interp
    );
    
    WallConditions wall{
        .temperature = wall_temp
    };
    
/*     // Compute derivatives (simplified for now)
    const double d_ue_dxi = edge.d_ue_dx / edge.d_xi_dx;
    edge.d_he_dxi = edge.d_he_dx / edge.d_xi_dx; */
    
    return BoundaryConditions{
        .edge = std::move(edge),
        .wall = std::move(wall),
        .beta = compute_beta(station, xi, sim_config, edge.velocity, d_ue_dxi),
        .station = station,
        .xi = xi
    };
}

} // namespace blast::boundary_layer::conditions