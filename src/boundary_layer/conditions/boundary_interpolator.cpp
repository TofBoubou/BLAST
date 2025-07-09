#include "blast/boundary_layer/conditions/boundary_interpolator.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include <algorithm>
#include <format>
#include <expected>

namespace blast::boundary_layer::conditions {

using blast::boundary_layer::grid::coordinate_transform::search_interval;
using blast::boundary_layer::grid::coordinate_transform::linear_interpolate;
using blast::boundary_layer::grid::coordinate_transform::hermite_interpolate;

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

// Overload with optional derivatives for Hermite interpolation
template<typename Container>
[[nodiscard]] auto interpolate_property(
    const Container& values,
    std::size_t i1,
    std::size_t i2,
    const std::vector<double>& x_grid,
    double x_target,
    const Container& derivatives
) -> double {
    
    if (i1 == i2) return values[i1];
    
    // Use Hermite interpolation when derivatives are available
    return hermite_interpolate(x_target,
                              x_grid[i1], x_grid[i2],
                              values[i1], values[i2],
                              derivatives[i1], derivatives[i2]);
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

// Compute derivatives reusing existing high-order function
[[nodiscard]] auto compute_property_derivatives(
    const std::vector<double>& values,
    const std::vector<double>& x_grid
) -> std::vector<double> {
    
    if (values.size() != x_grid.size() || values.size() < 2) {
        std::vector<double> derivatives(values.size(), 0.0);
        return derivatives;
    }
    
    // Check if grid is uniform (within tolerance)
    constexpr double tolerance = 1e-12;
    bool is_uniform = true;
    const double dx_first = x_grid[1] - x_grid[0];
    
    for (std::size_t i = 2; i < x_grid.size(); ++i) {
        const double dx_current = x_grid[i] - x_grid[i-1];
        if (std::abs(dx_current - dx_first) > tolerance) {
            is_uniform = false;
            break;
        }
    }
    
    if (is_uniform) {
        // Reuse existing high-order O(h^4) function for uniform grids
        return coefficients::derivatives::compute_eta_derivative(values, dx_first);
    } else {
        // Handle non-uniform grids with weighted central differences
        const auto n = values.size();
        std::vector<double> derivatives(n);
        
        // Forward difference at start
        derivatives[0] = (values[1] - values[0]) / (x_grid[1] - x_grid[0]);
        
        // Central differences in interior
        for (std::size_t i = 1; i < n - 1; ++i) {
            const double dx_left = x_grid[i] - x_grid[i-1];
            const double dx_right = x_grid[i+1] - x_grid[i];
            const double dx_total = x_grid[i+1] - x_grid[i-1];
            
            // Weighted central difference for non-uniform grids
            derivatives[i] = ((values[i+1] - values[i]) / dx_right * dx_left +
                             (values[i] - values[i-1]) / dx_left * dx_right) / dx_total;
        }
        
        // Backward difference at end
        derivatives[n-1] = (values[n-1] - values[n-2]) / (x_grid[n-1] - x_grid[n-2]);
        
        return derivatives;
    }
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
    
    // Lambda for standard linear interpolation
    auto interp_linear = [&](auto member_ptr) {
        auto values = extract_edge_property(edge_config.edge_points, member_ptr);
        return interpolate_property(values, ix1, ix2, x_grid, x_interp);
    };
    
    // Lambda for high-accuracy Hermite interpolation (for critical properties)
    auto interp_hermite = [&](auto member_ptr) {
        auto values = extract_edge_property(edge_config.edge_points, member_ptr);
        auto derivatives = compute_property_derivatives(values, x_grid);
        return interpolate_property(values, ix1, ix2, x_grid, x_interp, derivatives);
    };
    
    // Interpolate edge properties with appropriate method
    // Use Hermite for critical properties with strong gradients
    EdgeConditions edge{
        .pressure = interp_hermite(&io::OuterEdgeConfig::EdgePoint::pressure),    // Critical: strong gradients
        .viscosity = interp_linear(&io::OuterEdgeConfig::EdgePoint::viscosity),   // Linear sufficient
        .velocity = interp_hermite(&io::OuterEdgeConfig::EdgePoint::velocity),    // Critical: momentum boundary layer
        .enthalpy = interp_hermite(&io::OuterEdgeConfig::EdgePoint::enthalpy),    // Critical: energy boundary layer
        .density = interp_hermite(&io::OuterEdgeConfig::EdgePoint::density),      // Critical: compressible flow
        .species_fractions = {}, // Handle separately with Hermite below
        .d_xi_dx = 0.0,  // TODO: Compute from grid
        .d_ue_dx = 0.0,  // TODO: Compute derivative with Hermite
        .d_he_dx = 0.0,  // TODO: Compute derivative with Hermite
        .d_he_dxi = 0.0,
        .body_radius = interp_linear(&io::OuterEdgeConfig::EdgePoint::radius)     // Linear sufficient
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
            // Use Hermite for species fractions (critical for chemical accuracy)
            auto species_derivatives = compute_property_derivatives(species_values, x_grid);
            edge.species_fractions[s] = interpolate_property(
                species_values, ix1, ix2, x_grid, x_interp, species_derivatives
            );
        }
    }
    
    // Interpolate wall temperature using Hermite (critical for thermal boundary layer)
    auto wall_temp_derivatives = compute_property_derivatives(wall_config.wall_temperatures, x_grid);
    const double wall_temp = interpolate_property(
        wall_config.wall_temperatures, ix1, ix2, x_grid, x_interp, wall_temp_derivatives
    );
    
    WallConditions wall{
        .temperature = wall_temp
    };
    
    // Compute derivatives (simplified for now)
    const double d_ue_dxi = edge.d_ue_dx / edge.d_xi_dx;
    edge.d_he_dxi = edge.d_he_dx / edge.d_xi_dx;
    
    return BoundaryConditions{
        .edge = std::move(edge),
        .wall = std::move(wall),
        .beta = compute_beta(station, xi, sim_config, edge.velocity, d_ue_dxi),
        .station = station,
        .xi = xi
    };
}

} // namespace blast::boundary_layer::conditions