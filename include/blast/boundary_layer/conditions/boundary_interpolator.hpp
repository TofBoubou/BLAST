#pragma once
#include "boundary_conditions.hpp"
#include "../grid/coordinate_transform.hpp"
#include "../../io/config_types.hpp"
#include <span>

namespace blast::boundary_layer::conditions {

// Compute beta based on body type and station
[[nodiscard]] constexpr auto compute_beta(
    int station,
    double xi,
    const io::SimulationConfig& sim_config,
    double ue = 0.0,
    double d_ue_dxi = 0.0
) noexcept -> double;

// Main interpolation function
[[nodiscard]] auto interpolate_boundary_conditions(
    int station,
    double xi,
    std::span<const double> xi_grid,
    const io::OuterEdgeConfig& edge_config,
    const io::WallParametersConfig& wall_config,
    const io::SimulationConfig& sim_config
) -> std::expected<BoundaryConditions, BoundaryConditionError>;

// Specialized for stagnation point
[[nodiscard]] auto create_stagnation_conditions(
    const io::OuterEdgeConfig& edge_config,
    const io::WallParametersConfig& wall_config,
    const io::SimulationConfig& sim_config
) -> std::expected<BoundaryConditions, BoundaryConditionError>;


} // namespace blast::boundary_layer::conditions