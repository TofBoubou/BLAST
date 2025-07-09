#pragma once
#include "equation_types.hpp"
#include <vector>
#include <span>
#include <expected>

namespace blast::boundary_layer::equations {

// Simple integration for continuity equation
// Solves: dV/dη = -right_handed where V is the velocity-like variable
template<NumericRange VField>
[[nodiscard]] constexpr auto solve_continuity(
    VField&& d_velocity_field,
    PhysicalQuantity auto d_eta,
    PhysicalQuantity auto initial_value = 0.0
) noexcept -> std::vector<double>
requires std::ranges::sized_range<VField>;

} // namespace blast::boundary_layer::equations