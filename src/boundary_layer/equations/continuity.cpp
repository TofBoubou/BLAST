#include "blast/boundary_layer/equations/continuity.hpp"
#include "blast/boundary_layer/grid/coordinate_transform.hpp"
#include <algorithm>
#include <ranges>

namespace blast::boundary_layer::equations {

template<NumericRange VField>
constexpr auto solve_continuity(
    VField&& d_velocity_field,
    PhysicalQuantity auto d_eta,
    PhysicalQuantity auto initial_value
) noexcept -> std::vector<double>
requires std::ranges::sized_range<VField> {
    
    return grid::coordinate_transform::simpson_integrate(
        std::forward<VField>(d_velocity_field), d_eta, initial_value
    );
}

// Explicit instantiations for common use cases
template auto solve_continuity(std::span<const double>&&, double, double) noexcept -> std::vector<double>;
template auto solve_continuity(std::vector<double>&, double, double) noexcept -> std::vector<double>;

} // namespace blast::boundary_layer::equations