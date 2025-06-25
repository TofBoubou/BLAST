
#pragma once
#include "../../core/exceptions.hpp"
#include "../../io/config_types.hpp"
#include <vector>
#include <span>
#include <expected>
#include <concepts>
#include <ranges>
#include <algorithm>

namespace blast::boundary_layer::grid::coordinate_transform {

// Simpson integration constants (no more magic numbers!)
namespace simpson_constants {
    constexpr double coeff_5pt_1 = 17.0;
    constexpr double coeff_5pt_2 = 42.0;
    constexpr double coeff_5pt_3 = -16.0;
    constexpr double coeff_5pt_4 = 6.0;
    constexpr double coeff_5pt_5 = -1.0;
    constexpr double divisor_5pt = 48.0;
    
    constexpr double coeff_3pt_1 = 1.0;
    constexpr double coeff_3pt_2 = 4.0;
    constexpr double coeff_3pt_3 = 1.0;
    constexpr double divisor_3pt = 3.0;
}

// Modern concepts for type safety
template<typename T>
concept PhysicalQuantity = std::floating_point<T> && requires(T t) {
    { t >= T{0} } -> std::convertible_to<bool>;
};

template<typename Range>
concept PhysicalRange = std::ranges::sized_range<Range> && 
                       PhysicalQuantity<std::ranges::range_value_t<Range>>;


// Core transformation functions with concepts
template<PhysicalRange EdgeData>
[[nodiscard]] auto compute_xi_from_integration(
    EdgeData&& x_edge, // && can be used only in a template context. Outside a template && is only for rvalue. With a template it's an universal reference and && can be the 2
    EdgeData&& rho_edge,
    EdgeData&& mu_edge,
    EdgeData&& u_edge,
    EdgeData&& r_body
) -> std::expected<std::vector<double>, TransformError>;

template<PhysicalRange RhoData>
[[nodiscard]] auto compute_physical_y_from_eta(
    PhysicalQuantity auto eta,
    PhysicalQuantity auto xi,
    PhysicalQuantity auto rho_wall,
    PhysicalQuantity auto K_bl,
    RhoData&& rho_eta,
    PhysicalQuantity auto d_eta,
    const io::SimulationConfig& sim_config,
    PhysicalQuantity auto d_ue_dx = 0.0,
    PhysicalQuantity auto rho_e = 0.0,
    PhysicalQuantity auto mu_e = 0.0,
    PhysicalQuantity auto u_e = 0.0,
    PhysicalQuantity auto r_body = 1.0
) -> std::expected<double, TransformError>;

[[nodiscard]] constexpr auto compute_derivative_factor(
    int station,
    PhysicalQuantity auto xi,
    const io::SimulationConfig& sim_config,
    PhysicalQuantity auto d_ue_dx,
    PhysicalQuantity auto rho_e,
    PhysicalQuantity auto mu_e,
    PhysicalQuantity auto u_e = 0.0,
    PhysicalQuantity auto r_body = 1.0
) noexcept -> double;

// Interpolation utilities with constexpr
[[nodiscard]] constexpr auto linear_interpolate(
    std::floating_point auto x,
    std::floating_point auto x1, 
    std::floating_point auto x2,
    std::floating_point auto y1, 
    std::floating_point auto y2
) noexcept -> double;

template<std::ranges::random_access_range Grid>
[[nodiscard]] auto search_interval(Grid&& grid, std::floating_point auto target) noexcept
    -> std::expected<std::pair<int, int>, TransformError>;

// Numerical integration with concepts
template<PhysicalRange Values> 
[[nodiscard]] auto simpson_integrate(
    Values&& f_values,
    PhysicalQuantity auto dx,
    PhysicalQuantity auto f0 = 0.0
) -> std::vector<double>;

} // namespace blast::boundary_layer::grid::coordinate_transform