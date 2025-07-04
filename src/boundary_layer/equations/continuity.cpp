#include "blast/boundary_layer/equations/continuity.hpp"
#include <algorithm>
#include <ranges>

namespace blast::boundary_layer::equations {

auto solve_continuity(
    std::span<const double> y_coordinates,
    std::span<const double> velocity_field,
    double d_eta,
    double initial_value
) noexcept -> std::vector<double> {
    
    return detail::integrate_trapezoidal(velocity_field, d_eta, initial_value);
}

namespace detail {

template<NumericRange InputRange>
constexpr auto integrate_trapezoidal(
    InputRange&& values,
    PhysicalQuantity auto dx,
    PhysicalQuantity auto initial_value
) noexcept -> std::vector<double> {
    
    const auto n = std::ranges::size(values);
    std::vector<double> result;
    result.reserve(n);
    
    if (n == 0) return result;
    
    result.push_back(initial_value);
    
    if (n == 1) return result;
    
    // Trapezoidal integration: ∫f(x)dx ≈ Σ[(f[i] + f[i+1])/2 * dx]
    auto values_it = std::ranges::begin(values);
    double prev_value = *values_it;
    ++values_it;
    
    for (std::size_t i = 1; i < n; ++i) {
        const double current_value = *values_it;
        const double integral_step = (prev_value + current_value) * 0.5 * dx;
        result.push_back(result.back() + integral_step);
        
        prev_value = current_value;
        ++values_it;
    }
    
    return result;
}

// Explicit instantiation for common types
template auto integrate_trapezoidal(const std::vector<double>&, double, double) noexcept -> std::vector<double>;
template auto integrate_trapezoidal(std::span<const double>&&, double, double) noexcept -> std::vector<double>;

} // namespace detail

// Template instantiation for the main function
template<NumericRange YField, NumericRange VField>
constexpr auto solve_continuity(
    YField&& y_coordinates,
    VField&& velocity_field,
    PhysicalQuantity auto d_eta,
    PhysicalQuantity auto initial_value
) noexcept -> std::vector<double>
requires std::ranges::sized_range<YField> && std::ranges::sized_range<VField> {
    
    return detail::integrate_trapezoidal(
        std::forward<VField>(velocity_field), d_eta, initial_value
    );
}

// Explicit instantiations for common use cases
template auto solve_continuity(const std::vector<double>&, const std::vector<double>&, double, double) noexcept -> std::vector<double>;
template auto solve_continuity(std::span<const double>&&, std::span<const double>&&, double, double) noexcept -> std::vector<double>;

} // namespace blast::boundary_layer::equations