
#include "blast/boundary_layer/grid/coordinate_transform.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <ranges>
#include <format>

namespace blast::boundary_layer::grid::coordinate_transform {
using blast::core::TransformError;

template<PhysicalRange RhoData>
auto compute_physical_y_from_eta(
    PhysicalQuantity auto eta,
    PhysicalQuantity auto xi,
    PhysicalQuantity auto rho_wall,
    RhoData&& rho_eta,
    PhysicalQuantity auto d_eta,
    const io::SimulationConfig& sim_config,
    PhysicalQuantity auto d_ue_dx,
    PhysicalQuantity auto rho_e,
    PhysicalQuantity auto mu_e,
    PhysicalQuantity auto u_e,
    PhysicalQuantity auto r_body
) -> std::expected<double, TransformError> {
    
    const auto eta_index = static_cast<size_t>(eta / d_eta + 0.5);
    if (eta_index >= std::ranges::size(rho_eta)) {
        return std::unexpected(TransformError("Eta index out of bounds"));
    }
    
    auto rho_inv_range = rho_eta
                    | std::views::take(eta_index + 1)
                    | std::views::transform([](auto rho) { return 1.0 / rho; });

    std::vector<double> rho_inv_values;
    std::ranges::copy(rho_inv_range, std::back_inserter(rho_inv_values));
 
    auto integral_result = simpson_integrate(rho_inv_values, d_eta);
    if (integral_result.empty()) {
        return std::unexpected(TransformError("Failed to compute rho integral"));
    }
    
    const auto y_factor = (xi == 0.0) ? 
        compute_stagnation_y_factor(sim_config, d_ue_dx, rho_e, mu_e) :
        std::sqrt(2.0 * xi) / (u_e * r_body);
    
    return integral_result[eta_index] * y_factor;
}

constexpr auto compute_derivative_factor(
    int station,
    PhysicalQuantity auto xi,
    const io::SimulationConfig& sim_config,
    PhysicalQuantity auto d_ue_dx,
    PhysicalQuantity auto rho_e,
    PhysicalQuantity auto mu_e,
    PhysicalQuantity auto u_e,
    PhysicalQuantity auto r_body
) noexcept -> double {
    
    if (station == 0) {
        // Compile-time dispatch for stagnation point
        switch (sim_config.body_type) {
            case io::SimulationConfig::BodyType::Axisymmetric:
                return std::sqrt(2.0 * d_ue_dx / (rho_e * mu_e));
            case io::SimulationConfig::BodyType::TwoD:
                return std::sqrt(d_ue_dx / (rho_e * mu_e));
            case io::SimulationConfig::BodyType::Cone:
            case io::SimulationConfig::BodyType::FlatPlate:
                return 1.0 / std::sqrt(rho_e * mu_e);
        }
    }
    
    return u_e * r_body / std::sqrt(2.0 * xi);
}

// Helper function for stagnation point y-factor calculation
constexpr auto compute_stagnation_y_factor(
    const io::SimulationConfig& sim_config,
    PhysicalQuantity auto d_ue_dx,
    PhysicalQuantity auto rho_e,
    PhysicalQuantity auto mu_e
) noexcept -> double {
    
    switch (sim_config.body_type) {
        case io::SimulationConfig::BodyType::Axisymmetric:
            return std::sqrt(rho_e * mu_e / (2.0 * d_ue_dx));
        case io::SimulationConfig::BodyType::TwoD:
            return std::sqrt(rho_e * mu_e / d_ue_dx);
        case io::SimulationConfig::BodyType::Cone:
        case io::SimulationConfig::BodyType::FlatPlate:
            return 1.0;
    }
    return 1.0; // Default fallback
}


} // namespace blast::boundary_layer::grid::coordinate_transform