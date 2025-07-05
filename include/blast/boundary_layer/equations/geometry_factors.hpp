#pragma once
#include "equation_types.hpp"
#include "../conditions/boundary_conditions.hpp"
#include "../../io/config_types.hpp"
#include <cmath>

namespace blast::boundary_layer::equations {

// Compute geometry-dependent factors for equation coefficients
[[nodiscard]] constexpr auto compute_geometry_factors(
    int station,
    PhysicalQuantity auto xi,
    const conditions::BoundaryConditions& bc,
    const io::SimulationConfig& sim_config
) noexcept -> GeometryFactors {
    
    if (station == 0) {
        // Stagnation point - factors depend on body type
        switch (sim_config.body_type) {
            case io::SimulationConfig::BodyType::Axisymmetric: {
                const double sqrt_term = std::sqrt(2.0 * bc.rho_e() * bc.mu_e() * bc.d_ue_dx());
                return GeometryFactors(
                    1.0 / sqrt_term,                    // Diffusion flux factor
                    1.0 / (2.0 * bc.d_ue_dx()),        // Chemical production factor (missing /rho)
                    std::sqrt(bc.rho_e() * bc.mu_e() / (2.0 * bc.d_ue_dx()))  // BC factor
                );
            }
            case io::SimulationConfig::BodyType::TwoD: {
                const double sqrt_term = std::sqrt(bc.rho_e() * bc.mu_e() * bc.d_ue_dx());
                return GeometryFactors(
                    1.0 / sqrt_term,                    
                    1.0 / bc.d_ue_dx(),                // Missing /rho
                    std::sqrt(bc.rho_e() * bc.mu_e() / bc.d_ue_dx())
                );
            }
            case io::SimulationConfig::BodyType::Cone:
            case io::SimulationConfig::BodyType::FlatPlate:
                return GeometryFactors(
                    1.0,
                    0.0,   // No chemical production for these geometries
                    0.0   // No boundary condition factor
                );
        }
    } else {
        // Downstream station - general expressions
        return GeometryFactors(
            bc.r_body() * std::sqrt(2.0 * xi) / bc.d_xi_dx(),
            2.0 * xi / (bc.ue() * bc.d_xi_dx()),       // Missing /rho  
            std::sqrt(2.0 * xi) / (bc.ue() * bc.r_body())
        );
    }
    
    // Default fallback (should never reach here)
    return GeometryFactors(1.0, 0.0, 0.0);
}

// Compute J_fact for energy equation (slightly different from species)
[[nodiscard]] constexpr auto compute_energy_j_factor(
    int station,
    PhysicalQuantity auto xi,
    const conditions::BoundaryConditions& bc,
    const io::SimulationConfig& sim_config
) noexcept -> double {
    
    if (station == 0) {
        // Stagnation point
        switch (sim_config.body_type) {
            case io::SimulationConfig::BodyType::Axisymmetric:
                return 1.0 / std::sqrt(2.0 * bc.rho_e() * bc.mu_e() * bc.d_ue_dx());
            case io::SimulationConfig::BodyType::TwoD:
                return 1.0 / std::sqrt(bc.rho_e() * bc.mu_e() * bc.d_ue_dx());
            case io::SimulationConfig::BodyType::Cone:
            case io::SimulationConfig::BodyType::FlatPlate:
                return 1.0;
        }
    } else {
        // Downstream station
        return std::sqrt(2.0 * xi) * bc.r_body() / bc.d_xi_dx();
    }
    
    return 1.0; // Default fallback
}

} // namespace blast::boundary_layer::equations