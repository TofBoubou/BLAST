#pragma once
#include "equation_types.hpp"
#include "geometry_factors.hpp"
#include "../coefficients/coefficient_types.hpp"
#include "../coefficients/xi_derivatives.hpp"
#include "../conditions/boundary_conditions.hpp"
#include "../../io/config_types.hpp"
#include <span>
#include <expected>

namespace blast::boundary_layer::equations {

// Solve momentum equation: tridiagonal system for F (dimensionless stream function)
// Equation: d/dη[l0 * dF/dη] - V * dF/dη - 2ξλ₀FdF/dξ + d[i] = 0
[[nodiscard]] auto solve_momentum(
    std::span<const double> F_previous,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    const coefficients::XiDerivatives& xi_der,
    std::span<const double> V_field,
    PhysicalQuantity auto xi,
    PhysicalQuantity auto d_eta
) -> std::expected<std::vector<double>, EquationError>;

// Build tridiagonal coefficients for momentum equation
namespace detail {
    
    struct MomentumCoefficients {
        std::vector<double> a;  // Lower diagonal
        std::vector<double> b;  // Main diagonal  
        std::vector<double> c;  // Upper diagonal
        std::vector<double> d;  // Right-hand side
    };
    
    [[nodiscard]] auto build_momentum_coefficients(
        std::span<const double> F_previous,
        const coefficients::CoefficientSet& coeffs,
        const conditions::BoundaryConditions& bc,
        const coefficients::XiDerivatives& xi_der,
        std::span<const double> V_field,
        PhysicalQuantity auto xi,
        PhysicalQuantity auto d_eta
    ) noexcept -> MomentumCoefficients;
    
    // Momentum-specific boundary conditions
    struct MomentumBoundaryConditions {
        double f_bc = 0.0; 
        double g_bc = 1.0;  
        double h_bc = 0.0;
    };
    
    [[nodiscard]] constexpr auto get_momentum_boundary_conditions() noexcept 
        -> MomentumBoundaryConditions {
        return MomentumBoundaryConditions{};
    }
    
} // namespace detail

} // namespace blast::boundary_layer::equations