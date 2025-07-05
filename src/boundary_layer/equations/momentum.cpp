#include "blast/boundary_layer/equations/momentum.hpp"
#include "blast/boundary_layer/solvers/tridiagonal_solvers.hpp"
#include <algorithm>
#include <format>

namespace blast::boundary_layer::equations {

auto solve_momentum(
    std::span<const double> F_previous,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    const coefficients::XiDerivatives& xi_der,
    std::span<const double> V_field,
    PhysicalQuantity auto xi,
    PhysicalQuantity auto d_eta
) -> std::expected<std::vector<double>, EquationError> {
    
    const auto n_eta = F_previous.size();
    
    if (n_eta != V_field.size()) {
        return std::unexpected(EquationError(
            "Momentum: F_previous and V_field size mismatch: {} vs {}", 
            std::source_location::current(), n_eta, V_field.size()
        ));
    }
    
    if (n_eta != coeffs.transport.l0.size()) {
        return std::unexpected(EquationError(
            "Momentum: incompatible coefficient sizes"
        ));
    }
    
    // Build tridiagonal system coefficients
    auto momentum_coeffs = detail::build_momentum_coefficients(
        F_previous, coeffs, bc, xi_der, V_field, xi, d_eta
    );
    
    // Get boundary conditions
    auto boundary_conds = detail::get_momentum_boundary_conditions();
    
    // Solve tridiagonal system using existing solver
    auto solution_result = solvers::solve_momentum_energy_tridiagonal(
        F_previous,
        boundary_conds.f_bc,
        boundary_conds.g_bc, 
        boundary_conds.h_bc,
        momentum_coeffs.a,
        momentum_coeffs.b,
        momentum_coeffs.c,
        momentum_coeffs.d
    );
    
    if (!solution_result) {
        return std::unexpected(EquationError(
            "Momentum: tridiagonal solver failed: {}", 
            std::source_location::current(), solution_result.error().message()
        ));
    }
    
    return solution_result.value();
}

namespace detail {

auto build_momentum_coefficients(
    std::span<const double> F_previous,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    const coefficients::XiDerivatives& xi_der,
    std::span<const double> V_field,
    PhysicalQuantity auto xi,
    PhysicalQuantity auto d_eta
) noexcept -> MomentumCoefficients {
    
    const auto n_eta = F_previous.size();
    const double d_eta_sq = d_eta * d_eta;
    const double lambda0 = xi_der.lambda0();
    const auto F_derivatives = xi_der.F_derivative();
    
    MomentumCoefficients coeffs_out;
    coeffs_out.a.reserve(n_eta);
    coeffs_out.b.reserve(n_eta);
    coeffs_out.c.reserve(n_eta);
    coeffs_out.d.reserve(n_eta);
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        // a[i] = l0[i] / d_eta²
        coeffs_out.a.push_back(coeffs.transport.l0[i] / d_eta_sq);
        
        // b[i] = (dl0_deta[i] - V[i]) / d_eta 
        coeffs_out.b.push_back((coeffs.transport.dl0_deta[i] - V_field[i]) / d_eta);
        
        // c[i] = -2*xi*lambda0*F[i]
        coeffs_out.c.push_back(-2.0 * xi * lambda0 * F_previous[i]);
        
        // d[i] = -beta*(rho_e/rho[i] - F[i]^2) + 2*xi*F[i]*F_der[i]
        const double d_term = -bc.beta * (bc.rho_e() / coeffs.thermodynamic.rho[i] - F_previous[i] * F_previous[i]) + 
                             2.0 * xi * F_previous[i] * F_derivatives[i];
        coeffs_out.d.push_back(d_term);
    }
    
    return coeffs_out;
}

} // namespace detail

// Explicit instantiations for common use cases
template auto solve_momentum<double, double>(
    std::span<const double> F_previous,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    const coefficients::XiDerivatives& xi_der,
    std::span<const double> V_field,
    double xi,
    double d_eta
) -> std::expected<std::vector<double>, EquationError>;

} // namespace blast::boundary_layer::equations