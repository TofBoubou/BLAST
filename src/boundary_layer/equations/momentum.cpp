#include "blast/boundary_layer/equations/momentum.hpp"
#include "blast/boundary_layer/solvers/tridiagonal_solvers.hpp"
#include <algorithm>
#include <format>
#include <iostream>
#include <iomanip>

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

    // std::cout << std::scientific << "From solve_momentum xi = " << xi << std::endl;
    
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

    constexpr double R_universal = 8314.462618;

    const double T_edge = bc.P_e() * coeffs.thermodynamic.MW[n_eta-1] / 
                         (coeffs.thermodynamic.rho[n_eta-1] * R_universal);

    const double rho_e_actual = bc.P_e() * coeffs.thermodynamic.MW[n_eta-1] / 
                               (T_edge * R_universal);
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        // a[i] = l0[i] / d_eta²
        coeffs_out.a.push_back(coeffs.transport.l0[i] / d_eta_sq);
/*         std::cout << "---------------------------------------------------------------" << std::endl;
        std::cout << "a dans momentum : " << coeffs_out.a[i] << std::endl; */
        
        // b[i] = (dl0_deta[i] - V[i]) / d_eta 
        coeffs_out.b.push_back((coeffs.transport.dl0_deta[i] - V_field[i]) / d_eta);
/*         std::cout << "b dans momentum : " << coeffs_out.b[i] << std::endl;
        std::cout << "coeffs.transport.dl0_deta : " << coeffs.transport.dl0_deta[i] << std::endl;
        std::cout << "V_field : " << V_field[i] << std::endl; */

        // c[i] = -2*xi*lambda0*F[i]
        coeffs_out.c.push_back(-2.0 * xi * lambda0 * F_previous[i]);
        // std::cout << "c dans momentum : " << coeffs_out.c[i] << std::endl;
        
        // d[i] = -beta*(rho_e/rho[i] - F[i]^2) + 2*xi*F[i]*F_der[i]
/*         const double d_term = - bc.beta * (bc.rho_e() / coeffs.thermodynamic.rho[i] - F_previous[i] * F_previous[i]) + 
                             2.0 * xi * F_previous[i] * F_derivatives[i];
        coeffs_out.d.push_back(d_term); */

/*         const double d_term_si_linearisation = - bc.beta * (bc.rho_e() / coeffs.thermodynamic.rho[i]) + 
                             2.0 * xi * F_previous[i] * F_derivatives[i]; */

        const double d_term = - bc.beta * (rho_e_actual / coeffs.thermodynamic.rho[i] - F_previous[i] * F_previous[i]) + 
                             2.0 * xi * F_previous[i] * F_derivatives[i];
        coeffs_out.d.push_back(d_term);

        const double d_term_si_linearisation = - bc.beta * (rho_e_actual / coeffs.thermodynamic.rho[i]) + 
                             2.0 * xi * F_previous[i] * F_derivatives[i];
        


/*         std::cout << rho_e_actual << std::endl;
        std::cout << std::scientific << "d dans momentum si linéarisation : " << d_term_si_linearisation << std::endl;
        std::cout << "---------------------------------------------------------------" << std::endl; */
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