#include "blast/boundary_layer/equations/energy.hpp"
#include "blast/boundary_layer/solvers/tridiagonal_solvers.hpp"
#include <algorithm>
#include <cmath>
#include <format>

namespace blast::boundary_layer::equations {

auto solve_energy(
    std::span<const double> g_previous,
    const coefficients::CoefficientInputs& inputs,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    const coefficients::XiDerivatives& xi_der,
    const io::SimulationConfig& sim_config,
    std::span<const double> F_field,
    std::span<const double> dF_deta,
    std::span<const double> V_field,
    int station,
    PhysicalQuantity auto d_eta
) -> std::expected<std::vector<double>, EquationError> {
    
    const auto n_eta = g_previous.size();
    
    // Validation
    if (n_eta != F_field.size() || n_eta != V_field.size()) {
        return std::unexpected(EquationError(
            "Energy: field size mismatch"
        ));
    }
    
    // Build coefficients
    auto energy_coeffs = detail::build_energy_coefficients(
        g_previous, inputs, coeffs, bc, xi_der, sim_config,
        F_field, dF_deta, V_field, station, d_eta
    );
    
    // Build boundary conditions  
    auto boundary_conds = detail::build_energy_boundary_conditions(
        inputs, coeffs, bc, sim_config, station, d_eta
    );
    
    // Solve tridiagonal system
    auto solution_result = solvers::solve_momentum_energy_tridiagonal(
        g_previous,
        boundary_conds.f_bc,
        boundary_conds.g_bc,
        boundary_conds.h_bc,
        energy_coeffs.a,
        energy_coeffs.b,
        energy_coeffs.c,
        energy_coeffs.d
    );
    
    if (!solution_result) {
        return std::unexpected(EquationError(
            "Energy: tridiagonal solver failed: {}", 
            std::source_location::current(), solution_result.error().message()
        ));
    }
    
    return solution_result.value();
}

namespace detail {

auto build_energy_coefficients(
    std::span<const double> g_previous,
    const coefficients::CoefficientInputs& inputs,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    const coefficients::XiDerivatives& xi_der,
    const io::SimulationConfig& sim_config,
    std::span<const double> F_field,
    std::span<const double> dF_deta,
    std::span<const double> V_field,
    int station,
    PhysicalQuantity auto d_eta
) -> EnergyCoefficients {
    
    const auto n_eta = g_previous.size();
    const auto n_species = inputs.c.rows();
    const double d_eta_sq = d_eta * d_eta;
    const double xi = inputs.xi;
    const double lambda0 = xi_der.lambda0();
    const auto g_derivatives = xi_der.g_derivative();
    
    // Compute geometry factor for diffusion fluxes
    const double J_fact = compute_energy_j_factor(station, xi, bc, sim_config);
    
    EnergyCoefficients energy_coeffs;
    energy_coeffs.a.reserve(n_eta);
    energy_coeffs.b.reserve(n_eta);
    energy_coeffs.c.reserve(n_eta);
    energy_coeffs.d.reserve(n_eta);
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        
        // a[i] = l3[i] / d_eta²
        energy_coeffs.a.push_back(coeffs.transport.l3[i] / d_eta_sq);
        
        // b[i] = (dl3_deta[i] - V[i]) / d_eta
        energy_coeffs.b.push_back((coeffs.transport.dl3_deta[i] - V_field[i]) / d_eta);
        
        // c[i] = -2*xi*F[i]*d_he_dxi/he - 2*xi*F[i]*lambda0
        const double c_term = -2.0 * xi * F_field[i] * bc.d_he_dxi() / bc.he() - 
                             2.0 * xi * F_field[i] * lambda0;
        energy_coeffs.c.push_back(c_term);
        
        // Compute species enthalpy terms
        auto [tmp1, tmp2] = compute_species_enthalpy_terms(
            inputs, coeffs, bc, J_fact, i
        );
        
        // d[i] = -ue²/he * [l0[i]*dF_deta[i]² - beta*rho_e/rho[i]*F[i]] + 
        //        2*xi*F[i]*g_der[i] + tmp1 + tmp2
        const double d_term = 
            -bc.ue() * bc.ue() / bc.he() * 
                (coeffs.transport.l0[i] * dF_deta[i] * dF_deta[i] - 
                 bc.beta * bc.rho_e() / coeffs.thermodynamic.rho[i] * F_field[i]) +
            2.0 * xi * F_field[i] * g_derivatives[i] + 
            tmp1 + tmp2;
        
        energy_coeffs.d.push_back(d_term);
    }
    
    return energy_coeffs;
}

auto build_energy_boundary_conditions(
    const coefficients::CoefficientInputs& inputs,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    const io::SimulationConfig& sim_config,
    int station,
    PhysicalQuantity auto d_eta
) -> EnergyBoundaryConditions {
    
        return EnergyBoundaryConditions{
            .f_bc = 0.0,
            .g_bc = 1.0, 
            .h_bc = coeffs.thermodynamic.h_wall / bc.he()
        };
}

auto compute_species_enthalpy_terms(
    const coefficients::CoefficientInputs& inputs,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    double J_fact,
    std::size_t eta_index
) -> std::tuple<double, double> {
    
    const auto n_species = inputs.c.rows();
    double tmp1 = 0.0;
    double tmp2 = 0.0;
    
    for (std::size_t j = 0; j < n_species; ++j) {
        
        // tmp1: concentration and enthalpy derivative terms
        const double dc_deta_j = inputs.dc_deta(j, eta_index);
        const double dc_deta2_j = inputs.dc_deta2(j, eta_index);
        const double h_sp_j = coeffs.h_species(j, eta_index);
        const double dh_sp_deta_j = coeffs.dh_species_deta(j, eta_index);
        const double dl3_deta = coeffs.transport.dl3_deta[eta_index];
        const double l3 = coeffs.transport.l3[eta_index];
        
        tmp1 += dc_deta_j * h_sp_j / bc.he() * dl3_deta + 
                l3 * h_sp_j / bc.he() * dc_deta2_j +                    // Missing term from original BLAST
                l3 * dc_deta_j * dh_sp_deta_j / bc.he();
        
        // tmp2: diffusion flux terms
        const double J_j = coeffs.diffusion.J(j, eta_index);
        const double dJ_deta_j = coeffs.diffusion.dJ_deta(j, eta_index);
        
        tmp2 += J_j * dh_sp_deta_j / bc.he() + dJ_deta_j * h_sp_j / bc.he();
    }
    
    // Apply multiplication factors
    tmp2 *= J_fact;

    return {tmp1, tmp2};
}

} // namespace detail

// Explicit instantiations for common use cases
template auto solve_energy<double>(
    std::span<const double> g_previous,
    const coefficients::CoefficientInputs& inputs,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    const coefficients::XiDerivatives& xi_der,
    const io::SimulationConfig& sim_config,
    std::span<const double> F_field,
    std::span<const double> dF_deta,
    std::span<const double> V_field,
    int station,
    double d_eta
) -> std::expected<std::vector<double>, EquationError>;

} // namespace blast::boundary_layer::equations