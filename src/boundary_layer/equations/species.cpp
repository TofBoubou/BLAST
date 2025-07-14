#include "blast/boundary_layer/equations/species.hpp"
#include "blast/boundary_layer/solvers/tridiagonal_solvers.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <format>
#include <iostream>

namespace blast::boundary_layer::equations {

auto solve_species(
    const core::Matrix<double>& c_previous,
    const coefficients::CoefficientInputs& inputs,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    const coefficients::XiDerivatives& xi_der,
    const thermophysics::MixtureInterface& mixture,
    const io::SimulationConfig& sim_config,
    std::span<const double> F_field,
    std::span<const double> V_field,
    int station,
    PhysicalQuantity auto d_eta
) -> std::expected<core::Matrix<double>, EquationError> {
    
    const auto n_eta = c_previous.cols();
    const auto n_species = c_previous.rows();
    
    std::cout << "=== SPECIES SOLVER DEBUG ===" << std::endl;
    std::cout << "Matrix size: " << n_eta << "x" << n_species << std::endl;
    std::cout << "Station: " << station << ", d_eta: " << d_eta << std::endl;
    
    // Validation
    if (n_species != mixture.n_species()) {
        return std::unexpected(EquationError(
            "Species: mixture species count mismatch"
        ));
    }
    
    if (!sim_config.chemical_non_equilibrium) {
        std::cout << "Using equilibrium composition" << std::endl;
        return compute_equilibrium_composition(inputs.T, bc.P_e(), mixture);
    }
    
    std::cout << "Using non-equilibrium transport equations" << std::endl;
    
    // Show input concentrations
    std::cout << "Input concentrations:" << std::endl;
    for (std::size_t j = 0; j < n_species; ++j) {
        std::cout << "  species[" << j << "]: ";
        for (std::size_t i = 0; i < std::min(n_eta, static_cast<size_t>(3)); ++i) {
            std::cout << c_previous(j, i) << " ";
        }
        std::cout << std::endl;
    }
    
    // Build species coefficients
    std::cout << "Building species coefficients..." << std::endl;
    auto species_coeffs_result = detail::build_species_coefficients(
        c_previous, inputs, coeffs, bc, xi_der, mixture, sim_config,
        F_field, V_field, station, d_eta
    );
    
    if (!species_coeffs_result) {
        std::cout << "ERROR: Failed to build species coefficients" << std::endl;
        return std::unexpected(species_coeffs_result.error());
    }
    auto species_coeffs = species_coeffs_result.value();
    std::cout << "Species coefficients built successfully" << std::endl;
    
    // Build boundary conditions  
    std::cout << "Building boundary conditions..." << std::endl;
    auto boundary_result = detail::build_species_boundary_conditions(
        c_previous, coeffs, bc, mixture, sim_config
    );
    
    if (!boundary_result) {
        std::cout << "ERROR: Failed to build boundary conditions" << std::endl;
        return std::unexpected(boundary_result.error());
    }
    auto boundary_conds = boundary_result.value();
    std::cout << "Boundary conditions built successfully" << std::endl;
    
    // DÉBOGAGE: Vérifier les conditions aux limites
    std::cout << "\nBoundary conditions check:" << std::endl;
    for (std::size_t j = 0; j < n_species; ++j) {
        std::cout << "  species[" << j << "]: f_bc=" << boundary_conds.f_bc[j] 
                  << " g_bc=" << boundary_conds.g_bc[j] 
                  << " h_bc=" << boundary_conds.h_bc[j] << std::endl;
    }
    
    // DÉBOGAGE: Vérifier les coefficients de la matrice
    std::cout << "\nMatrix coefficients (first few points):" << std::endl;
    for (std::size_t i = 0; i < std::min(n_eta, static_cast<size_t>(3)); ++i) {
        std::cout << "  eta=" << i << ":" << std::endl;
        for (std::size_t j = 0; j < n_species; ++j) {
            std::cout << "    species[" << j << "]: a=" << species_coeffs.a(i, j)
                      << " b=" << species_coeffs.b(i, j) 
                      << " c=" << species_coeffs.c(i, j)
                      << " d=" << species_coeffs.d(i, j) << std::endl;
            
            // Vérifier les valeurs non-finies
            if (!std::isfinite(species_coeffs.a(i, j)) || !std::isfinite(species_coeffs.b(i, j)) || 
                !std::isfinite(species_coeffs.c(i, j)) || !std::isfinite(species_coeffs.d(i, j))) {
                std::cout << "ERROR: Non-finite coefficient detected!" << std::endl;
                return std::unexpected(EquationError("Non-finite coefficients detected"));
            }
            
            // Vérifier la dominance diagonale pour les points intérieurs
            if (i > 0 && i < n_eta-1) {
                double diag = std::abs(species_coeffs.b(i, j));
                double off_diag = std::abs(species_coeffs.a(i, j)) + std::abs(species_coeffs.c(i, j));
                double ratio = diag / (off_diag + 1e-15);  // Éviter division par zéro
                
                if (ratio < 1.0) {
                    std::cout << "      WARNING: Diagonal dominance ratio = " << ratio << " < 1.0" << std::endl;
                }
                
                // Vérifier les coefficients très grands
                if (diag > 1e12 || off_diag > 1e12) {
                    std::cout << "      WARNING: Very large coefficients detected!" << std::endl;
                    std::cout << "        |b| = " << diag << ", |a|+|c| = " << off_diag << std::endl;
                }
            }
        }
    }
    
    // DÉBOGAGE: Vérifier les champs d'entrée
    std::cout << "\nInput fields check:" << std::endl;
    std::cout << "  F_field: ";
    for (std::size_t i = 0; i < std::min(static_cast<size_t>(3), F_field.size()); ++i) {
        std::cout << F_field[i] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "  V_field: ";
    for (std::size_t i = 0; i < std::min(static_cast<size_t>(3), V_field.size()); ++i) {
        std::cout << V_field[i] << " ";
    }
    std::cout << std::endl;
    
    // Solve block tridiagonal system
    std::cout << "\nCalling block tridiagonal solver..." << std::endl;
    auto solution_result = solvers::solve_species_block_tridiagonal(
        c_previous,
        boundary_conds.f_bc,
        boundary_conds.g_bc,
        boundary_conds.h_bc,
        species_coeffs.a,
        species_coeffs.b,
        species_coeffs.c,
        species_coeffs.d,
        mixture.has_electrons()
    );
    
    if (!solution_result) {
        std::cout << "ERROR: Block tridiagonal solver failed!" << std::endl;
        std::cout << "Error: " << solution_result.error().message() << std::endl;
        return std::unexpected(EquationError(
            "Species: block tridiagonal solver failed: {}", 
            std::source_location::current(), solution_result.error().message()
        ));
    }
    
    auto result = solution_result.value();
    std::cout << "Block tridiagonal solver completed" << std::endl;
    
    // DÉBOGAGE: Vérifier le résultat détaillé
    std::cout << "\nSolver output check:" << std::endl;
    for (std::size_t i = 0; i < std::min(n_eta, static_cast<size_t>(5)); ++i) {
        double sum = 0.0;
        std::cout << "  eta=" << i << ": ";
        for (std::size_t j = 0; j < n_species; ++j) {
            std::cout << "species[" << j << "]=" << result(j, i) << " ";
            sum += result(j, i);
        }
        std::cout << " (sum=" << sum << ")" << std::endl;
        
        // Vérifier la divergence
        if (sum > 1000.0 || sum < -1000.0 || !std::isfinite(sum)) {
            std::cout << "ERROR: Divergence detected at eta=" << i << std::endl;
            std::cout << "Sum = " << sum << std::endl;
            
            // Montrer toutes les espèces pour ce point
            for (std::size_t j = 0; j < n_species; ++j) {
                std::cout << "  species[" << j << "] = " << result(j, i) << std::endl;
            }
            
            return std::unexpected(EquationError("Species solver diverged"));
        }
        
        // Vérifier les valeurs négatives importantes
        for (std::size_t j = 0; j < n_species; ++j) {
            if (result(j, i) < -0.1) {
                std::cout << "WARNING: Large negative concentration at eta=" << i 
                          << " species=" << j << " value=" << result(j, i) << std::endl;
            }
        }
    }
    
    // Apply charge neutrality for ionized mixtures
    if (mixture.has_electrons()) {
        std::cout << "Applying charge neutrality..." << std::endl;
        apply_charge_neutrality(result, mixture);
    }
    
    std::cout << "Species solver completed successfully" << std::endl;
    std::cout << "=== END SPECIES SOLVER DEBUG ===" << std::endl;
    
    return result;
}

auto compute_equilibrium_composition(
    std::span<const double> temperature_field,
    double pressure,
    const thermophysics::MixtureInterface& mixture
) -> std::expected<core::Matrix<double>, EquationError> {
    
    const auto n_eta = temperature_field.size();
    const auto n_species = mixture.n_species();
    
    core::Matrix<double> c_equilibrium(n_species, n_eta);
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        auto eq_result = mixture.equilibrium_composition(temperature_field[i], pressure);
        if (!eq_result) {
            return std::unexpected(EquationError(
                "Failed to compute equilibrium at eta={}: {}", 
                std::source_location::current(), i, eq_result.error().message()
            ));
        }
        
        const auto& eq_composition = eq_result.value();
        for (std::size_t j = 0; j < n_species; ++j) {
            c_equilibrium(j, i) = eq_composition[j];
        }
    }
    
    return c_equilibrium;
}

auto apply_charge_neutrality(
    core::Matrix<double>& species_matrix,
    const thermophysics::MixtureInterface& mixture
) -> void {
    
    if (!mixture.has_electrons()) return;
    
    const auto n_eta = species_matrix.cols();
    const auto n_species = species_matrix.rows();
    auto charges = mixture.species_charges();
    
    // Electrons are typically the first species (index 0)
    for (std::size_t i = 0; i < n_eta; ++i) {
        double charge_sum = 0.0;
        
        // Sum charges from all non-electron species
        for (std::size_t j = 1; j < n_species; ++j) {
            charge_sum += species_matrix(j, i) * charges[j];
        }
        
        // Set electron concentration to ensure neutrality
        if (std::abs(charges[0]) > 1e-15) {
            species_matrix(0, i) = -charge_sum / charges[0];
        }
    }
}

namespace detail {

auto build_species_coefficients(
    const core::Matrix<double>& c_previous,
    const coefficients::CoefficientInputs& inputs,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    const coefficients::XiDerivatives& xi_der,
    const thermophysics::MixtureInterface& mixture,
    const io::SimulationConfig& sim_config,
    std::span<const double> F_field,
    std::span<const double> V_field,
    int station,
    PhysicalQuantity auto d_eta
) -> std::expected<SpeciesCoefficients, EquationError> {
    
    const auto n_eta = c_previous.cols();
    const auto n_species = c_previous.rows();
    const double d_eta_sq = d_eta * d_eta;
    const double xi = inputs.xi;
    const double lambda0 = xi_der.lambda0();
    const auto c_derivatives = xi_der.c_derivative();
    
    // Compute geometry factors
    const auto factors_result = compute_geometry_factors(station, xi, bc, sim_config);
    if (!factors_result) {
        return std::unexpected(factors_result.error());
    }
    const auto factors = factors_result.value();
    
    // Parameters for fake fluxes
    constexpr double Le = 1.2;  // Lewis number
    constexpr double Pr = 0.72; // Prandtl number
    
    // Fix concentration derivatives for mass conservation
    auto dc_deta_fixed = fix_concentration_derivatives(c_previous, inputs.dc_deta);
    
    // Compute fake fluxes and their derivatives
    core::Matrix<double> J_fake, dJ_fake_deta;
    try {
        auto result = compute_fake_fluxes(dc_deta_fixed, coeffs, d_eta, Le, Pr);
        J_fake = std::move(result.first);
        dJ_fake_deta = std::move(result.second);
    } catch (const std::exception& e) {
        return std::unexpected(EquationError("Failed to compute fake fluxes"));
    }
    
    SpeciesCoefficients species_coeffs;
    species_coeffs.a = core::Matrix<double>(n_eta, n_species);
    species_coeffs.b = core::Matrix<double>(n_eta, n_species);
    species_coeffs.c = core::Matrix<double>(n_eta, n_species);
    species_coeffs.d = core::Matrix<double>(n_eta, n_species);
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        // std::cout << "l0 " << coeffs.transport.l0[i] << std::endl;
        for (std::size_t j = 0; j < n_species; ++j) {
            
            // a[i][j] = -l0[i]*Le/Pr / d_eta²
            species_coeffs.a(i, j) = -coeffs.transport.l0[i] * Le / Pr / d_eta_sq;
            
            // b[i][j] = (V[i] - Le/Pr*dl0_deta[i]) / d_eta
            species_coeffs.b(i, j) = (V_field[i] - Le / Pr * coeffs.transport.dl0_deta[i]) / d_eta;
            // std::cout << "V " << V_field[i] << std::endl;
            // std::cout << "dl0 " << coeffs.transport.dl0_deta[i] << std::endl;
            
            // c[i][j] = 2*xi*F[i]*lambda0
            species_coeffs.c(i, j) = 2.0 * xi * F_field[i] * lambda0;
            
            // d[i][j] = -dJ_fake_deta[j][i] - dJ_deta[j][i]*J_fact - 2*xi*c_der[j][i]*F[i]
            const double d_term = 
                -dJ_fake_deta(j, i) - 
                coeffs.diffusion.dJ_deta(j, i) * factors.J_fact - 
                2.0 * xi * c_derivatives(j, i) * F_field[i];

            const double wi_term = coeffs.chemical.wi(i, j) * factors.W_fact / coeffs.thermodynamic.rho[i];

/*             // DÉBOGAGE DU TERME CHIMIQUE
            if (i < 3 && j >= 3) {  // Vérifier les espèces qui explosent (N2=3, O2=4)
                std::cout << "DEBUG wi_term: eta=" << i << " species=" << j << std::endl;
                std::cout << "  wi=" << coeffs.chemical.wi(i, j) << std::endl;
                std::cout << "  W_fact=" << factors.W_fact << std::endl;
                std::cout << "  rho=" << coeffs.thermodynamic.rho[i] << std::endl;
                std::cout << "  wi_term=" << wi_term << std::endl;
                std::cout << "  d_term=" << d_term << std::endl;
                std::cout << "  total_d=" << d_term + wi_term << std::endl;
            } */

            species_coeffs.d(i, j) = d_term + wi_term;
            
            // Apply molecular weight normalization
/*             const double Mw_j = mixture.species_molecular_weight(j);
            species_coeffs.a(i, j) /= Mw_j;
            species_coeffs.b(i, j) /= Mw_j;
            species_coeffs.c(i, j) /= Mw_j;
            species_coeffs.d(i, j) /= Mw_j; */
        }
    }
    
    return species_coeffs;
}

auto build_species_boundary_conditions(
    const core::Matrix<double>& c_wall,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    const thermophysics::MixtureInterface& mixture,
    const io::SimulationConfig& sim_config
) -> std::expected<SpeciesBoundaryConditions, EquationError> {
    
    const auto n_species = mixture.n_species();
    
    SpeciesBoundaryConditions boundary_conds;
    boundary_conds.f_bc.resize(n_species);
    boundary_conds.g_bc.resize(n_species);
    boundary_conds.h_bc.resize(n_species);
    
    for (std::size_t i = 0; i < n_species; ++i) {
        boundary_conds.f_bc[i] = 0.0;
        boundary_conds.g_bc[i] = 1.0;
    }
    
    auto eq_wall_result = compute_equilibrium_wall(bc, mixture);
    if (!eq_wall_result) {
        return std::unexpected(eq_wall_result.error());
    }
    
    boundary_conds.h_bc = eq_wall_result.value();
        
    // Apply molecular weight normalization
/*     for (std::size_t i = 0; i < n_species; ++i) {
        const double Mw_i = mixture.species_molecular_weight(i);
        std::cout << Mw_i << std::endl;
        boundary_conds.f_bc[i] /= Mw_i;
        boundary_conds.g_bc[i] /= Mw_i;
        boundary_conds.h_bc[i] /= Mw_i;
    } */
    
    return boundary_conds;
}

auto compute_fake_fluxes(
    const core::Matrix<double>& dc_deta_fixed,
    const coefficients::CoefficientSet& coeffs,
    PhysicalQuantity auto d_eta,
    double Le,
    double Pr
) -> std::pair<core::Matrix<double>, core::Matrix<double>> {
    
    const auto n_species = dc_deta_fixed.rows();
    const auto n_eta = dc_deta_fixed.cols();
    
    core::Matrix<double> J_fake(n_species, n_eta);
    core::Matrix<double> dJ_fake_deta(n_species, n_eta);
    
    // Compute fake fluxes: J_fake[j][i] = Le/Pr * l0[i] * dc_deta_fix[j][i]
    for (std::size_t i = 0; i < n_eta; ++i) {
        for (std::size_t j = 0; j < n_species; ++j) {
            J_fake(j, i) = Le / Pr * coeffs.transport.l0[i] * dc_deta_fixed(j, i);
        }
    }
    
    // Compute derivatives of fake fluxes
    for (std::size_t j = 0; j < n_species; ++j) {
        auto J_row = J_fake.eigen().row(j);
        auto dJ_result = coefficients::derivatives::compute_eta_derivative(
            std::span(J_row.data(), n_eta), d_eta
        );
        if (!dJ_result) {
            throw std::runtime_error("Failed to compute diffusion flux derivative");
        }
        auto dJ = dJ_result.value();
        
        for (std::size_t i = 0; i < n_eta; ++i) {
            dJ_fake_deta(j, i) = dJ[i];
        }
    }
    
    return std::make_pair(std::move(J_fake), std::move(dJ_fake_deta));
}

auto fix_concentration_derivatives(
    const core::Matrix<double>& c_matrix,
    const core::Matrix<double>& dc_deta
) -> core::Matrix<double> {
    
    const auto n_species = c_matrix.rows();
    const auto n_eta = c_matrix.cols();
    
    core::Matrix<double> dc_deta_fixed(n_species, n_eta);
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        // Compute sums
        double sum_c = 0.0;
        double sum_dc = 0.0;
        
        for (std::size_t j = 0; j < n_species; ++j) {
            sum_c += c_matrix(j, i);
            sum_dc += dc_deta(j, i);
        }
        
        // Apply mass conservation correction
        for (std::size_t j = 0; j < n_species; ++j) {
            dc_deta_fixed(j, i) = dc_deta(j, i) - c_matrix(j, i) / sum_c * sum_dc;
        }
    }
    
    return dc_deta_fixed;
}

auto compute_equilibrium_wall(
    const conditions::BoundaryConditions& bc,
    const thermophysics::MixtureInterface& mixture
) -> std::expected<std::vector<double>, EquationError> {
    
    auto eq_result = mixture.equilibrium_composition(bc.Tw(), bc.P_e());
    if (!eq_result) {
        return std::unexpected(EquationError(
            "Failed to compute equilibrium wall composition: {}", 
            std::source_location::current(), eq_result.error().message()
        ));
    }
    
    return eq_result.value();
}

} // namespace detail

// Explicit instantiations for common use cases
template auto solve_species<double>(
    const core::Matrix<double>& c_previous,
    const coefficients::CoefficientInputs& inputs,
    const coefficients::CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    const coefficients::XiDerivatives& xi_der,
    const thermophysics::MixtureInterface& mixture,
    const io::SimulationConfig& sim_config,
    std::span<const double> F_field,
    std::span<const double> V_field,
    int station,
    double d_eta
) -> std::expected<core::Matrix<double>, EquationError>;

} // namespace blast::boundary_layer::equations