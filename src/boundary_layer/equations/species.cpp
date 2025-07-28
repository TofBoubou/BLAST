#include "blast/boundary_layer/equations/species.hpp"
#include "blast/boundary_layer/solvers/tridiagonal_solvers.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <format>
#include <iostream>
#include <iomanip>

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
    
    if (n_species != mixture.n_species()) {
        return std::unexpected(EquationError(
            "Species: mixture species count mismatch"
        ));
    }
    
    if (!sim_config.chemical_non_equilibrium) {
        return compute_equilibrium_composition(inputs.T, bc.P_e(), mixture);
    }
    
    auto species_coeffs_result = detail::build_species_coefficients(
        c_previous, inputs, coeffs, bc, xi_der, mixture, sim_config,
        F_field, V_field, station, d_eta
    );
    
    if (!species_coeffs_result) {
        return std::unexpected(species_coeffs_result.error());
    }
    auto species_coeffs = species_coeffs_result.value();
    
    auto boundary_result = detail::build_species_boundary_conditions(
        c_previous, coeffs, bc, mixture, sim_config
    );
    
    if (!boundary_result) {
        return std::unexpected(boundary_result.error());
    }
    auto boundary_conds = boundary_result.value();
    
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
        return std::unexpected(EquationError(
            "Species: block tridiagonal solver failed: {}", 
            std::source_location::current(), solution_result.error().message()
        ));
    }
    
    auto result = solution_result.value();
    
    if (mixture.has_electrons()) {
        apply_charge_neutrality(result, mixture);
    }
    
    // Debug: afficher l'état des fractions massiques à chaque eta
/*     std::cout << "\n=== DEBUG SOLVE_SPECIES - Station " << station << " ===\n";
    for (std::size_t i = 0; i < n_eta; ++i) {
        std::cout << "eta[" << i << "]: ";
        double sum = 0.0;
        for (std::size_t j = 0; j < n_species; ++j) {
            std::cout << std::format("c[{}]={:.6e} ", j, result(j, i));
            sum += result(j, i);
        }
        std::cout << std::format("(sum={:.6e})\n", sum);
    }
    std::cout << "=========================================\n\n"; */
    
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
    std::cout << "Pressure dans l'équilibre chimique : " << pressure << std::endl;
    std::cout << "Temperature dans l'équilibre chimique : " << temperature_field[19] << std::endl;
    
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

/*     std::cout << "Fractions dans equilibre chimique c_19_0: " << c_equilibrium(0,19) << std::endl;
    std::cout << "Fractions dans equilibre chimique c_19_1: " << c_equilibrium(1,19) << std::endl;
    std::cout << "Fractions dans equilibre chimique c_19_2: " << c_equilibrium(2,19) << std::endl;
    std::cout << "Fractions dans equilibre chimique c_19_3: " << c_equilibrium(3,19) << std::endl;
    std::cout << "Fractions dans equilibre chimique c_19_4: " << c_equilibrium(4,19) << std::endl; */
    
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
        
        // Get mass fractions for this eta point
        std::vector<double> mass_fractions(n_species);
        for (std::size_t j = 0; j < n_species; ++j) {
            mass_fractions[j] = species_matrix(j, i);
        }
        
        // Convert mass fractions to mole fractions
        auto mole_fractions_result = mixture.mass_fractions_to_mole_fractions(mass_fractions);
        if (!mole_fractions_result) {
            // If conversion fails, fall back to original (incorrect) method
            double charge_sum = 0.0;
            for (std::size_t j = 1; j < n_species; ++j) {
                charge_sum += species_matrix(j, i) * charges[j];
            }
            if (std::abs(charges[0]) > 1e-15) {
                species_matrix(0, i) = -charge_sum / charges[0];
            }
            continue;
        }
        auto mole_fractions = mole_fractions_result.value();
        
        // Calculate mixture molecular weight for molar concentration conversion
        auto mixture_mw_result = mixture.mixture_molecular_weight(mass_fractions);
        if (!mixture_mw_result) {
            // Fall back to original method if MW calculation fails
            double charge_sum = 0.0;
            for (std::size_t j = 1; j < n_species; ++j) {
                charge_sum += species_matrix(j, i) * charges[j];
            }
            if (std::abs(charges[0]) > 1e-15) {
                species_matrix(0, i) = -charge_sum / charges[0];
            }
            continue;
        }
        double mixture_mw = mixture_mw_result.value();
        
        // Convert mole fractions to molar concentrations
        // Assume total molar density = 1.0 mol/m³ (relative calculation)
        const double total_molar_density = 1.0;
        std::vector<double> molar_concentrations(n_species);
        for (std::size_t j = 0; j < n_species; ++j) {
            molar_concentrations[j] = mole_fractions[j] * total_molar_density;
        }
        
        // Apply charge neutrality: ∑(c_i × z_i) = 0
        // Calculate charge sum from all non-electron species
        double charge_sum = 0.0;
        for (std::size_t j = 1; j < n_species; ++j) {
            // Convert charge from C/kg to elementary charges
            const double species_mw = mixture.species_molecular_weight(j);
            const double elementary_charge_per_mol = charges[j] * species_mw / (1.602176634e-19 * 6.02214076e23);
            charge_sum += molar_concentrations[j] * elementary_charge_per_mol;
        }
        
        // Set electron molar concentration to ensure neutrality
        const double electron_mw = mixture.species_molecular_weight(0);
        const double electron_elementary_charge_per_mol = charges[0] * electron_mw / (1.602176634e-19 * 6.02214076e23);
        if (std::abs(electron_elementary_charge_per_mol) > 1e-15) {
            molar_concentrations[0] = -charge_sum / electron_elementary_charge_per_mol;
        }
        
        // Convert back to mole fractions
        double total_moles = 0.0;
        for (std::size_t j = 0; j < n_species; ++j) {
            total_moles += molar_concentrations[j];
        }
        if (total_moles > 1e-15) {
            for (std::size_t j = 0; j < n_species; ++j) {
                mole_fractions[j] = molar_concentrations[j] / total_moles;
            }
        }
        
        // Convert mole fractions back to mass fractions
        // Y_i = X_i * M_i / MW_mix
        double new_mixture_mw = 0.0;
        for (std::size_t j = 0; j < n_species; ++j) {
            new_mixture_mw += mole_fractions[j] * mixture.species_molecular_weight(j);
        }
        
        for (std::size_t j = 0; j < n_species; ++j) {
            mass_fractions[j] = mole_fractions[j] * mixture.species_molecular_weight(j) / new_mixture_mw;
        }
        
        // Renormalize to ensure mass conservation
        double sum = 0.0;
        for (std::size_t j = 0; j < n_species; ++j) {
            sum += mass_fractions[j];
        }
        if (sum > 1e-15) {
            for (std::size_t j = 0; j < n_species; ++j) {
                mass_fractions[j] /= sum;
            }
        }
        
        // Update species matrix
        for (std::size_t j = 0; j < n_species; ++j) {
            species_matrix(j, i) = mass_fractions[j];
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
        for (std::size_t j = 0; j < n_species; ++j) {
            
            // a[i][j] = -l0[i]*Le/Pr / d_eta²
/*             species_coeffs.a(i, j) = -coeffs.transport.l0[i] * Le / Pr / d_eta_sq;
            
            // b[i][j] = (V[i] - Le/Pr*dl0_deta[i]) / d_eta
            species_coeffs.b(i, j) = (V_field[i] - Le / Pr * coeffs.transport.dl0_deta[i]) / d_eta;
            
            // c[i][j] = 2*xi*F[i]*lambda0
            species_coeffs.c(i, j) = 2.0 * xi * F_field[i] * lambda0; */

/*             std::cout << "-------------------------------------------------------------------" << std::endl;
            std::cout << "COEF DANS SPECIES : " << std::endl; */

            species_coeffs.a(i, j) = -coeffs.transport.l0[i] * Le / Pr / d_eta_sq;
            // std::cout << "[i=" << i << "][j=" << j << "] a = " << std::scientific << species_coeffs.a(i, j) << std::endl;

            species_coeffs.b(i, j) = (V_field[i] - Le / Pr * coeffs.transport.dl0_deta[i]) / d_eta;
            // std::cout << "[i=" << i << "][j=" << j << "] b = " << std::scientific << species_coeffs.b(i, j) << std::endl;

            species_coeffs.c(i, j) = 2.0 * xi * F_field[i] * lambda0;
            // std::cout << "[i=" << i << "][j=" << j << "] c = " << std::scientific << species_coeffs.c(i, j) << std::endl;
            // std::cout << "-------------------------------------------------------------------" << std::endl;
            
            // d[i][j] = -dJ_fake_deta[j][i] - dJ_deta[j][i]*J_fact - 2*xi*c_der[j][i]*F[i]
            const double d_term = 
                -dJ_fake_deta(j, i) - 
                coeffs.diffusion.dJ_deta(j, i) * factors.J_fact - 
                2.0 * xi * c_derivatives(j, i) * F_field[i];

            const double wi_term = coeffs.chemical.wi(i, j) * factors.W_fact / coeffs.thermodynamic.rho[i];

            species_coeffs.d(i, j) = d_term + wi_term;

/*             std::cout << std::scientific << std::setprecision(8)
                    << "[d_coef] i=" << i << ", j=" << j << "\n"
                    << "  - dJ_fake_deta      = " << -dJ_fake_deta(j, i) << "\n"
                    << "  - J_term            = " << -coeffs.diffusion.dJ_deta(j, i) * factors.J_fact << "\n"
                    << "  - c_derivative_term = " << -2.0 * xi * c_derivatives(j, i) * F_field[i] << "\n"
                    << "  = d_term            = " << d_term << "\n"
                    << "  + wi_term           = " << wi_term << "\n"
                    << "  => species_coeffs.d = " << species_coeffs.d(i, j) << "\n"; */

            // std::cout << "[i=" << i << "][j=" << j << "] d = " << std::scientific << species_coeffs.d(i, j) << std::endl;
            // std::cout << "-------------------------------------------------------------------" << std::endl;
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
    
    return boundary_conds;
}

/* auto compute_fake_fluxes(
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
            // std::cout << std::scientific << std::setprecision(8)
                // << "[J_fake] j=" << j << ", i=" << i << "\n"
                // << "  Le         = " << Le << "\n"
                // << "  Pr         = " << Pr << "\n"
                // << "  l0[i]      = " << coeffs.transport.l0[i] << "\n"
                // << "  dc_deta    = " << dc_deta_fixed(j, i) << "\n"
                // << "  -> J_fake  = " << Le / Pr * coeffs.transport.l0[i] * dc_deta_fixed(j, i) << "\n";
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

    std::cout << std::scientific << std::setprecision(8);
    std::cout << "[dJ_fake_deta] tableau complet (n_species=" << n_species
            << ", n_eta=" << n_eta << ")\n";
    for (std::size_t j = 0; j < n_species; ++j) {
        std::cout << "species " << j << " : ";
        for (std::size_t i = 0; i < n_eta; ++i) {
            std::cout << dJ_fake_deta(j, i) << " ";
        }
        std::cout << "\n";
    }
    
    return std::make_pair(std::move(J_fake), std::move(dJ_fake_deta));
} */

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
        // Copie élément par élément pour garantir la contiguïté
        std::vector<double> J_row_data(n_eta);
        for (std::size_t i = 0; i < n_eta; ++i) {
            J_row_data[i] = J_fake(j, i);
        }
        
        auto dJ_result = coefficients::derivatives::compute_eta_derivative(
            std::span(J_row_data.data(), n_eta), d_eta
        );

        if (!dJ_result) {
            throw std::runtime_error("Failed to compute diffusion flux derivative");
        }
        auto dJ = dJ_result.value();
        
        for (std::size_t i = 0; i < n_eta; ++i) {
            dJ_fake_deta(j, i) = dJ[i];
        }
    }

    std::cout << std::scientific << std::setprecision(8);
    std::cout << "[dJ_fake_deta] tableau complet (n_species=" << n_species
            << ", n_eta=" << n_eta << ")\n";
    for (std::size_t j = 0; j < n_species; ++j) {
        std::cout << "species " << j << " : ";
        for (std::size_t i = 0; i < n_eta; ++i) {
            std::cout << dJ_fake_deta(j, i) << " ";
        }
        std::cout << "\n";
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