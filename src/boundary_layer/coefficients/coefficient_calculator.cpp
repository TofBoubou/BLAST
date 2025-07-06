#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include "blast/boundary_layer/coefficients/diffusion.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <ranges>
#include <iostream>

namespace blast::boundary_layer::coefficients {

using namespace blast::boundary_layer;

auto CoefficientCalculator::calculate(
    const CoefficientInputs& inputs,
    const conditions::BoundaryConditions& bc,
    const XiDerivatives& xi_der
) const -> std::expected<CoefficientSet, CoefficientError> {
    
    CoefficientSet coeffs;
    
    // 1. Calculate thermodynamic coefficients first (needed by others)
    auto thermo_result = calculate_thermodynamic_coefficients(inputs, bc);
    if (!thermo_result) {
        return std::unexpected(thermo_result.error());
    }
    coeffs.thermodynamic = std::move(thermo_result.value());
    
    // 3. Transport coefficients
    auto transport_result = calculate_transport_coefficients(inputs, coeffs.thermodynamic, bc);
    if (!transport_result) {
        return std::unexpected(transport_result.error());
    }
    coeffs.transport = std::move(transport_result.value());
    
    // 4. Diffusion coefficients
    auto diffusion_result = calculate_diffusion_coefficients(inputs, bc, xi_der);
    if (!diffusion_result) {
        return std::unexpected(diffusion_result.error());
    }
    coeffs.diffusion = std::move(diffusion_result.value());
    
    // 5. Chemical coefficients (if non-equilibrium)
    if (sim_config_.chemical_non_equilibrium) {
        auto chemical_result = calculate_chemical_coefficients(inputs, coeffs.thermodynamic, bc);
        if (!chemical_result) {
            return std::unexpected(chemical_result.error());
        }
        coeffs.chemical = std::move(chemical_result.value());
    } else {
        // Initialize with zeros
        const auto n_eta = inputs.T.size();
        const auto n_species = inputs.c.rows();
        coeffs.chemical.wi = core::Matrix<double>(n_eta, n_species);
        coeffs.chemical.wi.setZero();
        coeffs.chemical.d_wi_dc.resize(n_eta);
        for (auto& mat : coeffs.chemical.d_wi_dc) {
            mat = core::Matrix<double>(n_species, n_species);
            mat.setZero();
        }
    }
    
    // 6. Thermal diffusion (if enabled)
    if (sim_config_.consider_thermal_diffusion) {
        auto tdr_result = calculate_thermal_diffusion(inputs, coeffs.thermodynamic, bc);
        if (!tdr_result) {
            return std::unexpected(tdr_result.error());
        }
        coeffs.thermal_diffusion = std::move(tdr_result.value());
    }
    
    // 7. Wall properties
    auto wall_result = calculate_wall_properties(inputs, bc, coeffs.transport, coeffs.thermodynamic);
    if (!wall_result) {
        return std::unexpected(wall_result.error());
    }
    coeffs.wall = std::move(wall_result.value());
    
    // 8. Species enthalpies
    auto enthalpy_result = calculate_species_enthalpies(inputs);
    if (!enthalpy_result) {
        return std::unexpected(enthalpy_result.error());
    }
    auto [h_sp, dh_sp_deta] = std::move(enthalpy_result.value());
    coeffs.h_species = std::move(h_sp);
    coeffs.dh_species_deta = std::move(dh_sp_deta);

    // 9. Compute diffusion fluxes
    auto flux_result = diffusion::compute_stefan_maxwell_fluxes(
        inputs, coeffs, bc, xi_der, sim_config_, mixture_, d_eta_
    );
    if (!flux_result) {
        return std::unexpected(flux_result.error());
    }
    
    return coeffs;
}

auto CoefficientCalculator::calculate_thermodynamic_coefficients(
    const CoefficientInputs& inputs,
    const conditions::BoundaryConditions& bc
) const -> std::expected<ThermodynamicCoefficients, CoefficientError> {
    
    const auto n_eta = inputs.T.size();
    const auto n_species = inputs.c.rows();
    
    ThermodynamicCoefficients thermo;
    thermo.rho.reserve(n_eta);
    thermo.MW.reserve(n_eta);
    
    // Calculate molecular weights and densities
    for (std::size_t i = 0; i < n_eta; ++i) {
        // Extract species concentrations at this eta point
        std::vector<double> c_local(n_species);
        for (std::size_t j = 0; j < n_species; ++j) {
            c_local[j] = inputs.c(j, i);
        }
        
        // Get molecular weight
        auto mw_result = mixture_.mixture_molecular_weight(c_local);
        if (!mw_result) {
            return std::unexpected(CoefficientError(
                std::format("Failed to compute MW at eta={}: {}", i, mw_result.error().message())
            ));
        }
        thermo.MW.push_back(mw_result.value());
        
        // Compute density using ideal gas law: ρ = P·MW/(R·T)
        // Always use edge pressure to avoid circular dependency
        thermo.rho.push_back(bc.P_e() * thermo.MW[i] / 
                            (inputs.T[i] * thermophysics::constants::R_universal));
    }
    
    // Compute molecular weight derivative
    thermo.d_MW_deta.resize(n_eta);
    for (std::size_t i = 0; i < n_eta; ++i) {
        std::vector<double> dc_deta_local(n_species);
        for (std::size_t j = 0; j < n_species; ++j) {
            dc_deta_local[j] = inputs.dc_deta(j, i);
        }
        
        // d(MW)/dη = -MW² · Σ(dc_j/dη / Mw_j)
        double sum = 0.0;
        for (std::size_t j = 0; j < n_species; ++j) {
            sum += dc_deta_local[j] / mixture_.species_molecular_weight(j);
        }
        thermo.d_MW_deta[i] = -thermo.MW[i] * thermo.MW[i] * sum;
    }
    
    // Compute density derivative
    thermo.d_rho_deta = derivatives::compute_eta_derivative(thermo.rho, d_eta_);
    
    // Compute wall enthalpy
    std::vector<double> c_wall(n_species);
    for (std::size_t j = 0; j < n_species; ++j) {
        c_wall[j] = inputs.c(j, 0);
    }
    
    auto h_wall_result = mixture_.mixture_enthalpy(c_wall, bc.Tw(), bc.P_e());
    if (!h_wall_result) {
        return std::unexpected(CoefficientError("Failed to compute wall enthalpy"));
    }
    thermo.h_wall = h_wall_result.value();
    
    return thermo;
}

auto CoefficientCalculator::calculate_transport_coefficients(
    const CoefficientInputs& inputs,
    const ThermodynamicCoefficients& thermo,
    const conditions::BoundaryConditions& bc
) const -> std::expected<TransportCoefficients, CoefficientError> {
    
    const auto n_eta = inputs.T.size();
    const auto n_species = inputs.c.rows();
    
    TransportCoefficients transport;
    transport.l0.reserve(n_eta);
    transport.l3.reserve(n_eta);
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        // Extract local composition
        std::vector<double> c_local(n_species);
        for (std::size_t j = 0; j < n_species; ++j) {
            c_local[j] = inputs.c(j, i);
        }
        
        // Use edge pressure consistently to avoid circular dependency
        const double P_edge = bc.P_e();
        
        // Get transport properties
        auto mu_result = mixture_.viscosity(c_local, inputs.T[i], P_edge);
        auto cp_result = mixture_.frozen_cp(c_local, inputs.T[i], P_edge);
        auto k_result = mixture_.frozen_thermal_conductivity(c_local, inputs.T[i], P_edge);
        
        if (!mu_result || !cp_result || !k_result) {
            return std::unexpected(CoefficientError(
                std::format("Failed to compute transport properties at eta={}", i)
            ));
        }
        
        const double mu = mu_result.value();
        const double Cp = cp_result.value();
        const double k_fr = k_result.value();
        const double Pr = mu * Cp / k_fr;
        
        // Compute coefficients
        transport.l0.push_back(thermo.rho[i] * mu / (bc.rho_e() * bc.mu_e()));
        transport.l3.push_back(transport.l0[i] / Pr);
    }
    
    // Compute derivatives
    transport.dl0_deta = derivatives::compute_eta_derivative(transport.l0, d_eta_);
    transport.dl3_deta = derivatives::compute_eta_derivative(transport.l3, d_eta_);
    
    return transport;
}

auto CoefficientCalculator::calculate_diffusion_coefficients(
    const CoefficientInputs& inputs,
    const conditions::BoundaryConditions& bc,
    const XiDerivatives& xi_der
) const -> std::expected<DiffusionCoefficients, CoefficientError> {
    
    const auto n_eta = inputs.T.size();
    const auto n_species = inputs.c.rows();
    
    DiffusionCoefficients diff;
    diff.Dij_bin = core::Matrix<double>(n_eta * n_species, n_species);
    diff.y.reserve(n_eta);
    
    // Use edge pressure consistently
    const double P_edge = bc.P_e();
    
    // Calculate binary diffusion coefficients
    for (std::size_t i = 0; i < n_eta; ++i) {
        // Get binary diffusion coefficients using edge pressure
        auto dij_result = mixture_.binary_diffusion_coefficients(inputs.T[i], P_edge);
        if (!dij_result) {
            return std::unexpected(CoefficientError(
                std::format("Failed to compute Dij at eta={}", i)
            ));
        }
        
        // Store in matrix (row i*n_species to (i+1)*n_species)
        const auto& dij_local = dij_result.value();

        for (std::size_t j = 0; j < n_species; ++j) {
            for (std::size_t k = 0; k < n_species; ++k) {
                diff.Dij_bin(i * n_species + j, k) = dij_local(j, k);
            }
        }
        
        // Compute y coordinate
        diff.y.push_back(-(2.0 * inputs.xi * xi_der.lambda0() + 1.0) * inputs.F[i] - 
                         2.0 * inputs.xi * xi_der.F_derivative()[i]);
    }
    
    return diff;
}

auto CoefficientCalculator::calculate_chemical_coefficients(
    const CoefficientInputs& inputs,
    const ThermodynamicCoefficients& thermo,
    const conditions::BoundaryConditions& bc
) const -> std::expected<ChemicalCoefficients, CoefficientError> {
    
    const auto n_eta = inputs.T.size();
    const auto n_species = inputs.c.rows();
    
    ChemicalCoefficients chem;
    chem.wi = core::Matrix<double>(n_eta, n_species);
    chem.d_wi_dc.reserve(n_eta);
    
    // Calculate production rates for each eta station
    for (std::size_t i = 0; i < n_eta; ++i) {
        auto result = calculate_station_coefficients(i, inputs, thermo, bc);
        if (!result) {
            return std::unexpected(result.error());
        }
        
        // Store results
        const auto& [wi_local, dwi_dc] = result.value();
        for (std::size_t j = 0; j < n_species; ++j) {
            chem.wi(i, j) = wi_local[j];
        }
        chem.d_wi_dc.push_back(std::move(dwi_dc));
    }
    
    return chem;
}

auto CoefficientCalculator::calculate_station_coefficients(
    std::size_t station_index,
    const CoefficientInputs& inputs,
    const ThermodynamicCoefficients& thermo,
    const conditions::BoundaryConditions& bc
) const -> std::expected<std::pair<std::vector<double>, core::Matrix<double>>, CoefficientError> {
    
    const auto n_species = inputs.c.rows();
    
    // Extract and normalize composition
    auto [c_local, rho_species] = extract_local_composition(
        station_index, inputs.c, thermo.rho[station_index]
    );
    
    // Get production rates
    auto wi_result = mixture_.production_rates(rho_species, inputs.T[station_index]);
    if (!wi_result) {
        return std::unexpected(CoefficientError(
            std::format("Failed to compute production rates at eta={}: {}", 
                       station_index, wi_result.error().message())
        ));
    }
    const auto& wi_local = wi_result.value();
    
    // Get Jacobian d(wi)/d(rho_j)
    auto jac_result = mixture_.production_rate_jacobian(rho_species, inputs.T[station_index]);
    if (!jac_result) {
        return std::unexpected(CoefficientError(
            std::format("Failed to compute Jacobian at eta={}: {}", 
                       station_index, jac_result.error().message())
        ));
    }
    
    // Transform Jacobian from d(wi)/d(rho_j) to d(wi)/d(c_j)
    auto dwi_dc_result = transform_jacobian(
        station_index, inputs, thermo, bc, 
        c_local, rho_species, wi_local, jac_result.value()
    );
    
    if (!dwi_dc_result) {
        return std::unexpected(dwi_dc_result.error());
    }
    
    return std::make_pair(wi_local, dwi_dc_result.value());
}

auto CoefficientCalculator::extract_local_composition(
    std::size_t station_index,
    const core::Matrix<double>& c,
    double rho_total
) const -> std::pair<std::vector<double>, std::vector<double>> {
    
    const auto n_species = c.rows();
    std::vector<double> c_local(n_species);
    std::vector<double> rho_species(n_species);
    
    double sum_c = 0.0;
    for (std::size_t j = 0; j < n_species; ++j) {
        c_local[j] = c(j, station_index);
        sum_c += c_local[j];
    }
    
    // Normalize and compute partial densities
    for (std::size_t j = 0; j < n_species; ++j) {
        rho_species[j] = rho_total * c_local[j] / sum_c;
    }
    
    return {c_local, rho_species};
}

auto CoefficientCalculator::transform_jacobian(
    std::size_t station_index,
    const CoefficientInputs& inputs,
    const ThermodynamicCoefficients& thermo,
    const conditions::BoundaryConditions& bc,
    const std::vector<double>& c_local,
    const std::vector<double>& rho_species,
    const std::vector<double>& wi_local,
    const core::Matrix<double>& dwi_drho
) const -> std::expected<core::Matrix<double>, CoefficientError> {
    
    const auto n_species = c_local.size();
    const double T = inputs.T[station_index];
    const double P_edge = bc.P_e();
    
    // Get species enthalpies and Cp
    auto h_sp_result = mixture_.species_enthalpies(T);
    auto cp_result = mixture_.frozen_cp(c_local, T, P_edge);
    
    if (!h_sp_result || !cp_result) {
        return std::unexpected(CoefficientError("Failed to get enthalpies or Cp"));
    }
    
    const auto& h_species = h_sp_result.value();
    const double Cp = cp_result.value();
    
    // Compute dT/dc_j = -h_j/Cp
    std::vector<double> dT_dc(n_species);
    for (std::size_t j = 0; j < n_species; ++j) {
        dT_dc[j] = -h_species[j] / Cp;
    }
    
    // Compute d(rho_k)/d(c_j)
    auto drho_dc = compute_density_derivatives(
        rho_species, thermo.rho[station_index], thermo.MW[station_index]
    );
    
    // Compute d(wi)/dT
    auto dwi_dT_result = compute_temperature_derivatives(
        rho_species, T, wi_local
    );
    if (!dwi_dT_result) {
        return std::unexpected(CoefficientError("Failed to compute temperature derivatives"));
    }
    const auto& dwi_dT = dwi_dT_result.value();
    
    // Assemble complete Jacobian
    core::Matrix<double> dwi_dc(n_species, n_species);
    for (std::size_t i = 0; i < n_species; ++i) {
        for (std::size_t j = 0; j < n_species; ++j) {
            double sum = 0.0;
            for (std::size_t k = 0; k < n_species; ++k) {
                sum += dwi_drho(i, k) * drho_dc(k, j);
            }
            dwi_dc(i, j) = sum + dwi_dT[i] * dT_dc[j];
        }
    }
    
    return dwi_dc;
}

auto CoefficientCalculator::compute_density_derivatives(
    const std::vector<double>& rho_species,
    double rho_total,
    double MW_mixture
) const -> core::Matrix<double> {
    
    const auto n_species = rho_species.size();
    core::Matrix<double> drho_dc(n_species, n_species);
    
    for (std::size_t k = 0; k < n_species; ++k) {
        for (std::size_t j = 0; j < n_species; ++j) {
            drho_dc(k, j) = (k == j) ? rho_total : 0.0;
            drho_dc(k, j) -= rho_species[k] * MW_mixture / mixture_.species_molecular_weight(j);
        }
    }
    
    return drho_dc;
}

auto CoefficientCalculator::compute_temperature_derivatives(
    const std::vector<double>& rho_species,
    double T,
    const std::vector<double>& wi_base
) const -> std::expected<std::vector<double>, CoefficientError> {
    
    constexpr double eps = 1e-6;
    const double dT = eps * T;
    
    auto wi_pert_result = mixture_.production_rates(rho_species, T + dT);
    if (!wi_pert_result) {
        return std::unexpected(CoefficientError("Failed to compute perturbed rates"));
    }
    
    const auto& wi_pert = wi_pert_result.value();
    std::vector<double> dwi_dT(wi_base.size());
    
    for (std::size_t j = 0; j < wi_base.size(); ++j) {
        dwi_dT[j] = (wi_pert[j] - wi_base[j]) / dT;
    }
    
    return dwi_dT;
}

auto CoefficientCalculator::calculate_thermal_diffusion(
    const CoefficientInputs& inputs,
    const ThermodynamicCoefficients& thermo,
    const conditions::BoundaryConditions& bc
) const -> std::expected<ThermalDiffusionCoefficients, CoefficientError> {
    
    const auto n_eta = inputs.T.size();
    const auto n_species = inputs.c.rows();
    
    ThermalDiffusionCoefficients tdr;
    tdr.tdr = core::Matrix<double>(n_species, n_eta);
    
    // Use edge pressure consistently
    const double P_edge = bc.P_e();
    
    // Calculate thermal diffusion ratios
    for (std::size_t i = 0; i < n_eta; ++i) {
        std::vector<double> c_local(n_species);
        for (std::size_t j = 0; j < n_species; ++j) {
            c_local[j] = inputs.c(j, i);
        }
        
        // Get thermal diffusion ratios using edge pressure
        auto tdr_result = mixture_.thermal_diffusion_ratios(c_local, inputs.T[i], P_edge);
        if (!tdr_result) {
            return std::unexpected(CoefficientError(
                std::format("Failed to compute TDR at eta={}: {}", i, tdr_result.error().message())
            ));
        }
        
        const auto& tdr_local = tdr_result.value();
        for (std::size_t j = 0; j < n_species; ++j) {
            tdr.tdr(j, i) = tdr_local[j];
        }
    }
    
    // Compute TDR term for diffusion fluxes
    if (sim_config_.consider_thermal_diffusion) {
        // Compute dT/deta
        auto dT_deta = derivatives::compute_eta_derivative(inputs.T, d_eta_);
        
        // Compute TDR term
        tdr.tdr_term = core::Matrix<double>(n_eta, n_species);
        for (std::size_t i = 0; i < n_eta; ++i) {
            for (std::size_t j = 0; j < n_species; ++j) {
                tdr.tdr_term(i, j) = tdr.tdr(j, i) / inputs.T[i] * dT_deta[i];
            }
        }
    }
    
    return tdr;
}

auto CoefficientCalculator::calculate_species_enthalpies(
    const CoefficientInputs& inputs
) const -> std::expected<std::pair<core::Matrix<double>, core::Matrix<double>>, CoefficientError> {
    
    const auto n_eta = inputs.T.size();
    const auto n_species = mixture_.n_species();
    
    core::Matrix<double> h_species(n_species, n_eta);
    
    // Get species enthalpies at each temperature
    for (std::size_t i = 0; i < n_eta; ++i) {
        auto h_result = mixture_.species_enthalpies(inputs.T[i]);
        if (!h_result) {
            return std::unexpected(CoefficientError(
                std::format("Failed to get species enthalpies at eta={}: {}", i, h_result.error().message())
            ));
        }
        
        const auto& h_local = h_result.value();
        for (std::size_t j = 0; j < n_species; ++j) {
            h_species(j, i) = h_local[j];
        }
    }
    
    // Compute derivatives
    core::Matrix<double> dh_species_deta(n_species, n_eta);
    
    for (std::size_t j = 0; j < n_species; ++j) {
        std::vector<double> h_row(n_eta);
        for (std::size_t i = 0; i < n_eta; ++i) {
            h_row[i] = h_species(j, i);
        }
        
        auto dh = derivatives::compute_eta_derivative(h_row, d_eta_);
        for (std::size_t i = 0; i < n_eta; ++i) {
            dh_species_deta(j, i) = dh[i];
        }
    }
    
    return std::make_pair(std::move(h_species), std::move(dh_species_deta));
}

auto CoefficientCalculator::calculate_wall_properties(
    const CoefficientInputs& inputs,
    const conditions::BoundaryConditions& bc,
    const TransportCoefficients& transport,
    const ThermodynamicCoefficients& thermo
) const -> std::expected<WallProperties, CoefficientError> {
    
    const auto n_species = inputs.c.rows();
    
    // Extract wall composition
    std::vector<double> c_wall(n_species);
    for (std::size_t j = 0; j < n_species; ++j) {
        c_wall[j] = inputs.c(j, 0);
    }
    
    // Get wall transport properties using edge pressure
    auto k_result = mixture_.frozen_thermal_conductivity(c_wall, inputs.T[0], bc.P_e());
    if (!k_result) {
        return std::unexpected(CoefficientError("Failed to compute wall thermal conductivity"));
    }
    
    auto mu_result = mixture_.viscosity(c_wall, inputs.T[0], bc.P_e());
    if (!mu_result) {
        return std::unexpected(CoefficientError("Failed to compute wall viscosity"));
    }
    
    auto cp_result = mixture_.frozen_cp(c_wall, inputs.T[0], bc.P_e());
    if (!cp_result) {
        return std::unexpected(CoefficientError("Failed to compute wall specific heat"));
    }
    
    WallProperties wall;
    wall.k_wall = k_result.value();
    wall.mu_wall = mu_result.value();
    wall.Cp_wall = cp_result.value();
    wall.rho_wall = thermo.rho[0];
    wall.Pr_wall = (wall.k_wall > 0) ? wall.mu_wall * wall.Cp_wall / wall.k_wall : 0.72;
    
    return wall;
}

// Derivative computation implementations
namespace derivatives {

template<std::ranges::sized_range Range>
auto compute_eta_derivative(Range&& values, double d_eta) -> std::vector<double> {
    const auto n = std::ranges::size(values);
    std::vector<double> derivatives(n);
    
    if (n < 5) {
        // Simple finite differences for small arrays
        if (n >= 2) {
            derivatives[0] = (values[1] - values[0]) / d_eta;
            for (std::size_t i = 1; i < n - 1; ++i) {
                derivatives[i] = (values[i + 1] - values[i - 1]) / (2.0 * d_eta);
            }
            derivatives[n - 1] = (values[n - 1] - values[n - 2]) / d_eta;
        }
        return derivatives;
    }
    
    // 5-point stencil for higher accuracy
    const double dx12 = 12.0 * d_eta;
    
    // Forward difference at start
    derivatives[0] = (-25.0 * values[0] + 48.0 * values[1] - 36.0 * values[2] + 
                      16.0 * values[3] - 3.0 * values[4]) / dx12;
    derivatives[1] = (-3.0 * values[0] - 10.0 * values[1] + 18.0 * values[2] - 
                      6.0 * values[3] + 1.0 * values[4]) / dx12;
    
    // Central differences
    for (std::size_t i = 2; i < n - 2; ++i) {
        derivatives[i] = (values[i - 2] - 8.0 * values[i - 1] + 
                         8.0 * values[i + 1] - values[i + 2]) / dx12;
    }
    
    // Backward difference at end
    derivatives[n - 2] = (3.0 * values[n - 5] - 16.0 * values[n - 4] + 
                         36.0 * values[n - 3] - 48.0 * values[n - 2] + 
                         25.0 * values[n - 1]) / dx12;
    derivatives[n - 1] = derivatives[n - 2]; // Same as n-2 for stability
    
    return derivatives;
}

template<std::ranges::sized_range Range>
auto compute_eta_second_derivative(Range&& values, double d_eta) -> std::vector<double> {
    const auto n = std::ranges::size(values);
    std::vector<double> second_derivatives(n);
    
    if (n < 5) {
        // Simple finite differences for small arrays
        if (n >= 3) {
            second_derivatives[0] = (values[0] - 2.0 * values[1] + values[2]) / (d_eta * d_eta);
            for (std::size_t i = 1; i < n - 1; ++i) {
                second_derivatives[i] = (values[i-1] - 2.0 * values[i] + values[i+1]) / (d_eta * d_eta);
            }
            second_derivatives[n-1] = (values[n-3] - 2.0 * values[n-2] + values[n-1]) / (d_eta * d_eta);
        }
        return second_derivatives;
    }
    
    // 5-point stencil for higher accuracy: f''(x) = [-f(x-2h) + 16f(x-h) - 30f(x) + 16f(x+h) - f(x+2h)]/(12h²)
    const double dx2_12 = 12.0 * d_eta * d_eta;
    
    // Forward difference at boundaries
    second_derivatives[0] = (2.0 * values[0] - 5.0 * values[1] + 4.0 * values[2] - values[3]) / (d_eta * d_eta);
    second_derivatives[1] = (values[0] - 2.0 * values[1] + values[2]) / (d_eta * d_eta);
    
    // Central differences with 5-point stencil
    for (std::size_t i = 2; i < n - 2; ++i) {
        second_derivatives[i] = (-values[i-2] + 16.0 * values[i-1] - 30.0 * values[i] + 
                               16.0 * values[i+1] - values[i+2]) / dx2_12;
    }
    
    // Backward difference at boundaries
    second_derivatives[n-2] = (values[n-3] - 2.0 * values[n-2] + values[n-1]) / (d_eta * d_eta);
    second_derivatives[n-1] = (2.0 * values[n-1] - 5.0 * values[n-2] + 4.0 * values[n-3] - values[n-4]) / (d_eta * d_eta);
    
    return second_derivatives;
}

template<typename Matrix>
auto compute_matrix_eta_second_derivative(const Matrix& values, double d_eta) -> Matrix {
    const auto n_rows = values.rows();
    const auto n_cols = values.cols();
    Matrix result(n_rows, n_cols);
    
    for (std::size_t i = 0; i < n_rows; ++i) {
        std::vector<double> row_values(n_cols);
        for (std::size_t j = 0; j < n_cols; ++j) {
            row_values[j] = values(i, j);
        }
        auto row_derivatives = compute_eta_second_derivative(row_values, d_eta);
        for (std::size_t j = 0; j < n_cols; ++j) {
            result(i, j) = row_derivatives[j];
        }
    }
    
    return result;
}

// Explicit instantiations
template auto compute_eta_derivative(std::span<const double>&&, double) -> std::vector<double>;
template auto compute_eta_derivative(const std::vector<double>&, double) -> std::vector<double>;
template auto compute_eta_derivative(std::vector<double>&&, double) -> std::vector<double>;
template auto compute_eta_derivative(std::span<double>&&, double) -> std::vector<double>;

template auto compute_eta_second_derivative(std::span<const double>&&, double) -> std::vector<double>;
template auto compute_eta_second_derivative(const std::vector<double>&, double) -> std::vector<double>;
template auto compute_eta_second_derivative(std::vector<double>&&, double) -> std::vector<double>;
template auto compute_eta_second_derivative(std::vector<double>&, double) -> std::vector<double>;
template auto compute_eta_second_derivative(std::span<double>&&, double) -> std::vector<double>;

} // namespace derivatives

} // namespace blast::boundary_layer::coefficients