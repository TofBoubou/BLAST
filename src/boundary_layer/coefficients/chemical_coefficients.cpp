#include "blast/boundary_layer/coefficients/chemical_coefficients.hpp"
#include "mutation++.h"
#include <cmath>

extern Mutation::Mixture* mix_obj_ptr;
extern int n_species;

namespace blast::boundary_layer::coefficients {

namespace {

constexpr double EPSILON_FD = 1e-6;

auto compute_mass_production_rates(
    double T, 
    const std::vector<double>& c, 
    double rho,
    double MW,
    const std::vector<double>& h_sp,
    double Cp,
    std::vector<double>& wi,
    core::Matrix<double>& dwi_dc
) -> void {
    
    const auto n_sp = c.size();
    
    // Compute species densities
    std::vector<double> rho_sp(n_sp);
    double sum_c = std::accumulate(c.begin(), c.end(), 0.0);
    
    for (std::size_t i = 0; i < n_sp; ++i) {
        rho_sp[i] = rho * c[i] / sum_c;
    }
    
    // Set state and get production rates
    mix_obj_ptr->setState(rho_sp.data(), &T, 1);
    mix_obj_ptr->netProductionRates(wi.data());
    
    // Compute Jacobian
    std::vector<double> dT_dc_j(n_sp);
    core::Matrix<double> drho_k_dc_j(n_sp, n_sp);
    
    // dT/dc_j = -h_j/Cp
    for (std::size_t i = 0; i < n_sp; ++i) {
        dT_dc_j[i] = -h_sp[i] / Cp;
    }
    
    // drho_k/dc_j
    for (std::size_t i = 0; i < n_sp; ++i) {
        for (std::size_t j = 0; j < n_sp; ++j) {
            drho_k_dc_j(i,j) = (i == j ? 1.0 : 0.0) * rho - 
                               rho_sp[i] * MW / mix_obj_ptr->speciesMw(j);
        }
    }
    
    // Get Jacobian w.r.t. density
    std::vector<double> jacobian_rho(n_sp * n_sp);
    mix_obj_ptr->jacobianRho(jacobian_rho.data());
    
    core::Matrix<double> dwi_drho_j(n_sp, n_sp);
    for (std::size_t i = 0; i < n_sp; ++i) {
        for (std::size_t j = 0; j < n_sp; ++j) {
            dwi_drho_j(i,j) = jacobian_rho[i * n_sp + j];
        }
    }
    
    // Compute dwi/dT by finite difference
    std::vector<double> dwi_dT(n_sp);
    std::vector<double> wi_perturbed(n_sp);
    double T_new = T + EPSILON_FD * T;
    
    mix_obj_ptr->setState(rho_sp.data(), &T_new, 1);
    mix_obj_ptr->netProductionRates(wi_perturbed.data());
    
    for (std::size_t i = 0; i < n_sp; ++i) {
        dwi_dT[i] = (wi_perturbed[i] - wi[i]) / (EPSILON_FD * T);
    }
    
    // Combine to get full Jacobian
    for (std::size_t i = 0; i < n_sp; ++i) {
        for (std::size_t j = 0; j < n_sp; ++j) {
            double tmp_sum = 0.0;
            for (std::size_t k = 0; k < n_sp; ++k) {
                tmp_sum += dwi_drho_j(i,k) * drho_k_dc_j(k,j);
            }
            dwi_dc(i,j) = tmp_sum + dwi_dT[i] * dT_dc_j[j];
        }
    }
}

} // anonymous namespace

auto ChemicalCoefficients::compute(
    std::span<const double> temperature,
    const std::vector<std::vector<double>>& mass_fractions,
    std::span<const double> density,
    std::span<const double> molecular_weight,
    const std::vector<std::vector<double>>& species_enthalpies,
    std::span<const double> heat_capacity,
    bool chemical_non_equilibrium,
    bool frozen_boundary_layer
) -> std::expected<ChemicalCoefficients, CoefficientError> {
    
    const auto n_eta = temperature.size();
    const auto n_sp = mass_fractions.size();
    
    if (n_eta == 0 || n_sp == 0) {
        return std::unexpected(CoefficientError("Invalid dimensions"));
    }
    
    const bool is_frozen = frozen_boundary_layer || !chemical_non_equilibrium;
    
    ChemicalProperties props;
    props.mass_production_rates.resize(n_eta);
    props.production_jacobian.reserve(n_eta);
    
    if (is_frozen) {
        // Zero production rates
        for (std::size_t i = 0; i < n_eta; ++i) {
            props.mass_production_rates[i].resize(n_sp, 0.0);
            props.production_jacobian.emplace_back(n_sp, n_sp);
            props.production_jacobian.back().eigen().setZero();
        }
    } else {
        // Compute actual rates
        std::vector<double> c_tmp(n_sp);
        std::vector<double> h_sp_tmp(n_sp);
        std::vector<double> wi_tmp(n_sp);
        
        for (std::size_t i = 0; i < n_eta; ++i) {
            // Extract data at this eta
            for (std::size_t j = 0; j < n_sp; ++j) {
                c_tmp[j] = mass_fractions[j][i];
                h_sp_tmp[j] = species_enthalpies[j][i];
            }
            
            core::Matrix<double> dwi_dc_tmp(n_sp, n_sp);
            
            compute_mass_production_rates(
                temperature[i], c_tmp, density[i], molecular_weight[i],
                h_sp_tmp, heat_capacity[i], wi_tmp, dwi_dc_tmp
            );
            
            props.mass_production_rates[i] = wi_tmp;
            props.production_jacobian.emplace_back(std::move(dwi_dc_tmp));
        }
    }
    
    return ChemicalCoefficients(std::move(props), n_eta, n_sp, is_frozen);
}

auto ChemicalCoefficients::production_jacobian(std::size_t eta_idx) const 
    -> const core::Matrix<double>& {
    if (eta_idx >= properties_.production_jacobian.size()) {
        throw CoefficientError("Production jacobian index out of range");
    }
    return properties_.production_jacobian[eta_idx];
}

} // namespace blast::boundary_layer::coefficients