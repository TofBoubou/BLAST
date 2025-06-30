#include "blast/boundary_layer/coefficients/transport_coefficients.hpp"
#include "mutation++.h"
#include <algorithm>
#include <cmath>

// Access to global mixture object
extern Mutation::Mixture* mix_obj_ptr;
extern int n_species;

namespace blast::boundary_layer::coefficients {

namespace {

// Helper to compute derivatives
template<typename Container>
auto compute_derivative(const Container& f, double dx) -> std::vector<double> {
    const auto n = f.size();
    std::vector<double> df(n);
    
    if (n < 5) return df;
    
    const double dx12 = 12.0 * dx;
    
    // Fourth order at boundaries
    df[0] = (-25.0*f[0] + 48.0*f[1] - 36.0*f[2] + 16.0*f[3] - 3.0*f[4]) / dx12;
    df[1] = (-3.0*f[0] - 10.0*f[1] + 18.0*f[2] - 6.0*f[3] + f[4]) / dx12;
    
    // Central differences
    for (std::size_t i = 2; i < n-2; ++i) {
        df[i] = (f[i-2] - 8.0*f[i-1] + 8.0*f[i+1] - f[i+2]) / dx12;
    }
    
    df[n-2] = (3.0*f[n-5] - 16.0*f[n-4] + 36.0*f[n-3] - 48.0*f[n-2] + 25.0*f[n-1]) / dx12;
    df[n-1] = (3.0*f[n-5] - 16.0*f[n-4] + 36.0*f[n-3] - 48.0*f[n-2] + 25.0*f[n-1]) / dx12;
    
    return df;
}

// Get transport properties from Mutation++
auto get_viscosity(const std::vector<double>& c, double T, double P) -> double {
    double vars[2] = {P, T};
    mix_obj_ptr->setState(const_cast<double*>(c.data()), vars, 2);
    return mix_obj_ptr->viscosity();
}

auto get_thermal_conductivity(const std::vector<double>& c, double T, double P) -> double {
    double vars[2] = {P, T};
    mix_obj_ptr->setState(const_cast<double*>(c.data()), vars, 2);
    return mix_obj_ptr->frozenThermalConductivity();
}

auto get_heat_capacity(const std::vector<double>& c, double T, double P) -> double {
    double vars[2] = {P, T};
    mix_obj_ptr->setState(const_cast<double*>(c.data()), vars, 2);
    return mix_obj_ptr->mixtureFrozenCpMass();
}

auto get_binary_diffusion(double T, double P, int n_sp) -> core::Matrix<double> {
    double vars[2] = {P, T};
    mix_obj_ptr->setState(nullptr, vars, 2);
    
    core::Matrix<double> Dij(n_sp, n_sp);
    
    // Get from Mutation++
    Eigen::VectorXd Dij_vec = mix_obj_ptr->collisionDB().nDij();
    
    // Handle electrons if present
    int start_idx = mix_obj_ptr->hasElectrons() ? 1 : 0;
    
    if (mix_obj_ptr->hasElectrons()) {
        Eigen::VectorXd Dei_vec = mix_obj_ptr->collisionDB().nDei();
        for (int i = 0; i < n_sp; ++i) {
            Dij(0, i) = Dei_vec(i);
            Dij(i, 0) = Dei_vec(i);
        }
    }
    
    // Heavy species
    int cont = 0;
    for (int i = start_idx; i < n_sp; ++i) {
        for (int j = i; j < n_sp; ++j) {
            Dij(i, j) = Dij_vec(cont);
            Dij(j, i) = Dij_vec(cont);
            cont++;
        }
    }
    
    // Divide by number density
    const double n = P / (Mutation::KB * T);
    for (int i = 0; i < n_sp; ++i) {
        for (int j = 0; j < n_sp; ++j) {
            Dij(i, j) /= n;
        }
    }
    
    return Dij;
}

auto get_thermal_diffusion_ratios(const std::vector<double>& c, double T, double P) 
    -> std::vector<double> {
    std::vector<double> tdr(n_species);
    double vars[2] = {P, T};
    mix_obj_ptr->setState(const_cast<double*>(c.data()), vars, 2);
    mix_obj_ptr->heavyThermalDiffusionRatios(tdr.data());
    return tdr;
}

} // anonymous namespace

auto TransportCoefficients::compute(
    std::span<const double> temperature,
    std::span<const double> density,
    const std::vector<std::vector<double>>& mass_fractions,
    double pressure_edge,
    double density_edge,
    double viscosity_edge,
    double d_eta,
    bool consider_tdr
) -> std::expected<TransportCoefficients, CoefficientError> {
    
    const auto n_eta = temperature.size();
    const auto n_sp = mass_fractions.size();
    
    if (n_eta == 0 || n_sp == 0) {
        return std::unexpected(CoefficientError("Invalid dimensions"));
    }
    
    TransportProperties props;
    props.viscosity.reserve(n_eta);
    props.thermal_cond.reserve(n_eta);
    props.prandtl.reserve(n_eta);
    props.l0.reserve(n_eta);
    props.l3.reserve(n_eta);
    props.binary_diffusion.reserve(n_eta);
    
    // Extract mass fractions at each eta
    std::vector<double> c_tmp(n_sp);
    std::vector<double> Cpf_vec(n_eta);
    
    // Compute properties at each eta point
    for (std::size_t i = 0; i < n_eta; ++i) {
        // Extract species at this eta
        for (std::size_t j = 0; j < n_sp; ++j) {
            c_tmp[j] = mass_fractions[j][i];
        }
        
        // Local pressure accounting for composition
        double R_mix_local = 0.0;
        for (std::size_t j = 0; j < n_sp; ++j) {
            double R_j = Mutation::RU / mix_obj_ptr->speciesMw(j);
            R_mix_local += c_tmp[j] * R_j;
        }
        double P_local = density[i] * R_mix_local * temperature[i];
        
        // Get transport properties
        double mu = get_viscosity(c_tmp, temperature[i], P_local);
        double k_fr = get_thermal_conductivity(c_tmp, temperature[i], P_local);
        double Cpf = get_heat_capacity(c_tmp, temperature[i], P_local);
        double Pr = mu * Cpf / k_fr;
        
        props.viscosity.push_back(mu);
        props.thermal_cond.push_back(k_fr);
        props.prandtl.push_back(Pr);
        Cpf_vec[i] = Cpf;
        
        // Compute l0, l3
        props.l0.push_back(density[i] * mu / (density_edge * viscosity_edge));
        props.l3.push_back(props.l0.back() / Pr);
        
        // Binary diffusion
        props.binary_diffusion.emplace_back(get_binary_diffusion(temperature[i], P_local, n_sp));
    }
    
    // Store wall properties
    props.k_wall = props.thermal_cond[0];
    props.mu_wall = props.viscosity[0];
    props.Pr_wall = props.prandtl[0];
    
    // Compute derivatives
    props.dl0_deta = compute_derivative(props.l0, d_eta);
    props.dl3_deta = compute_derivative(props.l3, d_eta);
    
    // Thermal diffusion ratios if needed
    if (consider_tdr) {
        props.thermal_diffusion_ratios.resize(n_sp);
        props.d_tdr_deta.resize(n_sp);
        props.tdr_term.resize(n_eta);
        
        std::vector<std::vector<double>> tdr_transpose(n_eta);
        
        for (std::size_t i = 0; i < n_eta; ++i) {
            for (std::size_t j = 0; j < n_sp; ++j) {
                c_tmp[j] = mass_fractions[j][i];
            }
            
            // Get local pressure
            double R_mix_local = 0.0;
            for (std::size_t j = 0; j < n_sp; ++j) {
                double R_j = Mutation::RU / mix_obj_ptr->speciesMw(j);
                R_mix_local += c_tmp[j] * R_j;
            }
            double P_local = density[i] * R_mix_local * temperature[i];
            
            auto tdr = get_thermal_diffusion_ratios(c_tmp, temperature[i], P_local);
            tdr_transpose[i] = tdr;
            
            // Store in correct format
            for (std::size_t j = 0; j < n_sp; ++j) {
                if (i == 0) props.thermal_diffusion_ratios[j].resize(n_eta);
                props.thermal_diffusion_ratios[j][i] = tdr[j];
            }
        }
        
        // Compute temperature derivative
        auto dT_deta = compute_derivative(temperature, d_eta);
        
        // Compute tdr_term
        for (std::size_t i = 0; i < n_eta; ++i) {
            props.tdr_term[i].resize(n_sp);
            for (std::size_t j = 0; j < n_sp; ++j) {
                props.tdr_term[i][j] = props.thermal_diffusion_ratios[j][i] / 
                                      temperature[i] * dT_deta[i];
            }
        }
        
        // Compute derivatives of tdr
        for (std::size_t i = 0; i < n_sp; ++i) {
            props.d_tdr_deta[i] = compute_derivative(props.thermal_diffusion_ratios[i], d_eta);
        }
    }
    
    return TransportCoefficients(std::move(props), n_eta, n_sp);
}

auto TransportCoefficients::binary_diffusion(std::size_t eta_idx) const 
    -> const core::Matrix<double>& {
    if (eta_idx >= properties_.binary_diffusion.size()) {
        throw CoefficientError("Binary diffusion index out of range");
    }
    return properties_.binary_diffusion[eta_idx];
}

} // namespace blast::boundary_layer::coefficients