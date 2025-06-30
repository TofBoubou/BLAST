#include "blast/boundary_layer/coefficients/thermodynamic_coefficients.hpp"
#include "mutation++.h"
#include <algorithm>
#include <numeric>

extern Mutation::Mixture* mix_obj_ptr;
extern int n_species;

namespace blast::boundary_layer::coefficients {

namespace {

template<typename Container>
auto compute_derivative(const Container& f, double dx) -> std::vector<double> {
    const auto n = f.size();
    std::vector<double> df(n);
    
    if (n < 5) return df;
    
    const double dx12 = 12.0 * dx;
    
    df[0] = (-25.0*f[0] + 48.0*f[1] - 36.0*f[2] + 16.0*f[3] - 3.0*f[4]) / dx12;
    df[1] = (-3.0*f[0] - 10.0*f[1] + 18.0*f[2] - 6.0*f[3] + f[4]) / dx12;
    
    for (std::size_t i = 2; i < n-2; ++i) {
        df[i] = (f[i-2] - 8.0*f[i-1] + 8.0*f[i+1] - f[i+2]) / dx12;
    }
    
    df[n-2] = (3.0*f[n-5] - 16.0*f[n-4] + 36.0*f[n-3] - 48.0*f[n-2] + 25.0*f[n-1]) / dx12;
    df[n-1] = (3.0*f[n-5] - 16.0*f[n-4] + 36.0*f[n-3] - 48.0*f[n-2] + 25.0*f[n-1]) / dx12;
    
    return df;
}

auto get_molecular_weight(const std::vector<double>& c) -> double {
    double MW = 0.0;
    for (int i = 0; i < n_species; ++i) {
        MW += c[i] / mix_obj_ptr->speciesMw(i);
    }
    return 1.0 / MW;
}

auto get_dMW_deta(const std::vector<double>& dc_deta, double M) -> double {
    double dMW_deta = 0.0;
    for (int i = 0; i < n_species; ++i) {
        dMW_deta += dc_deta[i] / mix_obj_ptr->speciesMw(i);
    }
    return -dMW_deta * M * M;
}

auto get_enthalpy(const std::vector<double>& c, double T, double P) -> double {
    double vars[2] = {P, T};
    mix_obj_ptr->setState(const_cast<double*>(c.data()), vars, 2);
    return mix_obj_ptr->mixtureHMass();
}

auto get_species_enthalpies(double T, int n_sp) -> std::vector<double> {
    std::vector<double> h_sp(n_sp);
    mix_obj_ptr->speciesHOverRT(T, h_sp.data());
    
    // Convert to J/kg
    for (int i = 0; i < n_sp; ++i) {
        h_sp[i] *= Mutation::RU * T / mix_obj_ptr->speciesMw(i);
    }
    return h_sp;
}

auto get_heat_capacity(const std::vector<double>& c, double T, double P) -> double {
    double vars[2] = {P, T};
    mix_obj_ptr->setState(const_cast<double*>(c.data()), vars, 2);
    return mix_obj_ptr->mixtureFrozenCpMass();
}

} // anonymous namespace

auto ThermodynamicCoefficients::compute(
    std::span<const double> temperature,
    const std::vector<std::vector<double>>& mass_fractions,
    double pressure,
    double wall_temperature,
    std::span<const double> wall_mass_fractions,
    double d_eta
) -> std::expected<ThermodynamicCoefficients, CoefficientError> {
    
    const auto n_eta = temperature.size();
    const auto n_sp = mass_fractions.size();
    
    if (n_eta == 0 || n_sp == 0) {
        return std::unexpected(CoefficientError("Invalid dimensions"));
    }
    
    ThermodynamicProperties props;
    props.density.reserve(n_eta);
    props.molecular_weight.reserve(n_eta);
    props.heat_capacity.reserve(n_eta);
    props.species_enthalpies.resize(n_sp);
    props.dh_sp_deta.resize(n_sp);
    
    // Wall enthalpy
    std::vector<double> c_wall(wall_mass_fractions.begin(), wall_mass_fractions.end());
    props.h_wall = get_enthalpy(c_wall, wall_temperature, pressure);
    
    // Properties at each eta
    std::vector<double> c_tmp(n_sp);
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        // Extract species
        for (std::size_t j = 0; j < n_sp; ++j) {
            c_tmp[j] = mass_fractions[j][i];
        }
        
        // Molecular weight
        double MW = get_molecular_weight(c_tmp);
        props.molecular_weight.push_back(MW);
        
        // Density using ideal gas
        props.density.push_back(pressure * MW / temperature[i] / Mutation::RU);
        
        // Heat capacity
        props.heat_capacity.push_back(get_heat_capacity(c_tmp, temperature[i], pressure));
        
        // Species enthalpies
        auto h_sp = get_species_enthalpies(temperature[i], n_sp);
        for (std::size_t j = 0; j < n_sp; ++j) {
            if (i == 0) props.species_enthalpies[j].resize(n_eta);
            props.species_enthalpies[j][i] = h_sp[j];
        }
    }
    
    // Wall properties
    props.rho_wall = props.density[0];
    props.Cp_wall = props.heat_capacity[0];
    
    // Compute MW derivative
    props.d_MW_deta.resize(n_eta);
    std::vector<double> dc_deta_tmp(n_sp);
    std::vector<std::vector<double>> dc_deta(n_sp);
    
    // Compute dc_deta
    for (std::size_t j = 0; j < n_sp; ++j) {
        dc_deta[j] = compute_derivative(mass_fractions[j], d_eta);
    }
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        for (std::size_t j = 0; j < n_sp; ++j) {
            dc_deta_tmp[j] = dc_deta[j][i];
        }
        props.d_MW_deta[i] = get_dMW_deta(dc_deta_tmp, props.molecular_weight[i]);
    }
    
    // Density derivative
    props.d_rho_deta = compute_derivative(props.density, d_eta);
    
    // Species enthalpy derivatives
    for (std::size_t i = 0; i < n_sp; ++i) {
        props.dh_sp_deta[i] = compute_derivative(props.species_enthalpies[i], d_eta);
    }
    
    return ThermodynamicCoefficients(std::move(props), n_eta, n_sp);
}

} // namespace blast::boundary_layer::coefficients