#pragma once
#include "coefficient_concepts.hpp"
#include <vector>
#include <span>
#include <expected>

namespace blast::boundary_layer::coefficients {

struct ThermodynamicProperties {
    std::vector<double> density;           // rho at each eta
    std::vector<double> molecular_weight;  // MW at each eta
    std::vector<double> d_MW_deta;         // derivative of MW
    std::vector<double> d_rho_deta;        // derivative of density
    std::vector<double> heat_capacity;     // Cp frozen at each eta
    
    // Species enthalpies [species][eta]
    std::vector<std::vector<double>> species_enthalpies;
    std::vector<std::vector<double>> dh_sp_deta;
    
    // Wall values
    double h_wall{};
    double rho_wall{};
    double Cp_wall{};
};

class ThermodynamicCoefficients {
private:
    ThermodynamicProperties properties_;
    int n_eta_;
    int n_species_;
    
public:
    [[nodiscard]] static auto compute(
        std::span<const double> temperature,
        const std::vector<std::vector<double>>& mass_fractions,
        double pressure,
        double wall_temperature,
        std::span<const double> wall_mass_fractions,
        double d_eta
    ) -> std::expected<ThermodynamicCoefficients, CoefficientError>;
    
    // Accessors
    [[nodiscard]] constexpr auto density() const noexcept -> std::span<const double> { 
        return properties_.density; 
    }
    [[nodiscard]] constexpr auto molecular_weight() const noexcept -> std::span<const double> { 
        return properties_.molecular_weight; 
    }
    [[nodiscard]] constexpr auto d_MW_deta() const noexcept -> std::span<const double> { 
        return properties_.d_MW_deta; 
    }
    [[nodiscard]] constexpr auto d_rho_deta() const noexcept -> std::span<const double> { 
        return properties_.d_rho_deta; 
    }
    
    [[nodiscard]] auto species_enthalpies() const noexcept 
        -> const std::vector<std::vector<double>>& {
        return properties_.species_enthalpies;
    }
    
    [[nodiscard]] auto dh_sp_deta() const noexcept 
        -> const std::vector<std::vector<double>>& {
        return properties_.dh_sp_deta;
    }
    
    [[nodiscard]] constexpr auto wall_enthalpy() const noexcept -> double { 
        return properties_.h_wall; 
    }
    [[nodiscard]] constexpr auto wall_density() const noexcept -> double { 
        return properties_.rho_wall; 
    }
    [[nodiscard]] constexpr auto wall_heat_capacity() const noexcept -> double { 
        return properties_.Cp_wall; 
    }
    
private:
    ThermodynamicCoefficients(ThermodynamicProperties props, int n_eta, int n_species)
        : properties_(std::move(props)), n_eta_(n_eta), n_species_(n_species) {}
};

} // namespace blast::boundary_layer::coefficients