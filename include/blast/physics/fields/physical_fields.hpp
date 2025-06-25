#pragma once
#include "../../core/containers.hpp"
#include <concepts>

namespace blast::physics::fields {

// Field type concepts
template<typename T>
concept ScalarField = std::same_as<T, core::Vector<double>>;

template<typename T>  
concept SpeciesField = std::same_as<T, core::Matrix<double>>;

// Primary solution fields
using TemperatureField = core::Vector<double>;        // T(η)
using VelocityField = core::Vector<double>;           // F(η) - non-dimensional
using EnthalpyField = core::Vector<double>;           // g(η) - non-dimensional  
using DensityField = core::Vector<double>;            // ρ(η)
using ViscosityField = core::Vector<double>;          // μ(η)
using SpeciesField = core::Matrix<double>;            // c_i(η) - [n_species x n_eta]

// Complete solution state
struct SolutionFields {
    TemperatureField temperature;
    VelocityField velocity;
    EnthalpyField enthalpy;
    SpeciesField species_fractions;
    DensityField density;
    ViscosityField viscosity;
    
    explicit SolutionFields(std::size_t n_eta, std::size_t n_species) 
        : temperature(n_eta)
        , velocity(n_eta)
        , enthalpy(n_eta)
        , species_fractions(n_species, n_eta)
        , density(n_eta)
        , viscosity(n_eta) {}
        
    SolutionFields(const SolutionFields&) = default;
    SolutionFields(SolutionFields&&) = default;
    SolutionFields& operator=(const SolutionFields&) = default;
    SolutionFields& operator=(SolutionFields&&) = default;
    
    [[nodiscard]] constexpr auto n_eta() const noexcept -> std::size_t { 
        return temperature.size(); 
    }
    
    [[nodiscard]] constexpr auto n_species() const noexcept -> std::size_t { 
        return species_fractions.rows(); 
    }
};

// Auxiliary fields for derivatives and fluxes
struct DerivativeFields {
    core::Vector<double> dF_deta;                     // dF/dη
    core::Matrix<double> dc_deta;                     // dc_i/dη
    core::Matrix<double> dc_deta2;                    // d²c_i/dη²
    
    explicit DerivativeFields(std::size_t n_eta, std::size_t n_species)
        : dF_deta(n_eta)
        , dc_deta(n_species, n_eta) 
        , dc_deta2(n_species, n_eta) {}
};

struct FluxFields {
    core::Matrix<double> diffusion_fluxes;            // J_i(η)
    core::Matrix<double> dJ_deta;                     // dJ_i/dη
    core::Vector<double> normal_velocity;             // V(η)
    
    explicit FluxFields(std::size_t n_eta, std::size_t n_species)
        : diffusion_fluxes(n_species, n_eta)
        , dJ_deta(n_species, n_eta)
        , normal_velocity(n_eta) {}
};

} // namespace blast::physics::fields