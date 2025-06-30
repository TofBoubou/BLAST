#pragma once
#include "coefficient_concepts.hpp"
#include "../../core/containers.hpp"
#include <vector>
#include <span>
#include <memory>

namespace blast::boundary_layer::coefficients {

// Immutable transport property storage
struct TransportProperties {
    std::vector<double> viscosity;        // mu at each eta point
    std::vector<double> thermal_cond;     // k_fr at each eta point  
    std::vector<double> prandtl;          // Pr at each eta point
    std::vector<double> l0;               // momentum coefficient
    std::vector<double> l3;               // energy coefficient
    std::vector<double> dl0_deta;         // derivative of l0
    std::vector<double> dl3_deta;         // derivative of l3
    
    // Binary diffusion coefficients [eta][species][species]
    std::vector<core::Matrix<double>> binary_diffusion;
    
    // Thermal diffusion ratios [species][eta]
    std::vector<std::vector<double>> thermal_diffusion_ratios;
    std::vector<std::vector<double>> d_tdr_deta;
    std::vector<std::vector<double>> tdr_term;
    
    // Wall properties
    double k_wall{};
    double mu_wall{};
    double Pr_wall{};
};

class TransportCoefficients {
private:
    TransportProperties properties_;
    int n_eta_;
    int n_species_;
    
public:
    // Factory function for creating transport coefficients
    [[nodiscard]] static auto compute(
        std::span<const double> temperature,
        std::span<const double> density,
        const std::vector<std::vector<double>>& mass_fractions,
        double pressure_edge,
        double density_edge,
        double viscosity_edge,
        double d_eta,
        bool consider_tdr = false
    ) -> std::expected<TransportCoefficients, CoefficientError>;
    
    // Accessors
    [[nodiscard]] constexpr auto viscosity() const noexcept -> std::span<const double> { 
        return properties_.viscosity; 
    }
    [[nodiscard]] constexpr auto thermal_conductivity() const noexcept -> std::span<const double> { 
        return properties_.thermal_cond; 
    }
    [[nodiscard]] constexpr auto l0() const noexcept -> std::span<const double> { 
        return properties_.l0; 
    }
    [[nodiscard]] constexpr auto l3() const noexcept -> std::span<const double> { 
        return properties_.l3; 
    }
    [[nodiscard]] constexpr auto dl0_deta() const noexcept -> std::span<const double> { 
        return properties_.dl0_deta; 
    }
    [[nodiscard]] constexpr auto dl3_deta() const noexcept -> std::span<const double> { 
        return properties_.dl3_deta; 
    }
    
    [[nodiscard]] auto binary_diffusion(std::size_t eta_idx) const 
        -> const core::Matrix<double>&;
    
    [[nodiscard]] auto thermal_diffusion_ratios() const noexcept 
        -> std::span<const std::vector<double>> {
        return properties_.thermal_diffusion_ratios;
    }
    
    [[nodiscard]] auto tdr_term() const noexcept 
        -> const std::vector<std::vector<double>>& {
        return properties_.tdr_term;
    }
    
    [[nodiscard]] constexpr auto wall_properties() const noexcept {
        return std::tuple{properties_.k_wall, properties_.mu_wall, properties_.Pr_wall};
    }
    
private:
    TransportCoefficients(TransportProperties props, int n_eta, int n_species)
        : properties_(std::move(props)), n_eta_(n_eta), n_species_(n_species) {}
};

} // namespace blast::boundary_layer::coefficients