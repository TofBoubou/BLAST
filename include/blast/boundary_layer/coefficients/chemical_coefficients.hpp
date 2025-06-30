#pragma once
#include "coefficient_concepts.hpp"
#include "../../core/containers.hpp"
#include <vector>
#include <span>
#include <expected>

namespace blast::boundary_layer::coefficients {

struct ChemicalProperties {
    // Mass production rates [eta][species]
    std::vector<std::vector<double>> mass_production_rates;
    
    // Jacobian d_wi_dc [eta][species][species]
    std::vector<core::Matrix<double>> production_jacobian;
};

class ChemicalCoefficients {
private:
    ChemicalProperties properties_;
    int n_eta_;
    int n_species_;
    bool is_frozen_;
    
public:
    [[nodiscard]] static auto compute(
        std::span<const double> temperature,
        const std::vector<std::vector<double>>& mass_fractions,
        std::span<const double> density,
        std::span<const double> molecular_weight,
        const std::vector<std::vector<double>>& species_enthalpies,
        std::span<const double> heat_capacity,
        bool chemical_non_equilibrium,
        bool frozen_boundary_layer
    ) -> std::expected<ChemicalCoefficients, CoefficientError>;
    
    // Accessors
    [[nodiscard]] auto mass_production_rates() const noexcept 
        -> const std::vector<std::vector<double>>& {
        return properties_.mass_production_rates;
    }
    
    [[nodiscard]] auto production_jacobian(std::size_t eta_idx) const 
        -> const core::Matrix<double>&;
    
    [[nodiscard]] constexpr auto is_frozen() const noexcept -> bool { 
        return is_frozen_; 
    }
    
private:
    ChemicalCoefficients(ChemicalProperties props, int n_eta, int n_species, bool frozen)
        : properties_(std::move(props)), n_eta_(n_eta), n_species_(n_species), 
          is_frozen_(frozen) {}
};

} // namespace blast::boundary_layer::coefficients