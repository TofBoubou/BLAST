#pragma once
#include "transport_coefficients.hpp"
#include "thermodynamic_coefficients.hpp"
#include "chemical_coefficients.hpp"
#include "../../boundary_layer/conditions/boundary_conditions.hpp"
#include "../../io/config_types.hpp"
#include <memory>
#include <expected>

namespace blast::boundary_layer::coefficients {

// Main facade class that combines all coefficient types
class BoundaryLayerCoefficients {
private:
    std::unique_ptr<TransportCoefficients> transport_;
    std::unique_ptr<ThermodynamicCoefficients> thermodynamic_;
    std::unique_ptr<ChemicalCoefficients> chemical_;
    
    // Computed auxiliary values
    std::vector<double> y_continuity_;
    
    int n_eta_;
    int n_species_;
    double d_eta_;
    
public:
    // Main factory function
    [[nodiscard]] static auto compute(
        double xi,
        std::span<const double> F,
        const std::vector<std::vector<double>>& mass_fractions,
        const std::vector<std::vector<double>>& dc_deta,
        std::span<const double> temperature,
        const conditions::BoundaryConditions& boundary_conditions,
        const io::SimulationConfig& sim_config,
        double d_eta,
        // Xi derivatives
        double lambda0,
        std::span<const double> F_derivative
    ) -> std::expected<BoundaryLayerCoefficients, CoefficientError>;
    
    // Accessors for sub-coefficients
    [[nodiscard]] auto transport() const noexcept -> const TransportCoefficients& { 
        return *transport_; 
    }
    [[nodiscard]] auto thermodynamic() const noexcept -> const ThermodynamicCoefficients& { 
        return *thermodynamic_; 
    }
    [[nodiscard]] auto chemical() const noexcept -> const ChemicalCoefficients& { 
        return *chemical_; 
    }
    
    // Direct accessors for commonly used values
    [[nodiscard]] auto y() const noexcept -> std::span<const double> { 
        return y_continuity_; 
    }
    [[nodiscard]] auto rho() const noexcept -> std::span<const double> { 
        return thermodynamic_->density(); 
    }
    [[nodiscard]] auto MW() const noexcept -> std::span<const double> { 
        return thermodynamic_->molecular_weight(); 
    }
    
    // Legacy compatibility accessors
    [[nodiscard]] auto l0() const noexcept { return transport_->l0(); }
    [[nodiscard]] auto l3() const noexcept { return transport_->l3(); }
    [[nodiscard]] auto wi(std::size_t eta_idx) const -> std::span<const double>;
    [[nodiscard]] auto Dij_bin(std::size_t eta_idx) const { 
        return transport_->binary_diffusion(eta_idx); 
    }
    
private:
    BoundaryLayerCoefficients(
        std::unique_ptr<TransportCoefficients> transport,
        std::unique_ptr<ThermodynamicCoefficients> thermo,
        std::unique_ptr<ChemicalCoefficients> chem,
        std::vector<double> y_cont,
        int n_eta,
        int n_species,
        double d_eta
    ) : transport_(std::move(transport)),
        thermodynamic_(std::move(thermo)),
        chemical_(std::move(chem)),
        y_continuity_(std::move(y_cont)),
        n_eta_(n_eta),
        n_species_(n_species),
        d_eta_(d_eta) {}
};

} // namespace blast::boundary_layer::coefficients