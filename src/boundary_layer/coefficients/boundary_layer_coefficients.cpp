#include "blast/boundary_layer/coefficients/boundary_layer_coefficients.hpp"
#include <algorithm>

extern int n_eta;
extern int n_species;

namespace blast::boundary_layer::coefficients {

auto BoundaryLayerCoefficients::compute(
    double xi,
    std::span<const double> F,
    const std::vector<std::vector<double>>& mass_fractions,
    const std::vector<std::vector<double>>& dc_deta,
    std::span<const double> temperature,
    const conditions::BoundaryConditions& boundary_conditions,
    const io::SimulationConfig& sim_config,
    double d_eta,
    double lambda0,
    std::span<const double> F_derivative
) -> std::expected<BoundaryLayerCoefficients, CoefficientError> {
    
    const auto n_eta_val = temperature.size();
    const auto n_sp = mass_fractions.size();
    
    // Validate inputs
    if (n_eta_val == 0 || n_sp == 0) {
        return std::unexpected(CoefficientError("Invalid grid dimensions"));
    }
    
    if (F.size() != n_eta_val || F_derivative.size() != n_eta_val) {
        return std::unexpected(CoefficientError("Inconsistent array sizes"));
    }
    
    // 1. Compute thermodynamic coefficients
    std::vector<double> wall_c(n_sp);
    for (std::size_t i = 0; i < n_sp; ++i) {
        wall_c[i] = mass_fractions[i][0];
    }
    
    auto thermo_result = ThermodynamicCoefficients::compute(
        temperature,
        mass_fractions,
        boundary_conditions.P_e(),
        boundary_conditions.Tw(),
        wall_c,
        d_eta
    );
    
    if (!thermo_result) {
        return std::unexpected(CoefficientError(
            "Failed to compute thermodynamic coefficients: " + 
            thermo_result.error().what()
        ));
    }
    
    auto thermo = std::make_unique<ThermodynamicCoefficients>(
        std::move(thermo_result.value())
    );
    
    // 2. Compute transport coefficients
    auto transport_result = TransportCoefficients::compute(
        temperature,
        thermo->density(),
        mass_fractions,
        boundary_conditions.P_e(),
        boundary_conditions.rho_e(),
        boundary_conditions.mu_e(),
        d_eta,
        sim_config.consider_thermal_diffusion
    );
    
    if (!transport_result) {
        return std::unexpected(CoefficientError(
            "Failed to compute transport coefficients: " + 
            transport_result.error().what()
        ));
    }
    
    auto transport = std::make_unique<TransportCoefficients>(
        std::move(transport_result.value())
    );
    
    // 3. Compute chemical coefficients
    auto chemical_result = ChemicalCoefficients::compute(
        temperature,
        mass_fractions,
        thermo->density(),
        thermo->molecular_weight(),
        thermo->species_enthalpies(),
        std::span(thermo->density()), // Using density as placeholder for Cp
        sim_config.chemical_non_equilibrium,
        false // frozen_boundary_layer removed as requested
    );
    
    if (!chemical_result) {
        return std::unexpected(CoefficientError(
            "Failed to compute chemical coefficients: " + 
            chemical_result.error().what()
        ));
    }
    
    auto chemical = std::make_unique<ChemicalCoefficients>(
        std::move(chemical_result.value())
    );
    
    // 4. Compute y for continuity equation
    std::vector<double> y_continuity(n_eta_val);
    for (std::size_t i = 0; i < n_eta_val; ++i) {
        y_continuity[i] = -(2*xi*lambda0 + 1.0)*F[i] - 2*xi*F_derivative[i];
    }
    
    return BoundaryLayerCoefficients(
        std::move(transport),
        std::move(thermo),
        std::move(chemical),
        std::move(y_continuity),
        n_eta_val,
        n_sp,
        d_eta
    );
}

auto BoundaryLayerCoefficients::wi(std::size_t eta_idx) const 
    -> std::span<const double> {
    const auto& rates = chemical_->mass_production_rates();
    if (eta_idx >= rates.size()) {
        throw CoefficientError("Mass production rate index out of range");
    }
    return rates[eta_idx];
}

} // namespace blast::boundary_layer::coefficients