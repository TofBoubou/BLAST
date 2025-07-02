#pragma once
#include "coefficient_types.hpp"
#include "xi_derivatives.hpp"
#include "../conditions/boundary_conditions.hpp"
#include "../../thermophysics/mixture_interface.hpp"
#include "../../io/config_types.hpp"
#include <expected>
#include <concepts>

namespace blast::boundary_layer::coefficients {

// Derivative computation utilities
namespace derivatives {
    
    template<std::ranges::sized_range Range>
    [[nodiscard]] auto compute_eta_derivative(
        Range&& values,
        double d_eta
    ) -> std::vector<double>;
    
    template<typename Matrix>
    [[nodiscard]] auto compute_matrix_eta_derivative(
        const Matrix& values,
        double d_eta
    ) -> Matrix;
}

class CoefficientCalculator {
private:
    const thermophysics::MixtureInterface& mixture_;
    const io::SimulationConfig& sim_config_;
    const io::NumericalConfig& num_config_;
    double d_eta_;

    [[nodiscard]] auto calculate_transport_coefficients(
        const CoefficientInputs& inputs,
        const ThermodynamicCoefficients& thermo,
        const conditions::BoundaryConditions& bc
    ) const -> std::expected<TransportCoefficients, CoefficientError>;
    
    [[nodiscard]] auto calculate_thermodynamic_coefficients(
        const CoefficientInputs& inputs,
        const conditions::BoundaryConditions& bc
    ) const -> std::expected<ThermodynamicCoefficients, CoefficientError>;
    
    [[nodiscard]] auto calculate_diffusion_coefficients(
        const CoefficientInputs& inputs,
        const conditions::BoundaryConditions& bc,
        const XiDerivatives& xi_der
    ) const -> std::expected<DiffusionCoefficients, CoefficientError>;
    
    [[nodiscard]] auto calculate_chemical_coefficients(
        const CoefficientInputs& inputs,
        const ThermodynamicCoefficients& thermo,
        const conditions::BoundaryConditions& bc
    ) const -> std::expected<ChemicalCoefficients, CoefficientError>;
    
    [[nodiscard]] auto calculate_thermal_diffusion(
        const CoefficientInputs& inputs,
        const ThermodynamicCoefficients& thermo,
        const conditions::BoundaryConditions& bc
    ) const -> std::expected<ThermalDiffusionCoefficients, CoefficientError>;
    
    [[nodiscard]] auto calculate_wall_properties(
        const CoefficientInputs& inputs,
        const conditions::BoundaryConditions& bc,
        const TransportCoefficients& transport,
        const ThermodynamicCoefficients& thermo
    ) const -> std::expected<WallProperties, CoefficientError>;
    
    [[nodiscard]] auto calculate_species_enthalpies(
        const CoefficientInputs& inputs
    ) const -> std::expected<std::pair<core::Matrix<double>, core::Matrix<double>>, CoefficientError>;

public:
    CoefficientCalculator(
        const thermophysics::MixtureInterface& mixture,
        const io::SimulationConfig& sim_config,
        const io::NumericalConfig& num_config
    ) : mixture_(mixture), 
        sim_config_(sim_config), 
        num_config_(num_config),
        d_eta_(num_config.eta_max / static_cast<double>(num_config.n_eta - 1)) {}
    
    // Main calculation interface
    [[nodiscard]] auto calculate(
        const CoefficientInputs& inputs,
        const conditions::BoundaryConditions& bc,
        const XiDerivatives& xi_der
    ) const -> std::expected<CoefficientSet, CoefficientError>;
};

} // namespace blast::boundary_layer::coefficients