#pragma once
#include "../../core/containers.hpp"
#include "../../core/exceptions.hpp"
#include <vector>
#include <span>

namespace blast::boundary_layer::coefficients {

// Transport coefficients (l0, l3 and derivatives)
struct TransportCoefficients {
    std::vector<double> l0;        
    std::vector<double> l3;        
    std::vector<double> dl0_deta;  // Derivative of l0
    std::vector<double> dl3_deta;  // Derivative of l3
};

// Thermodynamic properties
struct ThermodynamicCoefficients {
    std::vector<double> rho;        // Density field
    std::vector<double> MW;         // Molecular weight field
    std::vector<double> d_MW_deta;  // MW derivative
    std::vector<double> d_rho_deta; // Density derivative
    double h_wall;                  // Wall enthalpy
};

// Diffusion related coefficients
struct DiffusionCoefficients {
    core::Matrix<double> Dij_bin;   // Binary diffusion coefficients [n_eta x n_species x n_species]
    std::vector<double> y;          // Transformed coordinate for continuity
    core::Matrix<double> J;         // Diffusion fluxes [n_species x n_eta]
    core::Matrix<double> dJ_deta;   // Flux derivatives [n_species x n_eta]
};

// Chemical production rates
struct ChemicalCoefficients {
    core::Matrix<double> wi;                      // Production rates [n_eta x n_species]
    std::vector<core::Matrix<double>> d_wi_dc;   // Jacobian [n_eta] of [n_species x n_species]
};

// Thermal diffusion ratios
struct ThermalDiffusionCoefficients {
    core::Matrix<double> tdr;       // Thermal diffusion ratios [n_species x n_eta]
    core::Matrix<double> tdr_term;   // TDR term for diffusion fluxes [n_eta x n_species]
};

// Wall transport properties
struct WallProperties {
    double k_wall;      // Wall thermal conductivity
    double rho_wall;    // Wall density
    double mu_wall;     // Wall viscosity
    double Cp_wall;     // Wall heat capacity
    double Pr_wall;     // Wall Prandtl number
};

// Complete coefficient set
struct CoefficientSet {
    TransportCoefficients transport;
    ThermodynamicCoefficients thermodynamic;
    DiffusionCoefficients diffusion;
    ChemicalCoefficients chemical;
    ThermalDiffusionCoefficients thermal_diffusion;
    WallProperties wall;
    
    // Species enthalpy data
    core::Matrix<double> h_species;       // [n_species x n_eta]
    core::Matrix<double> dh_species_deta; // Derivatives
};

// Input data structure for coefficient calculation
struct CoefficientInputs {
    double xi;
    std::span<const double> F;
    const core::Matrix<double>& c;       // Species concentrations [n_species x n_eta]
    const core::Matrix<double>& dc_deta; // Concentration derivatives
    std::span<const double> T;           // Temperature field
};

// Error type for coefficient calculations
class CoefficientError : public core::BlastException {
public:
    explicit CoefficientError(std::string_view message,
                             std::source_location location = std::source_location::current())
        : BlastException(std::format("Coefficient Error: {}", message), location) {}
};

} // namespace blast::boundary_layer::coefficients