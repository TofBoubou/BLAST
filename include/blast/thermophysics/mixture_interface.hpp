#pragma once
#include "../core/containers.hpp"
#include <expected>
#include <string_view>
#include <vector>
#include <span>

namespace blast::thermophysics {

namespace constants {
    inline constexpr double R_universal = 8314.46261815324; // J/(kmol·K)
    inline constexpr double k_Boltzmann = 1.380649e-23;     // J/K
    inline constexpr double N_Avogadro = 6.02214076e23;     // 1/mol
    inline constexpr double Pi = 3.14159265358979323846;
}

// Abstract interface for mixture properties
class MixtureInterface {
public:
    virtual ~MixtureInterface() = default;
    
    // Composition and molecular weight
    [[nodiscard]] virtual auto mixture_molecular_weight(
        std::span<const double> mass_fractions
    ) const -> std::expected<double, std::string_view> = 0;
    
    [[nodiscard]] virtual auto species_molecular_weight(
        std::size_t species_index
    ) const noexcept -> double = 0;
    
    [[nodiscard]] virtual auto species_charges() const noexcept 
        -> std::span<const double> = 0;
    
    // Thermodynamic properties
    [[nodiscard]] virtual auto mixture_enthalpy(
        std::span<const double> mass_fractions,
        double temperature,
        double pressure
    ) const -> std::expected<double, std::string_view> = 0;
    
    [[nodiscard]] virtual auto species_enthalpies(
        double temperature
    ) const -> std::expected<std::vector<double>, std::string_view> = 0;
    
    [[nodiscard]] virtual auto frozen_cp(
        std::span<const double> mass_fractions,
        double temperature
    ) const -> std::expected<double, std::string_view> = 0;
    
    // Transport properties
    [[nodiscard]] virtual auto viscosity(
        std::span<const double> mass_fractions,
        double temperature,
        double pressure
    ) const -> std::expected<double, std::string_view> = 0;
    
    [[nodiscard]] virtual auto frozen_thermal_conductivity(
        std::span<const double> mass_fractions,
        double temperature,
        double pressure
    ) const -> std::expected<double, std::string_view> = 0;
    
    [[nodiscard]] virtual auto binary_diffusion_coefficients(
        double temperature,
        double pressure
    ) const -> std::expected<core::Matrix<double>, std::string_view> = 0;
    
    // Chemical properties
    [[nodiscard]] virtual auto production_rates(
        std::span<const double> partial_densities,
        double temperature,
        double density
    ) const -> std::expected<std::vector<double>, std::string_view> = 0;
    
    [[nodiscard]] virtual auto production_rate_jacobian(
        std::span<const double> partial_densities,
        double temperature,
        double density
    ) const -> std::expected<core::Matrix<double>, std::string_view> = 0;
    
    // Thermal diffusion
    [[nodiscard]] virtual auto thermal_diffusion_ratios(
        std::span<const double> mass_fractions,
        double temperature,
        double pressure
    ) const -> std::expected<std::vector<double>, std::string_view> = 0;
    
    // Equilibrium composition
    [[nodiscard]] virtual auto equilibrium_composition(
        double temperature,
        double pressure
    ) const -> std::expected<std::vector<double>, std::string_view> = 0;
    
    // Number of species
    [[nodiscard]] virtual auto n_species() const noexcept -> std::size_t = 0;
    [[nodiscard]] virtual auto has_electrons() const noexcept -> bool = 0;
};

} // namespace blast::thermophysics