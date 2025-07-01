#pragma once
#include "mixture_interface.hpp"
#include <mutation++.h>
#include <memory>
#include <vector>

namespace blast::thermophysics {

// Concrete implementation using Mutation++
class MutationMixture final : public MixtureInterface { // final because can't be derived
private:
    std::unique_ptr<Mutation::Mixture> mixture_;
    const std::size_t n_species_;
    const bool has_electrons_;
    
    // Cached data
    mutable std::vector<double> species_charges_; // mutable = can be modified even in a const function
    mutable std::vector<double> species_mw_;
    mutable std::vector<std::string> species_names_;
    
    // Helper to validate composition array
    [[nodiscard]] auto validate_composition(std::span<const double> fractions) const 
        -> std::expected<void, ThermophysicsError>;

public:
    explicit MutationMixture(const io::MixtureConfig& config);
    ~MutationMixture() override = default;
    
    // Delete copy operations
    MutationMixture(const MutationMixture&) = delete;
    MutationMixture& operator=(const MutationMixture&) = delete;
    
    // Default move operations
    MutationMixture(MutationMixture&&) = default;
    MutationMixture& operator=(MutationMixture&&) = default;
    
    // Species information
    [[nodiscard]] auto n_species() const noexcept -> std::size_t override { 
        return n_species_; 
    }
    
    [[nodiscard]] auto has_electrons() const noexcept -> bool override { 
        return has_electrons_; 
    }
    
    // Implement all virtual methods
    [[nodiscard]] auto mixture_molecular_weight(
        std::span<const double> mass_fractions
    ) const -> std::expected<double, ThermophysicsError> override;
    
    [[nodiscard]] auto species_molecular_weight(
        std::size_t species_index
    ) const noexcept -> double override;
    
    [[nodiscard]] auto species_charges() const noexcept 
        -> std::span<const double> override;
    
    [[nodiscard]] auto set_state(
        std::span<const double> mass_fractions,
        double temperature,
        double pressure
    ) const -> std::expected<void, ThermophysicsError> override;
    
    [[nodiscard]] auto mixture_enthalpy(
        std::span<const double> mass_fractions,
        double temperature,
        double pressure
    ) const -> std::expected<double, ThermophysicsError> override;
    
    [[nodiscard]] auto species_enthalpies(
        double temperature
    ) const -> std::expected<std::vector<double>, ThermophysicsError> override;
    
    [[nodiscard]] auto frozen_cp(
        std::span<const double> mass_fractions,
        double temperature,
        double pressure
    ) const -> std::expected<double, ThermophysicsError> override;
    
    [[nodiscard]] auto viscosity(
        std::span<const double> mass_fractions,
        double temperature,
        double pressure
    ) const -> std::expected<double, ThermophysicsError> override;
    
    [[nodiscard]] auto frozen_thermal_conductivity(
        std::span<const double> mass_fractions,
        double temperature,
        double pressure
    ) const -> std::expected<double, ThermophysicsError> override;
    
    [[nodiscard]] auto binary_diffusion_coefficients(
        double temperature,
        double pressure
    ) const -> std::expected<core::Matrix<double>, ThermophysicsError> override;
    
    [[nodiscard]] auto production_rates(
        std::span<const double> partial_densities,
        double temperature
    ) const -> std::expected<std::vector<double>, ThermophysicsError> override;
    
    [[nodiscard]] auto production_rate_jacobian(
        std::span<const double> partial_densities,
        double temperature
    ) const -> std::expected<core::Matrix<double>, ThermophysicsError> override;
    
    [[nodiscard]] auto thermal_diffusion_ratios(
        std::span<const double> mass_fractions,
        double temperature,
        double pressure
    ) const -> std::expected<std::vector<double>, ThermophysicsError> override;
    
    [[nodiscard]] auto equilibrium_composition(
        double temperature,
        double pressure
    ) const -> std::expected<std::vector<double>, ThermophysicsError> override;
    
    [[nodiscard]] auto surface_reaction_rates(
        std::span<const double> partial_densities,
        double wall_temperature
    ) const -> std::expected<std::vector<double>, ThermophysicsError> override;
};

} // namespace blast::thermophysics