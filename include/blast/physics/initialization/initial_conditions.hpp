#pragma once
#include "../fields/physical_fields.hpp"
#include "../fields/boundary_conditions.hpp"
#include "../../io/config_types.hpp"
#include "../../core/exceptions.hpp"
#include <expected>

namespace blast::physics::initialization {

class FieldInitializer {
private:
    const io::InitialGuessConfig& initial_config_;
    const BoundaryConditionInterpolator& bc_interpolator_;
    
    [[nodiscard]] auto create_default_velocity_profile(
        std::size_t n_eta, double eta_max
    ) const noexcept -> fields::VelocityField;
    
    [[nodiscard]] auto create_default_temperature_profile(
        std::size_t n_eta, double wall_temperature
    ) const noexcept -> fields::TemperatureField;
    
    [[nodiscard]] auto create_default_enthalpy_profile(
        std::size_t n_eta
    ) const noexcept -> fields::EnthalpyField;
    
    [[nodiscard]] auto create_default_species_profile(
        std::size_t n_eta, std::size_t n_species, 
        const core::Vector<double>& wall_composition
    ) const noexcept -> fields::SpeciesField;

public:
    explicit FieldInitializer(
        const io::InitialGuessConfig& initial_config,
        const BoundaryConditionInterpolator& bc_interpolator
    ) noexcept;
    
    [[nodiscard]] auto initialize_solution_fields(
        std::size_t n_eta, std::size_t n_species, 
        double eta_max, double x_station = 0.0
    ) const -> std::expected<fields::SolutionFields, core::BlastException>;
    
    [[nodiscard]] auto validate_initial_guess(
        const fields::SolutionFields& fields
    ) const -> std::expected<void, core::ValidationError>;
};

} // namespace blast::physics::initialization