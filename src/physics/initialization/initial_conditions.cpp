#include "blast/physics/initialization/initial_conditions.hpp"
#include <algorithm>
#include <cmath>

namespace blast::physics::initialization {

FieldInitializer::FieldInitializer(
    const io::InitialGuessConfig& initial_config,
    const BoundaryConditionInterpolator& bc_interpolator
) noexcept : initial_config_(initial_config), bc_interpolator_(bc_interpolator) {}

auto FieldInitializer::create_default_velocity_profile(
    std::size_t n_eta, double eta_max
) const noexcept -> fields::VelocityField {
    
    fields::VelocityField velocity(n_eta);
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        const auto eta = static_cast<double>(i) * eta_max / static_cast<double>(n_eta - 1);
        const auto eta_norm = eta / eta_max;
        
        // Standard boundary layer profile: F = η/η_max * (2 - η/η_max)
        velocity[i] = eta_norm * (2.0 - eta_norm);
    }
    
    return velocity;
}

auto FieldInitializer::create_default_temperature_profile(
    std::size_t n_eta, double wall_temperature
) const noexcept -> fields::TemperatureField {
    
    fields::TemperatureField temperature(n_eta);
    std::ranges::fill(temperature, wall_temperature);
    return temperature;
}

auto FieldInitializer::create_default_enthalpy_profile(
    std::size_t n_eta
) const noexcept -> fields::EnthalpyField {
    
    fields::EnthalpyField enthalpy(n_eta);
    std::ranges::fill(enthalpy, 1.0); // Non-dimensional, g = h/h_e
    return enthalpy;
}

auto FieldInitializer::create_default_species_profile(
    std::size_t n_eta, std::size_t n_species,
    const core::Vector<double>& wall_composition
) const noexcept -> fields::SpeciesField {
    
    fields::SpeciesField species(n_species, n_eta);
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        for (std::size_t j = 0; j < n_species; ++j) {
            species(j, i) = (j < wall_composition.size()) ? wall_composition[j] : 0.0;
        }
    }
    
    return species;
}

auto FieldInitializer::initialize_solution_fields(
    std::size_t n_eta, std::size_t n_species,
    double eta_max, double x_station
) const -> std::expected<fields::SolutionFields, core::BlastException> {
    
    auto wall_props_result = bc_interpolator_.compute_wall_properties_at_x(x_station);
    if (!wall_props_result) {
        return std::unexpected(wall_props_result.error());
    }
    
    auto edge_props_result = bc_interpolator_.compute_edge_properties_at_x(x_station);
    if (!edge_props_result) {
        return std::unexpected(edge_props_result.error());
    }
    
    const auto& wall_props = wall_props_result.value();
    const auto& edge_props = edge_props_result.value();
    
    fields::SolutionFields fields(n_eta, n_species);
    
    if (initial_config_.use_initial_guess) {
        // Load from initial guess configuration
        if (initial_config_.temperature_profile && 
            initial_config_.temperature_profile->size() == n_eta) {
            
            std::ranges::copy(*initial_config_.temperature_profile, fields.temperature.begin());
        } else {
            fields.temperature = create_default_temperature_profile(n_eta, wall_props.temperature);
        }
        
        if (initial_config_.enthalpy_profile && 
            initial_config_.enthalpy_profile->size() == n_eta) {
            
            std::ranges::copy(*initial_config_.enthalpy_profile, fields.enthalpy.begin());
        } else {
            fields.enthalpy = create_default_enthalpy_profile(n_eta);
        }
        
        if (initial_config_.species_profiles && 
            initial_config_.species_profiles->size() == n_species) {
            
            for (std::size_t i = 0; i < n_species; ++i) {
                if ((*initial_config_.species_profiles)[i].size() == n_eta) {
                    for (std::size_t j = 0; j < n_eta; ++j) {
                        fields.species_fractions(i, j) = (*initial_config_.species_profiles)[i][j];
                    }
                }
            }
        } else {
            fields.species_fractions = create_default_species_profile(
                n_eta, n_species, edge_props.species_fractions
            );
        }
    } else {
        // Default initialization
        fields.temperature = create_default_temperature_profile(n_eta, wall_props.temperature);
        fields.enthalpy = create_default_enthalpy_profile(n_eta);
        fields.species_fractions = create_default_species_profile(
            n_eta, n_species, edge_props.species_fractions
        );
    }
    
    // Always use default velocity profile
    fields.velocity = create_default_velocity_profile(n_eta, eta_max);
    
    // Density and viscosity will be computed by thermodynamic interface later
    std::ranges::fill(fields.density, edge_props.density);
    std::ranges::fill(fields.viscosity, edge_props.viscosity);
    
    return fields;
}

auto FieldInitializer::validate_initial_guess(
    const fields::SolutionFields& fields
) const -> std::expected<void, core::ValidationError> {
    
    // Basic validation: non-negative values, species sum
    const auto n_eta = fields.n_eta();
    const auto n_species = fields.n_species();
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        if (fields.temperature[i] <= 0.0) {
            return std::unexpected(core::ValidationError(
                "temperature", "Temperature must be positive"
            ));
        }
        
        if (fields.density[i] <= 0.0) {
            return std::unexpected(core::ValidationError(
                "density", "Density must be positive"
            ));
        }
        
        // Check species mass fraction sum
        double species_sum = 0.0;
        for (std::size_t j = 0; j < n_species; ++j) {
            if (fields.species_fractions(j, i) < 0.0) {
                return std::unexpected(core::ValidationError(
                    "species_fractions", "Species fractions must be non-negative"
                ));
            }
            species_sum += fields.species_fractions(j, i);
        }
        
        constexpr auto tolerance = 1e-6;
        if (std::abs(species_sum - 1.0) > tolerance) {
            return std::unexpected(core::ValidationError(
                "species_fractions", "Species fractions must sum to 1.0"
            ));
        }
    }
    
    return {};
}

} // namespace blast::physics::initialization