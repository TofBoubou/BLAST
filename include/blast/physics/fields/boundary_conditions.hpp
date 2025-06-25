#pragma once
#include "../../io/config_types.hpp"
#include "../../core/exceptions.hpp"
#include <expected>
#include <concepts>

namespace blast::physics::fields {

// Edge properties at a given xi station
struct EdgeProperties {
    double pressure;
    double temperature; 
    double velocity;
    double enthalpy;
    double density;
    double viscosity;
    core::Vector<double> species_fractions;
    
    explicit EdgeProperties(std::size_t n_species) 
        : species_fractions(n_species) {}
};

// Wall properties at a given xi station  
struct WallProperties {
    double temperature;
    // Future: catalysis parameters, heat flux conditions
    
    explicit WallProperties(double temp) : temperature(temp) {}
};

class BoundaryConditionInterpolator {
private:
    const io::OuterEdgeConfig& edge_config_;
    const io::WallParametersConfig& wall_config_;
    core::Vector<double> x_coordinates_;
    
    [[nodiscard]] auto find_interval(double x_target) const noexcept
        -> std::expected<std::pair<std::size_t, std::size_t>, core::BlastException>;
        
    [[nodiscard]] constexpr auto linear_interpolate(
        double x, double x1, double x2, double y1, double y2
    ) const noexcept -> double;

public:
    explicit BoundaryConditionInterpolator(
        const io::OuterEdgeConfig& edge_config,
        const io::WallParametersConfig& wall_config
    );
    
    [[nodiscard]] auto compute_edge_properties_at_x(double x_station) const
        -> std::expected<EdgeProperties, core::BlastException>;
        
    [[nodiscard]] auto compute_wall_properties_at_x(double x_station) const  
        -> std::expected<WallProperties, core::BlastException>;
        
    [[nodiscard]] constexpr auto x_min() const noexcept -> double;
    [[nodiscard]] constexpr auto x_max() const noexcept -> double;
};

} // namespace blast::physics::fields