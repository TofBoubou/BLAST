#pragma once
#include "../../core/exceptions.hpp"
#include "../../io/config_types.hpp"
#include <vector>
#include <optional>
#include <expected>

namespace blast::boundary_layer::conditions {

// Immutable value types for boundary conditions
struct EdgeConditions {
    double pressure;
    double viscosity;
    double velocity;
    double enthalpy;
    double density;
    std::vector<double> species_fractions;
    
    // Derivatives
    double d_xi_dx;
    double d_ue_dx;
    double d_he_dx;
    double d_he_dxi;
    
    // Derived quantities
    double body_radius;
};

struct WallConditions {
    double temperature;
    // std::optional<std::vector<double>> catalycity;  // Only if catalytic
};

struct BoundaryConditions {
    EdgeConditions edge;
    WallConditions wall;
    double beta;
    int station;
    double xi;
    
    // Convenience getters with legacy names for compatibility
    [[nodiscard]] constexpr auto P_e() const noexcept { return edge.pressure; }
    [[nodiscard]] constexpr auto mu_e() const noexcept { return edge.viscosity; }
    [[nodiscard]] constexpr auto ue() const noexcept { return edge.velocity; }
    [[nodiscard]] constexpr auto he() const noexcept { return edge.enthalpy; }
    [[nodiscard]] constexpr auto rho_e() const noexcept { return edge.density; }
    [[nodiscard]] constexpr auto Tw() const noexcept { return wall.temperature; }
    [[nodiscard]] constexpr auto r_body() const noexcept { return edge.body_radius; }
    [[nodiscard]] constexpr auto d_xi_dx() const noexcept { return edge.d_xi_dx; }
    [[nodiscard]] constexpr auto d_ue_dx() const noexcept { return edge.d_ue_dx; }
    [[nodiscard]] constexpr auto d_he_dx() const noexcept { return edge.d_he_dx; }
    [[nodiscard]] constexpr auto d_he_dxi() const noexcept { return edge.d_he_dxi; }
    
    // Access species safely
    [[nodiscard]] auto c_e() const noexcept -> const std::vector<double>& { 
        return edge.species_fractions; 
    }
    
};

class BoundaryConditionError : public core::BlastException {
public:
    explicit BoundaryConditionError(std::string_view message,
                                   std::source_location location = std::source_location::current())
        : BlastException(std::format("Boundary Condition Error: {}", message), location) {}
};

} // namespace blast::boundary_layer::conditions