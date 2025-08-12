#pragma once
#include "../../core/exceptions.hpp"
#include "../../io/config_types.hpp"
#include "../../thermophysics/mixture_interface.hpp"
#include <expected>
#include <optional>
#include <vector>

namespace blast::boundary_layer::conditions {

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
  bool catalytic = false;
};

struct BoundaryConditions {
  EdgeConditions edge;
  WallConditions wall;
  double beta;
  int station;
  double xi;

  [[nodiscard]] constexpr auto P_e() const noexcept { return edge.pressure; }
  [[nodiscard]] constexpr auto mu_e() const noexcept { return edge.viscosity; }
  [[nodiscard]] constexpr auto ue() const noexcept { return edge.velocity; }
  [[nodiscard]] constexpr auto he() const noexcept { return edge.enthalpy; }
  [[nodiscard]] constexpr auto rho_e() const noexcept { return edge.density; }
  [[nodiscard]] constexpr auto Tw() const noexcept { return wall.temperature; }
  [[nodiscard]] constexpr auto catalytic() const noexcept { return wall.catalytic; }
  [[nodiscard]] constexpr auto r_body() const noexcept { return edge.body_radius; }
  [[nodiscard]] constexpr auto d_xi_dx() const noexcept { return edge.d_xi_dx; }
  [[nodiscard]] constexpr auto d_ue_dx() const noexcept { return edge.d_ue_dx; }
  [[nodiscard]] constexpr auto d_he_dx() const noexcept { return edge.d_he_dx; }
  [[nodiscard]] constexpr auto d_he_dxi() const noexcept { return edge.d_he_dxi; }

  [[nodiscard]] auto c_e() const noexcept -> const std::vector<double>& { return edge.species_fractions; }

  // Dynamic update methods for thermodynamic consistency
  void update_edge_density(double new_density) noexcept { edge.density = new_density; }

  void update_edge_viscosity(double new_viscosity) noexcept { edge.viscosity = new_viscosity; }
};

class BoundaryConditionError : public core::BlastException {
public:
  explicit BoundaryConditionError(std::string_view message,
                                  std::source_location location = std::source_location::current())
      : BlastException(std::format("Boundary Condition Error: {}", message), location) {}
};

} // namespace blast::boundary_layer::conditions