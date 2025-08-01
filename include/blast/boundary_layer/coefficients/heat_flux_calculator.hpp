#pragma once
#include "../../io/config_types.hpp"
#include "../../thermophysics/mixture_interface.hpp"
#include "../conditions/boundary_conditions.hpp"
#include "coefficient_types.hpp"
#include "xi_derivatives.hpp"
#include <concepts>
#include <expected>

namespace blast::boundary_layer::coefficients {

// Heat flux geometry factors
struct HeatFluxGeometryFactors {
  double der_fact;           // Derivative factor for conductive flux
  double dy_deta_factor;     // Coordinate transformation factor η → y
  bool valid_geometry;       // False for Cone/FlatPlate at stagnation
};

// Error type for heat flux calculations  
class HeatFluxError : public core::BlastException {
public:
  explicit HeatFluxError(std::string_view message, std::source_location location = std::source_location::current())
      : BlastException(std::format("Heat Flux Error: {}", message), location) {}
};

class HeatFluxCalculator {
private:
  const thermophysics::MixtureInterface& mixture_;
  const io::SimulationConfig& sim_config_;
  double d_eta_;

  [[nodiscard]] auto compute_heat_flux_geometry_factors(
      int station, double xi, const conditions::BoundaryConditions& bc,
      const WallProperties& wall_props) const 
      -> std::expected<HeatFluxGeometryFactors, HeatFluxError>;

  [[nodiscard]] auto compute_local_conductivities(
      const CoefficientInputs& inputs, const conditions::BoundaryConditions& bc) const
      -> std::expected<std::vector<double>, HeatFluxError>;

  [[nodiscard]] auto compute_reference_flux(
      const conditions::BoundaryConditions& bc, const CoefficientSet& coeffs) const -> double;

  [[nodiscard]] auto compute_conductive_flux_profile(
      const std::vector<double>& dT_deta, const std::vector<double>& k_local,
      const HeatFluxGeometryFactors& geo_factors) const -> std::vector<double>;

  [[nodiscard]] auto compute_diffusive_flux_profile(
      const CoefficientSet& coeffs) const 
      -> std::pair<std::vector<double>, core::Matrix<double>>;

  [[nodiscard]] auto compute_wall_heat_fluxes(
      const CoefficientInputs& inputs, const CoefficientSet& coeffs,
      const conditions::BoundaryConditions& bc, int station, double xi) const
      -> std::expected<std::tuple<double, double, double>, HeatFluxError>;

public:
  HeatFluxCalculator(const thermophysics::MixtureInterface& mixture, 
                     const io::SimulationConfig& sim_config,
                     const io::NumericalConfig& num_config)
      : mixture_(mixture), sim_config_(sim_config),
        d_eta_(num_config.eta_max / static_cast<double>(num_config.n_eta - 1)) {}

  [[nodiscard]] auto calculate(
      const CoefficientInputs& inputs, const CoefficientSet& coeffs,
      const conditions::BoundaryConditions& bc, const std::vector<double>& dT_deta,
      int station, double xi) const -> std::expected<HeatFluxCoefficients, HeatFluxError>;
};

} // namespace blast::boundary_layer::coefficients