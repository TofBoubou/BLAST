#pragma once
#include "../../io/config_types.hpp"
#include "../../thermophysics/mixture_interface.hpp"
#include "../conditions/boundary_conditions.hpp"
#include "coefficient_types.hpp"
#include "xi_derivatives.hpp"
#include <expected>

namespace blast::boundary_layer::coefficients::diffusion {

// Calculate Stefan-Maxwell diffusion fluxes and derivatives
[[nodiscard]] auto compute_stefan_maxwell_fluxes(const CoefficientInputs& inputs, CoefficientSet& coeffs,
                                                 const conditions::BoundaryConditions& bc, const XiDerivatives& xi_der,
                                                 const io::SimulationConfig& sim_config,
                                                 const thermophysics::MixtureInterface& mixture,
                                                 double d_eta) -> std::expected<void, CoefficientError>;

} // namespace blast::boundary_layer::coefficients::diffusion