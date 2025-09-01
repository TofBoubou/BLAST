#include "gasp2/catalysis/gamma/bose_flux.hpp"

#include <algorithm>

namespace gasp2::catalysis::gamma::bose {

[[nodiscard]] CatalysisFluxes
compute_fluxes(const SpeciesData &species, std::size_t ns, std::size_t idx_O,
               std::size_t idx_CO, std::size_t idx_O2, std::size_t idx_CO2,
               double gamma, double p2) {
  CatalysisFluxes fluxes;
  fluxes.cat_fluxes.assign(ns, 0.0);

  // Impinging number fluxes precomputed in SpeciesData.
  double imp_O = species.imp_flux[idx_O];
  double imp_CO = species.imp_flux[idx_CO];

  // Maximum usable oxygen recombination probability.
  double gamma_O_max = std::min(imp_CO / ((1.0 - p2) * imp_O), 1.0);

  // Number fluxes according to Bose's model.
  double Gamma_O = gamma * gamma_O_max * imp_O;
  double Gamma_CO = (1.0 - p2) * Gamma_O;
  double Gamma_CO2 = Gamma_CO;
  double Gamma_O2 = 0.5 * p2 * Gamma_O;

  // Convert number fluxes to mass fluxes.
  fluxes.cat_fluxes[idx_O] = Gamma_O * species.m[idx_O];
  fluxes.cat_fluxes[idx_CO] = Gamma_CO * species.m[idx_CO];
  fluxes.cat_fluxes[idx_CO2] = -Gamma_CO2 * species.m[idx_CO2];
  fluxes.cat_fluxes[idx_O2] = -Gamma_O2 * species.m[idx_O2];

  return fluxes;
}

} // namespace gasp2::catalysis::gamma::bose
