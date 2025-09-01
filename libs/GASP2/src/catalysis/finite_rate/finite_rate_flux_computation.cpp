#include "gasp2/catalysis/finite_rate/finite_rate_flux_computation.hpp"
#include "gasp2/catalysis/finite_rate/stoichiometry.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace gasp2::catalysis::finite_rate {
namespace details {

//---------------------------------------------------------------------------
// Flux computation
//---------------------------------------------------------------------------
CatalysisFluxes compute_reaction_fluxes(const ReactionProperties &props,
                                        std::span<const double> coverage) {
  CatalysisFluxes fluxes;

  // ========== STOICHIOMETRY ==========
  auto stoich = build_stoichiometry(props);
  const auto &nu_r = stoich.nu_reactants;
  const auto &nu_p = stoich.nu_products;
  const auto &species_order = stoich.species_order;
  const auto &surface_index = stoich.surface_index;

  const std::size_t n_species = species_order.size();
  const std::size_t n_reactions = props.reactions.size();
  const std::size_t n_ads = coverage.size() - 1; // last entry: empty site

  // Prepare flux containers
  fluxes.cat_fluxes.assign(ctx.species_order.size(), 0.0);
  fluxes.destruction_fractions.assign(ctx.species_order.size(),
                                      std::vector<double>(n_reactions, 0.0));

  // ========== DEBUG: REACTION CONSTANTS ==========
  if (ctx.debug) {
    for (std::size_t r = 0; r < n_reactions; ++r) {
      const auto &rx = props.reactions[r];
      std::cout << "Reaction " << r << ": kf=" << rx.kf << ", kb=" << rx.kb
                << ", Kc=" << rx.Kc << ", Ka=" << rx.Ka << '\n';
    }
  }

  // ========== CONCENTRATIONS ==========
  std::vector<double> conc(n_species, 0.0);
  for (std::size_t i = 0; i < ctx.species_order.size(); ++i) {
    const auto &name = ctx.species_order[i];
    conc[i] = props.concentration.at(name);
  }
  for (std::size_t j = 0; j < n_ads; ++j)
    conc[surface_index[j]] = coverage[j];
  conc[surface_index.back()] = coverage.back();

  // ========== REACTION RATES ==========
  std::vector<double> rate_f(n_reactions, 0.0);
  std::vector<double> rate_b(n_reactions, 0.0);
  for (std::size_t r = 0; r < n_reactions; ++r) {
    const auto &rx = props.reactions[r];
    rate_f[r] = rx.kf;
    rate_b[r] = rx.kb;
    for (std::size_t m = 0; m < n_species; ++m) {
      int nu_r_m = nu_r[m][r];
      if (nu_r_m > 0)
        rate_f[r] *= std::pow(conc[m], nu_r_m);
      int nu_p_m = nu_p[m][r];
      if (nu_p_m > 0)
        rate_b[r] *= std::pow(conc[m], nu_p_m);
    }
    if (ctx.debug)
      std::cout << "  rate_f[" << r << "]=" << rate_f[r] << ", rate_b[" << r
                << "]=" << rate_b[r] << '\n';
  }

  // ========== SPECIES FLUXES ==========
  for (std::size_t s = 0; s < ctx.species_order.size(); ++s) {
    double ws = 0.0;
    for (std::size_t r = 0; r < n_reactions; ++r) {
      int nu_p_s = nu_p[s][r];
      int nu_r_s = nu_r[s][r];
      double w_ir = static_cast<double>(nu_p_s - nu_r_s) *
                    (rate_f[r] - rate_b[r]) * ctx.molar_masses[s];
      ws += w_ir;
      fluxes.destruction_fractions[s][r] = w_ir;
    }
    fluxes.cat_fluxes[s] = ws;
    if (ctx.debug)
      std::cout << "Flux for species " << ctx.species_order[s] << " = "
                << fluxes.cat_fluxes[s] << '\n';
  }

  // Flip sign so that destruction>0
  for (auto &flux : fluxes.cat_fluxes)
    flux = -flux;
  for (auto &vec : fluxes.destruction_fractions)
    for (auto &val : vec)
      val = -val;

  // ========== FRACTIONAL CONTRIBUTIONS ==========
  for (std::size_t s = 0; s < fluxes.destruction_fractions.size(); ++s) {
    double wi = fluxes.cat_fluxes[s];
    if (wi != 0.0) {
      for (auto &val : fluxes.destruction_fractions[s])
        val /= wi;
    } else {
      std::fill(fluxes.destruction_fractions[s].begin(),
                fluxes.destruction_fractions[s].end(), 0.0);
    }
  }

  return fluxes;
}

} // namespace details
} // namespace gasp2::catalysis::finite_rate
