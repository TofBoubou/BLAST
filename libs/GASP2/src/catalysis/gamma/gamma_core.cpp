#include "gasp2/catalysis/gamma/gamma_core.hpp"

#include <array>
#include <iostream>
#include <stdexcept>

namespace gasp2::catalysis::gamma {
using gasp2::Error;
using gasp2::ErrorCode;

// Compute the species mass flux at the wall. The recombination probability must
// be in the physical range \f$0 \leq \gamma \leq 1\f$. For second-order
// formulations the effective probability becomes \f$2\gamma/(2-\gamma)\f$.
[[nodiscard]] Result<double> compute_flux(double gamma, double mass_imp_flux,
                                          bool first_order) {
  if (gamma < 0.0 || gamma > 1.0) {
    return std::unexpected(
        Error{ErrorCode::OutOfRange,
              "Recombination probability must be between 0 and 1"});
  }
  double gamma_eff = first_order ? gamma : (2.0 * gamma / (2.0 - gamma));
  return gamma_eff * mass_imp_flux;
}

[[nodiscard]] Result<CatalysisFluxes>
compute_all_fluxes(const SpeciesData &species, std::size_t ns, std::size_t nr,
                   const std::vector<std::vector<double>> &gammas,
                   const std::vector<std::vector<int>> &nu,
                   const std::vector<std::vector<int>> &nu_p,
                   const std::vector<std::vector<std::vector<int>>> &mu,
                   bool first_order, bool limiting_fluxes,
                   const std::vector<CatalysisModel> &types,
                   const std::vector<bool> &heterogeneous,
                   const std::vector<std::string> &species_order, bool debug) {
  CatalysisFluxes fluxes;
  fluxes.cat_fluxes.assign(ns, 0.0);

  for (std::size_t r = 0; r < nr; ++r) {
    // ========== LIMITING FLUX ADJUSTMENT ==========
    std::array<int, 2> react_idx{-1, -1};
    std::array<int, 2> stoich{0, 0};
    std::array<double, 2> imp_mod{0.0, 0.0};
    int nreact = 0;
    bool adjust = limiting_fluxes && heterogeneous[r] &&
                  (types[r] == CatalysisModel::GammaGiven ||
                   types[r] == CatalysisModel::GammaT);
    if (adjust) {
      for (std::size_t j = 0; j < ns; ++j) {
        if (r < nu.size() && j < nu[r].size() && nu[r][j] != 0 && nreact < 2) {
          react_idx[nreact] = static_cast<int>(j);
          stoich[nreact] =
              (r < nu_p.size() && j < nu_p[r].size()) ? nu_p[r][j] : 0;
          imp_mod[nreact] = species.imp_flux[j];
          ++nreact;
        }
      }
      if (nreact == 2) {
        double a = static_cast<double>(stoich[0]);
        double b = static_cast<double>(stoich[1]);
        std::size_t ia = static_cast<std::size_t>(react_idx[0]);
        std::size_t ib = static_cast<std::size_t>(react_idx[1]);
        double gA = gammas[r][ia];
        double gB = gammas[r][ib];
        double impA = imp_mod[0];
        double impB = imp_mod[1];
        double left = b * gA * impA;
        double right = a * gB * impB;
        if (left < right) {
          imp_mod[1] = (b / a) * (gA / gB) * impA;
        } else if (left > right) {
          imp_mod[0] = (a / b) * (gB / gA) * impB;
        }
      }
    }

    // ========== FLUX ACCUMULATION ==========
    for (std::size_t j = 0; j < ns; ++j) {
      if (r >= nu.size() || j >= nu[r].size() || nu[r][j] == 0) {
        continue; // species j not consumed in reaction r
      }
      double gamma = gammas[r][j];
      if (gamma == 0.0)
        continue;

      double imp = species.imp_flux[j];
      if (adjust) {
        if (react_idx[0] == static_cast<int>(j)) {
          imp = imp_mod[0];
        } else if (react_idx[1] == static_cast<int>(j)) {
          imp = imp_mod[1];
        }
      }

      double mass_imp_flux = species.m[j] * imp;
      auto mflux_res = compute_flux(gamma, mass_imp_flux, first_order);
      if (!mflux_res) {
        auto err = mflux_res.error();
        err.message = "Reaction " + std::to_string(r) + ", species " +
                      species_order[j] + ": " + err.message;
        return std::unexpected(err);
      }
      double mflux = *mflux_res;
      fluxes.cat_fluxes[j] += static_cast<double>(nu[r][j]) * mflux;

      if (debug) {
        std::cout << "r=" << r << ", sp=" << species_order[j]
                  << ": nu=" << nu[r][j] << ", gamma=" << gamma
                  << ", mass_imp_flux=" << mass_imp_flux << ", mflux=" << mflux
                  << '\n';
      }

      for (std::size_t i = 0; i < ns; ++i) {
        if (r < mu.size() && i < mu[r].size() && j < mu[r][i].size()) {
          int coeff = mu[r][i][j];
          if (coeff != 0) {
            fluxes.cat_fluxes[i] -= static_cast<double>(coeff) * mflux;
            if (debug) {
              std::cout << "  produce " << coeff << " of " << species_order[i]
                        << " from " << species_order[j] << " -> "
                        << -static_cast<double>(coeff) * mflux << '\n';
            }
          }
        }
      }
    }
  }

  if (debug) {
    std::cout << "Accumulated fluxes:" << '\n';
    for (std::size_t i = 0; i < ns; ++i) {
      std::cout << "  " << species_order[i] << ": " << fluxes.cat_fluxes[i]
                << '\n';
    }
  }

  return fluxes;
}

} // namespace gasp2::catalysis::gamma
