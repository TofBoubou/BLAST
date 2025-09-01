#include "gasp2/catalysis/gamma/gamma_dispatcher.hpp"
#include "gasp2/catalysis/gamma/bose_flux.hpp"
#include "gasp2/catalysis/gamma/gamma_core.hpp"
#include "gasp2/catalysis/gamma/gamma_build.hpp"
#include "gasp2/catalysis/gamma/rini_flux.hpp"

#include <cmath>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace gasp2::catalysis::gamma {
using gasp2::Error;
using gasp2::ErrorCode;

namespace {
// ========== DEBUG HELPERS ==========
// Print reaction definitions including reactants and products.
void print_reactions(const std::vector<ReactionInput> &reactions) {
  for (std::size_t r = 0; r < reactions.size(); ++r) {
    const auto &reaction = reactions[r];
    std::cout << "Reaction " << r << ": " << reaction.formula << '\n';
    std::cout << "  Type: ";
    switch (reaction.type) {
    case CatalysisModel::GammaGiven:
      std::cout << "GammaGiven";
      break;
    case CatalysisModel::GammaT:
      std::cout << "GammaT";
      break;
    case CatalysisModel::GammaConsistent:
      std::cout << "GammaConsistent";
      break;
    case CatalysisModel::GammaBose:
      std::cout << "GammaBose";
      break;
    case CatalysisModel::SuperCatalytic:
      std::cout << "SuperCatalytic";
      break;
    case CatalysisModel::RiniModel:
      std::cout << "Rini_model";
      break;
    }
    std::cout << '\n';

    std::cout << "  Reactants:";
    for (const auto &[sp, coeff] : reaction.reactants) {
      std::cout << ' ' << coeff << sp;
    }
    std::cout << '\n';

    std::cout << "  Products:";
    for (const auto &[sp, coeff] : reaction.products) {
      std::cout << ' ' << coeff << sp;
    }
    std::cout << '\n';
  }
}

// Print gamma values for each reaction and species.
void print_gammas(const std::vector<std::vector<double>> &gammas,
                  const std::vector<ReactionInput> &reactions,
                  const std::unordered_map<std::string, std::size_t> &index) {
  std::cout << "Computed gammas:\n";
  for (std::size_t r = 0; r < reactions.size(); ++r) {
    const auto &reaction = reactions[r];
    std::cout << "Reaction " << r << ": " << reaction.formula << '\n';
    for (const auto &[sp, _] : reaction.reactants) {
      std::size_t idx_r = index.at(sp);
      std::cout << "  " << sp << ": " << gammas[r][idx_r] << '\n';
    }
  }
}

// Print nu and mu matrices.
void print_nu_mu(const std::vector<std::vector<int>> &nu,
                 const std::vector<std::vector<std::vector<int>>> &mu,
                 const std::vector<std::string> &species_order) {
  std::size_t ns = species_order.size();
  std::cout << "nu matrix:\n";
  for (std::size_t r = 0; r < nu.size(); ++r) {
    std::cout << "Reaction " << r << ':';
    for (std::size_t s = 0; s < ns; ++s) {
      std::cout << ' ' << species_order[s] << '=' << nu[r][s];
    }
    std::cout << '\n';
  }

  std::cout << "mu matrices:\n";
  for (std::size_t r = 0; r < mu.size(); ++r) {
    std::cout << "Reaction " << r << "\n";
    for (std::size_t i = 0; i < ns; ++i) {
      std::cout << "  to " << species_order[i] << ':';
      for (std::size_t j = 0; j < ns; ++j) {
        std::cout << ' ' << species_order[j] << '=' << mu[r][i][j];
      }
      std::cout << '\n';
    }
  }
}
} // namespace

/// @brief Dispatches the gamma calculation for the given reactions.
/// @todo Make gammaT available for gamma consistent
[[nodiscard]] Result<CatalysisFluxes>
gamma_dispatcher(const SpeciesData &species, std::size_t ns, std::size_t nr,
                 const std::vector<ReactionInput> &reactions, bool first_order,
                 bool limiting_fluxes,
                 const std::unordered_map<std::string, std::size_t> &index,
                 const StoichiometricMatrices stoic_matrices, bool all_super,
                 bool all_rini, const std::vector<std::string> &species_order,
                 bool debug) {
  try {
    // Normalize numeric formatting for debug output
    auto old_flags = std::cout.flags();
    auto old_prec = std::cout.precision();
    std::cout.setf(static_cast<std::ios::fmtflags>(0), std::ios::floatfield);
    std::cout.precision(6);
    std::cout << "DEBUG GAMMA_DISP: gamma_dispatcher called with ns=" << ns << ", nr=" << nr << std::endl;
    
    // Debug: Check what reactions we received
    std::cout << "DEBUG GAMMA_DISP: Received " << reactions.size() << " reactions" << std::endl;
    for (size_t i = 0; i < reactions.size(); ++i) {
      const auto &reaction = reactions[i];
      std::cout << "  Reaction " << i << " gammas has_value: " << (reaction.gammas.has_value() ? "true" : "false") << std::endl;
      if (reaction.gammas) {
        std::cout << "    gammas size: " << reaction.gammas->size() << std::endl;
        for (const auto& [sp, val] : *reaction.gammas) {
          std::cout << "    " << sp << " -> " << val << std::endl;
        }
      }
    }
    
    // Display parsed reactions when debugging.
    if (debug) {
      print_reactions(reactions);
    }

    if (!reactions.empty() && reactions[0].type == CatalysisModel::GammaBose) {
      // Dedicated handling for the Bose recombination model.
      const auto &p = reactions[0].bose.value();
      std::size_t idx_O = index.at("O");
      std::size_t idx_CO = index.at("CO");
      std::size_t idx_O2 = index.at("O2");
      std::size_t idx_CO2 = index.at("CO2");
      return Result<CatalysisFluxes>{bose::compute_fluxes(
          species, ns, idx_O, idx_CO, idx_O2, idx_CO2, p.gamma, p.p2)};
    }
    // Extract the 3 stoichiometric matrices.
    const auto &nu = stoic_matrices.nu;
    const auto &mu = stoic_matrices.mu;
    const auto &nu_p = stoic_matrices.nu_p;
    const auto &nu_diff = stoic_matrices.nu_diff;

    // Compute gammas or handle special models
    std::cout << "DEBUG GAMMA_DISP: all_super=" << all_super << ", all_rini=" << all_rini << std::endl;
    if (all_super || all_rini) {
      std::cout << "DEBUG GAMMA_DISP: Taking all_super/all_rini early return path!" << std::endl;
      auto reactions_copy = reactions;
      if (all_super) {
        for (auto &r : reactions_copy)
          r.gamma_w = 1.0;
      }
      // Cache Rini lookup matrices once since they are independent of the
      // dynamic state (rho and Tw).
      static std::optional<rini::RiniMatrices> matrices;
      if (!matrices) {
        matrices = rini::precompute_matrices(ns, nr, index, reactions_copy);
      }
      return Result<CatalysisFluxes>{
          rini::compute_fluxes(species, ns, nr, nu_p, nu_diff, *matrices)};
    }

    std::vector<std::vector<double>> gammas(nr, std::vector<double>(ns, 0.0));
    std::vector<CatalysisModel> types(nr, CatalysisModel::GammaGiven);
    std::vector<bool> hetero(nr, false);
    for (std::size_t r = 0; r < nr; ++r) {
      const auto &reaction = reactions[r];
      types[r] = reaction.type;
      hetero[r] = reaction.heterogeneous;
      
      // For GammaGiven reactions, use gamma values from the computed_gamma map
      if (reaction.type == CatalysisModel::GammaGiven) {
        for (const auto &[sp, coeff] : reaction.reactants) {
          std::size_t idx_r = index.at(sp);
          // Use computed_gamma which was pre-populated during initialization
          auto it = reaction.computed_gamma.find(sp);
          double gamma_value = (it != reaction.computed_gamma.end()) ? it->second : 0.0;
          gammas[r][idx_r] = gamma_value;
        }
      } else {
        // For non-GammaGiven reactions, gamma values are computed dynamically
        for (const auto &[sp, coeff] : reaction.reactants) {
          std::size_t idx_r = index.at(sp);
          gammas[r][idx_r] = 0.0; // Default or computed value
        }
      }
    }
    
    
    if (debug) {
      print_gammas(gammas, reactions, index);
      print_nu_mu(nu, mu, species_order);
    }
    auto result = compute_all_fluxes(species, ns, nr, gammas, nu, nu_p, mu,
                              first_order, limiting_fluxes, types, hetero,
                              species_order, debug);
    // Restore previous formatting
    std::cout.flags(old_flags);
    std::cout.precision(old_prec);
    return result;
  } catch (const std::exception &e) {
    return std::unexpected(Error{ErrorCode::RuntimeError, e.what()});
  }
}

} // namespace gasp2::catalysis::gamma
