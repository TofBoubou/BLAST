#include "gasp2/catalysis/gamma/gamma_reaction_validation.hpp"
#include "gasp2/catalysis/gamma/gamma_build.hpp"
#include <algorithm>
#include <cctype>
#include <map>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <iostream>

namespace gasp2::catalysis::gamma {
namespace details {

// Comprehensive validation routines for gamma-model reactions. Validation is
// split into two stages:
//
// 1. static_validate_gamma_reactions: executed once at initialization. It
//    checks for species validity, model-specific constraints, and mass/charge
//    conservation using only the input file information.
// 2. verify_gamma_reactions: executed for each flux evaluation. It computes
//    gamma values for the current wall state and enforces physical bounds
//    (0 <= gamma <= 1) and cumulative sums.
//
// VALIDATION CATEGORIES:
//  - Species validation: all species names must exist in the allowed list.
//  - Gamma parameter validation: ensure 0 <= gamma <= 1 and that cumulative
//    gammas per species remain within physical bounds.
//  - Model-specific requirements: checks unique constraints of GammaGiven,
//    GammaT, GammaConsistent, GammaBose, SuperCatalytic, and Rini models.
//  - Temperature-dependent parameters: verifies Arrhenius coefficients when
//    specified via GammaT.
//  - Mass & charge conservation: decomposes species into elemental
//    compositions and verifies equality across reactants/products along with
//    ionic charge balance.

// The helper utilities below are shared by both validation stages and are kept
// in an anonymous namespace to restrict their linkage to this translation
// unit, following modern C++ best practices instead of using `static`.

//---------------------------------------------------------------------------
// Helper utilities used by both validation stages
//---------------------------------------------------------------------------

namespace {

// Break a species name into its elemental components. Example: "CO2" ->
// {C:1, O:2}. Handles multi-character elements and multi-digit counts.
//
// Parameters:
//   sp - species string view to decompose into element/count pairs.
// Returns: map from element symbol to the number of occurrences.
std::map<std::string, int> decompose(std::string_view sp) {
  std::map<std::string, int> elems;
  for (size_t i = 0; i < sp.size();) {
    if (!std::isalpha(static_cast<unsigned char>(sp[i]))) {
      ++i;
      continue;
    }
    size_t j = i + 1;
    while (j < sp.size() && std::islower(static_cast<unsigned char>(sp[j])))
      ++j;
    std::string el(sp.substr(i, j - i));
    size_t k = j;
    while (k < sp.size() && std::isdigit(static_cast<unsigned char>(sp[k])))
      ++k;
    int count = (k == j) ? 1 : std::stoi(std::string(sp.substr(j, k - j)));
    elems[el] += count;
    i = k;
  }
  return elems;
}

// Compute total elemental composition for one side of a reaction by summing the
// contributions of each species after removing any ionic charge markers.
std::map<std::string, int>
compute_mass(const std::map<std::string, int> &side) {
  std::map<std::string, int> mass;
  for (const auto &[sp, coeff] : side) {
    std::string base = sp;
    while (!base.empty() && (base.back() == '+' || base.back() == '-'))
      base.pop_back();
    if (base == "e")
      continue;
    for (const auto &[el, count] : decompose(sp))
      mass[el] += count * coeff;
  }
  return mass;
}

// Compute electric charge for a reaction side by counting trailing '+' and '-'
// symbols on each species name.
int compute_charge(const std::map<std::string, int> &side) {
  int charge = 0;
  for (const auto &[sp, coeff] : side) {
    int q = 0;
    for (auto it = sp.rbegin(); it != sp.rend(); ++it) {
      if (*it == '+')
        ++q;
      else if (*it == '-')
        --q;
      else
        break;
    }
    charge += q * coeff;
  }
  return charge;
}

} // namespace

//---------------------------------------------------------------------------
// Initialization-time validation
//---------------------------------------------------------------------------

void static_validate_gamma_reactions(std::vector<ReactionInput> &reactions,
                                     const std::vector<std::string> &species) {
  std::unordered_set<std::string> allowed(species.begin(), species.end());

  for (auto &reaction : reactions) {
    // ========== SPECIES VALIDATION ==========
    auto check_allowed = [&](const std::map<std::string, int> &side) {
      for (const auto &[sp, _] : side) {
        if (!allowed.count(sp)) {
          throw std::runtime_error("Species not allowed: " + sp);
        }
      }
    };
    check_allowed(reaction.reactants);
    check_allowed(reaction.products);

    // ========== MODEL-SPECIFIC VALIDATION ==========
    switch (reaction.type) {
    case CatalysisModel::GammaConsistent: {
      if (!reaction.heterogeneous) {
        throw std::runtime_error(
            "GammaConsistent requires heterogeneous reactants");
      }
      if (!reaction.gammas || reaction.gammas->size() != 1) {
        throw std::runtime_error(
            "GammaConsistent requires exactly one gamma value");
      }
      auto given = *reaction.gammas->begin();
      if (!reaction.reactants.count(given.first)) {
        throw std::runtime_error("Gamma specified for non-reactant: " +
                                 given.first);
      }
      if (given.second < 0.0 || given.second > 1.0) {
        throw std::runtime_error("Gamma for species " + given.first +
                                 " must be within [0,1]");
      }
      break;
    }
    case CatalysisModel::GammaGiven: {
      if (reaction.gammaT && !reaction.gammaT->empty()) {
        throw std::runtime_error(
            "GammaT parameters specified for GammaGiven reaction");
      }
        if (!reaction.gammas || reaction.gammas->empty()) {
          throw std::runtime_error("GammaGiven reaction must specify at least one gamma value");
        }
          // Ensure every reactant species has a gamma specified
          for (const auto &[sp, _] : reaction.reactants) {
            if (!reaction.gammas->count(sp)) {
              throw std::runtime_error("GammaGiven reaction missing gamma for reactant: " + sp);
            }
          }
        for (const auto &[sp, val] : *reaction.gammas) {
          if (!allowed.count(sp)) {
            throw std::runtime_error("Species not allowed: " + sp);
          }
          if (!reaction.reactants.count(sp)) {
            throw std::runtime_error("Gamma specified for non-reactant: " + sp);
          }
          if (val < 0.0 || val > 1.0) {
            throw std::runtime_error("Gamma for species " + sp +
                                     " must be within [0,1]");
          }
        }
      break;
    }
    case CatalysisModel::GammaT: {
      if (reaction.gammas && !reaction.gammas->empty()) {
        throw std::runtime_error("Gammas specified for GammaT reaction");
      }
      if (reaction.gammaT) {
        for (const auto &[sp, _] : *reaction.gammaT) {
          if (!allowed.count(sp)) {
            throw std::runtime_error("Species not allowed: " + sp);
          }
          if (!reaction.reactants.count(sp)) {
            throw std::runtime_error(
                "GammaT parameters specified for non-reactant: " + sp);
          }
        }
      }
      break;
    }
    case CatalysisModel::GammaBose: {
      if (!reaction.bose) {
        throw std::runtime_error("GammaBose parameters missing");
      }
      const auto &b = reaction.bose.value();
      if (b.gamma < 0.0 || b.gamma > 1.0) {
        throw std::runtime_error("GammaBose gamma must be within [0,1]");
      }
      if (b.p2 < 0.0 || b.p2 > 1.0) {
        throw std::runtime_error("GammaBose p2 must be within [0,1]");
      }
      break;
    }
    case CatalysisModel::SuperCatalytic: {
      if ((reaction.gammas && !reaction.gammas->empty()) ||
          (reaction.gammaT && !reaction.gammaT->empty())) {
        throw std::runtime_error(
            "SuperCatalytic reactions must not specify gammas");
      }
      break;
    }
    case CatalysisModel::RiniModel: {
      if ((reaction.gammas && !reaction.gammas->empty()) ||
          (reaction.gammaT && !reaction.gammaT->empty())) {
        throw std::runtime_error(
            "Rini_model reactions must not specify gammas");
      }
      if (!reaction.gamma_w) {
        throw std::runtime_error("Rini_model requires gamma_w parameter");
      }
      if (*reaction.gamma_w < 0.0 || *reaction.gamma_w > 1.0) {
        throw std::runtime_error("gamma_w must be within [0,1]");
      }
      break;
    }
    }

    // ========== TAUTOLOGICAL REACTION CHECK ==========
    // Reject reactions where reactant and product stoichiometry are identical.
    // Example: "O + O -> 2 O" yields no net change and should be considered an
    // input error.
    if (reaction.reactants == reaction.products) {
      throw std::runtime_error(
          "Reaction has identical reactants and products: " + reaction.formula);
    }

    // ========== MASS AND CHARGE CONSISTENCY ==========
    if (compute_mass(reaction.reactants) != compute_mass(reaction.products)) {
      throw std::runtime_error("Reaction not mass balanced: " +
                               reaction.formula);
    }
    if (compute_charge(reaction.reactants) !=
        compute_charge(reaction.products)) {
      throw std::runtime_error("Reaction not charge balanced: " +
                               reaction.formula);
    }
  }
}

//---------------------------------------------------------------------------
// Runtime verification: compute gammas for the current wall state
//---------------------------------------------------------------------------

void verify_gamma_reactions(
    std::vector<ReactionInput> &reactions, const SpeciesData &species,
    const std::unordered_map<std::string, std::size_t> &index) {
  std::unordered_map<std::string, double> gamma_sums;

  for (auto &reaction : reactions) {
    if (reaction.type != CatalysisModel::GammaGiven &&
        reaction.type != CatalysisModel::GammaT &&
        reaction.type != CatalysisModel::GammaConsistent) {
      continue; // Other models handled separately
    }

    for (const auto &[sp, _] : reaction.reactants) {
      std::cout << "DEBUG verify_gamma_reactions - calling compute_gamma for species: " << sp << std::endl;
      std::cout.flush();
      double gamma =
          compute_gamma(reaction, sp, species.T_wall, species, index);
      std::cout << "DEBUG verify_gamma_reactions - compute_gamma returned: " << gamma << " for species: " << sp << std::endl;
      std::cout.flush();
      reaction.computed_gamma[sp] = gamma;
      if (gamma < 0.0 || gamma > 1.0) {
        throw std::runtime_error("Gamma for species " + sp +
                                 " must be within [0,1]");
      }
      gamma_sums[sp] += gamma;
    }
  }

  for (const auto &[sp, total] : gamma_sums) {
    if (total < 0.0 || total > 1.0) {
      throw std::runtime_error("Sum of gammas for species " + sp +
                               " must be within [0,1]");
    }
  }
}

} // namespace details
} // namespace gasp2::catalysis::gamma
