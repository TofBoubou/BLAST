#include "gasp2/catalysis/gamma/gamma_model.hpp"
#include "gasp2/catalysis/gamma/gamma_dispatcher.hpp"
#include "gasp2/catalysis/gamma/gamma_input_parser.hpp"
#include "gasp2/catalysis/gamma/gamma_reaction_validation.hpp"
#include "gasp2/catalysis/gamma/types.hpp"
#include <algorithm>
#include <exception>
#include <iostream>
#include <unordered_map>

//===----------------------------------------------------------------------===//
// Gamma-model solver
//
// Provides an initialization step that parses reaction input, validates
// mass/charge conservation, and builds stoichiometric matrices. A subsequent
// lightweight flux computation uses the cached data to evaluate catalysis
// fluxes for arbitrary wall states without re-reading the input file.
//===----------------------------------------------------------------------===//

namespace gasp2::catalysis::gamma {
using gasp2::Error;
using gasp2::ErrorCode;

namespace {
//---------------------------------------------------------------------------
// Persistent storage for the gamma-model data and diagnostic matrices
//--------------------------------------------------------------------------

// Build stoichiometric matrices and production/consumption tensors.
//
// @param reactions Parsed reaction list.
// @param species Ordering of gas-phase species.
// @param index Mapping from species name to its index in \p species.
// @returns Populated stoichiometric matrices and tensors.
[[nodiscard]] StoichiometricMatrices compute_stoichiometric_matrices(
    const std::vector<ReactionInput> &reactions,
    const std::vector<std::string> &species,
    const std::unordered_map<std::string, std::size_t> &index) {
  StoichiometricMatrices M;
  std::size_t n_reac = reactions.size();
  std::size_t n_species = species.size();
  M.nu_p.assign(n_reac, std::vector<int>(n_species, 0));
  M.nu_pp.assign(n_reac, std::vector<int>(n_species, 0));
  M.nu_diff.assign(n_reac, std::vector<int>(n_species, 0));
  M.nu.assign(n_reac, std::vector<int>(n_species, 0));
  M.mu.assign(n_reac, std::vector<std::vector<int>>(
                          n_species, std::vector<int>(n_species, 0)));
  for (std::size_t i = 0; i < n_reac; ++i) {
    for (const auto &[sp, coeff] : reactions[i].reactants) {
      std::size_t idx_r = index.at(sp);
      M.nu_p[i][idx_r] = coeff;
      M.nu[i][idx_r] = 1; // species consumed
      for (const auto &[prod, coeff_p] : reactions[i].products) {
        std::size_t idx_p = index.at(prod);
        M.mu[i][idx_p][idx_r] = coeff_p > 0 ? 1 : 0;
      }
    }
    for (const auto &[sp, coeff] : reactions[i].products) {
      std::size_t idx_p = index.at(sp);
      M.nu_pp[i][idx_p] = coeff;
    }
    for (std::size_t j = 0; j < n_species; ++j)
      M.nu_diff[i][j] = M.nu_pp[i][j] - M.nu_p[i][j];
  }
  return M;
}

// Utility to print a stoichiometric matrix with species headers.
void print_matrix(const std::vector<std::vector<int>> &M,
                  const std::vector<std::string> &species, const char *name) {
  std::cout << name << " matrix:\n  ";
  for (const auto &sp : species)
    std::cout << sp << ' ';
  std::cout << '\n';
  for (std::size_t i = 0; i < M.size(); ++i) {
    std::cout << "  r" << i << ": ";
    for (int v : M[i])
      std::cout << v << ' ';
    std::cout << '\n';
  }
}

struct GammaContext {
  std::vector<std::string> species_order; ///< Species ordering.
  std::vector<double> molar_masses;       ///< Molar masses (kg/mol).
  std::vector<ReactionInput> reactions;   ///< Parsed reactions.
  bool first_order{false};     ///< True if reactions are first order.
  bool limiting_fluxes{false}; ///< Enforce reactant compatibility for
                               ///< heterogeneous reactions.
  std::unordered_map<std::string, std::size_t>
      index;               ///< Species name -> index map.
  bool debug{true};        ///< Debug flag.
  bool initialized{false}; ///< Initialization state.
  std::size_t n_species{0};
  std::size_t n_reactions{0};
  bool all_super{false};     ///< True if all reactions are super-catalytic.
  bool all_rini{false};      ///< True if all reactions use the Rini model.
  StoichiometricMatrices nu; ///< Cached stoichiometric and gamma matrices.
} ctx;
} // namespace

//---------------------------------------------------------------------------
// Initialization
//---------------------------------------------------------------------------
/// Initialize the gamma-model context by parsing reactions. All species
/// descriptors (names, sizes, and molar masses) must be validated by the caller
/// prior to invoking this function.
///
/// @param species_order Ordering of gas-phase species.
/// @param molar_masses Molar mass for each species (kg/mol).
/// @param input_filename XML file describing surface reactions.
/// @param debug Enable verbose diagnostic output.
///
/// @returns An empty result on success or an error describing the failure.
Result<void> initialize(const std::vector<std::string> &species_order,
                        std::span<const double> molar_masses,
                        const std::filesystem::path &input_filename,
                        bool debug) {
  try {
    // Normalize cout formatting for all debug prints in this routine
    auto cout_flags_guard = std::cout.flags();
    auto cout_prec_guard = std::cout.precision();
    std::cout.setf(static_cast<std::ios::fmtflags>(0), std::ios::floatfield);
    std::cout.precision(6);
    ctx = GammaContext{}; // Reset previous state
    ctx.debug = debug;

    // Species data are assumed to have been validated upstream by
    // gasp2::validate_init_inputs, so no further checks are performed here.
    ctx.species_order = species_order;
    ctx.molar_masses.assign(molar_masses.begin(), molar_masses.end());

    auto parsed =
        details::read_gamma_input_file(input_filename, ctx.species_order);
    ctx.reactions = std::move(parsed.reactions);
    ctx.first_order = parsed.first_order;
    ctx.limiting_fluxes = parsed.limiting_fluxes;
    ctx.n_species = ctx.species_order.size();
    ctx.n_reactions = ctx.reactions.size();

    ctx.index.clear();
    for (std::size_t i = 0; i < ctx.species_order.size(); ++i)
      ctx.index[ctx.species_order[i]] = i;

    ctx.all_super =
        !ctx.reactions.empty() &&
        std::all_of(ctx.reactions.begin(), ctx.reactions.end(),
                    [](const auto &r) {
                      return r.type == CatalysisModel::SuperCatalytic;
                    });
    ctx.all_rini = !ctx.reactions.empty() &&
                   std::all_of(ctx.reactions.begin(), ctx.reactions.end(),
                               [](const auto &r) {
                                 return r.type == CatalysisModel::RiniModel;
                               });

    details::static_validate_gamma_reactions(ctx.reactions, ctx.species_order);

    // Build stoichiometric matrices and nu/mu tensors once. These are useful
    // for diagnostics and for downstream numerical routines that require quick
    // access to stoichiometric coefficients.
    ctx.nu = compute_stoichiometric_matrices(ctx.reactions, ctx.species_order,
                                             ctx.index);
    
    // Debug: Print gammas from parsed reactions before pre-population
    std::cout << "DEBUG BEFORE PRE-POPULATION: Checking gammas in ctx.reactions:" << std::endl;
    for (size_t i = 0; i < ctx.reactions.size(); ++i) {
      const auto &reaction = ctx.reactions[i];
      std::cout << "  Reaction " << i << " gammas has_value: " << (reaction.gammas.has_value() ? "true" : "false") << std::endl;
      if (reaction.gammas) {
        std::cout << "    gammas size: " << reaction.gammas->size() << std::endl;
        for (const auto& [sp, val] : *reaction.gammas) {
          std::cout << "    " << sp << " -> " << val << std::endl;
        }
      }
    }
    
    std::cout << "DEBUG INIT MARKER: About to start pre-population" << std::endl;
    // Pre-populate computed_gamma values during initialization
    // This avoids the need to call verify_gamma_reactions during each flux computation
    std::cout << "DEBUG INIT: Starting pre-population with " << ctx.reactions.size() << " reactions" << std::endl;
    for (size_t i = 0; i < ctx.reactions.size(); ++i) {
      auto &reaction = ctx.reactions[i];  // This is a reference to modify the actual element
      std::cout << "DEBUG INIT: Reaction " << i << " type: " << static_cast<int>(reaction.type) << std::endl;
      std::cout << "DEBUG INIT: Reaction " << i << " gammas has_value: " << (reaction.gammas.has_value() ? "true" : "false") << std::endl;
      
      if (reaction.type == CatalysisModel::GammaGiven) {
        std::cout << "DEBUG INIT: Processing GammaGiven reaction " << i << std::endl;
        for (const auto &[sp_ref, _] : reaction.reactants) {
          std::string sp = sp_ref;  // Create a copy to avoid reference issues
          std::cout << "DEBUG INIT: Checking species " << sp << std::endl;
          if (reaction.gammas && reaction.gammas->find(sp) != reaction.gammas->end()) {
            double gamma_value = reaction.gammas->at(sp);
            std::cout << "DEBUG INIT: About to assign gamma_value " << gamma_value << " to computed_gamma[" << sp << "]" << std::endl;
            reaction.computed_gamma[sp] = gamma_value;
            std::cout << "DEBUG INIT: gamma_value from gammas map = " << gamma_value << std::endl;
            std::cout << "DEBUG INIT: Pre-populated computed_gamma[" << sp << "] = " << reaction.computed_gamma[sp] << std::endl;
            // Double check
            std::cout << "DEBUG INIT: Verification - computed_gamma[" << sp << "] = " << reaction.computed_gamma.at(sp) << std::endl;
            // Triple check with direct access to ctx
            std::cout << "DEBUG INIT: Direct check - ctx.reactions[" << i << "].computed_gamma[" << sp << "] = " << ctx.reactions[i].computed_gamma[sp] << std::endl;
          } else {
            std::cout << "DEBUG INIT: ERROR - No gamma found for species " << sp << std::endl;
            if (!reaction.gammas) {
              std::cout << "DEBUG INIT: reaction.gammas is null!" << std::endl;
            } else {
              std::cout << "DEBUG INIT: reaction.gammas size: " << reaction.gammas->size() << std::endl;
            }
          }
        }
      }
    }
    
    // Debug: Verify pre-population was successful
    std::cout << "DEBUG AFTER PRE-POPULATION: Verifying computed_gamma values:" << std::endl;
    for (size_t i = 0; i < ctx.reactions.size(); ++i) {
      const auto &reaction = ctx.reactions[i];
      std::cout << "  Reaction " << i << " computed_gamma size: " << reaction.computed_gamma.size() << std::endl;
      for (const auto& [sp, val] : reaction.computed_gamma) {
        std::cout << "    " << sp << " -> " << val << std::endl;
      }
    }

    if (ctx.debug) {
      print_matrix(ctx.nu.nu_p, ctx.species_order, "nu_p");
      print_matrix(ctx.nu.nu_pp, ctx.species_order, "nu_pp");
      print_matrix(ctx.nu.nu_diff, ctx.species_order, "nu_diff");
    }
    
    std::cout << "DEBUG INIT PRE-CHECK: Right after debug matrices print" << std::endl;
    std::cout << "DEBUG INIT END: ctx.reactions pointer = " << &ctx.reactions << std::endl;
    std::cout << "DEBUG INIT END: ctx.reactions.size() = " << ctx.reactions.size() << std::endl;
    std::cout << "DEBUG INIT END: Final verification of computed_gamma in ctx.reactions[0]: " << ctx.reactions[0].computed_gamma.size() << std::endl;
    for (const auto& [sp, val] : ctx.reactions[0].computed_gamma) {
      std::cout << "DEBUG INIT END: " << sp << " -> " << val << std::endl;
    }

    ctx.initialized = true;
    // Restore cout formatting
    std::cout.flags(cout_flags_guard);
    std::cout.precision(cout_prec_guard);
    return {};
  } catch (const std::exception &e) {
    return std::unexpected(Error{ErrorCode::RuntimeError, e.what()});
  }
}

//---------------------------------------------------------------------------
// Flux computation
//---------------------------------------------------------------------------
/// Compute catalysis fluxes using the gamma-model for a given wall state. Wall
/// temperature and species densities must be validated by the caller prior to
/// invocation.
///
/// @param T_wall Wall temperature (K).
/// @param rho_wall Species mass densities at the wall (kg/m^3) in the ordering
///                 specified during initialization.
///
/// @returns CatalysisFluxes on success or an error if the model has not been
///          initialized.
Result<CatalysisFluxes> compute_fluxes(double T_wall,
                                       std::span<const double> rho_wall) {
  std::cout << "DEBUG COMPUTE_FLUXES: Function called with T_wall=" << T_wall << std::endl;
  std::cout << "DEBUG COMPUTE_FLUXES: ctx.reactions pointer = " << &ctx.reactions << std::endl;
  std::cout << "DEBUG COMPUTE_FLUXES: ctx.reactions.size() = " << ctx.reactions.size() << std::endl;
  if (!ctx.initialized) {
    return std::unexpected(
        Error{ErrorCode::RuntimeError, "initialize_catalysis not called"});
  }
  SpeciesData species{rho_wall, ctx.molar_masses, T_wall};

  try {
    // Note: computed_gamma values are pre-populated during initialization
    // so we don't need to call verify_gamma_reactions here
    // Ensure debug numbers print with decimals (reset any external stream formatting)
    auto old_flags = std::cout.flags();
    auto old_prec = std::cout.precision();
    std::cout.setf(static_cast<std::ios::fmtflags>(0), std::ios::floatfield);
    std::cout.precision(6);
    std::cout << "DEBUG COMPUTE_FLUXES: About to call gamma_dispatcher" << std::endl;
    // Debug: Check computed_gamma values before dispatcher call
    std::cout << "DEBUG COMPUTE_FLUXES: Checking ctx.reactions computed_gamma values:" << std::endl;
    for (size_t i = 0; i < ctx.reactions.size(); ++i) {
      std::cout << "  Reaction " << i << " computed_gamma size: " << ctx.reactions[i].computed_gamma.size() << std::endl;
      for (const auto& [sp, val] : ctx.reactions[i].computed_gamma) {
        std::cout << "    " << sp << " -> " << val << std::endl;
      }
    }
    
    // Debug: Force re-check just before dispatcher call
    std::cout << "DEBUG COMPUTE_FLUXES: Final check just before dispatcher:" << std::endl;
    for (size_t i = 0; i < ctx.reactions.size(); ++i) {
      const auto &reaction = ctx.reactions[i];
      std::cout << "  Reaction " << i << " computed_gamma size: " << reaction.computed_gamma.size() << std::endl;
      for (const auto& [sp, val] : reaction.computed_gamma) {
        std::cout << "    computed_gamma[" << sp << "] = " << val << std::endl;
      }
      std::cout << "  Reaction " << i << " gammas has_value: " << (reaction.gammas.has_value() ? "true" : "false") << std::endl;
      if (reaction.gammas) {
        std::cout << "    gammas size: " << reaction.gammas->size() << std::endl;
        for (const auto& [sp, val] : *reaction.gammas) {
          std::cout << "    " << sp << " -> " << val << std::endl;
        }
      }
    }
    
    auto fluxes = gamma_dispatcher(
        species, ctx.n_species, ctx.n_reactions, ctx.reactions, ctx.first_order,
        ctx.limiting_fluxes, ctx.index, ctx.nu, ctx.all_super, ctx.all_rini,
        ctx.species_order, ctx.debug);
    if (!fluxes)
      return fluxes;
    if (ctx.debug) {
      std::cout << "Computed catalysis fluxes:\n";
      for (std::size_t i = 0; i < fluxes->cat_fluxes.size(); ++i) {
        std::cout << "  " << ctx.species_order[i] << ": "
                  << fluxes->cat_fluxes[i] << '\n';
      }
    }
    // Restore previous formatting
    std::cout.flags(old_flags);
    std::cout.precision(old_prec);
    return fluxes;
  } catch (const std::exception &e) {
    return std::unexpected(Error{ErrorCode::RuntimeError, e.what()});
  }
}

} // namespace gasp2::catalysis::gamma
