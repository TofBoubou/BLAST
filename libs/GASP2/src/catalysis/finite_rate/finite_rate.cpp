#include "gasp2/catalysis/finite_rate/finite_rate.hpp"
#include "gasp2/catalysis/finite_rate/finite_rate_flux_computation.hpp"
#include "gasp2/catalysis/finite_rate/finite_rate_input_parser.hpp"
#include "gasp2/catalysis/finite_rate/finite_rate_reaction_properties.hpp"
#include "gasp2/catalysis/finite_rate/finite_rate_reaction_validation.hpp"
#include "gasp2/catalysis/finite_rate/steady_surface_coverage.hpp"
#include "gasp2/catalysis/finite_rate/types.hpp"
#include <iostream>
#include <map>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace gasp2::catalysis::finite_rate {
using details::FiniteRateContext;
using gasp2::Error;
using gasp2::ErrorCode;

namespace details {
// Definition of the global finite-rate context declared in the header.
FiniteRateContext ctx;
} // namespace details

namespace {
// Bring the global context into the anonymous namespace for convenient access.
using details::ctx;

// Convert a ReactionType enum to a human-readable string for debugging output.
static std::string to_string(ReactionType t) {
  switch (t) {
  case ReactionType::Adsorption:
    return "Adsorption";
  case ReactionType::Desorption:
    return "Desorption";
  case ReactionType::AdsorptionDesorption:
    return "Adsorption/Desorption";
  case ReactionType::EleyRideal:
    return "Eley-Rideal";
  case ReactionType::LangmuirHinshelwood:
    return "Langmuir-Hinshelwood";
  }
  return "Unknown";
}

} // namespace

//----------------------------------------------------------------------------//
// Initialization
//----------------------------------------------------------------------------//
Result<void> initialize(const std::vector<std::string> &species_order,
                        std::span<const double> molar_masses,
                        const std::filesystem::path &input_filename,
                        bool debug) {
  try {
    ctx = FiniteRateContext{}; // Reset previous state
    ctx.species_order = species_order;
    ctx.molar_masses.assign(molar_masses.begin(), molar_masses.end());
    ctx.debug = debug;
    // Parse the input file, extracting the reactions and site density
    auto parsed = details::read_finite_rate_input_file(input_filename);
    ctx.site_density = parsed.site_density;
    ctx.reactions = std::move(parsed.reactions);
    // Validate the reactions (static validation)
    details::static_validate_finite_rate_reactions(ctx.reactions,
                                                   ctx.species_order);

    // ========== SPECIES LISTS ==========
    // Save all surface species (empty site excluded)
    std::set<std::string> surf_set;
    for (const auto &r : ctx.reactions) {
      auto gather = [&](const std::map<std::string, int> &side) {
        for (const auto &[sp, coeff] : side) {
          if (sp != "(s)" && sp.size() > 3 && sp.substr(sp.size() - 3) == "(s)")
            surf_set.insert(sp);
        }
      };
      gather(r.reactants);
      gather(r.products);
    }
    ctx.surface_species.assign(surf_set.begin(), surf_set.end());
    // Save all species (including surface species)
    ctx.all_species = ctx.species_order;
    ctx.all_species.insert(ctx.all_species.end(), ctx.surface_species.begin(),
                           ctx.surface_species.end());
    // Create an index mapping species names to their indices
    ctx.species_index.clear();
    for (std::size_t i = 0; i < ctx.all_species.size(); ++i)
      ctx.species_index[ctx.all_species[i]] = i;
    // ========== DEBUG ==========
    if (ctx.debug) {
      std::cout << "Parsed " << ctx.reactions.size()
                << " finite-rate reactions:\n";
      for (const auto &r : ctx.reactions) {
        std::cout << "Reaction: " << r.formula << '\n';
        std::cout << "  Type: " << to_string(r.type) << '\n';
        std::cout << "  Reactants:\n";
        for (const auto &[sp, coeff] : r.reactants)
          std::cout << "    " << coeff << ' ' << sp << '\n';
        std::cout << "  Products:\n";
        for (const auto &[sp, coeff] : r.products)
          std::cout << "    " << coeff << ' ' << sp << '\n';
        std::cout << "  Parameters:\n";
        for (const auto &[key, val] : r.parameters)
          std::cout << "    " << key << ": " << val << '\n';
        if (!r.kc.empty()) {
          std::cout << "  Kc:\n";
          for (const auto &[key, val] : r.kc)
            std::cout << "    " << key << ": " << val << '\n';
        }
      }
      if (!ctx.surface_species.empty()) {
        std::cout << "Surface species:";
        for (const auto &sp : ctx.surface_species)
          std::cout << ' ' << sp;
        std::cout << '\n';
      }
      // Print site density
      std::cout << "Site density: " << ctx.site_density << '\n';
    }
    // Set initialization flag
    ctx.initialized = true;
    return {};
  } catch (const std::exception &e) {
    return std::unexpected(Error{ErrorCode::RuntimeError, e.what()});
  }
}

//----------------------------------------------------------------------------//
// Flux computation (placeholder)
//----------------------------------------------------------------------------//
Result<CatalysisFluxes>
compute_fluxes(double T_wall, std::span<const double> rho_wall,
               std::span<const double> gibbs_free_energy) {
  if (!ctx.initialized) {
    return std::unexpected(
        Error{ErrorCode::RuntimeError,
              "catalysis::finite_rate::initialize not called"});
  }
  SpeciesData species{rho_wall, ctx.molar_masses, T_wall};
  (void)gibbs_free_energy;

  // ========== REACTION PROPERTIES ==========
  auto props =
      details::compute_reactions_properties(species, gibbs_free_energy);
  // ========== SURFACE COVERAGE ==========
  auto coverage = details::solve_surface_coverage(props);
  // ========== FLUXES ==========
  auto fluxes = details::compute_reaction_fluxes(props, coverage);

  // ========== SURFACE DENSITIES AND COVERAGES ==========
  fluxes.surface_densities.clear();
  fluxes.surface_coverages.clear();
  for (std::size_t j = 0; j < ctx.surface_species.size(); ++j) {
    double density = coverage[j];
    fluxes.surface_densities[ctx.surface_species[j]] = density;
    fluxes.surface_coverages[ctx.surface_species[j]] =
        density / ctx.site_density;
  }
  // Store the vacant surface site density and coverage using the generic
  // "(s)" key.
  fluxes.surface_densities["(s)"] = coverage.back();
  fluxes.surface_coverages["(s)"] = coverage.back() / ctx.site_density;

  // ========== LOSS FACTORS ==========
  fluxes.loss_factors.assign(ctx.species_order.size(), 0.0);
  for (std::size_t i = 0; i < ctx.species_order.size(); ++i) {
    double imp_mass_flux = species.imp_flux[i] * species.m[i];
    if (imp_mass_flux > 0.0) {
      fluxes.loss_factors[i] = fluxes.cat_fluxes[i] / imp_mass_flux;
    }
  }

  if (ctx.debug) {
    std::cout << "Finite-rate fluxes:\n";
    for (std::size_t i = 0; i < fluxes.cat_fluxes.size(); ++i) {
      std::cout << "  " << ctx.species_order[i] << ": " << fluxes.cat_fluxes[i]
                << ", gamma=" << fluxes.loss_factors[i] << '\n';
    }
    std::cout << "Surface densities (mol/m^2):\n";
    for (const auto &[name, theta] : fluxes.surface_densities) {
      std::cout << "  " << name << ": " << theta << '\n';
    }
    std::cout << "Surface coverages (fraction):\n";
    for (const auto &[name, theta] : fluxes.surface_coverages) {
      std::cout << "  " << name << ": " << theta << '\n';
    }
  }

  return fluxes;
}

} // namespace gasp2::catalysis::finite_rate
