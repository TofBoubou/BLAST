#include "gasp2/catalysis/finite_rate/finite_rate_reaction_properties.hpp"
#include <cmath>
#include <iostream>
#include <numbers>
#include <set>

namespace gasp2::catalysis::finite_rate {
namespace details {

//---------------------------------------------------------------------------
// Helper utilities
//---------------------------------------------------------------------------
namespace {

constexpr double Pr = 101325.0;        ///< Reference pressure (Pa) == 1 ATM
constexpr double R = 8.31446261815324; ///< Universal gas constant (J/mol/K)
constexpr double Na = 6.02214076e23;   ///< Avogadro's number (1/mol)

/// @brief Check whether a species string denotes a surface species.
/// @param sp The species string to check
/// @return True if the species is a surface species, false otherwise
/// @note Surface species are identified by either being exactly "(s)"
///           or ending with "(s)" when the string is longer than 3 characters
static bool is_surface(const std::string &sp) {
  return sp == "(s)" || (sp.size() > 3 && sp.substr(sp.size() - 3) == "(s)");
}

/// @brief Retrieve the gas species participating in an adsorption reaction.
/// @param r The reaction input containing reactants and products
/// @return The name of the gas species found, or empty string if none found
static std::string find_gas_species(const ReactionInput &r) {
  // Searches through reactants first, then products, returning the
  // first non-surface species encountered. Assumes only one gas
  // species per adsorption reaction.
  for (const auto &[sp, _] : r.reactants)
    if (!is_surface(sp))
      return sp;
  for (const auto &[sp, _] : r.products)
    if (!is_surface(sp))
      return sp;
  return {};
}

/// @brief Identify the surface species assumed mobile in a L-H reaction.
/// @param r The reaction input containing reactants and products
/// @return Name of the surface species used for the 2D thermal velocity
///
/// The Langmuir-Hinshelwood mechanism involves two surface species.  In the
/// current model we do not explicitly specify which species is mobile.  If the
/// reaction features a single unique surface reactant (e.g., `2 O(s) -> ...`)
/// its name is returned.  When multiple distinct surface reactants are present
/// (e.g., `A(s) + B(s) -> ...`) the first one encountered is selected.  This is
/// an approximation and should be revised once the input format allows the
/// mobile surface species to be specified explicitly.
/// In order to match the molar mass, the suffix "(s)" is removed
static std::string find_surface_mobile_species(const ReactionInput &r) {
  // Collect all surface reactant species
  std::vector<std::string> surface_reactants;
  for (const auto &[sp, _] : r.reactants)
    if (is_surface(sp))
      surface_reactants.push_back(sp);

  if (surface_reactants.empty())
    throw std::runtime_error(
        "Langmuir-Hinshelwood reaction without surface reactants");

  // If only one surface reactant is present there is no ambiguity.  For
  // multiple reactants we default to the first one.
  std::string result = surface_reactants.front();
  // Remove the "(s)" suffix to get the base species name
  if (result.size() > 3 && result.substr(result.size() - 3) == "(s)")
    result = result.substr(0, result.size() - 3);
  else
    throw std::runtime_error("Invalid surface species format: " + result + " (expected format: species(s))");
  return result;
}

/// @brief Compute nu_g for a reaction (sum over gas species of nu_p - nu_r).
/// @param r The reaction input containing reactants and products with
/// stoichiometric coefficients
/// @return The net stoichiometric change in gas species (products - reactants)
static double compute_nu_g(const ReactionInput &r) {
  // Collect all unique gas species from both reactants and products
  std::set<std::string> gas_species;
  auto gather = [&](const std::map<std::string, int> &side) {
    for (const auto &[sp, _] : side)
      if (!is_surface(sp))
        gas_species.insert(sp);
  };
  gather(r.reactants);
  gather(r.products);

  // Calculate nu_g = sum over all gas species of (nu_products - nu_reactants)
  double nu_g = 0.0;
  for (const auto &sp : gas_species) {
    int nu_r = 0; // stoichiometric coefficient in reactants
    int nu_p = 0; // stoichiometric coefficient in products
    if (auto it = r.reactants.find(sp); it != r.reactants.end())
      nu_r = it->second;
    if (auto it = r.products.find(sp); it != r.products.end())
      nu_p = it->second;
    nu_g += static_cast<double>(nu_p - nu_r);
  }
  return nu_g;
}

/// @brief Retrieve the thermal speed (3D) of the specified gas species.
/// @param species The species data containing thermal speed information for all
/// species
/// @param sp The name of the gas species to look up
/// @return The 3D thermal speed of the species in m/s, or 0.0 if not found
static double thermal_speed(const gasp2::SpeciesData &species,
                            const std::string &sp) {
  // Look up the species index in the global species mapping
  auto it = ctx.species_index.find(sp);
  // Verify the species exists and the index is within bounds
  if (it != ctx.species_index.end() &&
      it->second < species.thermal_speed.size())
    return species.thermal_speed[it->second];
  // Throw error if species not found or index out of bounds - should not happen
  throw std::runtime_error("Thermal speed not found for species: " + sp);
}

/// @brief Retrieve the molar mass (kg/mol) of the given gas species.
/// @param sp The name of the gas species to look up
/// @return The molar mass of the species in kg/mol, or 0.0 if not found
static double molar_mass(const std::string &sp) {
  // Look up the species index in the global species mapping
  auto it = ctx.species_index.find(sp);
  // Verify the species exists and the index is within bounds
  if (it != ctx.species_index.end() && it->second < ctx.molar_masses.size())
    return ctx.molar_masses[it->second];
  // Throw error if species not found or index out of bounds - should not happen
  throw std::runtime_error("Molar mass not found for species: " + sp);
}

/// @brief Calculate the site-density exponent from surface reactant
/// stoichiometry.
/// @param r The reaction input containing reactants with stoichiometric
/// coefficients
/// @return The sum of stoichiometric coefficients for all surface reactants
static double compute_surface_reactant_exponent(const ReactionInput &r) {
  double rs = 0.0; // Site-density exponent from surface reactant stoichiometry
  for (const auto &[sp, coeff] : r.reactants)
    if (is_surface(sp))
      rs += coeff;
  return rs;
}

/// @brief Compute the forward rate constant for a pure adsorption reaction.
/// @param r The reaction input containing parameters and species information
/// @param species The species data containing thermal speeds and wall
/// temperature
/// @return The forward rate constant for adsorption in appropriate units
static double adsorption_rate_constant(const ReactionInput &r,
                                       const gasp2::SpeciesData &species) {
  // Extract the gas species involved in the adsorption reaction
  std::string gas_sp = find_gas_species(r);
  // Get the thermal speed of the gas species
  double v = thermal_speed(species, gas_sp);
  // Extract reaction parameters from input
  double S0 = r.parameters.at("S0");     // Sticking coefficient
  double beta = r.parameters.at("beta"); // Temperature exponent
  double Ead = r.parameters.at("Ead");   // Activation energy for adsorption
  double rs = compute_surface_reactant_exponent(r);
  // Calculate rate constant using kinetic theory formula:
  // k = (v / (4 * site_density^rs)) * S0 * T^beta * exp(-Ead / (R * T))
  return (v / (4.0 * std::pow(ctx.site_density, rs))) * S0 *
         std::pow(species.T_wall, beta) * std::exp(-Ead / (R * species.T_wall));
}

/// @brief Compute the desorption rate constant for a pure desorption reaction.
/// @param r The reaction input containing desorption parameters
/// @param T The temperature at which to evaluate the rate constant
/// @return The desorption rate constant in appropriate units
static double desorption_rate_constant(const ReactionInput &r, double T) {
  // Extract desorption parameters from reaction input
  double Ades =
      r.parameters.at("Ades"); // Pre-exponential factor for desorption
  double beta = r.parameters.at("beta"); // Temperature exponent
  double nu = r.parameters.at("v");      // Vibrational frequency
  double Edes = r.parameters.at("Edes"); // Activation energy for desorption
  // Calculate rate constant using transition state theory:
  // k = Ades * T^beta * nu * exp(-Edes / (R * T))
  return Ades * std::pow(T, beta) * nu * std::exp(-Edes / (R * T));
}

/// @brief Compute the equilibrium constant for adsorption/desorption when
/// explicitly specified in the input.
/// @param r The reaction input containing equilibrium constant parameters
/// @param T The temperature at which to evaluate the equilibrium constant
/// @return The equilibrium constant, or 0.0 if no parameters are specified
static double adsorption_equilibrium_constant(const ReactionInput &r,
                                              double T) {
  // Check if equilibrium constant parameters are provided
  if (r.kc.empty())
    return 0.0;
  // Extract equilibrium constant parameters from reaction input
  double Aeq = r.kc.at("A_eq");  // Pre-exponential factor for equilibrium
  double beta = r.kc.at("beta"); // Temperature exponent
  double K0 = r.kc.at("K0");     // Reference equilibrium constant
  double dE = r.kc.at("dE");     // Energy difference for equilibrium
  // Calculate equilibrium constant using Arrhenius-type expression:
  // Kc = Aeq * T^beta * K0 * exp(dE / (R * T))
  return Aeq * std::pow(T, beta) * K0 * std::exp(dE / (R * T));
}

/// @brief Compute the forward rate constant for the Eley-Rideal mechanism.
/// @param r The reaction input containing Eley-Rideal parameters
/// @param species The species data containing thermal speeds and wall
/// temperature
/// @return The forward rate constant for the Eley-Rideal reaction in
/// appropriate units
static double eley_rideal_rate_constant(const ReactionInput &r,
                                        const gasp2::SpeciesData &species) {
  // Extract the gas species involved in the Eley-Rideal reaction
  std::string gas_sp = find_gas_species(r);
  // Get the thermal speed of the gas species
  double v = thermal_speed(species, gas_sp);
  // Extract Eley-Rideal reaction parameters from input
  double gamma0 =
      r.parameters.at("gamma_er");       // Reaction probability coefficient
  double beta = r.parameters.at("beta"); // Temperature exponent
  double Eer = r.parameters.at("Eer");   // Activation energy for Eley-Rideal
  double rs = compute_surface_reactant_exponent(r);
  // Calculate rate constant using kinetic theory formula for Eley-Rideal
  // mechanism: k = (v / (4 * site_density^rs)) * gamma0 * T^beta * exp(-Eer /
  // (R * T))
  return (v / (4.0 * std::pow(ctx.site_density, rs))) * gamma0 *
         std::pow(species.T_wall, beta) * std::exp(-Eer / (R * species.T_wall));
}

/// @brief Compute the forward rate constant for the Langmuir-Hinshelwood
/// mechanism.
/// @param r The reaction input containing Langmuir-Hinshelwood parameters
/// @param species The species data containing thermal speeds and wall
/// temperature
/// @return The forward rate constant for the Langmuir-Hinshelwood reaction in
/// appropriate units
static double
langmuir_hinshelwood_rate_constant(const ReactionInput &r,
                                   const gasp2::SpeciesData &species) {
  // Identify which surface species is assumed to be mobile.  When multiple
  // surface reactants are present the first one is used as a heuristic.
  std::string mobile_sp = find_surface_mobile_species(r);
  // Use the molar mass of the mobile surface species to compute the 2D thermal
  // velocity. This assumes the surface species mass is provided in the input.
  double M = molar_mass(mobile_sp);
  // Calculate 2D thermal velocity: v2d = sqrt(pi * R * T / (2 * M))
  double v2d = std::sqrt(std::numbers::pi * R * species.T_wall / (2.0 * M));

  // Calculate site-density exponent from surface reactant stoichiometry
  double rs = compute_surface_reactant_exponent(r);

  // Extract Langmuir-Hinshelwood reaction parameters from input
  double Clh = r.parameters.at("C_lh");  // Pre-exponential factor
  double beta = r.parameters.at("beta"); // Temperature exponent
  double Elh = r.parameters.at("Elh");   // Activation energy for L-H mechanism

  // Calculate rate constant using Langmuir-Hinshelwood formula:
  // k = v2d * site_density^(1.5-rs) * sqrt(Na) * Clh * T^beta * exp(-Elh / (R *
  // T))
  return v2d * std::pow(ctx.site_density, 1.5 - rs) * std::sqrt(Na) * Clh *
         std::pow(species.T_wall, beta) * std::exp(-Elh / (R * species.T_wall));
}

} // namespace

//---------------------------------------------------------------------------
// Main processing routine
//---------------------------------------------------------------------------
/// @brief Compute reaction properties including rate constants and equilibrium
/// constants for all reactions.
/// @param species The species data containing densities, temperatures, and
/// thermal speeds
/// @param gibbs_free_energy Array of Gibbs free energies for all gas species
/// @return Complete reaction properties structure with processed reactions and
/// concentrations
ReactionProperties
compute_reactions_properties(const gasp2::SpeciesData &species,
                             std::span<const double> gibbs_free_energy) {
  // Initialize reaction properties structure
  ReactionProperties props;
  props.site_density = ctx.site_density;
  props.all_species = ctx.all_species;
  // ========== CONCENTRATIONS ==========
  // Convert mass densities to molar concentrations for all species
  for (std::size_t i = 0; i < ctx.species_order.size(); ++i) {
    // Density (kg/m^3) divided by molar mass (kg/mol) yields molar
    // concentration (mol/m^3)
    props.concentration[ctx.species_order[i]] =
        species.rho[i] / ctx.molar_masses[i];
  }

  // Map storing Gibbs energy jumps for surface species: G_ads - G_site.
  std::map<std::string, double> deltaG_ads;

  // ========== ADSORPTION/DESORPTION PAIRS ==========
  // Create lookup maps to pair adsorption and desorption reactions
  std::map<std::string, const ReactionInput *> ads_map;
  std::map<std::string, const ReactionInput *> des_map;
  std::map<std::string, const ReactionInput *> ad_des_map;
  // Categorize reactions by type and surface species involved
  for (const auto &r : ctx.reactions) {
    std::string surf;
    if (r.type == ReactionType::Adsorption ||
        r.type == ReactionType::AdsorptionDesorption) {
      // Find the surface species in products (excluding bare sites)
      for (const auto &[sp, _] : r.products)
        if (sp != "(s)" && is_surface(sp))
          surf = sp;
      if (r.type == ReactionType::Adsorption)
        ads_map[surf] = &r;
      else
        ad_des_map[surf] = &r; // Combined adsorption/desorption reaction
    }
    if (r.type == ReactionType::Desorption) {
      // Find the surface species in reactants (excluding bare sites)
      for (const auto &[sp, _] : r.reactants)
        if (sp != "(s)" && is_surface(sp))
          surf = sp;
      des_map[surf] = &r;
    }
  }

  // Process each surface species to compute adsorption/desorption kinetics
  for (const auto &sp : ctx.surface_species) {
    const ReactionInput *ads = nullptr;
    const ReactionInput *des = nullptr;
    // Find corresponding adsorption and desorption reactions
    if (ad_des_map.count(sp)) {
      ads = ad_des_map[sp];
      des = ad_des_map[sp]; // same reaction describes both directions
    } else {
      if (ads_map.count(sp))
        ads = ads_map[sp];
      if (des_map.count(sp))
        des = des_map[sp];
    }
    if (!ads || !des)
      throw std::runtime_error(
          "Missing adsorption or desorption reaction for surface species: " +
          sp);

    // Compute forward (adsorption) and backward (desorption) rate constants
    // and derive the equilibrium constant. Treat the pair as a single
    // reversible reaction A + (s) <-> A(s).
    double kf_ads = adsorption_rate_constant(*ads, species);
    double kf_des = 0.0;
    double Kc = 0.0;
    // Handle different ways of specifying equilibrium
    if (ads == des) {
      // Single reversible reaction with explicit equilibrium constant
      Kc = adsorption_equilibrium_constant(*ads, species.T_wall);
      if (Kc != 0.0)
        kf_des = kf_ads / Kc;
    } else {
      // Separate adsorption and desorption reactions
      kf_des = desorption_rate_constant(*des, species.T_wall);
      if (kf_des != 0.0)
        Kc = kf_ads / kf_des;
    }
    // Convert concentration-based to activity-based equilibrium constant
    if (Kc == 0.0)
      throw std::runtime_error(
          "Zero equilibrium constant for surface species: " + sp);
    double nu_g = compute_nu_g(*ads);
    double Ka = Kc * std::pow(Pr / (R * species.T_wall), -nu_g);

    // Store the processed reaction
    props.reactions.push_back(ProcessedReaction{ads->reactants, ads->products,
                                                kf_ads, kf_des, Kc, Ka});

    // Calculate Gibbs energy difference for surface species
    std::string gas_sp = find_gas_species(*ads);
    auto it = ctx.species_index.find(gas_sp);
    if (it == ctx.species_index.end() || it->second >= gibbs_free_energy.size())
      throw std::runtime_error(
          "Cannot compute Gibbs energy difference for surface species " + sp +
          ": gas species " + gas_sp + " not found or missing energy data");
    double Ga = gibbs_free_energy[it->second];
    deltaG_ads[sp] = Ga - R * species.T_wall * std::log(Ka);
    // If debug, print them
    if (ctx.debug) {
      std::cout << "Debug: Surface species " << sp << ", Gas species " << gas_sp
                << ", Ga = " << Ga << ", Ka = " << Ka
                << ", DeltaG_ads = " << deltaG_ads[sp] << " J/mol" << std::endl;
    }
  }

  // ========== OTHER REACTIONS ==========
  // Process Eley-Rideal and Langmuir-Hinshelwood surface reactions
  for (const auto &r : ctx.reactions) {
    // Skip adsorption/desorption reactions (already processed above)
    if (r.type == ReactionType::Adsorption ||
        r.type == ReactionType::Desorption ||
        r.type == ReactionType::AdsorptionDesorption)
      continue;

    // Calculate reaction Gibbs energy change from species energies
    double deltaG = 0.0;
    auto accumulate = [&](const std::map<std::string, int> &side, int sign) {
      for (const auto &[sp, coeff] : side) {
        if (sp == "(s)") // Skip bare surface sites
          continue;
        if (is_surface(sp)) {
          // Use previously calculated surface species energies
          if (auto it = deltaG_ads.find(sp); it != deltaG_ads.end())
            deltaG += sign * coeff * it->second;
        } else {
          // Use gas-phase Gibbs energies
          auto idx = ctx.species_index.find(sp);
          if (idx != ctx.species_index.end() &&
              idx->second < gibbs_free_energy.size())
            deltaG += sign * coeff * gibbs_free_energy[idx->second];
        }
      }
    };
    accumulate(r.products, +1);  // Products contribute positively
    accumulate(r.reactants, -1); // Reactants contribute negatively

    // Calculate equilibrium constants from thermodynamics
    double Ka = std::exp(-deltaG / (R * species.T_wall));
    double nu_g = compute_nu_g(r);
    double Kc = Ka * std::pow(Pr / (R * species.T_wall), nu_g);

    // Calculate forward rate constant based on reaction mechanism
    double kf = 0.0;
    if (r.type == ReactionType::EleyRideal)
      kf = eley_rideal_rate_constant(r, species);
    else if (r.type == ReactionType::LangmuirHinshelwood)
      kf = langmuir_hinshelwood_rate_constant(r, species);
    // Calculate backward rate constant from equilibrium
    if (Kc == 0.0)
      throw std::runtime_error("Zero equilibrium constant for "
                               "Eley-Rideal/Langmuir-Hinshelwood reaction");
    double kb = kf / Kc;

    // Store the processed reaction
    props.reactions.push_back(
        ProcessedReaction{r.reactants, r.products, kf, kb, Kc, Ka});
  }

  return props;
}

} // namespace details
} // namespace gasp2::catalysis::finite_rate
