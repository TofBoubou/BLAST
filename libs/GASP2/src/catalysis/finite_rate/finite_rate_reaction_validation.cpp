#include "gasp2/catalysis/finite_rate/finite_rate_reaction_validation.hpp"
#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace gasp2::catalysis::finite_rate {
namespace details {

/**
 * @brief Removes surface site markers "(s)" from species names while preserving
 * ionic charges
 *
 * This function extracts the base species name by removing the surface marker
 * "(s)" suffix commonly used in heterogeneous catalysis to denote surface-bound
 * species. The function preserves any ionic charge markers (+ or -) to maintain
 * proper charge balance validation.
 *
 * @param sp The species name string to process (e.g., "O(s)", "H+(s)", "(s)")
 * @return The base species name without surface markers:
 *         - "O(s)" -> "O"
 *         - "H+(s)" -> "H+"
 *         - "(s)" -> "" (empty string for bare surface sites)
 *         - "CO2" -> "CO2" (unchanged if no surface marker)
 */
static std::string strip_surface(const std::string &sp) {
  // Handle bare surface site case - return empty string to denote surface site
  if (sp == "(s)")
    return {};

  // Check if species has surface marker suffix and strip it
  // Requires minimum length check to avoid substr out-of-bounds
  if (sp.size() > 3 && sp.substr(sp.size() - 3) == "(s)")
    return sp.substr(0, sp.size() - 3);

  // Return unchanged if no surface marker found
  return sp;
}

/**
 * @brief Parses a chemical species string into its elemental composition
 *
 * This function decomposes a chemical formula into a map of element symbols
 * and their stoichiometric counts. It first strips surface markers using
 * strip_surface(), then parses the remaining string for element-count pairs.
 * The parser handles standard chemical notation where elements start with
 * uppercase letters followed by optional lowercase letters, and counts are
 * represented by digits immediately following the element symbol.
 *
 * @param sp_in The input species string (e.g., "CO2(s)", "H2O", "CaCl2")
 * @return A map where keys are element symbols and values are their counts:
 *         - "CO2(s)" -> {{"C", 1}, {"O", 2}}
 *         - "H2O" -> {{"H", 2}, {"O", 1}}
 *         - "CaCl2" -> {{"Ca", 1}, {"Cl", 2}}
 *         - "O" -> {{"O", 1}} (implicit count of 1)
 */
static std::map<std::string, int> decompose(const std::string &sp_in) {
  // Remove surface markers to get clean chemical formula
  std::string sp = strip_surface(sp_in);
  std::map<std::string, int> elems;

  // Parse the formula character by character
  for (size_t i = 0; i < sp.size();) {
    // Skip non-alphabetic characters (charges, punctuation, etc.)
    if (!std::isalpha(static_cast<unsigned char>(sp[i]))) {
      ++i;
      continue;
    }

    // Find element symbol: uppercase letter followed by optional lowercase
    // letters
    size_t j = i + 1;
    while (j < sp.size() && std::islower(static_cast<unsigned char>(sp[j])))
      ++j;
    std::string el =
        sp.substr(i, j - i); // Extract element symbol (e.g., "Ca", "Cl")

    // Find the count: digits immediately following the element symbol
    size_t k = j;
    while (k < sp.size() && std::isdigit(static_cast<unsigned char>(sp[k])))
      ++k;

    // Parse count (default to 1 if no digits found)
    int count = (k == j) ? 1 : std::stoi(sp.substr(j, k - j));

    // Add to elemental composition (accumulate if element appears multiple
    // times)
    elems[el] += count;

    // Move to next element
    i = k;
  }
  return elems;
}

/**
 * @brief Computes the total elemental composition for one side of a chemical
 * reaction
 *
 * This function calculates the mass balance by determining the total number of
 * each element present on one side of a reaction, accounting for stoichiometric
 * coefficients. It excludes pure surface sites (empty base species) and
 * electrons ("e") from the mass balance calculation, as these don't contribute
 * to elemental mass conservation. This is essential for validating that
 * reactions conserve mass.
 *
 * @param side A map representing one side of a reaction, where keys are species
 * names and values are their stoichiometric coefficients
 * @return A map where keys are element symbols and values are the total count
 * of each element on this reaction side:
 *         - For side {"CO2(s)": 2, "H2O": 3}: returns {{"C": 2}, {"O": 7},
 * {"H": 6}}
 *         - Pure surface sites "(s)" and electrons "e" are ignored
 */
static std::map<std::string, int>
compute_mass(const std::map<std::string, int> &side) {
  std::map<std::string, int> mass;

  // Iterate through each species on this reaction side
  for (const auto &[sp, coeff] : side) {
    // Strip surface markers to get the base chemical formula
    std::string base = strip_surface(sp);

    // Skip pure surface sites (empty string) and electrons - they don't
    // contribute to mass
    if (base.empty() || base == "e")
      continue;

    // Decompose the species into elements and accumulate with stoichiometric
    // coefficient
    for (const auto &[el, count] : decompose(sp))
      mass[el] +=
          count * coeff; // Multiply element count by stoichiometric coefficient
  }
  return mass;
}

/**
 * @brief Calculates the total electric charge for one side of a chemical
 * reaction
 *
 * This function computes the net electric charge by parsing ionic charge
 * markers
 * (+ and -) from species names and accounting for stoichiometric coefficients.
 * It uses reverse iteration to parse charge markers from the end of species
 * names, where ionic charges are typically located (e.g., "Ca2+", "SO4-2").
 * Pure surface sites are excluded as they carry no charge. This is essential
 * for validating that reactions maintain charge neutrality.
 *
 * @param side A map representing one side of a reaction, where keys are species
 * names and values are their stoichiometric coefficients
 * @return The total electric charge on this reaction side:
 *         - For side {"Ca2+": 1, "Cl-": 2}: returns 1*2 + 2*(-1) = 0
 *         - For side {"H+": 3, "PO4-3": 1}: returns 3*1 + 1*(-3) = 0
 *         - Neutral species contribute 0 to the total charge
 */
static int compute_charge(const std::map<std::string, int> &side) {
  int charge = 0;

  // Iterate through each species on this reaction side
  for (const auto &[sp, coeff] : side) {
    // Strip surface markers to get the base species with charge markers
    std::string base = strip_surface(sp);

    // Skip pure surface sites (empty string) - they carry no charge
    if (base.empty())
      continue;

    // Parse ionic charge from the end of the species name
    int q = 0;
    // Start from the end and work backwards to find consecutive charge markers
    for (auto it = base.rbegin(); it != base.rend(); ++it) {
      if (*it == '+')
        ++q; // Positive charge marker
      else if (*it == '-')
        --q; // Negative charge marker
      else
        break; // Stop at first non-charge character
    }

    // Accumulate total charge considering stoichiometric coefficient
    charge += q * coeff;
  }
  return charge;
}

/**
 * @brief Counts the total number of surface sites consumed or produced on one
 * side of a reaction
 *
 * This function calculates the surface site balance by identifying all surface
 * species (those with "(s)" markers) and summing their stoichiometric
 * coefficients. Surface sites are fundamental resources in heterogeneous
 * catalysis - they represent locations on the catalyst surface where
 * adsorption, desorption, and surface reactions can occur. Both bare surface
 * sites "(s)" and occupied sites like "O(s)" are counted. This is essential for
 * validating that reactions conserve surface sites.
 *
 * @param side A map representing one side of a reaction, where keys are species
 * names and values are their stoichiometric coefficients
 * @return The total number of surface sites on this reaction side:
 *         - For side {"O(s)": 2, "(s)": 1, "CO2": 1}: returns 2 + 1 = 3
 *         - For side {"H2O": 1}: returns 0 (no surface species)
 *         - Both bare "(s)" and occupied "X(s)" sites are counted
 */
static int compute_sites(const std::map<std::string, int> &side) {
  int sites = 0;

  // Iterate through each species on this reaction side
  for (const auto &[sp, coeff] : side) {
    // Check if species is a surface species (has "(s)" marker)
    // Include both bare surface sites "(s)" and occupied surface sites "X(s)"
    if (sp == "(s)" || (sp.size() > 3 && sp.substr(sp.size() - 3) == "(s)"))
      sites += coeff; // Add stoichiometric coefficient to site count
  }
  return sites;
}

/**
 * @brief Determines whether a reaction side contains at least one surface
 * species
 *
 * This function checks if any species on one side of a reaction involves the
 * catalyst surface, identified by the "(s)" marker. In heterogeneous catalysis,
 * valid reactions must involve the surface on both sides to ensure proper
 * interaction with the catalyst. This includes both bare surface sites "(s)"
 * and occupied surface sites like "O(s)". This is a critical validation step to
 * ensure reactions are appropriate for surface catalysis modeling.
 *
 * @param side A map representing one side of a reaction, where keys are species
 * names and values are their stoichiometric coefficients
 * @return true if at least one surface species is found, false otherwise:
 *         - For side {"O(s)": 1, "CO2": 1}: returns true (has O(s))
 *         - For side {"(s)": 2, "H2O": 1}: returns true (has bare sites)
 *         - For side {"H2O": 1, "CO2": 1}: returns false (no surface species)
 */
static bool has_surface_species(const std::map<std::string, int> &side) {
  // Scan through all species looking for surface markers
  for (const auto &[sp, _] : side) {
    // Check for surface species: bare sites "(s)" or occupied sites "X(s)"
    // Use size check to avoid substr out-of-bounds for short strings
    if (sp == "(s)" || (sp.size() > 3 && sp.substr(sp.size() - 3) == "(s)"))
      return true; // Found at least one surface species - early return
  }
  return false; // No surface species found on this side
}

/**
 * @brief Counts the number of distinct gas-phase species on one side of a
 * reaction
 *
 * This function determines how many different gas-phase (non-surface) species
 * are present on one side of a reaction by excluding all surface species (those
 * with
 * "(s)" markers). In heterogeneous catalysis modeling, reactions are typically
 * limited to at most one gas species per side to maintain kinetic simplicity
 * and physical realism. This count is used to enforce reaction complexity
 * constraints.
 *
 * @param side A map representing one side of a reaction, where keys are species
 * names and values are their stoichiometric coefficients
 * @return The number of distinct gas-phase species on this reaction side:
 *         - For side {"CO2": 1, "O(s)": 2, "(s)": 1}: returns 1 (only CO2 is
 * gas)
 *         - For side {"H2": 1, "O2": 1, "(s)": 2}: returns 2 (H2 and O2 are
 * gas)
 *         - For side {"O(s)": 1, "(s)": 1}: returns 0 (no gas species)
 */
static int count_gas_species(const std::map<std::string, int> &side) {
  int count = 0;

  // Iterate through all species on this reaction side
  for (const auto &[sp, _] : side) {
    // Count species that are NOT surface species (inverse of surface check)
    // Exclude bare surface sites "(s)" and occupied surface sites "X(s)"
    if (!(sp == "(s)" || (sp.size() > 3 && sp.substr(sp.size() - 3) == "(s)")))
      ++count; // This is a gas-phase species - increment counter
  }
  return count;
}

/**
 * @brief Creates a canonical string representation for one side of a chemical
 * reaction
 *
 * This function generates a normalized, deterministic string key for a reaction
 * side by separating gas-phase and surface species, then sorting each group
 * alphabetically by species name. The canonical format ensures that equivalent
 * reaction sides (with the same species and coefficients) produce identical
 * keys regardless of input order. This is essential for duplicate reaction
 * detection and reaction signature comparison. The output format is:
 * gas_species+surface_species with coefficients prefixed.
 *
 * @param side A map representing one side of a reaction, where keys are species
 * names and values are their stoichiometric coefficients
 * @return A canonical string key for this reaction side:
 *         - For side {"CO2": 1, "O(s)": 2, "(s)": 1}: returns "1CO2+1(s)+2O(s)"
 *         - For side {"(s)": 2, "H2O": 1}: returns "1H2O+2(s)"
 *         - Gas species appear first, then surface species, both sorted
 * alphabetically
 */
static std::string make_key(const std::map<std::string, int> &side) {
  // Separate species into gas-phase and surface categories
  std::vector<std::pair<std::string, int>> gas;
  std::vector<std::pair<std::string, int>> surface;

  for (const auto &kv : side) {
    // Classify each species based on presence of "(s)" marker
    if (kv.first == "(s)" ||
        (kv.first.size() > 3 && kv.first.substr(kv.first.size() - 3) == "(s)"))
      surface.push_back(kv); // Surface species (bare or occupied sites)
    else
      gas.push_back(kv); // Gas-phase species
  }

  // Sort both categories alphabetically by species name for consistent ordering
  auto cmp = [](const auto &a, const auto &b) { return a.first < b.first; };
  std::sort(gas.begin(), gas.end(), cmp);
  std::sort(surface.begin(), surface.end(), cmp);

  // Build the canonical key string
  std::string key;
  auto append = [&](const std::vector<std::pair<std::string, int>> &vec) {
    for (const auto &kv : vec) {
      if (!key.empty())
        key += "+"; // Add separator between species
      // Format: coefficient + species_name (e.g., "2CO2", "1(s)")
      key += std::to_string(kv.second) + kv.first;
    }
  };

  // Append gas species first, then surface species
  append(gas);
  append(surface);
  return key;
}

/**
 * @brief Validates a complete set of finite rate heterogeneous catalysis
 * reactions
 *
 * This is the main validation function that performs comprehensive checks on a
 * set of reactions for heterogeneous catalysis modeling. It validates species
 * names, stoichiometric coefficients, conservation laws (mass, charge, surface
 * sites), reaction complexity constraints, required kinetic parameters, and
 * ensures proper adsorption/desorption pairing. The function enforces all
 * constraints necessary for physically meaningful and computationally tractable
 * finite rate kinetics.
 *
 * @param reactions Vector of reaction inputs to validate, each containing
 * reactants, products, reaction type, and kinetic parameters
 * @param species_order List of allowed species names for validation - only
 * species in this list (after surface marker removal) are permitted
 * @throws std::runtime_error for any validation failure with descriptive
 * message
 */
void static_validate_finite_rate_reactions(
    const std::vector<ReactionInput> &reactions,
    const std::vector<std::string> &species_order) {
  // Initialize validation state containers: allowed gas species
  std::unordered_set<std::string> allowed(species_order.begin(),
                                          species_order.end());
  std::set<std::string> seen_signatures;     // Track duplicate reactions
  std::set<std::string> adsorption_sigs;     // Forward reaction signatures
  std::set<std::string> desorption_sigs;     // Reverse reaction signatures
  std::set<std::string> all_surface_species; // All surface species found
  std::set<std::string>
      adsorption_species; // Surface species with adsorption reactions

  // Main validation loop - process each reaction
  for (const auto &r : reactions) {
    // Lambda for validating one side of a reaction
    auto check_side = [&](const std::map<std::string, int> &side) {
      for (const auto &[sp, coeff] : side) {
        // ========== BASIC SPECIES VALIDATION ==========
        std::string base = strip_surface(sp);
        // Check if species is in allowed list (skip empty surface sites)
        if (!base.empty() && !allowed.count(base))
          throw std::runtime_error("Species not allowed: " + base);
        // Check for unsupported ionic species (future extension possibility)
        if (base.find('+') != std::string::npos ||
            base.find('-') != std::string::npos) {
          throw std::runtime_error("Ionic species not supported: " + sp +
                                   " (maybe extend to ions in the future...)");
        }
        // ========== STOICHIOMETRIC COEFFICIENT CHECK ==========
        // This should be redundant, already done in parser
        if (coeff < 1) {
          throw std::runtime_error(
              "Stoichiometric coefficients must be >= 1: " + r.formula);
        }
      }
    };
    // Validate both sides of the reaction
    check_side(r.reactants);
    check_side(r.products);

    // ========== HETEROGENEOUS CATALYSIS CONSTRAINTS ==========
    // Require surface species on both sides for heterogeneous reactions
    if (!has_surface_species(r.reactants) || !has_surface_species(r.products)) {
      throw std::runtime_error(
          "At least one surface species required on each side: " + r.formula);
    }
    // Limit gas species complexity (max 1 per side for kinetic tractability)
    if (count_gas_species(r.reactants) > 1 ||
        count_gas_species(r.products) > 1) {
      throw std::runtime_error("Too many gas species in reaction: " +
                               r.formula);
    }
    // Limit overall reaction complexity (max 2 species per side)
    if (r.reactants.size() > 2)
      throw std::runtime_error("Too many reactants in reaction: " + r.formula);
    if (r.products.size() > 2)
      throw std::runtime_error("Too many products in reaction: " + r.formula);

    // ========== CONSERVATION LAW VALIDATION ==========
    // Ensure mass balance (elemental composition conservation)
    if (compute_mass(r.reactants) != compute_mass(r.products)) {
      throw std::runtime_error("Reaction not mass balanced: " + r.formula);
    }
    // Ensure charge balance (electric charge conservation)
    if (compute_charge(r.reactants) != compute_charge(r.products)) {
      throw std::runtime_error("Reaction not charge balanced: " + r.formula);
    }
    // Ensure surface site balance (catalyst site conservation)
    if (compute_sites(r.reactants) != compute_sites(r.products)) {
      throw std::runtime_error("Reaction not site balanced: " + r.formula);
    }

    // ========== PARAMETER VALIDATION ==========
    // Lambda to verify that all required kinetic parameters are present and
    // that no extraneous parameters were supplied by the user.
    auto validate_params = [&](std::initializer_list<const char *> names,
                               const std::map<std::string, double> &provided,
                               const std::string &ctx) {
      // Track allowed parameter names for quick membership checks
      std::unordered_set<std::string> allowed(names.begin(), names.end());

      // Ensure every required parameter is supplied
      for (const char *n : names) {
        if (!provided.count(n))
          throw std::runtime_error("Missing parameter '" + std::string(n) +
                                   "' for " + ctx + " reaction: " + r.formula);
      }

      // Reject any user-specified parameters that are not explicitly supported
      for (const auto &kv : provided) {
        const auto &key = kv.first;
        if (!allowed.count(key))
          throw std::runtime_error("Unsupported parameter '" + key + "' for " +
                                   ctx + " reaction: " + r.formula);
      }
    };
    // Check required parameters based on reaction mechanism type
    switch (r.type) {
    case ReactionType::Adsorption:
      validate_params({"S0", "beta", "Ead"}, r.parameters, "adsorption");
      validate_params({}, r.kc, "adsorption equilibrium constants");
      break;
    case ReactionType::Desorption:
      validate_params({"Ades", "beta", "v", "Edes"}, r.parameters,
                      "desorption");
      validate_params({}, r.kc, "desorption equilibrium constants");
      break;
    case ReactionType::AdsorptionDesorption:
      validate_params({"S0", "beta", "Ead"}, r.parameters,
                      "adsorption/desorption");
      validate_params({"A_eq", "beta", "K0", "dE"}, r.kc,
                      "adsorption/desorption equilibrium");
      break;
    case ReactionType::EleyRideal:
      validate_params({"gamma_er", "beta", "Eer"}, r.parameters, "Eley-Rideal");
      validate_params({}, r.kc, "Eley-Rideal equilibrium constants");
      break;
    case ReactionType::LangmuirHinshelwood:
      validate_params({"C_lh", "beta", "Elh"}, r.parameters,
                      "Langmuir-Hinshelwood");
      validate_params({}, r.kc, "Langmuir-Hinshelwood equilibrium constants");
      break;
    }

    // ========== SURFACE SPECIES COLLECTION ==========
    // Lambda to collect surface species (exclude bare surface sites)
    auto collect_surface = [&](const std::map<std::string, int> &side) {
      for (const auto &[sp, coeff] : side) {
        // Only collect occupied surface sites, not bare "(s)" sites
        if (sp != "(s)" && sp.size() > 3 && sp.substr(sp.size() - 3) == "(s)")
          all_surface_species.insert(sp);
      }
    };
    collect_surface(r.reactants);
    collect_surface(r.products);

    // Track surface species involved in adsorption reactions
    if (r.type == ReactionType::Adsorption ||
        r.type == ReactionType::AdsorptionDesorption) {
      auto add_ads = [&](const std::map<std::string, int> &side) {
        for (const auto &[sp, coeff] : side) {
          // Only collect occupied surface sites for adsorption tracking
          if (sp != "(s)" && sp.size() > 3 && sp.substr(sp.size() - 3) == "(s)")
            adsorption_species.insert(sp);
        }
      };
      add_ads(r.reactants);
      add_ads(r.products);
    }

    // ========== DUPLICATE REACTION DETECTION ==========
    // Check if reaction is reversible based on arrow notation
    bool reversible = r.formula.find("<->") != std::string::npos;
    // Create canonical reaction signature
    std::string forward = make_key(r.reactants) + "->" + make_key(r.products);

    // Check for duplicate forward reaction
    // std::set::insert() returns pair<iterator, bool> where .second indicates
    // success
    if (!seen_signatures.insert(forward).second) {
      throw std::runtime_error("Duplicate reaction after normalization: " +
                               forward);
    }

    // Handle reversible reactions
    if (reversible) {
      std::string reverse = make_key(r.products) + "->" + make_key(r.reactants);
      // Check for duplicate reverse reaction (e.g., same reaction entered in
      // opposite orientation)
      if (!seen_signatures.insert(reverse).second) {
        throw std::runtime_error("Duplicate reaction after normalization: " +
                                 reverse);
      }
      // Only adsorption/desorption mechanisms participate in pairing checks.
      // Eley-Rideal and Langmuir-Hinshelwood reactions are always reversible,
      // but they do not represent adsorption pathways and are therefore
      // excluded from adsorption/desorption bookkeeping.
      if (r.type == ReactionType::AdsorptionDesorption) {
        adsorption_sigs.insert(forward);
        desorption_sigs.insert(reverse);
      }
    } else {
      // Track irreversible adsorption/desorption reactions for pairing
      // validation
      if (r.type == ReactionType::Adsorption)
        adsorption_sigs.insert(forward);
      else if (r.type == ReactionType::Desorption)
        desorption_sigs.insert(forward);
    }

    // ========== MECHANISM-SPECIFIC STOICHIOMETRY VALIDATION ==========
    // Count gas-phase species and retrieve their total stoichiometric
    // coefficients on each side. For all reaction forms we require exact
    // counts, so coefficients greater than one are considered invalid.
    auto gas_info = [](const std::map<std::string, int> &side) {
      int count = 0;
      int coeff = 0;
      for (const auto &[sp, c] : side) {
        // Gas species are those without the "(s)" surface marker.
        if (!(sp == "(s)" ||
              (sp.size() > 3 && sp.substr(sp.size() - 3) == "(s)"))) {
          ++count;
          coeff = c; // Only one gas species allowed per side
        }
      }
      return std::pair<int, int>{count, coeff};
    };

    // Helper to count bare surface sites with species name "(s)".
    auto bare_sites = [](const std::map<std::string, int> &side) {
      auto it = side.find("(s)");
      return it == side.end() ? 0 : it->second;
    };

    // Extract counts for reactants and products.
    auto [gas_react, gas_react_coeff] = gas_info(r.reactants);
    auto [gas_prod, gas_prod_coeff] = gas_info(r.products);
    int sites_react = compute_sites(r.reactants);
    int sites_prod = compute_sites(r.products);
    int bare_react = bare_sites(r.reactants);
    int bare_prod = bare_sites(r.products);

    // Validate stoichiometric patterns for each mechanism type.
    switch (r.type) {
    case ReactionType::Adsorption:
      // Expect: A + (s) -> A(s)
      if (!(gas_react == 1 && gas_react_coeff == 1 && bare_react == 1 &&
            sites_react == bare_react && // only bare sites on reactant side
            gas_prod == 0 && bare_prod == 0 && sites_prod == 1 &&
            r.products.size() == 1)) {
        throw std::runtime_error("Adsorption must be A + (s) -> A(s): " +
                                 r.formula);
      }
      break;

    case ReactionType::Desorption:
      // Expect: A(s) -> A + (s)
      if (!(gas_react == 0 && bare_react == 0 && sites_react == 1 &&
            r.reactants.size() == 1 && gas_prod == 1 && gas_prod_coeff == 1 &&
            bare_prod == 1 && sites_prod == 1 && r.products.size() == 2)) {
        throw std::runtime_error("Desorption must be A(s) -> A + (s): " +
                                 r.formula);
      }
      break;

    case ReactionType::AdsorptionDesorption:
      // Expect: A + (s) <-> A(s)
      if (!(gas_react == 1 && gas_react_coeff == 1 && bare_react == 1 &&
            sites_react == bare_react && // reactant side has only bare sites
            gas_prod == 0 && bare_prod == 0 && sites_prod == 1 &&
            r.products.size() == 1)) {
        throw std::runtime_error(
            "Adsorption/Desorption must be A + (s) <-> A(s): " + r.formula);
      }
      break;

    case ReactionType::EleyRideal: {
      // Expect either:
      // 1) A + B(s) <-> AB + (s)
      // 2) A + (s) <-> B + C(s)
      //
      // The second form handles dissociative collisions of a gas species with
      // a bare surface site. An example is:
      //   N2 + (s) <-> N + N(s)
      bool standard_er =
          (gas_react == 1 && gas_react_coeff == 1 && bare_react == 0 &&
           sites_react == 1 && gas_prod == 1 && gas_prod_coeff == 1 &&
           bare_prod == 1 && sites_prod == 1);
      bool bare_site_er =
          (gas_react == 1 && gas_react_coeff == 1 && bare_react == 1 &&
           sites_react == 1 && gas_prod == 1 && gas_prod_coeff == 1 &&
           bare_prod == 0 && sites_prod == 1);
      if (!(r.reactants.size() == 2 && r.products.size() == 2 &&
            (standard_er || bare_site_er))) {
        throw std::runtime_error("Eley-Rideal must be A + B(s) <-> AB + (s) or "
                                 "A + (s) <-> B + C(s): " +
                                 r.formula);
      }
      break;
    }

    case ReactionType::LangmuirHinshelwood:
      // Expect: A(s) + B(s) <-> AB + 2(s)
      if (!(gas_react == 0 && bare_react == 0 && sites_react == 2 &&
            gas_prod == 1 && gas_prod_coeff == 1 && bare_prod == 2 &&
            sites_prod == 2 && r.products.size() == 2)) {
        throw std::runtime_error(
            "Langmuir-Hinshelwood must be A(s) + B(s) <-> AB + 2(s): " +
            r.formula);
      }
      break;
    }
  }

  // ========== POST-PROCESSING VALIDATION ==========
  // Ensure every adsorption reaction has a corresponding desorption reaction
  for (const auto &sig : adsorption_sigs) {
    // Create reverse signature by swapping reactants and products
    std::string rev =
        sig.substr(sig.find("->") + 2) + "->" + sig.substr(0, sig.find("->"));
    if (!desorption_sigs.count(rev)) {
      throw std::runtime_error("Missing desorption reaction for " + sig);
    }
  }

  // Ensure every desorption reaction has a corresponding adsorption reaction
  for (const auto &sig : desorption_sigs) {
    // Create reverse signature by swapping reactants and products
    std::string rev =
        sig.substr(sig.find("->") + 2) + "->" + sig.substr(0, sig.find("->"));
    if (!adsorption_sigs.count(rev)) {
      throw std::runtime_error("Missing adsorption reaction for " + sig);
    }
  }

  // Ensure every surface species has an adsorption pathway
  for (const auto &sp : all_surface_species) {
    if (!adsorption_species.count(sp))
      throw std::runtime_error(
          "Missing adsorption reaction for surface species: " + sp);
  }
}

} // namespace details
} // namespace gasp2::catalysis::finite_rate
