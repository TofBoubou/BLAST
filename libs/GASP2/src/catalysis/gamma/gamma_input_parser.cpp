#include "gasp2/catalysis/gamma/gamma_input_parser.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <optional>
#include <locale>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace gasp2::catalysis::gamma {
namespace details {

// Locale-independent double parsing (enforces '.' as decimal separator)
static double parse_double_c_locale(const std::string& s) {
  std::istringstream iss(s);
  iss.imbue(std::locale::classic());
  double v = 0.0;
  iss >> v;
  if (!iss) {
    throw std::runtime_error("Invalid numeric value: " + s);
  }
  return v;
}

// Fully parse the XML input file and extract the list of reactions and their
// associated parameters. This routine performs light validation and preserves
// the order of the species as given by the caller.
[[nodiscard]] ParsedInput
read_gamma_input_file(const std::filesystem::path &input_filename,
                      const std::vector<std::string> &species_order) {
  // Normalize cout formatting for the entire parsing routine (debug readability)
  auto cout_flags_guard = std::cout.flags();
  auto cout_prec_guard = std::cout.precision();
  std::cout.setf(static_cast<std::ios::fmtflags>(0), std::ios::floatfield);
  std::cout.precision(6);
  if (!std::filesystem::exists(input_filename)) {
    throw std::runtime_error("Input file does not exist: " +
                             input_filename.string());
  }

  // Read the entire XML input into a string to simplify parsing with regex
  std::ifstream in(input_filename);
  // Read entire file into string using iterator-based constructor
  // First iterator reads from file, second is end-of-stream sentinel
  std::string xml((std::istreambuf_iterator<char>(in)),
                  std::istreambuf_iterator<char>());

  std::string reactions_block;
  // Find start position of <reactions> tag (returns npos if not found)
  std::size_t reactions_start = xml.find("<reactions");
  if (reactions_start != std::string::npos) {
    std::size_t reactions_end = xml.find("</reactions>", reactions_start);
    if (reactions_end == std::string::npos) {
      throw std::runtime_error("Missing </reactions> closing tag");
    }
    if (xml.find("<reactions", reactions_start + 1) != std::string::npos) {
      throw std::runtime_error("Multiple <reactions> blocks found");
    }
    // Validate that ALL <reaction> tags are inside the <reactions> block
    std::regex reaction_tag_re{R"(<\s*reaction\b)", std::regex::icase};
    // Use regex iterator to find every <reaction> tag in the entire XML
    for (auto it =
             std::sregex_iterator(xml.begin(), xml.end(), reaction_tag_re);
         it != std::sregex_iterator(); ++it) {
      // Get character position of this <reaction> tag
      std::size_t pos = static_cast<std::size_t>(it->position());
      // Error if any <reaction> tag is outside the valid <reactions> block
      if (pos < reactions_start || pos > reactions_end) {
        throw std::runtime_error("<reaction> outside <reactions> block");
      }
    }
    // Extract just the <reactions> block content for further processing
    // substr(start_pos, length) where length = end_pos - start_pos
    reactions_block =
        xml.substr(reactions_start, reactions_end - reactions_start);
  } else {
    std::regex reaction_tag_re{R"(<\s*reaction\b)", std::regex::icase};
    if (std::regex_search(xml, reaction_tag_re)) {
      throw std::runtime_error("<reaction> outside <reactions> block");
    }
  }

  ParsedInput parsed;

  bool bose_found = false;
  bool super_found = false;
  bool rini_found = false;
  bool other_found = false;
  // Optional value: either contains a gamma_w parameter from XML, or is empty
  // Used for Rini model reactions that may have a global gamma_w setting
  std::optional<double> rini_gamma_w;

  // Detect Bose reaction block - a special self-closing reaction tag
  // Pattern matches: <reaction type="GammaBose" gamma="value" p2="value" />
  // Captures gamma and p2 parameter values in parentheses groups
  std::regex bose_re{
      R"(<reaction\s+type=\"GammaBose\"[^>]*gamma=\"([^\"]+)\"[^>]*p2=\"([^\"]+)\"[^>]*/>)",
      std::regex::icase};
  // Create iterator range to find all matches of the Bose pattern
  auto bose_begin = std::sregex_iterator(reactions_block.begin(),
                                         reactions_block.end(), bose_re);
  // Default constructor creates "end iterator" - signals no more matches
  auto bose_end = std::sregex_iterator();
  if (bose_begin != bose_end) {
    bose_found = true;
    // Ensure only one GammaBose reaction exists
    // std::next advances iterator to check if second match exists
    if (std::next(bose_begin) != bose_end) {
      throw std::runtime_error("Multiple GammaBose blocks found");
    }
    // Extract parameter values from regex capture groups and convert to double
    // (*bose_begin)[1] gets first capture group, .str() converts to string,
    // std::stod converts to double
    double gamma = parse_double_c_locale((*bose_begin)[1].str());
    double p2 = parse_double_c_locale((*bose_begin)[2].str());
    // GammaBose model requires these specific chemical species to be available
    // It implements: O + O -> O2 and CO + O -> CO2 reactions
    std::set<std::string> required{"O", "O2", "CO", "CO2"};
    for (const auto &sp : required) {
      // Check if each required species exists in the user-provided species list
      if (std::find(species_order.begin(), species_order.end(), sp) ==
          species_order.end()) {
        throw std::runtime_error("Species not found in order: " + sp);
      }
    }
    ReactionInput::BoseParams params{gamma, p2};
    ReactionInput r1;
    r1.formula = "O + O -> O2";
    r1.type = CatalysisModel::GammaBose;
    r1.bose = params;
    r1.reactants["O"] = 2;
    r1.products["O2"] = 1;
    parsed.reactions.push_back(r1);
    ReactionInput r2;
    r2.formula = "CO + O -> CO2";
    r2.type = CatalysisModel::GammaBose;
    r2.bose = params;
    r2.reactants["CO"] = 1;
    r2.reactants["O"] = 1;
    r2.products["CO2"] = 1;
    parsed.reactions.push_back(r2);
  }
  // Create a regex match result
  std::smatch fm;
  // Optional flag to choose first-order theory
  std::regex first_order_re{R"(<first_order>([^<]*)</first_order>)",
                            std::regex::icase};
  if (std::regex_search(xml, fm, first_order_re)) {
    std::string val = fm[1].str();
    std::transform(val.begin(), val.end(), val.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (val == "true" || val == "1") {
      parsed.first_order = true;
    } else if (val == "false" || val == "0") {
      parsed.first_order = false;
    } else {
      throw std::runtime_error("Invalid value in <first_order> block");
    }
  }

  // Regex to capture individual reaction blocks with attributes
  // [\s\S]*? matches any character (including newlines) non-greedily
  // \s = whitespace, \S = non-whitespace, so [\s\S] = any character
  // *? = zero or more (non-greedy), stops at first </reaction> it encounters
  // Without ?, would greedily match to the LAST </reaction> in the entire text
  std::regex reaction_re{
      R"(<reaction\s+type=\"([^\"]+)\"\s+formula=\"([^\"]+)\"\s*>([\s\S]*?)</reaction>)",
      std::regex::icase};
  // Regex to find the optional list of recombination probabilities within
  // <gammas>
  std::regex gamma_re{R"(<gammas>([^<]*)</gammas>)", std::regex::icase};
  std::regex gamma_w_re{R"(<gamma_w>([^<]*)</gamma_w>)", std::regex::icase};
  // Optional flag to activate limiting-flux adjustments
  std::regex limiting_fluxes_re{R"(<limiting_fluxes>([^<]*)</limiting_fluxes>)",
                                std::regex::icase};
  // Regex to find optional GammaT parameters A:E:Tmin:Tmax within <gammaT>
  std::regex gammaT_re{R"(<gammaT>([^<]*)</gammaT>)", std::regex::icase};
  // Regex to parse species tokens with optional stoichiometric coefficients and
  // optional ionic charge. Pattern: [optional digits][species][+/-].
  // Examples:
  //   "O"   -> coeff="", species="O"
  //   "2 N2"-> coeff="2", species="N2"
  //   "NO+" -> coeff="", species="NO+"
  //   "e-"  -> coeff="", species="e-"
  // The final non-capturing group consumes either a '+' separator or the end
  // of the string to allow ion species containing '+' characters.
  std::regex species_re{R"(\s*(\d*)\s*([A-Za-z0-9_]+(?:[+-]+)?)\s*(?:\+|$))"};

  std::set<std::string> seen_reactions;
  // Detect a single global gamma_w for the Rini model. This global parameter
  // is applied to every Rini_model reaction encountered below.
  auto gw_begin = std::sregex_iterator(reactions_block.begin(),
                                       reactions_block.end(), gamma_w_re);
  auto gw_end = std::sregex_iterator();
  if (gw_begin != gw_end) {
    if (std::next(gw_begin) != gw_end) {
      throw std::runtime_error("Multiple gamma_w entries found");
    }
    rini_gamma_w = parse_double_c_locale((*gw_begin)[1].str());
    // Remove <gamma_w> tags from XML after extracting value. This prevents the
    // tag from being misinterpreted as part of an individual <reaction> block
    // in the subsequent parsing loop.
    reactions_block =
        std::regex_replace(reactions_block, gamma_w_re, std::string());
  }

  // Extract optional <limiting_fluxes> flag before processing reactions.
  if (std::regex_search(reactions_block, fm, limiting_fluxes_re)) {
    std::string val = fm[1].str();
    std::transform(val.begin(), val.end(), val.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (val == "true" || val == "1") {
      parsed.limiting_fluxes = true;
    } else if (val == "false" || val == "0") {
      parsed.limiting_fluxes = false;
    } else {
      throw std::runtime_error("Invalid value in <limiting_fluxes> block");
    }
    reactions_block =
        std::regex_replace(reactions_block, limiting_fluxes_re, std::string());
  }
  if (parsed.limiting_fluxes && !parsed.first_order) {
    throw std::runtime_error(
        "<limiting_fluxes> requires <first_order>true</first_order>");
  }

  auto begin = std::sregex_iterator(reactions_block.begin(),
                                    reactions_block.end(), reaction_re);
  auto end = std::sregex_iterator();
  // Iterate over all found reactions
  for (auto it = begin; it != end; ++it) {
    if (bose_found) {
      throw std::runtime_error(
          "GammaBose cannot be combined with other reactions");
    }
    // Create a ReactionInput object
    ReactionInput r;
    // Extract the reaction type attribute (e.g., "GammaGiven",
    // "SuperCatalytic")
    std::string type_attr = (*it)[1].str();
    std::string type_lc = type_attr;
    // Convert to lowercase for case-insensitive comparison
    std::transform(type_lc.begin(), type_lc.end(), type_lc.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    // Map reaction type strings to enums and enforce compatibility rules
    if (type_lc == "gammagiven") {
      r.type = CatalysisModel::GammaGiven;
      other_found = true; // Track "other" model types for validation
    } else if (type_lc == "gammat") {
      r.type = CatalysisModel::GammaT;
      other_found = true;
    } else if (type_lc == "gammaconsistent") {
      r.type = CatalysisModel::GammaConsistent;
      other_found = true;
      // GammaConsistent REQUIRES first-order theory to be enabled
      if (!parsed.first_order) {
        throw std::runtime_error(
            "GammaConsistent requires <first_order>true</first_order>");
      }
    } else if (type_lc == "supercatalytic") {
      r.type = CatalysisModel::SuperCatalytic;
      super_found = true; // Track SuperCatalytic separately
      // SuperCatalytic CANNOT be used with first-order theory
      if (parsed.first_order) {
        throw std::runtime_error("SuperCatalytic cannot be used with "
                                 "<first_order>true</first_order>");
      }
    } else if (type_lc == "rini_model") {
      r.type = CatalysisModel::RiniModel;
      rini_found = true; // Track Rini model separately
      // Rini model CANNOT be used with first-order theory
      if (parsed.first_order) {
        throw std::runtime_error(
            "Rini_model cannot be used with <first_order>true</first_order>");
      }
    } else {
      // Unknown reaction type - fail with helpful error message
      throw std::runtime_error("Unknown reaction type: " + type_attr);
    }
    // Extract reaction formula
    r.formula = (*it)[2].str();

    // ========== ARROW VALIDATION ==========
    // Only irreversible reactions using the "->" arrow are supported for the
    // gamma model. Reject reversible arrows such as "<-", "<->" or the "="
    // symbol sometimes used for equilibrium reactions.
    if (r.formula.find("<->") != std::string::npos ||
        r.formula.find("<-") != std::string::npos ||
        r.formula.find('=') != std::string::npos) {
      throw std::runtime_error(
          "Only irreversible reactions with '->' are allowed: " + r.formula);
    }
    auto arrow = r.formula.find("->");
    if (arrow == std::string::npos) {
      throw std::runtime_error("Missing '->' in reaction: " + r.formula);
    }
    std::string block = (*it)[3].str();
    // Create a match variable for the gammas
    std::smatch gm;
    // Extract species-specific recombination probabilities from <gammas> tag
    // Expected XML format: <gammas>N:0.5, O:0.3; CO:0.8</gammas>
    if (std::regex_search(block, gm, gamma_re)) {
      std::string gammas_str =
          gm[1].str(); // Get content between <gammas>...</gammas>
      std::cout << "DEBUG XML Parsing - Found gammas string: '" << gammas_str << "'" << std::endl;
      // Normalize separators: replace commas with semicolons for consistent
      // parsing
      std::replace(gammas_str.begin(), gammas_str.end(), ',', ';');
      std::cout << "DEBUG XML Parsing - After comma replacement: '" << gammas_str << "'" << std::endl;
      std::stringstream ss(gammas_str);
      std::string token;
      // Split on semicolons to get individual species:value pairs
      while (std::getline(ss, token, ';')) {
        // Trim whitespace from both ends of token
        // Remove leading whitespace
        token.erase(0, token.find_first_not_of(" \t\n\r"));
        // Remove trailing whitespace
        token.erase(token.find_last_not_of(" \t\n\r") + 1);
        if (token.empty())
          continue; // Skip empty tokens

        // Look for colon separator between species name and value
        std::cout << "DEBUG XML Parsing - Processing token: '" << token << "'" << std::endl;
        auto pos = token.find(':');
        if (pos != std::string::npos) {
          std::string species = token.substr(0, pos);      // Before ':'
          std::string value_str = token.substr(pos + 1);   // After ':'
          double value = parse_double_c_locale(value_str);
          // Force sane numeric formatting for debug (avoid global stream state side effects)
          auto old_flags = std::cout.flags();
          auto old_prec = std::cout.precision();
          std::cout.setf(static_cast<std::ios::fmtflags>(0), std::ios::floatfield);
          std::cout.precision(6);
          std::cout << "DEBUG XML Parsing - Parsed species: '" << species << "', value: " << value << std::endl;
          std::cout.flags(old_flags);
          std::cout.precision(old_prec);
          // Create gammas map if it doesn't exist yet (lazy initialization)
          if (!r.gammas)
            r.gammas.emplace();
          // Store species -> gamma_value mapping
          (*r.gammas)[species] = value;
          std::cout << "DEBUG XML Parsing - Stored in map: " << species << " -> " << value << std::endl;
          std::cout << "DEBUG XML Parsing - Map size after storing: " << r.gammas->size() << std::endl;
        } else {
          std::cout << "DEBUG XML Parsing - No colon found in token: '" << token << "'" << std::endl;
        }
        // Note: Silently ignore tokens without ':' (malformed entries)
      }
    }

    // Extract temperature-dependent gamma parameters from <gammaT> tag
    // Expected XML format: <gammaT>species:A:E:Tmin:Tmax,
    // species2:A2:E2:Tmin2:Tmax2</gammaT> Where A = pre-exponential factor, E =
    // activation energy, Tmin/Tmax = temperature range
    if (std::regex_search(block, gm, gammaT_re)) {
      std::string params_str =
          gm[1].str(); // Extract content between <gammaT>...</gammaT>
      // Normalize separators: replace commas with semicolons for consistent
      // parsing
      std::replace(params_str.begin(), params_str.end(), ',', ';');
      std::stringstream ss(params_str);
      std::string token;
      // Split on semicolons to get individual species parameter sets
      while (std::getline(ss, token, ';')) {
        // Trim whitespace from both ends of token
        token.erase(0, token.find_first_not_of(" \t\n\r")); // Remove leading
        token.erase(token.find_last_not_of(" \t\n\r") + 1); // Remove trailing
        if (token.empty())
          continue; // Skip empty tokens

        // Split token by colons to extract individual parameters
        // Expected format: "species:A:E:Tmin:Tmax" (5 parts minimum)
        std::vector<std::string> parts;
        std::stringstream tss(token);
        std::string part;
        while (std::getline(tss, part, ':'))
          parts.push_back(part);

        // Validate we have all required parameters: species + 4 numeric values
        if (parts.size() >= 5) {
          // Create temperature-dependent parameter structure
          ReactionInput::GammaTParams params;
          params.A = parse_double_c_locale(parts[1]);    // Pre-exponential factor
          params.E = parse_double_c_locale(parts[2]);    // Activation energy
          params.tmin = parse_double_c_locale(parts[3]); // Minimum temperature
          params.tmax = parse_double_c_locale(parts[4]); // Maximum temperature
          // Initialize gammaT map if this is the first temperature-dependent
          // parameter
          if (!r.gammaT)
            r.gammaT.emplace(); // Create empty map in-place for efficiency
          // Store species-specific temperature parameters: species_name ->
          // GammaTParams
          (*r.gammaT)[parts[0]] = params; // parts[0] is species name
        }
        // Note: Silently ignore tokens with insufficient parameters (< 5 parts)
      }
    }
    // Assign the global Rini_model gamma_w value to each parsed reaction.  All
    // Rini reactions share the same recombination probability.
    if (r.type == CatalysisModel::RiniModel && rini_gamma_w) {
      r.gamma_w = *rini_gamma_w;
    }

    // Parse reaction formula to extract stoichiometric coefficients for
    // reactants and products. Expected format: "2N + O -> N2O" or
    // "NO+ + e- -> N + O".
    // Lambda function to parse one side of the chemical equation (reactants or
    // products). Uses regex iteration to correctly handle species names that
    // contain '+' or '-' characters for ionic charges.
    auto parse_side = [&](const std::string &side,
                          std::map<std::string, int> &out) {
      auto begin = std::sregex_iterator(side.begin(), side.end(), species_re);
      auto end = std::sregex_iterator();
      for (auto it = begin; it != end; ++it) {
        int coeff = (*it)[1].str().empty() ? 1 : std::stoi((*it)[1].str());
        std::string sp = (*it)[2].str();
        out[sp] += coeff;
      }
    };
    // Parse left side of arrow (reactants)
    parse_side(r.formula.substr(0, arrow), r.reactants);
    // Parse right side of arrow (products)
    parse_side(r.formula.substr(arrow + 2), r.products);

    // ========== STRUCTURAL VALIDATION ==========
    // Ensure reactions define a reasonable chemical transformation:
    // - at least one reactant (no empty left side)
    // - no more than two reactants (model limitation)
    // - exactly one product (model assumption)
    if (r.reactants.empty()) {
      throw std::runtime_error("Reaction must have at least one reactant: " +
                               r.formula);
    }
    if (r.reactants.size() > 2) {
      throw std::runtime_error("Reaction cannot have more than 2 reactants: " +
                               r.formula);
    }
    if (r.products.size() != 1) {
      throw std::runtime_error("Reaction must have exactly 1 product: " +
                               r.formula);
    }

    // Determine whether this reaction mixes different reactant species. A
    // size greater than one indicates a heterogeneous reaction, while a single
    // unique reactant implies a homogeneous process.
    r.heterogeneous = r.reactants.size() > 1;

    // Normalize stoichiometric coefficients to simplify reaction representation
    // Goal: Reduce coefficients by their greatest common divisor when possible
    // Example transformations: "2N + 2N -> 2N2" becomes "N + N -> N2"
    //                         "4CO + 2O2 -> 4CO2" becomes "2CO + O2 -> 2CO2"
    // Get the stoichiometric coefficient of the single product
    // r.products.begin()->second accesses the coefficient value of the first
    // (only) product
    int prod_coeff = r.products.begin()->second;
    // Only normalize if product coefficient > 1 (avoid division by 1)
    if (prod_coeff > 1) {
      // Check if ALL reactant coefficients are divisible by the product
      // coefficient This ensures the reaction remains chemically balanced
      // after division
      bool divisible = true;
      for (const auto &rc : r.reactants) { // rc = reactant coefficient pair
        // rc.second is the stoichiometric coefficient of this reactant
        if (rc.second % prod_coeff != 0) { // Check if evenly divisible
          divisible = false; // Found a reactant that can't be divided evenly
          break;             // Stop checking - normalization not possible
        }
      }
      // Perform the normalization only if all coefficients are divisible
      if (divisible) {
        // Divide all reactant coefficients by the product coefficient
        for (auto &rc : r.reactants)
          rc.second /= prod_coeff; // Example: 4CO becomes 2CO when prod_coeff=2
        // Divide all product coefficients by the product coefficient
        for (auto &pc : r.products)
          pc.second /= prod_coeff; // Example: 2N2 becomes 1N2 (written as N2)
        // Result: "4N + 2N2 -> 4N2" becomes "2N + N2 -> 2N2" if prod_coeff=2
      }
      // If not divisible, keep original coefficients to maintain chemical
      // balance
      // After normalization, update product coefficient for later checks
      prod_coeff = r.products.begin()->second;
    }
    // If prod_coeff !=1, throw an error for gammagiven, gammaT and
    // gammaConsistent models
    if (prod_coeff != 1 && (r.type == CatalysisModel::GammaGiven ||
                            r.type == CatalysisModel::GammaT ||
                            r.type == CatalysisModel::GammaConsistent)) {
      throw std::runtime_error(
          "Error: For GammaGiven,GammaT,GammaConsistent, the formulation is "
          "written for "
          "a product with unitary stoichiometric coefficient " +
          r.formula);
    }
    // Create a canonical string representation of a reaction
    // side (reactants or products) This lambda generates a deterministic key
    // from a map of species->coefficient pairs The key format concatenates
    // coefficient+species for each entry, joined by '+'
    auto make_key = [](const std::map<std::string, int> &side) {
      std::string key;
      // Iterate through the map (automatically sorted by species name due to
      // std::map) This ensures consistent ordering: "CO + N" and "N + CO" both
      // become "1CO+1N"
      for (auto it = side.begin(); it != side.end(); ++it) {
        if (it !=
            side.begin()) // Add '+' separator between species (skip for first)
          key += "+";
        // Concatenate coefficient and species name: "2N2" for coefficient=2,
        // species="N2" std::to_string converts integer coefficient to string
        key += std::to_string(it->second) +
               it->first; // it->second=coeff, it->first=species
      }
      return key; // Example results: "1CO+1O" or "2N+1O2"
    };

    // Generate unique signature for this reaction to detect duplicates
    // Format: "reactant_key->product_key" (e.g., "1CO+1O->1CO2")
    // The arrow "->" separates reactants from products in the signature
    std::string signature = make_key(r.reactants) + "->" + make_key(r.products);

    // Attempt to insert signature into seen_reactions set
    // std::set::insert() returns pair<iterator, bool> where bool indicates
    // success If insertion fails (.second == false), the signature already
    // exists = duplicate reaction
    if (!seen_reactions.insert(signature).second) {
      // Throw error with the normalized signature to help user identify the
      // duplicate
      throw std::runtime_error("Duplicate reaction after normalization: " +
                               signature);
    }

    // Move the completed reaction object into the parsed results vector
    // std::move transfers ownership, avoiding unnecessary copying of the
    // ReactionInput object
    parsed.reactions.push_back(std::move(r));
  }

  // Debug: summarize parsed gamma values before returning
  std::cout << "DEBUG XML Parsing - FINAL SUMMARY of gammas:" << std::endl;
  for (size_t i = 0; i < parsed.reactions.size(); ++i) {
    const auto &rx = parsed.reactions[i];
    std::cout << "  Reaction " << i << ": " << rx.formula << std::endl;
    if (rx.gammas) {
      std::cout << "    gammas size: " << rx.gammas->size() << std::endl;
      for (const auto &kv : *rx.gammas) {
        std::cout << "    " << kv.first << " -> " << kv.second << std::endl;
      }
    } else {
      std::cout << "    gammas: <none>" << std::endl;
    }
  }

  // Final validation: Enforce mutual exclusivity rules for catalysis models
  // These models have incompatible mathematical formulations and cannot coexist
  // SuperCatalytic: Uses infinite reaction rates (instantaneous equilibrium)
  // Rini_model: Uses specific kinetic equations with gamma_w parameter
  // Other models: Standard gamma-based approaches (GammaGiven, GammaT,
  // GammaConsistent)
  if ((super_found && other_found) || (super_found && rini_found) ||
      (rini_found && other_found)) {
    throw std::runtime_error(
        "SuperCatalytic or Rini_model cannot be combined with other reactions");
  }

  // Validate GammaBose model requirements
  // GammaBose is a specialized model for oxygen/carbon chemistry that REQUIRES
  // first-order theory First-order theory simplifies the mathematical treatment
  // of surface reactions
  if (bose_found && (!parsed.first_order)) {
    throw std::runtime_error(
        "GammaBose requires <first_order>true</first_order>");
  }

  // Restore cout formatting
  std::cout.flags(cout_flags_guard);
  std::cout.precision(cout_prec_guard);
  return parsed;
}

} // namespace details
} // namespace gasp2::catalysis::gamma
