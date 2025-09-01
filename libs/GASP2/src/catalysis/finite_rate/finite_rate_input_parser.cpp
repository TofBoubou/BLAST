#include "gasp2/catalysis/finite_rate/finite_rate_input_parser.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <locale>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace gasp2::catalysis::finite_rate {
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

/**
 * @brief Parse comma/semicolon-separated key:value pairs into a map.
 * @param text Input string with format "key1:value1,key2:value2" or using semicolons
 * @return Map of parameter names to their numeric values
 */
static std::map<std::string, double> parse_params(std::string text) {
  std::map<std::string, double> params;
  
  // Normalize delimiters: replace all commas with semicolons for consistent parsing
  std::replace(text.begin(), text.end(), ',', ';');
  
  // Use stringstream to process the normalized string
  std::stringstream ss(text);
  std::string token;
  
  // Split by semicolon delimiter and process each token
  while (std::getline(ss, token, ';')) {
    // Trim leading whitespace (space, tab, newline, carriage return)
    token.erase(0, token.find_first_not_of(" \t\n\r"));
    // Trim trailing whitespace
    token.erase(token.find_last_not_of(" \t\n\r") + 1);
    
    // Skip empty tokens (e.g., from consecutive delimiters like "a:1;;b:2")
    if (token.empty())
      continue;
    
    // Look for the colon separator between key and value
    auto pos = token.find(':');
    if (pos != std::string::npos) {
      // Extract key (everything before the colon)
      std::string key = token.substr(0, pos);
      // Extract and convert value (everything after the colon) to double
      double value = parse_double_c_locale(token.substr(pos + 1));
      // Store the key-value pair in the map
      params[key] = value;
    }
    // Silently ignore tokens without colon separator
  }
  return params;
}

/**
 * @brief Build a canonical string representation of a reaction side.
 * @param side Map of species names to their stoichiometric coefficients
 * @return Formatted string with gas-phase species first, then surface species
 */
static std::string side_to_string(const std::map<std::string, int> &side) {
  // Separate species into gas-phase and surface-phase categories
  std::vector<std::pair<std::string, int>> gas;
  std::vector<std::pair<std::string, int>> surface;
  // Scan the map (kv is a key-value pair)
  for (const auto &kv : side) {
    // Check if species is surface-bound: either "(s)" or ends with "(s)"
    if (kv.first == "(s)" ||
        (kv.first.size() > 3 && kv.first.substr(kv.first.size() - 3) == "(s)"))
      surface.push_back(kv);
    else
      gas.push_back(kv);
  }
  
  // Sort both categories alphabetically by species name for consistent output
  auto cmp = [](const auto &a, const auto &b) { return a.first < b.first; };
  std::sort(gas.begin(), gas.end(), cmp);
  std::sort(surface.begin(), surface.end(), cmp);
  
  std::string result;
  // Helper lambda to append species to the result string
  auto append = [&](const std::vector<std::pair<std::string, int>> &vec) {
    for (const auto &kv : vec) {
      // Add '+' separator between species (but not before the first one)
      if (!result.empty())
        result += "+";
      // Only include stoichiometric coefficient if it's not 1
      if (kv.second != 1)
        result += std::to_string(kv.second);
      // Append the species name
      result += kv.first;
    }
  };
  
  // Append gas-phase species first, then surface species
  append(gas);
  append(surface);
  return result;
  
  // Example: Input map {"O2": 2, "H2": 1, "CO(s)": 1, "O(s)": 2} 
  // produces output string "H2+2O2+CO(s)+2O(s)"
  // (gas species sorted alphabetically first, then surface species)
}

/**
 * @brief Parse finite rate reaction input file and extract reaction data.
 * @param input_filename Path to the XML input file containing reaction definitions
 * @return ParsedInput structure containing site density and reaction definitions
 */
ParsedInput read_finite_rate_input_file(
    const std::filesystem::path &input_filename) {
  // ========== FILE VALIDATION ==========
  // Check if the input file exists before attempting to read it
  if (!std::filesystem::exists(input_filename)) {
    throw std::runtime_error("Input file does not exist: " +
                             input_filename.string());
  }

  // Read entire XML into memory for simpler regex-based parsing
  // Using istreambuf_iterator to read the complete file content at once
  std::ifstream in(input_filename);
  std::string xml((std::istreambuf_iterator<char>(in)),
                  std::istreambuf_iterator<char>());

  // Extract <reactions> block and ensure no stray <reaction> tags exist
  std::string reactions_block;
  std::size_t reactions_start = xml.find("<reactions");
  if (reactions_start != std::string::npos) {
    // Find the corresponding closing tag
    std::size_t reactions_end = xml.find("</reactions>", reactions_start);
    if (reactions_end == std::string::npos)
      throw std::runtime_error("Missing </reactions> closing tag");
    // Ensure there's only one <reactions> block in the file
    if (xml.find("<reactions", reactions_start + 1) != std::string::npos)
      throw std::runtime_error("Multiple <reactions> blocks found");
    // Validate site_density tags are within <reactions> block and ensure only one exists
    std::size_t site_density_first = xml.find("<site_density>");
    if (site_density_first != std::string::npos) {
      // Check if site_density is within the reactions block
      if (site_density_first < reactions_start || site_density_first > reactions_end)
        throw std::runtime_error("<site_density> outside <reactions> block");
      // Check for multiple site_density blocks
      if (xml.find("<site_density>", site_density_first + 1) != std::string::npos)
        throw std::runtime_error("Multiple <site_density> blocks found");
    }
    // Regex to find any <reaction> tags and validate they're within <reactions> block:
    // - <\s* - opening bracket with optional whitespace  
    // - reaction\b - literal "reaction" followed by word boundary
    // This ensures we only match actual <reaction> tags, not partial matches
    std::regex reaction_tag_re{R"(<\s*reaction\b)", std::regex::icase};
    for (auto it =
             std::sregex_iterator(xml.begin(), xml.end(), reaction_tag_re);
         it != std::sregex_iterator(); ++it) {
      std::size_t pos = static_cast<std::size_t>(it->position());
      // Verify each <reaction> tag is within the <reactions> block
      if (pos < reactions_start || pos > reactions_end)
        throw std::runtime_error("<reaction> outside <reactions> block");
    }
    // Extract just the <reactions> block content for further processing
    reactions_block =
        xml.substr(reactions_start, reactions_end - reactions_start);
  } else {
    // No <reactions> block found - this is required for valid input
    throw std::runtime_error("Missing <reactions> block in input file");
  }

  // Initialize the output data structure
  ParsedInput data;

  // ========== SITE DENSITY ==========
  // Regex to capture site density value from XML:
  // - Matches opening tag: <site_density>
  // - Captures group 1: ([^<]+) - any characters except '<' (the density value)
  // - Matches closing tag: </site_density>
  std::regex sd_re{R"(<site_density>([^<]+)</site_density>)",
                   std::regex::icase};
  std::smatch fm;
  if (std::regex_search(reactions_block, fm, sd_re)) {
    // Parse the captured site density value and store it
    data.site_density = parse_double_c_locale(fm[1].str());
    // Remove the site_density tag from the reactions block to avoid reprocessing
    reactions_block = std::regex_replace(reactions_block, sd_re, std::string());
  } else {
    throw std::runtime_error("Missing <site_density> tag");
  }

  // ========== REACTION PARSING SETUP ==========
  // Regex to match complete reaction blocks:
  // - <reaction - opening tag start
  // - \s+type=\"([^\"]+)\" - captures group 1: reaction type attribute value
  // - \s+formula=\"([^\"]+)\" - captures group 2: formula attribute value
  // - \s*> - optional whitespace and closing >
  // - ([\s\S]*?) - captures group 3: reaction content (non-greedy, including newlines)
  // - </reaction> - closing tag
  std::regex reaction_re{
      R"(<reaction\s+type=\"([^\"]+)\"\s+formula=\"([^\"]+)\"\s*>([\s\S]*?)</reaction>)",
      std::regex::icase};
  
  // Regex to extract parameters block content:
  // - <parameters> - opening tag
  // - ([^<]*) - captures group 1: any characters except '<' (parameter content)
  // - </parameters> - closing tag
  std::regex params_re{R"(<parameters>([^<]*)</parameters>)",
                       std::regex::icase};
  
  // Regex to extract equilibrium constant (Kc) block content:
  // - <Kc> - opening tag
  // - ([^<]*) - captures group 1: any characters except '<' (Kc content)
  // - </Kc> - closing tag
  std::regex kc_re{R"(<Kc>([^<]*)</Kc>)", std::regex::icase};
  
  // Regex to parse individual species in reaction formula:
  // - \s* - optional leading whitespace
  // - ([0-9]+(?:\.[0-9]+)?)? - captures group 1: optional stoichiometric coefficient
  //   - [0-9]+ - one or more digits
  //   - (?:\.[0-9]+)? - optional decimal part (non-capturing group)
  // - \s* - optional whitespace after coefficient
  // - ([A-Za-z0-9_()]+(?:[+-]+)?) - captures group 2: species name
  //   - [A-Za-z0-9_()]+ - alphanumeric chars, underscores, parentheses
  //   - (?:[+-]+)? - optional charge notation (e.g., +, -, ++, --)
  // - \s*(?:\+|$) - optional whitespace followed by '+' separator or end of string
  // Allow optional decimal coefficients to detect and reject fractional values.
  std::regex species_re{
      R"(\s*([0-9]+(?:\.[0-9]+)?)?\s*([A-Za-z0-9_()]+(?:[+-]+)?)\s*(?:\+|$))"};

  // Create iterators to process all reaction matches in the reactions block
  auto begin = std::sregex_iterator(reactions_block.begin(),
                                    reactions_block.end(), reaction_re);
  auto end = std::sregex_iterator();

  // Process each reaction found in the XML
  for (auto it = begin; it != end; ++it) {
    ReactionInput r;
    // Extract reaction type attribute and normalize for comparison
    std::string type_attr = (*it)[1].str();
    std::string type_lc = type_attr; //lowercase string
    // Convert to lowercase for case-insensitive comparison
    std::transform(type_lc.begin(), type_lc.end(), type_lc.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    std::string type_clean;
    // Remove hyphens and spaces to handle variations like "Eley-Rideal" vs "EleyRideal"
    std::remove_copy_if(type_lc.begin(), type_lc.end(),
                        std::back_inserter(type_clean),
                        [](unsigned char c) { return c == '-' || c == ' '; });

    // Map cleaned reaction type string to enum value
    if (type_clean == "adsorption") {
      r.type = ReactionType::Adsorption;
    } else if (type_clean == "desorption") {
      r.type = ReactionType::Desorption;
    } else if (type_clean == "adsorption/desorption" ||
               type_clean == "adsorptiondesorption") {
      r.type = ReactionType::AdsorptionDesorption;
    } else if (type_clean == "eleyrideal" || type_clean == "er") {
      r.type = ReactionType::EleyRideal;
    } else if (type_clean == "langmuirhinshelwood" || type_clean == "lh") {
      r.type = ReactionType::LangmuirHinshelwood;
    } else {
      throw std::runtime_error("Unknown reaction type: " + type_attr);
    }

    // Extract formula and reaction block content
    r.formula = (*it)[2].str();
    std::string block = (*it)[3].str();

    // Determine if reaction is reversible by checking arrow type
    bool reversible = false;
    std::size_t arrow = r.formula.find("<->");
    std::size_t arrow_len = 3;
    if (arrow != std::string::npos) {
      reversible = true;
    } else {
      // Look for irreversible arrow
      arrow = r.formula.find("->");
      arrow_len = 2;
    }
    if (arrow == std::string::npos)
      throw std::runtime_error("Missing reaction arrow in: " + r.formula);

    // Validate arrow type matches reaction type constraints
    if ((r.type == ReactionType::Adsorption ||
         r.type == ReactionType::Desorption) &&
        reversible) {
      throw std::runtime_error(
          "Adsorption/Desorption reactions must use '->': " + r.formula);
    }
    if ((r.type == ReactionType::AdsorptionDesorption ||
         r.type == ReactionType::EleyRideal ||
         r.type == ReactionType::LangmuirHinshelwood) &&
        !reversible) {
      throw std::runtime_error("Reaction must use '<->': " + r.formula);
    }

    // Parse kinetic parameters from the reaction block
    std::smatch pm;
    if (std::regex_search(block, pm, params_re)) {
      r.parameters = parse_params(pm[1].str());
    } else {
      throw std::runtime_error("Missing <parameters> in reaction: " +
                               r.formula);
    }
    // For adsorption/desorption reactions, also parse equilibrium constants
    if (r.type == ReactionType::AdsorptionDesorption) {
      if (std::regex_search(block, pm, kc_re)) {
        r.kc = parse_params(pm[1].str());
      } else {
        throw std::runtime_error("Adsorption/Desorption requires <Kc>: " +
                                 r.formula);
      }
    }

    // Lambda function to parse one side of the reaction (reactants or products)
    auto parse_side = [&](const std::string &side,
                          std::map<std::string, int> &out) {
      auto sbegin = std::sregex_iterator(side.begin(), side.end(), species_re);
      auto send = std::sregex_iterator();
      for (auto sit = sbegin; sit != send; ++sit) {
        // ================= STOICHIOMETRIC VALIDATION =================
        // Capture the raw coefficient text, defaulting to one when absent.
        std::string coeff_str = (*sit)[1].str();
        double coeff_val = coeff_str.empty() ? 1.0 : parse_double_c_locale(coeff_str);
        // Reject coefficients below unity (e.g., 0.5O2) to enforce integer
        // stoichiometry. Users should multiply the entire reaction to eliminate
        // fractional values.
        if (coeff_val < 1.0)
          throw std::runtime_error(
              "Stoichiometric coefficients <1 are not allowed: " + r.formula);
        // Coefficients must be integers; fractional values like 1.5 trigger an
        // error to preserve discrete site counts.
        if (std::floor(coeff_val) != coeff_val)
          throw std::runtime_error(
              "Non-integer stoichiometric coefficient in reaction: " +
              r.formula);
        int coeff = static_cast<int>(coeff_val);
        std::string sp = (*sit)[2].str();
        // Add species to the output map (accumulate if already present)
        out[sp] += coeff;
      }
    };
    // Parse reactants (left side of arrow) and products (right side of arrow)
    parse_side(r.formula.substr(0, arrow), r.reactants);
    parse_side(r.formula.substr(arrow + arrow_len), r.products);

    // Simplify reaction by dividing all coefficients by their GCD
    int gcd_coeff = 0;
    for (const auto &kv : r.reactants)
      gcd_coeff = std::gcd(gcd_coeff, kv.second);
    for (const auto &kv : r.products)
      gcd_coeff = std::gcd(gcd_coeff, kv.second);
    if (gcd_coeff > 1) {
      // Reduce all coefficients by the common factor
      for (auto &kv : r.reactants)
        kv.second /= gcd_coeff;
      for (auto &kv : r.products)
        kv.second /= gcd_coeff;
    }

    // Rebuild the formula in canonical form using the parsed species data
    std::string arrow_str = reversible ? "<->" : "->";
    r.formula =
        side_to_string(r.reactants) + arrow_str + side_to_string(r.products);

    // Add the completed reaction to the output data structure
    data.reactions.push_back(std::move(r));
  }
  // All reactions have been processed, return the parsed data
  return data;
}

} // namespace details
} // namespace gasp2::catalysis::finite_rate
