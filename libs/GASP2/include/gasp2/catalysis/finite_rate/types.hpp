#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace gasp2::catalysis::finite_rate {

//---------------------------------------------------------------------------
// Core reaction and stoichiometry types
//---------------------------------------------------------------------------
/**
 * @brief Supported finite-rate reaction categories.
 * @param Adsorption Irreversible gas-phase species adsorption onto surface sites
 * @param Desorption Irreversible surface species desorption to gas phase
 * @param AdsorptionDesorption Reversible adsorption/desorption equilibrium
 * @param EleyRideal Gas-phase species reacts with adsorbed surface species
 * @param LangmuirHinshelwood Two adsorbed surface species react on the surface
 */
enum class ReactionType {
  Adsorption,           ///< Gas + site -> adsorbed species
  Desorption,           ///< Adsorbed species -> Gas + site
  AdsorptionDesorption, ///< Gas + site <-> adsorbed species
  EleyRideal,           ///< Eley-Rideal mechanism
  LangmuirHinshelwood   ///< Langmuir-Hinshelwood mechanism
};

/**
 * @brief Parsed information for a single finite-rate reaction.
 * @param formula Original reaction formula string as parsed from input
 * @param type Mechanism category (adsorption, desorption, Eley-Rideal, etc.)
 * @param reactants Map of reactant species names to their stoichiometric coefficients
 * @param products Map of product species names to their stoichiometric coefficients
 * @param parameters Map of kinetic model coefficients (rate constants, activation energies, etc.)
 * @param kc Map of equilibrium constants for reversible reactions
 */
struct ReactionInput {
  std::string formula;                      ///< Original reaction formula
  ReactionType type;                        ///< Mechanism category
  std::map<std::string, int> reactants;     ///< Reactant stoichiometry
  std::map<std::string, int> products;      ///< Product stoichiometry
  std::map<std::string, double> parameters; ///< Model coefficients
  std::map<std::string, double> kc;         ///< Equilibrium constants
};

namespace details {

//---------------------------------------------------------------------------
// Reaction property containers
//---------------------------------------------------------------------------
/**
 * @brief Data for a single finite-rate reaction after processing.
 * 
 * @param reactants Reactant stoichiometry mapping species names to their stoichiometric coefficients
 * @param products Product stoichiometry mapping species names to their stoichiometric coefficients
 * @param kf Forward rate constant for the reaction
 * @param kb Backward rate constant for the reaction
 * @param Kc Concentration-based equilibrium constant
 * @param Ka Activity-based equilibrium constant
 */
struct ProcessedReaction {
  std::map<std::string, int> reactants; ///< Reactant stoichiometry
  std::map<std::string, int> products;  ///< Product stoichiometry
  double kf{0.0};                       ///< Forward rate constant
  double kb{0.0};                       ///< Backward rate constant
  double Kc{0.0}; ///< Concentration-based equilibrium constant
  double Ka{0.0}; ///< Activity-based equilibrium constant
};

/**
 * @brief Bundle of precomputed properties for all finite-rate reactions.
 * 
 * @param reactions Vector containing all processed reactions with their rate constants and stoichiometry
 * @param site_density Surface site density in mol/m^2
 * @param all_species Ordered list of all species (gas + surface species)
 * @param concentration Gas-phase concentrations in mol/m^3, keyed by species name
 */
struct ReactionProperties {
  std::vector<ProcessedReaction> reactions; ///< All processed reactions
  double site_density{0.0};                 ///< Surface site density (mol/m^2)
  std::vector<std::string> all_species;     ///< Gas + surface species ordering
  std::map<std::string, double>
      concentration; ///< Gas-phase concentrations (mol/m^3) keyed by species
                     ///< name
};

//---------------------------------------------------------------------------
// Stoichiometric matrices
//---------------------------------------------------------------------------
/**
 * @brief Precomputed stoichiometric data for finite-rate reactions.
 * 
 * @param nu_reactants Reactant stoichiometric matrix ν' with dimensions (n_species x n_reactions)
 * @param nu_products Product stoichiometric matrix ν'' with dimensions (n_species x n_reactions)
 * @param species_order Species ordering: gas species followed by surface species and the empty site "(s)"
 * @param surface_index Indices in species_order that correspond to surface species and the empty site
 */
struct StoichiometryData {
  /// Reactant stoichiometric matrix \f$\nu'\f$ with dimensions
  /// (n\_species x n\_reactions).
  std::vector<std::vector<int>> nu_reactants;
  /// Product stoichiometric matrix \f$\nu''\f$ with dimensions
  /// (n\_species x n\_reactions).
  std::vector<std::vector<int>> nu_products;
  /// Species ordering: gas species followed by surface species and the empty
  /// site "(s)".
  std::vector<std::string> species_order;
  /// Indices in species\_order that correspond to surface species and the empty
  /// site.
  std::vector<std::size_t> surface_index;
};

//---------------------------------------------------------------------------
// Finite-rate global context
//---------------------------------------------------------------------------
/**
 * @brief Stores data required by the finite-rate module between initialization
 * and flux computations.
 * 
 * @param species_order Ordering of gas-phase species supplied by the user
 * @param molar_masses Molar masses in kg/mol for each species
 * @param site_density Surface site density in mol/m^2
 * @param reactions Vector of parsed reactions from input
 * @param surface_species Unique surface species detected in reactions
 * @param all_species Gas species followed by surface species
 * @param species_index Map from species name to index in all_species
 * @param debug Debug flag for additional output
 * @param initialized Initialization state flag
 */
struct FiniteRateContext {
  std::vector<std::string>
      species_order; ///< Ordering of gas-phase species supplied by the user
  std::vector<double> molar_masses;     ///< Molar masses (kg/mol)
  double site_density{0.0};             ///< Surface site density (mol/m^2)
  std::vector<ReactionInput> reactions; ///< Parsed reactions
  std::vector<std::string>
      surface_species; ///< Unique surface species detected in reactions
  std::vector<std::string>
      all_species; ///< Gas species followed by surface species
  std::map<std::string, std::size_t>
      species_index;       ///< Map from species name to index in all_species
  bool debug{true};        ///< Debug flag
  bool initialized{false}; ///< Initialization state
};

/// Global instance of the finite-rate context.
extern FiniteRateContext ctx;

} // namespace details
} // namespace gasp2::catalysis::finite_rate
