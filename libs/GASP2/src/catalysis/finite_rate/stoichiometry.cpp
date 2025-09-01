#include "gasp2/catalysis/finite_rate/stoichiometry.hpp"
#include <map>

namespace gasp2::catalysis::finite_rate {
namespace details {

/**
 * @brief Build stoichiometric matrices and species ordering for reaction system.
 *
 * This function extracts the stoichiometric information from the reaction
 * properties and organizes it into matrices suitable for numerical computation.
 * It creates:
 * 1. A unified species ordering (gas species first, then surface species)
 * 2. Reactant and product stoichiometric coefficient matrices
 * 3. Index mappings for surface species
 *
 * The stoichiometric matrices are indexed as [species][reaction], where each
 * element contains the stoichiometric coefficient for that species in that
 * reaction.
 *
 * @param props Reaction properties containing:
 *              - reactions: Vector of reaction data with reactants/products
 * @return StoichiometryData structure containing:
 *         - species_order: Ordered list of all species names
 *         - nu_reactants: Reactant stoichiometric matrix [species][reaction]
 *         - nu_products: Product stoichiometric matrix [species][reaction]  
 *         - surface_index: Indices of surface species in the global ordering
 */
StoichiometryData build_stoichiometry(const ReactionProperties &props) {
  StoichiometryData data;
  
  // Step 1: Build unified species ordering
  // Start with all known species (gas + surface), then add empty site marker
  data.species_order = ctx.all_species;
  data.species_order.push_back("(s)");  // Empty site representation

  const std::size_t n_species = data.species_order.size();
  const std::size_t n_reactions = props.reactions.size();

  // Step 2: Initialize stoichiometric matrices with zeros
  // nu_reactants[i][j] = stoichiometric coefficient of species i as reactant in reaction j
  // nu_products[i][j] = stoichiometric coefficient of species i as product in reaction j
  data.nu_reactants.assign(n_species, std::vector<int>(n_reactions, 0));
  data.nu_products.assign(n_species, std::vector<int>(n_reactions, 0));

  // Step 3: Create species name-to-index mapping for efficient lookup
  // This avoids O(n) searches when processing reaction species
  std::map<std::string, std::size_t> index_map;
  for (std::size_t i = 0; i < data.species_order.size(); ++i)
    index_map[data.species_order[i]] = i;

  // Step 4: Identify surface species indices
  // Surface species come after gas species in the ordering
  // This includes adsorbed species and the empty site "(s)"
  data.surface_index.clear();
  for (std::size_t i = ctx.species_order.size(); i < data.species_order.size();
       ++i)
    data.surface_index.push_back(i);

  // Step 5: Populate stoichiometric matrices from reaction data
  // Process each reaction and extract stoichiometric coefficients
  for (std::size_t r = 0; r < n_reactions; ++r) {
    const auto &reaction = props.reactions[r];
    
    // Fill reactant coefficients: A + B -> products becomes nu_reactants[A][r] = 1, nu_reactants[B][r] = 1
    for (const auto &[sp, coeff] : reaction.reactants) {
      data.nu_reactants[index_map[sp]][r] = coeff;
    }
    
    // Fill product coefficients: reactants -> C + D becomes nu_products[C][r] = 1, nu_products[D][r] = 1  
    for (const auto &[sp, coeff] : reaction.products) {
      data.nu_products[index_map[sp]][r] = coeff;
    }
  }

  return data;
}

} // namespace details
} // namespace gasp2::catalysis::finite_rate
