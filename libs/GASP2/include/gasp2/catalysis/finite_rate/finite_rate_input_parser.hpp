#pragma once
#include "gasp2/catalysis/finite_rate/finite_rate.hpp"
#include <filesystem>
#include <vector>

namespace gasp2::catalysis::finite_rate {
namespace details {

/**
 * @brief Structure holding information parsed from the finite-rate XML file.
 * @param site_density Surface site density in mol/m^2 extracted from <site_density> tag
 * @param reactions Vector of ReactionInput parsed reaction definitions with their parameters and stoichiometry
 */
struct ParsedInput {
  double site_density{0.0};             ///< Surface site density (mol/m^2)
  std::vector<ReactionInput> reactions; ///< List of parsed reactions
};

/// Read and minimally parse a finite-rate XML input file. This routine
/// extracts the surface site density and reaction descriptions but performs no
/// detailed validation, leaving that to dedicated checks.
[[nodiscard]] ParsedInput
read_finite_rate_input_file(const std::filesystem::path &input_filename);

} // namespace details
} // namespace gasp2::catalysis::finite_rate
