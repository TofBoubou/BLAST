#pragma once
#include "gasp2/types.hpp"
#include <filesystem>
#include <gasp2/catalysis/gamma/types.hpp>
#include <string>
#include <vector>

namespace gasp2::catalysis::gamma {
namespace details {

/// Structure holding the information extracted from the XML file.
struct ParsedInput {
  std::vector<ReactionInput> reactions; ///< List of reactions.
  bool first_order{false};              ///< Use first-order theory for fluxes.
  bool limiting_fluxes{false}; ///< Adjust impinging fluxes for heterogeneous
                               ///< reactions when true.
};

/// Fully parse a GASP2 XML input file for the gamma model and extract the
/// list of reactions and their associated parameters while preserving the
/// provided species order.
[[nodiscard]] ParsedInput
read_gamma_input_file(const std::filesystem::path &input_filename,
                      const std::vector<std::string> &species_order);

} // namespace details
} // namespace gasp2::catalysis::gamma
