#pragma once

#include <string>
#include <string_view>
#include <unordered_map>

#include <gasp2/catalysis/gamma/types.hpp>
#include <gasp2/types.hpp>

namespace gasp2::catalysis::gamma {

using gasp2::SpeciesData;

/// Compute the recombination probability for the given species in a reaction.
/// Supports GammaGiven, GammaT, and GammaConsistent reactions.
[[nodiscard]] double
compute_gamma(const ReactionInput &reaction, std::string_view species,
              double T_wall, const SpeciesData &data,
              const std::unordered_map<std::string, std::size_t> &index);

} // namespace gasp2::catalysis::gamma
