#include "blast/catalysis/mutation_catalysis.hpp"
#include <format>

namespace blast::catalysis {

MutationCatalysis::MutationCatalysis(const thermophysics::MixtureInterface& mixture)
    : mixture_(mixture) {
}

auto MutationCatalysis::compute_surface_fluxes(std::span<const double> partial_densities, double wall_temperature) const
    -> std::expected<std::vector<double>, CatalysisError> {
  
  // Delegate to the underlying Mutation++ mixture
  auto result = mixture_.surface_reaction_rates(partial_densities, wall_temperature);
  if (!result) {
    return std::unexpected(CatalysisError(std::format("Mutation++ catalysis failed: {}", result.error().message())));
  }
  
  return result.value();
}

auto MutationCatalysis::species_names() const -> std::vector<std::string> {
  std::vector<std::string> names;
  names.reserve(mixture_.n_species());
  
  for (std::size_t i = 0; i < mixture_.n_species(); ++i) {
    names.emplace_back(mixture_.species_name(i));
  }
  
  return names;
}

auto MutationCatalysis::n_species() const noexcept -> std::size_t {
  return mixture_.n_species();
}

} // namespace blast::catalysis