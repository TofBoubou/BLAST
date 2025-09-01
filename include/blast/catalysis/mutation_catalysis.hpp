#pragma once
#include "catalysis_interface.hpp"
#include "../thermophysics/mixture_interface.hpp"
#include <memory>

namespace blast::catalysis {

// Mutation++-based catalysis implementation (wrapper)
class MutationCatalysis : public CatalysisInterface {
public:
  // Constructor with existing MixtureInterface
  explicit MutationCatalysis(const thermophysics::MixtureInterface& mixture);

  // Destructor
  ~MutationCatalysis() override = default;

  // Implementation of CatalysisInterface
  [[nodiscard]] auto compute_surface_fluxes(std::span<const double> partial_densities, double wall_temperature) const
      -> std::expected<std::vector<double>, CatalysisError> override;

  [[nodiscard]] auto species_names() const -> std::vector<std::string> override;

  [[nodiscard]] auto n_species() const noexcept -> std::size_t override;

private:
  const thermophysics::MixtureInterface& mixture_;
};

} // namespace blast::catalysis