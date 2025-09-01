#pragma once
#include "../core/exceptions.hpp"
#include <expected>
#include <span>
#include <string_view>
#include <vector>

namespace blast::catalysis {

// Error type for catalysis operations
class CatalysisError : public core::BlastException {
public:
  explicit CatalysisError(std::string_view message, std::source_location location = std::source_location::current())
      : BlastException(std::format("Catalysis Error: {}", message), location) {}
};

// Abstract interface for catalysis providers
class CatalysisInterface {
public:
  virtual ~CatalysisInterface() = default;

  // Compute surface reaction rates (species destruction rates in kg/(m²·s))
  [[nodiscard]] virtual auto compute_surface_fluxes(std::span<const double> partial_densities, double wall_temperature) const
      -> std::expected<std::vector<double>, CatalysisError> = 0;

  // Get species names in the same order as flux computation
  [[nodiscard]] virtual auto species_names() const -> std::vector<std::string> = 0;

  // Get number of species
  [[nodiscard]] virtual auto n_species() const noexcept -> std::size_t = 0;
};

} // namespace blast::catalysis