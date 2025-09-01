#pragma once
#include "catalysis_interface.hpp"
#include <memory>
#include <string>

namespace blast::catalysis {

// GASP2-based catalysis implementation
class Gasp2Catalysis : public CatalysisInterface {
public:
  // Constructor with XML input file and species order
  explicit Gasp2Catalysis(const std::string& xml_file, const std::vector<std::string>& species_order, 
                          const std::vector<double>& molar_masses);

  // Destructor
  ~Gasp2Catalysis() override = default;

  // Disable copy/move for simplicity
  Gasp2Catalysis(const Gasp2Catalysis&) = delete;
  Gasp2Catalysis& operator=(const Gasp2Catalysis&) = delete;
  Gasp2Catalysis(Gasp2Catalysis&&) = delete;
  Gasp2Catalysis& operator=(Gasp2Catalysis&&) = delete;

  // Implementation of CatalysisInterface
  [[nodiscard]] auto compute_surface_fluxes(std::span<const double> partial_densities, double wall_temperature) const
      -> std::expected<std::vector<double>, CatalysisError> override;

  [[nodiscard]] auto species_names() const -> std::vector<std::string> override;

  [[nodiscard]] auto n_species() const noexcept -> std::size_t override;

private:
  std::vector<std::string> species_order_;
  std::vector<double> molar_masses_;
  bool initialized_;
};

} // namespace blast::catalysis