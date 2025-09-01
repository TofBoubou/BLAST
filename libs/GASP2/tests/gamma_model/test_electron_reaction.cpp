#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

// Validate that reactions involving electrons satisfy both mass and charge
// conservation when parsed and solved by the gamma model.
int main() {
  // Wall number densities for NO+, electron, and NO respectively. Electrons
  // have near-zero wall density.
  std::vector<double> rho_wall{1.0, 1e-12, 1.0};
  // Species order must match reaction participants: NO+, e-, NO.
  std::vector<std::string> species_order{"NO+", "e-", "NO"};
  // Approximate molar masses (kg/mol). Electron mass is ~5.49e-7 kg/mol.
  std::vector<double> molar_masses{30e-3, 5.49e-7, 30e-3};

  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses,
          "gamma_model/input_electron_reaction.xml");
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }
  // The electron consumes the charge of NO+, producing neutral NO.
  auto flux = gasp2::compute_catalysis_fluxes(300.0, rho_wall);
  if (!flux) {
    std::cerr << flux.error() << "\n";
    return 1;
  }
  return 0;
}
