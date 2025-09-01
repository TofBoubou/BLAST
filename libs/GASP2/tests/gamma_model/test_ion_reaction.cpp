#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

// Validate that reactions involving ionic species with balanced charge pass the
// parser and solver without errors.
int main() {
  std::vector<double> rho_wall{1.0, 1.0, 1e-12};
  std::vector<std::string> species_order{"O-", "O", "O2-"};
  std::vector<double> molar_masses{16e-3, 16e-3, 32e-3};

  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses, "gamma_model/input_ion_balanced.xml");
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }
  auto flux = gasp2::compute_catalysis_fluxes(300.0, rho_wall);
  if (!flux) {
    std::cerr << flux.error() << "\n";
    return 1;
  }
  return 0;
}
