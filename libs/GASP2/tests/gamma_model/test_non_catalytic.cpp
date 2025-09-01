#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  std::vector<double> rho_wall{1.0, 2.0};
  std::vector<std::string> species_order{"O", "N"};
  std::vector<double> molar_masses{16e-3, 14e-3};

  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses, "gamma_model/input_non_catalytic.xml");
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }

  auto fluxes_r = gasp2::compute_catalysis_fluxes(300.0, rho_wall);
  if (!fluxes_r) {
    std::cerr << fluxes_r.error() << "\n";
    return 1;
  }
  auto fluxes = *fluxes_r;
  if (fluxes.cat_fluxes.size() != rho_wall.size()) {
    std::cerr << "Unexpected number of fluxes\n";
    return 1;
  }
  for (double f : fluxes.cat_fluxes) {
    if (f != 0.0) {
      std::cerr << "Flux should be zero for non_catalytic surface\n";
      return 1;
    }
  }
  std::cout << "OK\n";
  return 0;
}
