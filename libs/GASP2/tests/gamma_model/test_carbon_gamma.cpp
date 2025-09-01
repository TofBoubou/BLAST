#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  std::vector<double> rho_wall{0.1, 0.2, 0.3, 1e-12, 1e-12};
  std::vector<std::string> species_order{"C", "O", "CO", "CO2", "O2"};
  std::vector<double> molar_masses{12e-3, 16e-3, 28e-3, 44e-3, 32e-3};

  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses, "gamma_model/input_carbon_gamma.xml");
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
  if (fluxes.cat_fluxes.size() != species_order.size()) {
    std::cerr << "Unexpected number of fluxes\n";
    return 1;
  }
  // O and CO should be consumed (positive), CO2 and O2 produced (negative)
  if (fluxes.cat_fluxes[1] <= 0 || fluxes.cat_fluxes[3] >= 0 ||
      fluxes.cat_fluxes[4] >= 0) {
    std::cerr << "Flux signs incorrect\n";
    return 1;
  }
  std::cout << "OK\n";
  return 0;
}
