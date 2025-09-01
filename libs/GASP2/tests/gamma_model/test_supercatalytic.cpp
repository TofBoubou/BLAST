#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  std::vector<double> rho_wall{1.2, 1.4, 1.4, 1.0, 1e-12};
  std::vector<std::string> species_order{"N2", "O2", "O", "N", "NO"};
  std::vector<double> molar_masses{28e-3, 32e-3, 16e-3, 14e-3, 30e-3};

  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses, "gamma_model/input_supercatalytic.xml");
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
  if (fluxes.cat_fluxes[2] <= 0 || fluxes.cat_fluxes[3] <= 0 ||
      fluxes.cat_fluxes[0] >= 0 || fluxes.cat_fluxes[1] >= 0) {
    std::cerr << "Supercatalytic flux signs incorrect\n";
    return 1;
  }
  auto res = gasp2::initialize_catalysis(
      species_order, molar_masses,
      "gamma_model/input_supercatalytic_first_order.xml");
  if (res) {
    std::cerr << "Expected failure for first_order=true\n";
    return 1;
  }
  std::cout << "OK\n";
  return 0;
}
