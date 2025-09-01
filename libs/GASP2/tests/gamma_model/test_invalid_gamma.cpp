#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  std::vector<double> rho_wall{1.2, 1.4, 1.4, 0.0, 0.0, 0.0};
  std::vector<std::string> species_order{"N2", "O2", "O", "N", "NO", "C"};
  std::vector<double> molar_masses{28e-3, 32e-3, 16e-3, 14e-3, 30e-3, 12e-3};

  auto init = gasp2::initialize_catalysis(
      species_order, molar_masses, "gamma_model/input_invalid_gamma.xml");
  if (init)
    std::cerr << "Expected failure for invalid gamma\n";
  return init ? 1 : 0;
}
