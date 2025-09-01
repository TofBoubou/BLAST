#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

// Ensure that reactions with unbalanced electric charge are rejected during
// validation.
int main() {
  std::vector<double> rho_wall{1.0, 0.0, 0.0};
  std::vector<std::string> species_order{"NO+", "N", "O"};
  std::vector<double> molar_masses{30e-3, 14e-3, 16e-3};

  auto init = gasp2::initialize_catalysis(
      species_order, molar_masses, "gamma_model/input_charge_imbalance.xml");
  if (init) {
    std::cerr << "Expected failure for charge imbalance\n";
    return 1;
  }
  return 0;
}
