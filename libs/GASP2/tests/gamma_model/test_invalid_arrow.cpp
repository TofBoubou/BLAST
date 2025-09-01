#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

// Verify that reversible reactions using invalid arrow symbols are rejected.
int main() {
  std::vector<double> rho_wall{1.0, 1.0};
  std::vector<std::string> species_order{"O", "O2"};
  std::vector<double> molar_masses{16e-3, 32e-3};

  auto init = gasp2::initialize_catalysis(
      species_order, molar_masses, "gamma_model/input_invalid_arrow.xml");
  if (init) {
    std::cerr << "Expected failure for reversible reaction\n";
    return 1;
  }
  return 0;
}
