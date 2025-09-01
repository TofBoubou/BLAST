#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

// Verifies that every surface species requires an adsorption reaction.
int main() {
  std::vector<std::string> species_order{"O", "O2"};
  std::vector<double> molar_masses{16e-3, 32e-3};
  auto res = gasp2::initialize_catalysis(
      species_order, molar_masses, "finite_rate/input_missing_adsorption.xml");
  if (res) {
    std::cerr << "Expected failure for missing adsorption reaction" << '\n';
    return 1;
  }
  std::cout << "OK\n";
  return 0;
}
