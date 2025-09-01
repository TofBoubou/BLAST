#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

// Verifies that missing kinetic parameters trigger initialization failure.
int main() {
  std::vector<std::string> species_order{"O"};
  std::vector<double> molar_masses{16e-3};
  auto res = gasp2::initialize_catalysis(
      species_order, molar_masses, "finite_rate/input_missing_parameter.xml");
  if (res) {
    std::cerr << "Expected failure for missing kinetic parameter\n";
    return 1;
  }
  std::cout << "OK\n";
  return 0;
}
