#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

// Tests finite-rate input parser validation by supplying invalid
// species data. Initialization should fail in both cases.
int main() {
  // Empty species name should trigger initialization failure.
  std::vector<std::string> bad_species{""};
  std::vector<double> bad_mass{28e-3};
  auto res = gasp2::initialize_catalysis(
      bad_species, bad_mass,
      "finite_rate/input_typical_simple_finite_rate.xml");
  if (res) {
    std::cerr << "Expected failure for empty species name\n";
    return 1;
  }

  // Negative molar mass should also fail.
  bad_species = {"O"};
  bad_mass = {-16e-3};
  res = gasp2::initialize_catalysis(
      bad_species, bad_mass,
      "finite_rate/input_typical_simple_finite_rate.xml");
  if (res) {
    std::cerr << "Expected failure for negative molar mass\n";
    return 1;
  }

  std::cout << "OK\n";
  return 0;
}
