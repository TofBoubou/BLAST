#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

// Validates that fractional stoichiometric coefficients below unity are
// rejected.
//
// Returns 0 when the parser correctly rejects the invalid input. A return value
// of 1 indicates that the parser accepted the malformed reaction, which should
// never happen.
int main() {
  std::vector<std::string> species_order{"O", "O2"};
  std::vector<double> molar_masses{16e-3, 32e-3};
  auto res = gasp2::initialize_catalysis(
      species_order, molar_masses,
      "finite_rate/input_invalid_fractional_coeff.xml");
  if (res) {
    std::cerr << "Expected failure for fractional stoichiometric coefficient"
              << '\n';
    return 1;
  }
  std::cout << "OK\n";
  return 0;
}
