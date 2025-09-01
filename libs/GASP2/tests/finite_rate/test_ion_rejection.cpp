#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

// Ensures ionic species are rejected by the finite-rate parser.
int main() {
  std::vector<std::string> species_order{"e-"};
  std::vector<double> molar_masses{1e-3};
  auto res = gasp2::initialize_catalysis(species_order, molar_masses,
                                         "finite_rate/input_invalid_ion.xml");
  if (res) {
    std::cerr << "Expected failure for ionic species" << '\n';
    return 1;
  }
  std::cout << "OK\n";
  return 0;
}
