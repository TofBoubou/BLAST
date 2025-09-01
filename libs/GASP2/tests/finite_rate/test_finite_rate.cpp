#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

// Ensures invalid finite-rate input file is rejected by the parser.
int main() {
  std::vector<double> rho_wall{1.0};
  std::vector<std::string> species_order{"O"};
  std::vector<double> molar_masses{16e-3};

  auto init = gasp2::initialize_catalysis(
      species_order, molar_masses,
      "finite_rate/input_typical_simple_finite_rate.xml");
  if (init)
    std::cerr << "Expected failure for invalid finite-rate input\n";
  return init ? 1 : 0;
}
