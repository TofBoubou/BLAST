#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  std::vector<double> rho_wall{1.0};
  std::vector<std::string> species_order{"O"};
  std::vector<double> molar_masses{16e-3};

  auto init = gasp2::initialize_catalysis(
      species_order, molar_masses, "gamma_model/input_wrong_surface.xml");
  if (init)
    std::cerr << "Expected failure for invalid surface model\n";
  return init ? 1 : 0;
}
