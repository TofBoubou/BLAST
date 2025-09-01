#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

// Ensures that reactions not matching their mechanism-specific stoichiometric
// patterns are rejected during initialization.
int main() {
  std::vector<std::string> species_order{"O", "O2"};
  std::vector<double> molar_masses{16e-3, 32e-3};
  std::vector<std::string> files{
      "finite_rate/input_invalid_adsorption_form.xml",
      "finite_rate/input_invalid_desorption_form.xml",
      "finite_rate/input_invalid_adsdes_form.xml",
      "finite_rate/input_invalid_er_form.xml",
      "finite_rate/input_invalid_lh_form.xml"};

  for (const auto &f : files) {
    auto res = gasp2::initialize_catalysis(species_order, molar_masses, f);
    if (res) {
      std::cerr << "Expected failure for invalid reaction form in " << f
                << '\n';
      return 1;
    }
  }
  std::cout << "OK\n";
  return 0;
}
