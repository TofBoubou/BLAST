#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

// Validates that a simple adsorption reaction is accepted by the parser.
/**
 * @brief Entry point for validating simple adsorption parsing.
 *
 * Creates a minimal species set and attempts to initialize the catalysis
 * module using an input file that contains the reaction O + (s) -> O(s).
 * Successful initialization confirms that the adsorption reaction passes
 * stoichiometric validation.
 *
 * @return 0 on success, 1 if initialization fails.
 */
int main() {
  std::vector<std::string> species_order{"O", "O2"};
  std::vector<double> molar_masses{16e-3, 32e-3};

  if (auto init =
          gasp2::initialize_catalysis(species_order, molar_masses, "input.xml");
      !init) {
    std::cerr << init.error() << '\n';
    return 1;
  }

  std::cout << "OK\n";
  return 0;
}
