#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

// Validates that an Eley-Rideal reaction with a bare-site reactant is accepted.
/**
 * @brief Ensures that ER reactions of the form A + (s) <-> B + C(s) are valid.
 *
 * Initializes the catalysis module using an input file containing the
 * reaction N2 + (s) <-> N + N(s). Successful initialization confirms that
 * bare surface sites can participate directly in Eley-Rideal reactions.
 *
 * @return 0 on success, 1 if initialization fails.
 */
int main() {
  std::vector<std::string> species_order{"N", "N2"};
  std::vector<double> molar_masses{14e-3, 28e-3};

  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses, "finite_rate/input_er_bare_site.xml");
      !init) {
    std::cerr << init.error() << '\n';
    return 1;
  }

  std::cout << "OK\n";
  return 0;
}
