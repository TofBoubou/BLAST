#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  // Empty species name should trigger initialization failure.
  std::vector<std::string> bad_species{""};
  std::vector<double> bad_mass{28e-3};
  auto res = gasp2::initialize_catalysis(bad_species, bad_mass,
                                         "gamma_model/input_non_catalytic.xml");
  if (res) {
    std::cerr << "Expected failure for empty species name\n";
    return 1;
  }

  // Negative molar mass should also fail.
  bad_species = {"N2"};
  bad_mass = {-28e-3};
  res = gasp2::initialize_catalysis(bad_species, bad_mass,
                                    "gamma_model/input_non_catalytic.xml");
  if (res) {
    std::cerr << "Expected failure for negative molar mass\n";
    return 1;
  }

  // Valid initialization for subsequent flux validation.
  std::vector<std::string> species{"O", "N"};
  std::vector<double> mass{16e-3, 14e-3};
  if (auto init = gasp2::initialize_catalysis(
          species, mass, "gamma_model/input_non_catalytic.xml");
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }

  // Negative density should be rejected.
  std::vector<double> rho{1.0, -1.0};
  auto flux = gasp2::compute_catalysis_fluxes(300.0, rho);
  if (flux) {
    std::cerr << "Expected failure for negative density\n";
    return 1;
  }

  // Non-positive temperature should be rejected.
  rho = {1.0, 1.0};
  flux = gasp2::compute_catalysis_fluxes(0.0, rho);
  if (flux) {
    std::cerr << "Expected failure for non-positive temperature\n";
    return 1;
  }

  std::cout << "OK\n";
  return 0;
}
