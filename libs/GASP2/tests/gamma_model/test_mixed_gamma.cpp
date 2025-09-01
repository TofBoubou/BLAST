#include <cmath>
#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  std::vector<double> rho_wall = {1.2, 1.4, 1.4, 1e-12, 1e-12};
  std::vector<std::string> species_order{"N2", "O2", "O", "N", "NO"};
  std::vector<double> molar_masses{28e-3, 32e-3, 16e-3, 14e-3, 30e-3};

  // Valid mixed GammaGiven/GammaT input
  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses, "gamma_model/input_mixed_gamma.xml");
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }
  auto flux_mixed_r = gasp2::compute_catalysis_fluxes(300.0, rho_wall);
  if (!flux_mixed_r) {
    std::cerr << flux_mixed_r.error() << "\n";
    return 1;
  }
  auto flux_mixed = *flux_mixed_r;

  // Reference all-GammaGiven input with equivalent coefficients
  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses, "gamma_model/input_mixed_given.xml");
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }
  auto flux_ref_r = gasp2::compute_catalysis_fluxes(300.0, rho_wall);
  if (!flux_ref_r) {
    std::cerr << flux_ref_r.error() << "\n";
    return 1;
  }
  auto flux_ref = *flux_ref_r;

  for (std::size_t i = 0; i < flux_mixed.cat_fluxes.size(); ++i) {
    if (std::abs(flux_mixed.cat_fluxes[i] - flux_ref.cat_fluxes[i]) > 1e-12) {
      std::cerr << "Mixed gamma flux mismatch at index " << i << "\n";
      return 1;
    }
  }

  // Invalid gamma outside [0,1]
  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses,
          "gamma_model/input_mixed_invalid_range.xml");
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }
  auto flux_bad = gasp2::compute_catalysis_fluxes(300.0, rho_wall);
  if (flux_bad) {
    std::cerr << "Expected failure for gamma outside [0,1]\n";
    return 1;
  }

  // Invalid sum of gammas >1 for shared species
  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses,
          "gamma_model/input_mixed_invalid_sum.xml");
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }
  auto flux_bad_sum_r = gasp2::compute_catalysis_fluxes(300.0, rho_wall);
  if (flux_bad_sum_r) {
    std::cerr << "Expected failure for gamma sum >1\n";
    return 1;
  }

  std::cout << "OK\n";
}
