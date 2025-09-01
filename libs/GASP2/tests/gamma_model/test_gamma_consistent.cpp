#include <cmath>
#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  std::vector<double> rho_wall{1.2, 1.4, 1.4, 1.0, 1e-12};
  std::vector<std::string> species_order{"N2", "O2", "O", "N", "NO"};
  std::vector<double> molar_masses{28e-3, 32e-3, 16e-3, 14e-3, 30e-3};

  if (auto init =
          gasp2::initialize_catalysis(species_order, molar_masses,
                                      "gamma_model/input_gamma_consistent.xml");
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }
  auto flux_consistent_r = gasp2::compute_catalysis_fluxes(300.0, rho_wall);
  if (!flux_consistent_r) {
    std::cerr << flux_consistent_r.error() << "\n";
    return 1;
  }
  auto flux_consistent = *flux_consistent_r;

  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses,
          "gamma_model/input_gamma_consistent_ref.xml");
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

  for (std::size_t i = 0; i < flux_consistent.cat_fluxes.size(); ++i) {
    if (std::abs(flux_consistent.cat_fluxes[i] - flux_ref.cat_fluxes[i]) >
        1e-12) {
      std::cerr << "GammaConsistent flux mismatch at index " << i << "\n";
      return 1;
    }
  }

  auto flux_bad = gasp2::initialize_catalysis(
      species_order, molar_masses,
      "gamma_model/input_gamma_consistent_homogeneous.xml");
  if (flux_bad) {
    std::cerr << "Expected failure for homogeneous GammaConsistent\n";
    return 1;
  }

  auto flux_missing = gasp2::initialize_catalysis(
      species_order, molar_masses,
      "gamma_model/input_gamma_consistent_no_first_order.xml");
  if (flux_missing) {
    std::cerr << "Expected failure for GammaConsistent without first_order\n";
    return 1;
  }

  std::cout << "OK\n";
}
