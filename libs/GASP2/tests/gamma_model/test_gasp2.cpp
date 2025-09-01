#include <cmath>
#include <filesystem>
#include <gasp2/gasp2.hpp>
#include <iostream>

int main() {
  std::vector<double> rho_wall = {1.2, 1.4, 1.4, 0.1, 0.1};
  std::vector<std::string> species_order{"O2", "N2", "NO", "O", "N"};
  std::vector<double> molar_masses{32e-3, 28e-3, 30e-3, 16e-3, 14e-3};

  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses, "gamma_model/input_mixed_given.xml");
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }
  auto f_norm_r = gasp2::compute_catalysis_fluxes(300.0, rho_wall);
  if (!f_norm_r) {
    std::cerr << f_norm_r.error() << "\n";
    return 1;
  }
  auto f_norm = *f_norm_r;

  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses, "gamma_model/input_mixed_given.xml");
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }
  auto f_unnorm_r = gasp2::compute_catalysis_fluxes(300.0, rho_wall);
  if (!f_unnorm_r) {
    std::cerr << f_unnorm_r.error() << "\n";
    return 1;
  }
  auto f_unnorm = *f_unnorm_r;

  for (std::size_t i = 0; i < species_order.size(); ++i) {
    if (std::fabs(f_norm.cat_fluxes[i] - f_unnorm.cat_fluxes[i]) > 1e-12) {
      std::cerr << "Flux mismatch\n";
      return 1;
    }
  }

  std::cout << "Catalysis fluxes: ";
  for (const auto &v : f_norm.cat_fluxes)
    std::cout << v << " ";
  std::cout << "\nOK\n";
}
