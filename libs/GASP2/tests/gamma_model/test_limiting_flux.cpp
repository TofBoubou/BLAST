#include <cmath>
#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  std::vector<double> rho_wall = {1e-12, 1e-12, 1e-6, 5e-6, 1e-12};
  std::vector<std::string> species_order{"N2", "O2", "O", "N", "NO"};
  std::vector<double> molar_masses{28e-3, 32e-3, 16e-3, 14e-3, 30e-3};

  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses, "gamma_model/input_limiting_flux.xml");
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }

  auto flux_r = gasp2::compute_catalysis_fluxes(300.0, rho_wall);
  if (!flux_r) {
    std::cerr << flux_r.error() << "\n";
    return 1;
  }
  auto flux = *flux_r;

  const double Na = 6.02214076e23;
  std::size_t idx_O = 2;
  std::size_t idx_N = 3;
  double gamma_O = 0.2;
  double gamma_N = 0.8;
  double m_O = molar_masses[idx_O] / Na;
  double m_N = molar_masses[idx_N] / Na;
  double imp_O = flux.cat_fluxes[idx_O] / (gamma_O * m_O);
  double imp_N = flux.cat_fluxes[idx_N] / (gamma_N * m_N);

  double lhs = gamma_O * imp_O;
  double rhs = gamma_N * imp_N;
  if (std::abs(lhs - rhs) / lhs > 1e-6) {
    std::cerr << "Compatibility equation not satisfied\n";
    return 1;
  }

  std::cout << "OK\n";
}
