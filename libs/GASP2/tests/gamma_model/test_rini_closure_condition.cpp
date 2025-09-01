#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  std::vector<std::string> species_order{"O", "CO", "X", "O2", "COX"};
  std::vector<double> M{16e-3, 28e-3, 1e-3, 32e-3, 45e-3};
  if (auto init = gasp2::initialize_catalysis(
          species_order, M, "gamma_model/input_rini_closure_fail.xml");
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }
  std::vector<double> rho{1.0, 1.0, 1e-12, 1e-12, 1e-12};
  auto fluxes = gasp2::compute_catalysis_fluxes(300.0, rho);
  if (fluxes) {
    std::cerr << "Expected failure due to closure condition\n";
    return 1;
  }
  if (fluxes.error().message.find("N + M") == std::string::npos) {
    std::cerr << "Unexpected error: " << fluxes.error() << "\n";
    return 1;
  }
  std::cout << "OK\n";
  return 0;
}
