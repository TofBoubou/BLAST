#include <cmath>
#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  std::vector<double> rho_wall{1.0, 1.0, 1e-12, 1e-12};
  std::vector<std::string> species_order{"O", "CO", "O2", "CO2"};
  std::vector<double> M{16e-3, 28e-3, 32e-3, 44e-3};
  const double Na = 6.02214076e23;
  std::vector<double> m{16e-3 / Na, 28e-3 / Na, 32e-3 / Na, 44e-3 / Na};

  if (auto init = gasp2::initialize_catalysis(
          species_order, M, "gamma_model/input_rini_gamma.xml");
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }
  auto fluxes_r = gasp2::compute_catalysis_fluxes(300.0, rho_wall);
  if (!fluxes_r) {
    std::cerr << fluxes_r.error() << "\n";
    return 1;
  }
  auto fluxes = *fluxes_r;
  const double tol = 1e-8;

  // Recover omega (production positive)
  std::vector<double> omega(fluxes.cat_fluxes.size());
  for (std::size_t i = 0; i < omega.size(); ++i)
    omega[i] = -fluxes.cat_fluxes[i];

  // Mass balance
  double mass_sum = 0.0;
  double max_mass = 0.0;
  for (double w : omega) {
    mass_sum += w;
    max_mass = std::max(max_mass, std::abs(w));
  }
  if (std::abs(mass_sum) > tol * std::max(1.0, max_mass)) {
    std::cerr << "Mass not conserved\n";
    return 1;
  }

  // Element balance (O and C)
  double o_balance = omega[0] / m[0] + omega[1] / m[1] + 2.0 * omega[2] / m[2] +
                     2.0 * omega[3] / m[3];
  double max_O = 0.0;
  max_O = std::max(max_O, std::abs(omega[0] / m[0]));
  max_O = std::max(max_O, std::abs(omega[1] / m[1]));
  max_O = std::max(max_O, std::abs(2.0 * omega[2] / m[2]));
  max_O = std::max(max_O, std::abs(2.0 * omega[3] / m[3]));
  double c_balance = omega[1] / m[1] + omega[3] / m[3];
  double max_C = 0.0;
  max_C = std::max(max_C, std::abs(omega[1] / m[1]));
  max_C = std::max(max_C, std::abs(omega[3] / m[3]));
  if (std::abs(o_balance) > tol * std::max(1.0, max_O) ||
      std::abs(c_balance) > tol * std::max(1.0, max_C)) {
    std::cerr << "Element balance failed\n";
    return 1;
  }

  // Reaction frequencies and gammas
  double chi1 = omega[2] / m[2]; // O2
  double chi2 = omega[3] / m[3]; // CO2
  if (chi1 < -tol || chi2 < -tol) {
    std::cerr << "Negative chi\n";
    return 1;
  }
  if (chi1 < 0)
    chi1 = 0;
  if (chi2 < 0)
    chi2 = 0;

  // Compute impinging fluxes
  const double kB = 1.380649e-23;
  double nO = rho_wall[0] / m[0];
  double nCO = rho_wall[1] / m[1];
  double Mdown_O =
      (2.0 / (2.0 - 0.2)) * nO * std::sqrt(kB * 300.0 / (2.0 * M_PI * m[0]));
  double Mdown_CO =
      (2.0 / (2.0 - 0.2)) * nCO * std::sqrt(kB * 300.0 / (2.0 * M_PI * m[1]));

  double gamma_CO2 = chi2 / Mdown_CO;
  double gamma_O2 = chi2 / Mdown_O;
  double gamma_O1 = 2.0 * chi1 / Mdown_O;

  if (gamma_CO2 < -tol || gamma_O2 < -tol || gamma_O1 < -tol) {
    std::cerr << "Negative gamma_r\n";
    return 1;
  }
  double gw = 0.2;
  if (gamma_CO2 > gw + tol || gamma_O2 > gw + tol || gamma_O1 > gw + tol) {
    std::cerr << "Gamma exceeds bound\n";
    return 1;
  }

  double resid_O = (2.0 * chi1 + chi2) / Mdown_O - gw;
  double resid_CO = chi2 / Mdown_CO - gw;
  if (std::abs(resid_O) > tol || std::abs(resid_CO) > tol) {
    std::cerr << "Constraint residual\n";
    return 1;
  }

  auto res = gasp2::initialize_catalysis(
      species_order, M, "gamma_model/input_rini_gamma_first_order.xml");
  if (res) {
    std::cerr << "Expected failure for first_order=true\n";
    return 1;
  }

  std::cout << "OK\n";
  return 0;
}
