#include <cmath>
#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  std::vector<double> rho_wall{1.0, 1.0, 0.05, 1e-12, 1e-12};
  std::vector<std::string> species_order{"O", "CO", "C", "O2", "CO2"};
  std::vector<double> M{16e-3, 28e-3, 12e-3, 32e-3, 44e-3};
  const double Na = 6.02214076e23;
  std::vector<double> m{16e-3 / Na, 28e-3 / Na, 12e-3 / Na, 32e-3 / Na,
                        44e-3 / Na};

  if (auto init = gasp2::initialize_catalysis(
          species_order, M, "gamma_model/input_rini_gamma_three.xml");
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

  std::vector<double> omega(fluxes.cat_fluxes.size());
  for (std::size_t i = 0; i < omega.size(); ++i)
    omega[i] = -fluxes.cat_fluxes[i];

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

  double o_balance = omega[0] / m[0] + omega[1] / m[1] + 2.0 * omega[3] / m[3] +
                     2.0 * omega[4] / m[4];
  double max_O = 0.0;
  max_O = std::max(max_O, std::abs(omega[0] / m[0]));
  max_O = std::max(max_O, std::abs(omega[1] / m[1]));
  max_O = std::max(max_O, std::abs(2.0 * omega[3] / m[3]));
  max_O = std::max(max_O, std::abs(2.0 * omega[4] / m[4]));
  double c_balance = omega[1] / m[1] + omega[2] / m[2] + omega[4] / m[4];
  double max_C = 0.0;
  max_C = std::max(max_C, std::abs(omega[1] / m[1]));
  max_C = std::max(max_C, std::abs(omega[2] / m[2]));
  max_C = std::max(max_C, std::abs(omega[4] / m[4]));
  if (std::abs(o_balance) > tol * std::max(1.0, max_O) ||
      std::abs(c_balance) > tol * std::max(1.0, max_C)) {
    std::cerr << "Element balance failed\n";
    return 1;
  }

  double chi1 = -fluxes.cat_fluxes[3] / m[3]; // O2 product
  double chi2 = fluxes.cat_fluxes[1] / m[1];  // CO reactant
  double chi3 = fluxes.cat_fluxes[2] / m[2];  // C reactant
  if (chi1 < -tol || chi2 < -tol || chi3 < -tol) {
    std::cerr << "Negative chi\n";
    return 1;
  }
  if (chi1 < 0)
    chi1 = 0;
  if (chi2 < 0)
    chi2 = 0;
  if (chi3 < 0)
    chi3 = 0;

  double gamma_w = 0.2;
  const double kB = 1.380649e-23;
  double nO = rho_wall[0] / m[0];
  double nCO = rho_wall[1] / m[1];
  double nC = rho_wall[2] / m[2];
  double Mdown_O = (2.0 / (2.0 - gamma_w)) * nO *
                   std::sqrt(kB * 300.0 / (2.0 * M_PI * m[0]));
  double Mdown_CO = (2.0 / (2.0 - gamma_w)) * nCO *
                    std::sqrt(kB * 300.0 / (2.0 * M_PI * m[1]));
  double Mdown_C = (2.0 / (2.0 - gamma_w)) * nC *
                   std::sqrt(kB * 300.0 / (2.0 * M_PI * m[2]));
  double gamma_CO2 = chi2 / Mdown_CO;
  double gamma_C3 = chi3 / Mdown_C;
  double gamma_O1 = 2.0 * chi1 / Mdown_O;
  double gamma_O2 = chi2 / Mdown_O;
  double gamma_O3 = 2.0 * chi3 / Mdown_O;
  if (gamma_CO2 < -tol || gamma_C3 < -tol || gamma_O1 < -tol ||
      gamma_O2 < -tol || gamma_O3 < -tol) {
    std::cerr << "Negative gamma_r\n";
    return 1;
  }
  if (gamma_CO2 > gamma_w + tol || gamma_C3 > gamma_w + tol ||
      gamma_O1 > gamma_w + tol || gamma_O2 > gamma_w + tol ||
      gamma_O3 > gamma_w + tol) {
    std::cerr << "Gamma exceeds bound\n";
    return 1;
  }

  double resid_O = (2.0 * chi1 + chi2 + 2.0 * chi3) / Mdown_O - gamma_w;
  double resid_CO = chi2 / Mdown_CO - gamma_w;
  double resid_C = chi3 / Mdown_C - gamma_w;
  if (std::abs(resid_O) > tol || std::abs(resid_CO) > tol ||
      std::abs(resid_C) > tol) {
    std::cerr << "Constraint residual\n";
    return 1;
  }

  std::cout << "OK\n";
  return 0;
}
