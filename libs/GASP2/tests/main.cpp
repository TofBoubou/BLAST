#include <filesystem>
#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  // std::vector<double> rho_wall{1.2, 1.4, 1.4, 0.1, 0.1};
  // // Switch between "air" and "co2" compositions
  // std::string gas_type = "air"; // Change to "co2" for CO2 composition

  // std::vector<std::string> species_order;
  // std::vector<double> molar_masses;
  // std::vector<double> gibbs_free_energy; // Gibbs free energies (J/mol)
  // if (gas_type == "air") {
  //   species_order = {"O2", "N2", "NO", "O", "N"};
  //   molar_masses = {32e-3, 28e-3, 30e-3, 16e-3, 14e-3};
  // } else if (gas_type == "co2") {
  //   species_order = {"CO2", "CO", "O2", "O", "C"};
  //   molar_masses = {44e-3, 28e-3, 32e-3, 16e-3, 12e-3};
  // }
  double Tw = 2000.0; //[K]
// Density of species O: 0.000192428 kg/m^3
// Density of species O2: 0.00346371 kg/m^3
// Gibbs free energy of species O: -117604 J/mol
// Gibbs free energy of species O2: -478313 J/mol
  std::vector<double> gibbs_free_energy;
  gibbs_free_energy = {-117604, -478313}; // Example values
  std::vector<double> rho_wall{0.000192428, 0.00346371};
  std::vector<std::string> species_order = {"O", "O2"};
  std::vector<double> molar_masses = {16e-3, 32e-3};
  const auto input_file = "input.xml";
  if (auto init =
          gasp2::initialize_catalysis(species_order, molar_masses, input_file);
      !init) {
    std::cerr << init.error() << "\n";
    return 1;
  }
  // gibbs_free_energy = {-0.5, 0.0, 0.1, -0.2, -0.3}; // Example values
  auto fluxes_r = gasp2::compute_catalysis_fluxes(Tw, rho_wall, gibbs_free_energy);
  if (!fluxes_r) {
    std::cerr << fluxes_r.error() << "\n";
    return 1;
  }

  std::cout << "Catalysis fluxes:";
  for (double v : fluxes_r->cat_fluxes)
    std::cout << ' ' << v;
  std::cout << '\n';
  // Compute loss factor: since w=m*gamma*imp_flux, gamma=w/m/imp_flux
  auto species = gasp2::SpeciesData(rho_wall,molar_masses,Tw);
  double gamma = fluxes_r->cat_fluxes[0] / species.imp_flux[0] / species.m[0];
  std::cout << "Loss factor (gamma): " << std::scientific << gamma << '\n';
}
