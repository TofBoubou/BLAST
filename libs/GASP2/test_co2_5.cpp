#include <gasp2/gasp2.hpp>
#include <iostream>
#include <string>
#include <vector>

int main() {
  // Use BLAST's species order: {"CO2", "CO", "O2", "O", "C"}
  std::vector<std::string> species_order{"CO2", "CO", "O2", "O", "C"};
  std::vector<double> molar_masses{44e-3, 28e-3, 32e-3, 16e-3, 12e-3};
  std::vector<double> rho_wall{0.1, 0.2, 0.3, 0.1, 1e-12}; // Non-zero for CO2, CO, O2, O

  std::cout << "Testing GASP2 with BLAST configuration..." << std::endl;
  std::cout << "Species order: ";
  for (const auto& sp : species_order) {
    std::cout << sp << " ";
  }
  std::cout << std::endl;

  // Test with our XML file
  if (auto init = gasp2::initialize_catalysis(
          species_order, molar_masses, "CO2_5_gamma.xml");
      !init) {
    std::cerr << "Error: " << init.error() << "\n";
    return 1;
  }
  
  std::cout << "GASP2 initialization successful!" << std::endl;
  
  auto fluxes_r = gasp2::compute_catalysis_fluxes(3000.0, rho_wall);
  if (!fluxes_r) {
    std::cerr << "Flux computation error: " << fluxes_r.error() << "\n";
    return 1;
  }
  
  auto fluxes = *fluxes_r;
  std::cout << "Catalysis fluxes computed successfully!" << std::endl;
  for (size_t i = 0; i < species_order.size(); ++i) {
    std::cout << "  " << species_order[i] << ": " << fluxes.cat_fluxes[i] << std::endl;
  }
  
  return 0;
}