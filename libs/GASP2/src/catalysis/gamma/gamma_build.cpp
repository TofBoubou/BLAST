#include "gasp2/catalysis/gamma/gamma_build.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace gasp2::catalysis::gamma {

[[nodiscard]] double
compute_gamma(const ReactionInput &reaction, std::string_view species,
              double T_wall, const SpeciesData &data,
              const std::unordered_map<std::string, std::size_t> &index) {
  std::cout << "DEBUG compute_gamma - species: " << species << ", reaction.type: " << static_cast<int>(reaction.type) << std::endl;
  std::cout.flush();
  std::cout << "DEBUG compute_gamma - reaction.gammas has_value: " << (reaction.gammas.has_value() ? "true" : "false") << std::endl;
  std::cout.flush();
  if (reaction.gammas.has_value()) {
    std::cout << "DEBUG compute_gamma - reaction.gammas size: " << reaction.gammas->size() << std::endl;
    for (const auto& [sp, val] : *reaction.gammas) {
      std::cout << "DEBUG compute_gamma - gammas map contains: " << sp << " -> " << val << std::endl;
    }
  }
  
  switch (reaction.type) {
  case CatalysisModel::GammaGiven: {
    if (!reaction.gammas.has_value()) {
      std::cout << "DEBUG compute_gamma - ERROR: reaction.gammas is null!" << std::endl;
      return 0.0;
    }
    if (reaction.gammas->find(std::string(species)) == reaction.gammas->end()) {
      std::cout << "DEBUG compute_gamma - ERROR: species " << species << " not found in gammas map!" << std::endl;
      return 0.0;
    }
    double gamma_val = reaction.gammas->at(std::string(species));
    std::cout << "DEBUG compute_gamma - returning gamma: " << gamma_val << " for species: " << species << std::endl;
    return gamma_val;
  }
  case CatalysisModel::GammaT: {
    const auto &p = reaction.gammaT->at(std::string(species));
    if ((p.tmin != 0.0 && T_wall < p.tmin) ||
        (p.tmax != -1.0 && T_wall > p.tmax)) {
      throw std::runtime_error("GammaT temperature out of range for species: " +
                               std::string(species));
    }
    return p.A * std::exp(p.E / T_wall);
  }
  case CatalysisModel::GammaConsistent: {
    auto given = *reaction.gammas->begin();
    std::size_t idx_g = index.at(given.first);
    double imp_g = data.imp_flux[idx_g];
    std::size_t idx_s = index.at(std::string(species));
    double imp_s = data.imp_flux[idx_s];
    if (species == given.first) {
      return given.second;
    }
    double gamma = given.second * imp_g / imp_s;
    return gamma;
  }
  default:
    break;
  }
  throw std::runtime_error("Unsupported catalysis model");
}

} // namespace gasp2::catalysis::gamma
