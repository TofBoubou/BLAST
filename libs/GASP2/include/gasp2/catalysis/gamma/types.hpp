#pragma once

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace gasp2::catalysis::gamma {

/// Available catalysis models for gamma-model reactions.
enum class CatalysisModel {
  GammaGiven,
  GammaT,
  GammaConsistent,
  GammaBose,
  SuperCatalytic,
  RiniModel
};

/// Stoichiometric matrices and tensors for surface reactions.
struct StoichiometricMatrices {
  std::vector<std::vector<int>> nu_p;    ///< Reactant coefficients.
  std::vector<std::vector<int>> nu_pp;   ///< Product coefficients.
  std::vector<std::vector<int>> nu_diff; ///< Net change (nu_pp - nu_p).
  std::vector<std::vector<int>> nu;      ///< Indicator of species consumption.
  /// mu_l,i,j: for reaction l, production of species i from reactant j
  std::vector<std::vector<std::vector<int>>> mu;
};

/// Input data for a single gamma-model surface reaction.
struct ReactionInput {
  /// Reaction formula, e.g. "O + O -> O2".
  std::string formula;
  /// Reaction model (e.g. GammaGiven).
  CatalysisModel type;
  /// Species-specific recombination probabilities.
  std::optional<std::map<std::string, double>> gammas;
  /// Global recombination probability (\(\gamma_w\)) for the Rini model.
  std::optional<double> gamma_w;
  /// Parameters for temperature-dependent gamma.
  struct GammaTParams {
    double A{0.0};     ///< Pre-exponential factor.
    double E{0.0};     ///< Activation energy (K).
    double tmin{0.0};  ///< Minimum valid temperature.
    double tmax{-1.0}; ///< Maximum valid temperature (-1 for no limit).
  };
  std::optional<std::map<std::string, GammaTParams>> gammaT;
  /// Parameters for the Bose model.
  struct BoseParams {
    double gamma{0.0}; ///< Base recombination probability.
    double p2{0.0};    ///< Probability of CO oxidation.
  };
  std::optional<BoseParams> bose; ///< Optional Bose parameters.
  /// Precomputed recombination probabilities for each reactant species.
  /// Populated during validation to avoid redundant gamma calculations later.
  std::map<std::string, double> computed_gamma;
  /// Indicates whether the reaction involves multiple distinct reactant
  /// species (heterogeneous) or a single species (homogeneous). Set during
  /// input parsing to avoid repeated size checks downstream.
  bool heterogeneous{false};
  /// Reactant stoichiometry.
  std::map<std::string, int> reactants;
  /// Product stoichiometry.
  std::map<std::string, int> products;
};

} // namespace gasp2::catalysis::gamma
