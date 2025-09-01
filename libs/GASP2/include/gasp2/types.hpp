#pragma once
#include <cmath>
#include <expected>
#include <map>
#include <numbers>
#include <ostream>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace gasp2 {

// ---------------- Catalysis (public types) ----------------
// This file contains only types and enums used in the public library.
// Reaction-specific data structures (e.g., catalysis::gamma::ReactionInput)
// are declared in their respective module headers.

/// Possible top-level surface interaction models.
enum class SurfaceModel { NonCatalytic, Catalysis, Ablation };

/// Available reaction-rate modeling approaches.
enum class ReactionsModel { GammaModel, FiniteRates };

/// Lightweight view of the catalysis subsystem state.
///
/// This struct exposes a minimal subset of the internal catalysis state so
/// applications can introspect the library configuration without accessing
/// private implementation details.
///
/// @param initialized Has initialization completed successfully?
/// @param surface_model Active surface model (NonCatalytic, Catalysis, or
/// Ablation).
/// @param model Active reactions model (GammaModel or FiniteRates).
struct CatalysisStatus {
  bool initialized{false}; ///< Has initialization completed successfully?
  SurfaceModel surface_model{
      SurfaceModel::NonCatalytic};                  ///< Active surface model.
  ReactionsModel model{ReactionsModel::GammaModel}; ///< Active reactions model.
};

enum class ErrorCode {
  InvalidInput,
  UnsupportedModel,
  NotImplemented,
  OutOfRange,
  SpeciesNotFound,
  RuntimeError
};

struct Error {
  ErrorCode code;
  std::string message;
};

inline std::ostream &operator<<(std::ostream &os, const Error &err) {
  return os << err.message;
}

template <typename T> using Result = std::expected<T, Error>;

/// Catalysis fluxes output.
struct CatalysisFluxes {
  /// Flux for each species in the order specified during initialization.
  std::vector<double> cat_fluxes;

  /// Loss factor (gamma) for each species in the same order as `cat_fluxes`.
  ///
  /// The value is defined as the ratio of the net mass flux to the
  /// impinging mass flux for that species. A positive value represents a net
  /// loss of the gas-phase species at the surface.
  std::vector<double> loss_factors;

  /// Steady-state surface densities (mol/m^2) for adsorbed species and vacant
  /// sites. The keys correspond to the species names, e.g. "O(s)" or "(s)" for
  /// empty sites.
  std::map<std::string, double> surface_densities;

  /// Dimensionless surface coverages computed by dividing the surface density
  /// by the total site density provided during initialization.
  std::map<std::string, double> surface_coverages;

  /// Fractional contribution of each reaction to the destruction flux of each
  /// species. Organized as `[species][reaction]` and reporting
  /// \f$w_{ir}/w_i\f$.
  std::vector<std::vector<double>> destruction_fractions;
};
/// Bundle of species-related quantities for convenient access.
struct SpeciesData {
  double T_wall;                      ///< Wall temperature (K).
  std::span<const double> rho;        ///< Mass densities (kg/m^3).
  std::span<const double> molar_mass; ///< Molar masses (kg/mol).
  std::vector<double> number_density; ///< Number densities (1/m^3).
  std::vector<double> m;              ///< Molecular masses (kg).
  std::vector<double> thermal_speed;  ///< Thermal speeds (m/s).
  std::vector<double> u;              ///< One quarter of thermal speeds (m/s).
  std::vector<double> imp_flux;       ///< Impinging fluxes (1/m^2/s).

  /// Construct species-related thermodynamic data from wall conditions.
  ///
  /// @param rho_in Mass densities (kg/m^3); all values must be positive.
  /// @param molar_mass_in Molar masses (kg/mol); all values must be positive.
  /// @param T_wall Wall temperature (K); must be greater than zero.
  ///
  /// @throws std::invalid_argument if sizes mismatch or any value is
  /// non-physical.
  SpeciesData(std::span<const double> rho_in,
              std::span<const double> molar_mass_in, double T_wall)
      : T_wall(T_wall), rho(rho_in), molar_mass(molar_mass_in) {
    // ========== VALIDATION ==========
    if (rho.size() != molar_mass.size()) {
      throw std::invalid_argument("SpeciesData: size mismatch");
    }
    if (T_wall <= 0.0) {
      throw std::invalid_argument("SpeciesData: non-positive temperature");
    }
    for (std::size_t i = 0; i < rho.size(); ++i) {
      if (rho[i] <= 0.0) {
        throw std::invalid_argument(
            "SpeciesData: non-positive density at index " + std::to_string(i));
      }
      if (molar_mass[i] <= 0.0) {
        throw std::invalid_argument(
            "SpeciesData: non-positive molar mass at index " +
            std::to_string(i));
      }
    }

    // ========== PRECOMPUTATION ==========
    constexpr double Na = 6.02214076e23;
    constexpr double R = 8.31446261815324;
    number_density.resize(rho.size());
    m.resize(molar_mass.size());
    thermal_speed.resize(rho.size());
    u.resize(rho.size());
    imp_flux.resize(rho.size());
    for (std::size_t i = 0; i < rho.size(); ++i) {
      m[i] = molar_mass[i] / Na;
      number_density[i] = rho[i] / m[i];
      thermal_speed[i] =
          std::sqrt(8.0 * R * T_wall / (std::numbers::pi * molar_mass[i]));
      u[i] = 0.25 * thermal_speed[i];
      imp_flux[i] = u[i] * number_density[i];
    }
  }
};

} // namespace gasp2
