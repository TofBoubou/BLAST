#include "gasp2/gasp2.hpp"
#include "gasp2/catalysis/finite_rate/finite_rate.hpp"
#include "gasp2/catalysis/gamma/gamma_model.hpp"
#include "utils/input_models.hpp"
#include <filesystem>
#include <regex>
#include <stdexcept>
#include <iostream>

namespace gasp2 {

namespace {
/**
 * @brief Global state container for catalysis computations.
 *
 * This struct maintains persistent state shared between initialization
 * and flux computation phases. It stores configuration parameters,
 * species information, and model selection that remain constant
 * throughout the catalysis lifecycle.
 *
 * @param initialized Has initialization been run?
 * @param surface_model Active surface model type (NonCatalytic, Catalysis,
 * Ablation)
 * @param model Active reactions model type (GammaModel, FiniteRates)
 * @param species_order Species order used throughout computations
 * @param molar_masses Corresponding molar masses (kg/mol)
 * @param debug Debug flag for verbose output
 */
struct CatalysisState {
  bool initialized{false}; ///< Has initialization been run?
  SurfaceModel surface_model{
      SurfaceModel::NonCatalytic}; ///< Active surface model type
  ReactionsModel model{
      ReactionsModel::GammaModel};        ///< Active reactions model type
  std::vector<std::string> species_order; ///< Species order used throughout.
  std::vector<double> molar_masses;       ///< Corresponding molar masses.
  bool debug{true};                       ///< Debug flag.
} state;

/**
 * @brief Validate species descriptors provided at initialization.
 *
 * Ensures matching sizes, positive molar masses and valid species names.
 * Chemical species names must start with a letter, followed by letters/digits,
 * and optionally end with charge notation (e.g., H2O, Ca2+, OH-).
 *
 * @param molar_masses Array of molar masses (kg/mol); must be positive values
 * @param species Vector of chemical species names; must follow naming
 * conventions
 *
 * @return Result<void> Success or error describing validation failure
 *
 * @throws None (returns Result type for error handling)
 */
Result<void> validate_init_inputs(std::span<const double> molar_masses,
                                  const std::vector<std::string> &species) {
  if (molar_masses.empty() || species.empty()) {
    return std::unexpected(
        Error{ErrorCode::InvalidInput, "Catalysis: empty species description"});
  }
  if (molar_masses.size() != species.size()) {
    return std::unexpected(Error{ErrorCode::InvalidInput,
                                 "Catalysis: species_order size mismatch"});
  }
  // Regex to validate chemical species names: must start with a letter,
  // followed by letters/digits, optionally ending with charge (e.g., H2O, Ca2+,
  // OH-)
  std::regex sp_regex("^[A-Za-z][A-Za-z0-9]*(?:[+-][0-9]*)?$");
  for (std::size_t i = 0; i < species.size(); ++i) {
    const auto &sp = species[i];
    if (sp.empty()) {
      return std::unexpected(
          Error{ErrorCode::InvalidInput, "Catalysis: empty species name"});
    }
    if (!std::regex_match(sp, sp_regex)) {
      return std::unexpected(Error{ErrorCode::InvalidInput,
                                   "Catalysis: invalid species name " + sp});
    }
    if (molar_masses[i] <= 0.0) {
      return std::unexpected(
          Error{ErrorCode::InvalidInput,
                "Catalysis: non-positive molar mass for " + sp});
    }
  }
  return {};
}

/**
 * @brief Validate flux computation inputs for physical consistency.
 *
 * Ensures temperature is positive, mass densities are positive, and array
 * sizes match the initialized species order.
 *
 * @param T_wall Wall temperature (K); must be positive
 * @param rho_wall Mass densities at wall (kg/m^3); must be positive values
 *
 * @return Result<void> Success or error describing validation failure
 */
Result<void> validate_flux_inputs(double T_wall,
                                  std::span<const double> rho_wall) {
  if (rho_wall.size() != state.species_order.size()) {
    return std::unexpected(
        Error{ErrorCode::InvalidInput, "rho_wall size mismatch"});
  }

  if (T_wall <= 0.0) {
    return std::unexpected(
        Error{ErrorCode::InvalidInput, "T_wall must be positive"});
  }

  for (std::size_t i = 0; i < rho_wall.size(); ++i) {
    if (rho_wall[i] <= 0.0) {
      return std::unexpected(
          Error{ErrorCode::InvalidInput,
                "rho_wall[" + state.species_order[i] + "] must be positive"});
    }
  }

  return {};
}

/**
 * @brief Validate Gibbs free energy inputs for finite-rate reactions.
 *
 * Ensures the Gibbs free energy array is provided and has the correct size
 * when finite-rate kinetics are active.
 *
 * @param gibbs_free_energy Gibbs free energies per mole (J/mol); required for
 * finite-rate
 * @param expected_size Expected array size to match species count
 *
 * @return Result<void> Success or error describing validation failure
 */
Result<void> check_gibbs_energies(std::span<const double> gibbs_free_energy,
                                  std::size_t expected_size) {
  if (gibbs_free_energy.data() == nullptr || gibbs_free_energy.size() == 0) {
    return std::unexpected(
        Error{ErrorCode::InvalidInput,
              "Finite rate model requires Gibbs free energy data"});
  }
  if (gibbs_free_energy.size() != expected_size) {
    return std::unexpected(
        Error{ErrorCode::InvalidInput, "gibbs_free_energy size mismatch"});
  }
  return {};
}
} // namespace

//---------------------------------------------------------------------------
// Initialization entry point
//---------------------------------------------------------------------------

/**
 * @brief Initialize the catalysis system with species and reaction models.
 *
 * Sets up global catalysis state by validating inputs, reading configuration,
 * and initializing the appropriate reaction model. Must be called before flux
 * computations.
 *
 * @param species_order Ordered list of chemical species names for consistent
 * indexing
 * @param molar_masses Array of molar masses (kg/mol) corresponding to
 * species_order
 * @param input_filename Name of XML input file containing reaction definitions
 * @param debug Enable verbose logging and diagnostic output
 *
 * @return Result<void> Success or error describing initialization failure
 *
 * @note Supports NonCatalytic (no reactions), Catalysis (Gamma/FiniteRate),
 * Ablation not allowed
 * @note species_order and molar_masses must have matching sizes with positive
 * masses
 */
Result<void> initialize_catalysis(const std::vector<std::string> &species_order,
                                  std::span<const double> molar_masses,
                                  const std::string &input_filename,
                                  bool debug) {
  try {
    // Validate the inputs passed by the user
    if (auto valid = validate_init_inputs(molar_masses, species_order);
        !valid) {
      return std::unexpected(valid.error());
    }
    // Initialize the Catalysis State
    state = CatalysisState{};
    state.species_order = species_order;
    state.molar_masses.assign(molar_masses.begin(), molar_masses.end());
    state.debug = debug;
    // Read the input model
    std::filesystem::path input_path =
        std::filesystem::path(INPUTS_DIR) / input_filename;
    if (debug) {
      std::cout << "[GASP2] Using input XML: " << input_path.string() << std::endl;
    }
    auto models = gasp2::details::read_input_models(input_path);

    // ========== HANDLE NO-REACTION CASE ==========
    if (!models.has_reactions) {
      // No reactions should only occur for non-catalytic surfaces. Any other
      // surface coupled with empty reactions is ill defined and rejected.
      if (models.surface_model == SurfaceModel::NonCatalytic) {
        state.surface_model = SurfaceModel::NonCatalytic;
        state.initialized = true;
        return {};
      }
      return std::unexpected(
          Error{ErrorCode::InvalidInput, "No reactions specified"});
    }

    // ========== DISPATCH SURFACE MODEL ==========
    if (models.surface_model == SurfaceModel::NonCatalytic) {
      // Reaching this block implies reactions are present, which is
      // incompatible with a non-catalytic surface.
      return std::unexpected(
          Error{ErrorCode::UnsupportedModel,
                "Reactions are not allowed for non_catalytic surfaces"});
    }

    if (models.surface_model == SurfaceModel::Ablation) {
      // Initialization is currently focused on catalysis; ablation requests
      // are rejected explicitly to avoid silent fall-through.
      return std::unexpected(Error{
          ErrorCode::InvalidInput,
          "The catalysis function has been called. Ablation is not permitted"});
    }
    // ========== CATALYSIS ==========
    if (models.surface_model == SurfaceModel::Catalysis) {
      // Dispatch between available reactions models for catalytic surfaces.
      if (models.model == ReactionsModel::GammaModel) {
        if (auto res = catalysis::gamma::initialize(
                state.species_order, state.molar_masses, input_path, debug);
            !res) {
          return res;
        }
        state.surface_model = SurfaceModel::Catalysis;
        state.model = ReactionsModel::GammaModel;
        state.initialized = true;
        return {};
      }
      if (models.model == ReactionsModel::FiniteRates) {
        if (auto res = catalysis::finite_rate::initialize(
                state.species_order, state.molar_masses, input_path, debug);
            !res) {
          return res;
        }
        state.surface_model = SurfaceModel::Catalysis;
        state.model = ReactionsModel::FiniteRates;
        state.initialized = true;
        return {};
      }
      return std::unexpected(
          Error{ErrorCode::NotImplemented, "Reactions model not implemented"});
    }

    // Surface model not recognized. Even if parsing validated the string, a
    // new enumerator may be missing a dispatch above.
    return std::unexpected(
        Error{ErrorCode::NotImplemented, "Surface model not implemented"});
  } catch (const std::exception &e) {
    return std::unexpected(Error{ErrorCode::RuntimeError, e.what()});
  }
}

//---------------------------------------------------------------------------
// Flux computation
//---------------------------------------------------------------------------

/// Compute surface mass fluxes due to CATALYSIS for the current wall state.
///
/// @param T_wall Wall temperature (K); must be positive.
/// @param rho_wall Mass densities at the wall (kg/m^3); each entry must be
///                 positive and the vector size must match the species order.
///
/// @param gibbs_free_energy Gibbs free energies per mole for each species
///                          (J/mol). Required only for finite-rate reactions
///                          and ignored for other models.
///
/// @return CatalysisFluxes on success or an error describing invalid inputs or
///         unsupported configurations.
Result<CatalysisFluxes>
compute_catalysis_fluxes(double T_wall, std::span<const double> rho_wall,
                         std::span<const double> gibbs_free_energy) {
  // Check if the catalysis state has been initialized
  if (!state.initialized) {
    return std::unexpected(
        Error{ErrorCode::RuntimeError, "initialize_catalysis not called"});
  }

  // ========== VALIDATION ==========
  // Check if the inputs are consistent and valid
  if (auto valid = validate_flux_inputs(T_wall, rho_wall); !valid) {
    return std::unexpected(valid.error());
  }
  // ========== NON CATALYTIC CASE ==========
  if (state.surface_model == SurfaceModel::NonCatalytic) {
    CatalysisFluxes fluxes;
    fluxes.cat_fluxes.assign(rho_wall.size(), 0.0);
    return fluxes;
  }
  // ========== CATALYTIC CASE ==========
  // GAMMA MODEL:
  if (state.surface_model == SurfaceModel::Catalysis &&
      state.model == ReactionsModel::GammaModel) {
    // Compute the fluxes
    return catalysis::gamma::compute_fluxes(T_wall, rho_wall);
  }
  // FINITE RATE
  if (state.surface_model == SurfaceModel::Catalysis &&
      state.model == ReactionsModel::FiniteRates) {
    // Validate Gibbs' free energies
    if (auto valid = check_gibbs_energies(gibbs_free_energy, rho_wall.size());
        !valid) {
      return std::unexpected(valid.error());
    }
    // Compute the fluxes
    return catalysis::finite_rate::compute_fluxes(T_wall, rho_wall,
                                                  gibbs_free_energy);
  }

  return std::unexpected(
      Error{ErrorCode::RuntimeError, "Unsupported catalysis configuration"});
}

//---------------------------------------------------------------------------
// State Introspection
//---------------------------------------------------------------------------

/**
 * @brief Return a snapshot of the current catalysis state.
 *
 * Exposes a minimal, read-only view of the global catalysis state so external
 * callers can confirm whether initialization ran and which models are active.
 *
 * @return CatalysisStatus containing initialization flag, surface model and
 *         reactions model.
 */
[[nodiscard]] CatalysisStatus get_catalysis_status() {
  // ========== STATE REPORTING ==========
  return CatalysisStatus{state.initialized, state.surface_model, state.model};
}

} // namespace gasp2
