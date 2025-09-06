#include "blast/boundary_layer/solver/initial_guess_factory.hpp"
#include "blast/boundary_layer/solver/expected_utils.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <format>

namespace blast::boundary_layer::solver {

InitialGuessFactory::InitialGuessFactory(const grid::BoundaryLayerGrid& grid,
                                        const thermophysics::MixtureInterface& mixture) noexcept
    : grid_(grid), mixture_(mixture) {}

auto InitialGuessFactory::create_initial_guess(
    int station,
    double xi,
    const conditions::BoundaryConditions& bc,
    double T_edge
) const -> std::expected<equations::SolutionState, SolverError> {

    const auto n_species = mixture_.n_species();

    if (n_species == 1) {
        return create_single_species_guess(bc, T_edge);
    } else {
        return create_multi_species_guess(bc, T_edge);
    }
}

auto InitialGuessFactory::create_single_species_guess(
    const conditions::BoundaryConditions& bc,
    double T_edge
) const -> std::expected<equations::SolutionState, SolverError> {

    const auto n_eta = grid_.n_eta();
    const auto n_species = mixture_.n_species();
    const double eta_max = grid_.eta_max();

    equations::SolutionState guess(n_eta, n_species);

    // Initialize V field to zero
    std::fill(guess.V.begin(), guess.V.end(), 0.0);

    // Set unit composition for single species
    guess.c.eigen().setOnes();

    // Compute wall equilibrium enthalpy
    std::array<double, 1> c_wall{{1.0}};
    auto h_wall_eq_result = mixture_.mixture_enthalpy(c_wall, bc.Tw(), bc.P_e());
    if (!h_wall_eq_result) {
        return std::unexpected(NumericError(
            std::format("Failed to compute wall equilibrium enthalpy: {}", h_wall_eq_result.error().message())));
    }
    double g_wall = h_wall_eq_result.value() / bc.he();

    // Create profiles
    for (std::size_t i = 0; i < n_eta; ++i) {
        const double eta = static_cast<double>(i) * eta_max / (n_eta - 1);
        const double eta_norm = static_cast<double>(i) / (n_eta - 1);

        guess.F[i] = 1.0 - 1.034 * std::exp(-5.628 * eta / eta_max);
        guess.g[i] = g_wall + eta_norm * (1.0 - g_wall);
        guess.T[i] = bc.Tw() + guess.F[i] * (T_edge - bc.Tw());
    }

    return guess;
}

// ORIGINAL VERSION - COMMENTED OUT
/*
auto InitialGuessFactory::create_multi_species_guess(
    const conditions::BoundaryConditions& bc,
    double T_edge
) const -> std::expected<equations::SolutionState, SolverError> {

    const auto n_eta = grid_.n_eta();
    const auto n_species = mixture_.n_species();
    const double eta_max = grid_.eta_max();

    equations::SolutionState guess(n_eta, n_species);

    // Initialize V field to zero
    std::fill(guess.V.begin(), guess.V.end(), 0.0);

    // Step 1: Compute equilibrium composition at the wall
    std::vector<double> equilibrium_result;
    BLAST_TRY_ASSIGN_CTX(
        equilibrium_result,
        mixture_.equilibrium_composition(bc.Tw(), bc.P_e()),
        "Failed to compute equilibrium composition at wall conditions"
    );

    // Step 2: Extract edge composition from boundary conditions
    const auto& c_edge = bc.c_e();
    if (c_edge.size() != n_species) {
        return std::unexpected(InitializationError(
            std::format("Edge composition size mismatch: expected {}, got {}", n_species, c_edge.size())));
    }

    // Step 3: Compute wall equilibrium enthalpy
    double h_wall_equilibrium;
    BLAST_TRY_ASSIGN_CTX(
        h_wall_equilibrium,
        mixture_.mixture_enthalpy(equilibrium_result, bc.Tw(), bc.P_e()),
        "Failed to compute wall equilibrium enthalpy"
    );
    double g_wall = h_wall_equilibrium / bc.he();

    // Step 4: Create profiles for multi-species case
    for (std::size_t i = 0; i < n_eta; ++i) {
        const double eta = static_cast<double>(i) * eta_max / (n_eta - 1);
        const double eta_norm = static_cast<double>(i) / (n_eta - 1);

        // Momentum profile (analytical)
        guess.F[i] = 1.0 - 1.034 * std::exp(-5.628 * eta / eta_max);
        
        // Energy profile (linear from wall equilibrium to edge)
        guess.g[i] = g_wall + eta_norm * (1.0 - g_wall);
        
        // Temperature profile
        guess.T[i] = bc.Tw() + guess.F[i] * (T_edge - bc.Tw());

        // Species profiles with smooth transition
        const double transition_factor = transition_function(eta_norm);
        double sum_interpolated = 0.0;

        for (std::size_t s = 0; s < n_species; ++s) {
            guess.c(s, i) = equilibrium_result[s] * (1.0 - transition_factor) + c_edge[s] * transition_factor;
            sum_interpolated += guess.c(s, i);
        }

        // Enforce mass conservation: normalize so that sum_j c_j = 1
        if (sum_interpolated > 1e-15) {
            for (std::size_t s = 0; s < n_species; ++s) {
                guess.c(s, i) /= sum_interpolated;
            }
        } else {
            return std::unexpected(InitializationError(
                std::format("Mass conservation violation: sum of species concentrations is too small ({})", sum_interpolated)));
        }
    }

    return guess;
}
*/

// MODIFIED VERSION
auto InitialGuessFactory::create_multi_species_guess(
    const conditions::BoundaryConditions& bc,
    double T_edge
) const -> std::expected<equations::SolutionState, SolverError> {

    const auto n_eta = grid_.n_eta();
    const auto n_species = mixture_.n_species();
    const double eta_max = grid_.eta_max();

    equations::SolutionState guess(n_eta, n_species);

    // Initialize V field to zero
    std::fill(guess.V.begin(), guess.V.end(), 0.0);

    // Step 1: Compute equilibrium composition at the wall
    std::vector<double> equilibrium_result;
    BLAST_TRY_ASSIGN_CTX(
        equilibrium_result,
        mixture_.equilibrium_composition(bc.Tw(), bc.P_e()),
        "Failed to compute equilibrium composition at wall conditions"
    );

    // Step 2: Extract edge composition from boundary conditions
    const auto& c_edge = bc.c_e();
    if (c_edge.size() != n_species) {
        return std::unexpected(InitializationError(
            std::format("Edge composition size mismatch: expected {}, got {}", n_species, c_edge.size())));
    }

    // Step 3: Compute wall equilibrium enthalpy
    double h_wall_equilibrium;
    BLAST_TRY_ASSIGN_CTX(
        h_wall_equilibrium,
        mixture_.mixture_enthalpy(equilibrium_result, bc.Tw(), bc.P_e()),
        "Failed to compute wall equilibrium enthalpy"
    );
    double g_wall = h_wall_equilibrium / bc.he();

    // Step 4: Create profiles for multi-species case
    for (std::size_t i = 0; i < n_eta; ++i) {
        const double eta = static_cast<double>(i) * eta_max / (n_eta - 1);
        const double eta_norm = static_cast<double>(i) / (n_eta - 1);

        // Momentum profile (analytical)
        guess.F[i] = (1.0 - 0.99 * std::exp(-5.628 * eta / eta_max)) / (1.0 - 0.99 * std::exp(-5.628));
        
        // MODIFIED: Energy profile now follows F profile instead of linear
        guess.g[i] = g_wall + guess.F[i] * (1.0 - g_wall);
        
        // Temperature profile
        guess.T[i] = bc.Tw() + guess.F[i] * (T_edge - bc.Tw());

        // MODIFIED: Species profiles with simple linear interpolation instead of tanh
        double sum_interpolated = 0.0;

        for (std::size_t s = 0; s < n_species; ++s) {
            // Simple linear interpolation
            guess.c(s, i) = equilibrium_result[s] * (1.0 - eta_norm) + c_edge[s] * eta_norm;
            sum_interpolated += guess.c(s, i);
        }

        // Enforce mass conservation: normalize so that sum_j c_j = 1
        if (sum_interpolated > 1e-15) {
            for (std::size_t s = 0; s < n_species; ++s) {
                guess.c(s, i) /= sum_interpolated;
            }
        } else {
            return std::unexpected(InitializationError(
                std::format("Mass conservation violation: sum of species concentrations is too small ({})", sum_interpolated)));
        }
    }

    return guess;
}

auto InitialGuessFactory::transition_function(double eta_norm, double eta_center, double sharpness) noexcept -> double {
    return 0.5 * (1.0 + std::tanh(sharpness * (eta_norm - eta_center)));
}

} // namespace blast::boundary_layer::solver
