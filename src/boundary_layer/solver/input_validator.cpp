#include "blast/boundary_layer/solver/input_validator.hpp"
#include <algorithm>
#include <cmath>
#include <format>

namespace blast::boundary_layer::solver {

InputValidator::InputValidator(const grid::BoundaryLayerGrid& grid,
                              const thermophysics::MixtureInterface& mixture) noexcept
    : grid_(grid), mixture_(mixture) {}

auto InputValidator::validate_station_parameters(int station, double xi) const noexcept
    -> std::expected<void, SolverError> {

    // Validate station number
    if (station < 0) {
        return std::unexpected(InitializationError(
            std::format("Invalid station number: {} (must be non-negative)", station)));
    }

    // Validate xi coordinate
    if (!std::isfinite(xi) || xi < 0.0) {
        return std::unexpected(InitializationError(
            std::format("Invalid xi coordinate: {} (must be finite and non-negative)", xi)));
    }

    return {};
}

auto InputValidator::validate_solution_state(const equations::SolutionState& solution) const noexcept
    -> std::expected<void, SolverError> {

    const auto expected_n_eta = grid_.n_eta();
    const auto expected_n_species = mixture_.n_species();

    // Validate field dimensions
    if (solution.F.size() != expected_n_eta || 
        solution.T.size() != expected_n_eta ||
        solution.g.size() != expected_n_eta ||
        solution.V.size() != expected_n_eta) {
        return std::unexpected(InitializationError(
            std::format("Solution field dimensions mismatch: expected {} eta points", expected_n_eta)));
    }

    // Validate species matrix dimensions
    if (solution.c.rows() != expected_n_species || solution.c.cols() != expected_n_eta) {
        return std::unexpected(InitializationError(
            std::format("Solution species matrix dimensions mismatch: expected {}x{}", 
                       expected_n_species, expected_n_eta)));
    }

    return {};
}

auto InputValidator::validate_xi_coordinate(int station, double xi) const noexcept
    -> std::expected<void, SolverError> {

    // For downstream stations, check consistency with grid
    if (station > 0) {
        const auto& xi_coords = grid_.xi_coordinates();
        if (station >= static_cast<int>(xi_coords.size())) {
            return std::unexpected(InitializationError(
                std::format("Station {} exceeds available xi coordinates (max: {})", 
                           station, xi_coords.size() - 1)));
        }

        // Allow some tolerance for floating point comparison
        const double expected_xi = xi_coords[station];
        if (std::abs(xi - expected_xi) > 1e-10) {
            return std::unexpected(InitializationError(
                std::format("Xi coordinate mismatch for station {}: provided {}, expected {}", 
                           station, xi, expected_xi)));
        }
    }

    return {};
}

auto InputValidator::validate_radiative_equilibrium_inputs(
    double emissivity,
    double T_infinity,
    double q_wall
) noexcept -> std::expected<void, SolverError> {

    if (emissivity <= 0.0) {
        return std::unexpected(InitializationError(
            std::format("Invalid emissivity: {} (must be > 0)", emissivity)));
    }

    if (T_infinity <= 0.0) {
        return std::unexpected(InitializationError(
            std::format("Invalid environment temperature: {} K (must be > 0)", T_infinity)));
    }

    if (!std::isfinite(q_wall)) {
        return std::unexpected(InitializationError(
            std::format("Invalid wall heat flux: {} (not finite)", q_wall)));
    }

    return {};
}

auto InputValidator::validate_solution_finite(const equations::SolutionState& solution) noexcept
    -> std::expected<void, SolverError> {

    // Check F field
    auto check_field_finite = [](const auto& field, const std::string& name) -> std::expected<void, SolverError> {
        for (std::size_t i = 0; i < field.size(); ++i) {
            if (!std::isfinite(field[i])) {
                return std::unexpected(NumericError(
                    std::format("Non-finite value detected in {} field at index {}: {}", name, i, field[i])));
            }
        }
        return {};
    };

    if (auto result = check_field_finite(solution.F, "F"); !result) return result;
    if (auto result = check_field_finite(solution.g, "g"); !result) return result;
    if (auto result = check_field_finite(solution.T, "T"); !result) return result;
    if (auto result = check_field_finite(solution.V, "V"); !result) return result;

    // Check species matrix
    for (std::size_t i = 0; i < solution.c.rows(); ++i) {
        for (std::size_t j = 0; j < solution.c.cols(); ++j) {
            if (!std::isfinite(solution.c(i, j))) {
                return std::unexpected(NumericError(
                    std::format("Non-finite value detected in species matrix at ({}, {}): {}", 
                               i, j, solution.c(i, j))));
            }
        }
    }

    return {};
}

} // namespace blast::boundary_layer::solver