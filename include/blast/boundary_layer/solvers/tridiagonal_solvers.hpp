#pragma once
#include "../../core/containers.hpp"
#include "../../core/exceptions.hpp"
#include <vector>
#include <span>
#include <expected>

namespace blast::boundary_layer::solvers {

// Error type for solver operations
class SolverError : public core::BlastException {
public:
    explicit SolverError(std::string_view message,
                        std::source_location location = std::source_location::current())
        : BlastException(std::format("Solver Error: {}", message), location) {}
};

// Solve scalar tridiagonal system for momentum/energy equations
[[nodiscard]] auto solve_momentum_energy_tridiagonal(
    std::span<const double> prev_solution,
    double f_bc,  // Wall boundary coefficient
    double g_bc,  // Wall boundary coefficient  
    double h_bc,  // Wall boundary RHS
    std::span<const double> a_coeffs,
    std::span<const double> b_coeffs,
    std::span<const double> c_coeffs,
    std::span<const double> d_rhs
) -> std::expected<std::vector<double>, SolverError>;

// Solve block tridiagonal system for species equations
[[nodiscard]] auto solve_species_block_tridiagonal(
    const core::Matrix<double>& prev_solution,  // [n_species x n_eta]
    std::span<const double> f_bc,  // Wall BC coefficients per species
    std::span<const double> g_bc,  // Wall BC coefficients per species
    std::span<const double> h_bc,  // Wall BC RHS per species
    const core::Matrix<double>& a_coeffs,  // [n_eta x n_species]
    const core::Matrix<double>& b_coeffs,  // [n_eta x n_species]
    const core::Matrix<double>& c_coeffs,  // [n_eta x n_species]
    const core::Matrix<double>& d_rhs,     // [n_eta x n_species]
    bool has_electrons = false
) -> std::expected<core::Matrix<double>, SolverError>;

// Internal Thomas algorithm implementation
namespace detail {
    void thomas_algorithm(
        std::span<double> lower_diag,
        std::span<double> main_diag,
        std::span<double> upper_diag,
        std::span<double> rhs,
        std::span<double> solution
    );
    
    void block_thomas_algorithm(
        std::vector<core::Matrix<double>>& lower_blocks,
        std::vector<core::Matrix<double>>& main_blocks,
        std::vector<core::Matrix<double>>& upper_blocks,
        core::Matrix<double>& rhs,
        core::Matrix<double>& solution,
        int start_species
    );

}

} // namespace blast::boundary_layer::solvers