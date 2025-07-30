#include "blast/boundary_layer/coefficients/derivative_calculator.hpp"
#include <algorithm>
#include <ranges>

namespace blast::boundary_layer::coefficients {

auto DerivativeCalculator::eta_derivative_O4(std::span<const double> values) const 
    -> std::expected<std::vector<double>, CoefficientError> {
    
    const auto n = values.size();
    std::vector<double> derivatives(n);

    if (d_eta_ <= 0.0) {
        return std::unexpected(CoefficientError("Invalid grid spacing: d_eta must be positive"));
    }

    if (n < 2) {
        return std::unexpected(CoefficientError("Need at least 2 points for derivative"));
    }

    if (n < 5) {
        // Simple finite differences for small arrays
        derivatives[0] = (values[1] - values[0]) / d_eta_;
        for (std::size_t i = 1; i < n - 1; ++i) {
            derivatives[i] = (values[i + 1] - values[i - 1]) / (2.0 * d_eta_);
        }
        derivatives[n - 1] = (values[n - 1] - values[n - 2]) / d_eta_;
        return derivatives;
    }

    const double dx12 = 12.0 * d_eta_;

    // Forward 5-point stencil O(h^4) for first 2 points
    derivatives[0] = (-25.0 * values[0] + 48.0 * values[1] - 36.0 * values[2] + 
                      16.0 * values[3] - 3.0 * values[4]) / dx12;

    derivatives[1] = (-25.0 * values[1] + 48.0 * values[2] - 36.0 * values[3] + 
                      16.0 * values[4] - 3.0 * values[5]) / dx12;

    // Central 5-point stencil O(h^4) for interior points
    for (std::size_t i = 2; i < n - 2; ++i) {
        derivatives[i] = (values[i - 2] - 8.0 * values[i - 1] + 8.0 * values[i + 1] - values[i + 2]) / dx12;
    }

    // Backward 5-point stencil O(h^4) for last 2 points
    derivatives[n - 2] = (-1.0 * values[n - 5] + 6.0 * values[n - 4] - 18.0 * values[n - 3] + 
                          10.0 * values[n - 2] + 3.0 * values[n - 1]) / dx12;

    derivatives[n - 1] = (3.0 * values[n - 5] - 16.0 * values[n - 4] + 36.0 * values[n - 3] - 
                          48.0 * values[n - 2] + 25.0 * values[n - 1]) / dx12;

    return derivatives;
}

auto DerivativeCalculator::eta_second_derivative_O4(std::span<const double> values) const
    -> std::expected<std::vector<double>, CoefficientError> {
    
    const auto n = values.size();
    std::vector<double> second_derivatives(n);

    if (d_eta_ <= 0.0) {
        return std::unexpected(CoefficientError("Invalid grid spacing: d_eta must be positive"));
    }

    if (n < 3) {
        return std::unexpected(CoefficientError("Need at least 3 points for second derivative"));
    }

    if (n < 5) {
        // Simple 3-point stencil for small arrays
        const double d_eta_sq = d_eta_ * d_eta_;
        second_derivatives[0] = (values[0] - 2.0 * values[1] + values[2]) / d_eta_sq;
        for (std::size_t i = 1; i < n - 1; ++i) {
            second_derivatives[i] = (values[i - 1] - 2.0 * values[i] + values[i + 1]) / d_eta_sq;
        }
        second_derivatives[n - 1] = (values[n - 3] - 2.0 * values[n - 2] + values[n - 1]) / d_eta_sq;
        return second_derivatives;
    }

    const double dx2_12 = 12.0 * d_eta_ * d_eta_;

    // Forward 5-point stencil O(h^4) for first 2 points
    second_derivatives[0] = (35.0 * values[0] - 104.0 * values[1] + 114.0 * values[2] - 
                             56.0 * values[3] + 11.0 * values[4]) / dx2_12;

    second_derivatives[1] = (35.0 * values[1] - 104.0 * values[2] + 114.0 * values[3] - 
                             56.0 * values[4] + 11.0 * values[5]) / dx2_12;

    // Central 5-point stencil O(h^4) for interior points
    for (std::size_t i = 2; i < n - 2; ++i) {
        second_derivatives[i] = (-values[i - 2] + 16.0 * values[i - 1] - 30.0 * values[i] + 
                                 16.0 * values[i + 1] - values[i + 2]) / dx2_12;
    }

    // Backward 5-point stencils O(h^4) for last 2 points
    second_derivatives[n - 2] = (-1.0 * values[n - 5] + 4.0 * values[n - 4] + 6.0 * values[n - 3] - 
                                 20.0 * values[n - 2] + 11.0 * values[n - 1]) / dx2_12;

    second_derivatives[n - 1] = (11.0 * values[n - 5] - 56.0 * values[n - 4] + 114.0 * values[n - 3] - 
                                 104.0 * values[n - 2] + 35.0 * values[n - 1]) / dx2_12;

    return second_derivatives;
}

auto DerivativeCalculator::compute_all_derivatives(const equations::SolutionState& solution) const
    -> std::expected<UnifiedDerivativeState, CoefficientError> {
    
    const auto n_eta = solution.F.size();
    const auto n_species = solution.c.rows();
    
    UnifiedDerivativeState derivatives(n_eta, n_species);
    
    // ===== DÉRIVÉES SCALAIRES =====
    
    // dF/deta
    auto dF_result = eta_derivative_O4(solution.F);
    if (!dF_result) {
        return std::unexpected(CoefficientError("Failed to compute dF/deta: " + dF_result.error().message()));
    }
    derivatives.dF_deta = std::move(dF_result.value());
    
    // dg/deta
    auto dg_result = eta_derivative_O4(solution.g);
    if (!dg_result) {
        return std::unexpected(CoefficientError("Failed to compute dg/deta: " + dg_result.error().message()));
    }
    derivatives.dg_deta = std::move(dg_result.value());
    
    // dV/deta
    auto dV_result = eta_derivative_O4(solution.V);
    if (!dV_result) {
        return std::unexpected(CoefficientError("Failed to compute dV/deta: " + dV_result.error().message()));
    }
    derivatives.dV_deta = std::move(dV_result.value());
    
    // dT/deta
    auto dT_result = eta_derivative_O4(solution.T);
    if (!dT_result) {
        return std::unexpected(CoefficientError("Failed to compute dT/deta: " + dT_result.error().message()));
    }
    derivatives.dT_deta = std::move(dT_result.value());
    
    // ===== DÉRIVÉES DES ESPÈCES =====
    
    for (std::size_t j = 0; j < n_species; ++j) {
        // Extraire la ligne j de la matrice c
        std::vector<double> c_row(n_eta);
        for (std::size_t i = 0; i < n_eta; ++i) {
            c_row[i] = solution.c(j, i);
        }
        
        // Première dérivée dc/deta
        auto dc_deta_result = eta_derivative_O4(c_row);
        if (!dc_deta_result) {
            return std::unexpected(CoefficientError("Failed to compute dc/deta for species " + std::to_string(j)));
        }
        auto dc_deta = dc_deta_result.value();
        for (std::size_t i = 0; i < n_eta; ++i) {
            derivatives.dc_deta(j, i) = dc_deta[i];
        }
        
        // Seconde dérivée dc/deta2
        auto dc_deta2_result = eta_second_derivative_O4(c_row);
        if (!dc_deta2_result) {
            return std::unexpected(CoefficientError("Failed to compute dc/deta2 for species " + std::to_string(j)));
        }
        auto dc_deta2 = dc_deta2_result.value();
        for (std::size_t i = 0; i < n_eta; ++i) {
            derivatives.dc_deta2(j, i) = dc_deta2[i];
        }
    }
    
    return derivatives;
}

auto DerivativeCalculator::compute_single_derivative(std::span<const double> field) const
    -> std::expected<std::vector<double>, CoefficientError> {
    return eta_derivative_O4(field);
}

auto DerivativeCalculator::compute_single_second_derivative(std::span<const double> field) const  
    -> std::expected<std::vector<double>, CoefficientError> {
    return eta_second_derivative_O4(field);
}

} // namespace blast::boundary_layer::coefficients