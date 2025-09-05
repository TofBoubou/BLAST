#include "blast/boundary_layer/coefficients/derivative_calculator.hpp"
#include <algorithm>
#include <ranges>

namespace blast::boundary_layer::coefficients {

/* auto DerivativeCalculator::eta_derivative_O4(std::span<const double> values) const
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
  derivatives[0] =
      (-25.0 * values[0] + 48.0 * values[1] - 36.0 * values[2] + 16.0 * values[3] - 3.0 * values[4]) / dx12;

  derivatives[1] =
      (-25.0 * values[1] + 48.0 * values[2] - 36.0 * values[3] + 16.0 * values[4] - 3.0 * values[5]) / dx12;

  // Central 5-point stencil O(h^4) for interior points
  for (std::size_t i = 2; i < n - 2; ++i) {
    derivatives[i] = (values[i - 2] - 8.0 * values[i - 1] + 8.0 * values[i + 1] - values[i + 2]) / dx12;
  }

  // Backward 5-point stencil O(h^4) for last 2 points
  derivatives[n - 2] =
      (-1.0 * values[n - 5] + 6.0 * values[n - 4] - 18.0 * values[n - 3] + 10.0 * values[n - 2] + 3.0 * values[n - 1]) /
      dx12;

  derivatives[n - 1] = (3.0 * values[n - 5] - 16.0 * values[n - 4] + 36.0 * values[n - 3] - 48.0 * values[n - 2] +
                        25.0 * values[n - 1]) /
                       dx12;

  return derivatives;
} */

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

  if (n == 2) {
    // Différence finie simple pour 2 points seulement
    derivatives[0] = (values[1] - values[0]) / d_eta_;
    derivatives[1] = derivatives[0];
    return derivatives;
  }

  // Schéma d'ordre 2 (O(h²))
  
  // Premier point : différence finie avant (forward) d'ordre 2
  // f'(x0) = (-3*f(x0) + 4*f(x1) - f(x2)) / (2*h) + O(h²)
  derivatives[0] = (-3.0 * values[0] + 4.0 * values[1] - values[2]) / (2.0 * d_eta_);

  // Points intérieurs : différence finie centrée d'ordre 2
  // f'(xi) = (f(xi+1) - f(xi-1)) / (2*h) + O(h²)
  for (std::size_t i = 1; i < n - 1; ++i) {
    derivatives[i] = (values[i + 1] - values[i - 1]) / (2.0 * d_eta_);
  }

  // Dernier point : différence finie arrière (backward) d'ordre 2
  // f'(xn) = (f(xn-2) - 4*f(xn-1) + 3*f(xn)) / (2*h) + O(h²)
  derivatives[n - 1] = (values[n - 3] - 4.0 * values[n - 2] + 3.0 * values[n - 1]) / (2.0 * d_eta_);

  return derivatives;
}

/* auto DerivativeCalculator::eta_second_derivative_O4(std::span<const double> values) const
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
  second_derivatives[0] =
      (35.0 * values[0] - 104.0 * values[1] + 114.0 * values[2] - 56.0 * values[3] + 11.0 * values[4]) / dx2_12;

  second_derivatives[1] =
      (35.0 * values[1] - 104.0 * values[2] + 114.0 * values[3] - 56.0 * values[4] + 11.0 * values[5]) / dx2_12;

  // Central 5-point stencil O(h^4) for interior points
  for (std::size_t i = 2; i < n - 2; ++i) {
    second_derivatives[i] =
        (-values[i - 2] + 16.0 * values[i - 1] - 30.0 * values[i] + 16.0 * values[i + 1] - values[i + 2]) / dx2_12;
  }

  // Backward 5-point stencils O(h^4) for last 2 points
  second_derivatives[n - 2] =
      (-1.0 * values[n - 5] + 4.0 * values[n - 4] + 6.0 * values[n - 3] - 20.0 * values[n - 2] + 11.0 * values[n - 1]) /
      dx2_12;

  second_derivatives[n - 1] = (11.0 * values[n - 5] - 56.0 * values[n - 4] + 114.0 * values[n - 3] -
                               104.0 * values[n - 2] + 35.0 * values[n - 1]) /
                              dx2_12;

  return second_derivatives;
} */

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

  const double d_eta_sq = d_eta_ * d_eta_;

  if (n == 3) {
    // Cas spécial : seulement 3 points, utilisation du schéma centré partout
    for (std::size_t i = 0; i < n; ++i) {
      second_derivatives[i] = (values[0] - 2.0 * values[1] + values[2]) / d_eta_sq;
    }
    return second_derivatives;
  }

  // Schéma d'ordre 2 (O(h²))

  // Premier point : différence finie avant (forward) d'ordre 2 avec 4 points
  // f''(x0) = (2*f(x0) - 5*f(x1) + 4*f(x2) - f(x3)) / h² + O(h²)
  if (n >= 4) {
    second_derivatives[0] = (2.0 * values[0] - 5.0 * values[1] + 4.0 * values[2] - values[3]) / d_eta_sq;
  } else {
    // Si on n'a que 3 points, utiliser le schéma centré
    second_derivatives[0] = (values[0] - 2.0 * values[1] + values[2]) / d_eta_sq;
  }

  // Points intérieurs : différence finie centrée d'ordre 2 (schéma classique)
  // f''(xi) = (f(xi-1) - 2*f(xi) + f(xi+1)) / h² + O(h²)
  for (std::size_t i = 1; i < n - 1; ++i) {
    second_derivatives[i] = (values[i - 1] - 2.0 * values[i] + values[i + 1]) / d_eta_sq;
  }

  // Dernier point : différence finie arrière (backward) d'ordre 2 avec 4 points
  // f''(xn) = (-f(xn-3) + 4*f(xn-2) - 5*f(xn-1) + 2*f(xn)) / h² + O(h²)
  if (n >= 4) {
    second_derivatives[n - 1] = (-values[n - 4] + 4.0 * values[n - 3] - 5.0 * values[n - 2] + 2.0 * values[n - 1]) / d_eta_sq;
  } else {
    // Si on n'a que 3 points, utiliser le schéma centré
    second_derivatives[n - 1] = (values[n - 3] - 2.0 * values[n - 2] + values[n - 1]) / d_eta_sq;
  }

  return second_derivatives;
}

auto DerivativeCalculator::compute_all_derivatives(const equations::SolutionState& solution) const
    -> std::expected<UnifiedDerivativeState, CoefficientError> {

  const auto n_eta = solution.F.size();
  const auto n_species = solution.c.rows();

  UnifiedDerivativeState derivatives(n_eta, n_species);


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


  for (std::size_t j = 0; j < n_species; ++j) {
    // Map direct de la ligne j (RowMajor) en span sans copie
    auto c_row = solution.c.eigen().row(static_cast<Eigen::Index>(j));
    auto dc_deta_result = eta_derivative_O4(std::span<const double>(c_row.data(), n_eta));
    if (!dc_deta_result) {
      return std::unexpected(CoefficientError("Failed to compute dc/deta for species " + std::to_string(j)));
    }
    const auto& dc_deta = dc_deta_result.value();
    Eigen::Map<const Eigen::RowVectorXd> dc_map(dc_deta.data(), static_cast<Eigen::Index>(n_eta));
    derivatives.dc_deta.eigen().row(static_cast<Eigen::Index>(j)) = dc_map;

    auto dc_deta2_result = eta_second_derivative_O4(std::span<const double>(c_row.data(), n_eta));
    if (!dc_deta2_result) {
      return std::unexpected(CoefficientError("Failed to compute dc/deta2 for species " + std::to_string(j)));
    }
    const auto& dc_deta2 = dc_deta2_result.value();
    Eigen::Map<const Eigen::RowVectorXd> dc2_map(dc_deta2.data(), static_cast<Eigen::Index>(n_eta));
    derivatives.dc_deta2.eigen().row(static_cast<Eigen::Index>(j)) = dc2_map;
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
