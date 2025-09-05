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

  const std::size_t n = values.size();
  std::vector<double> derivatives(n);

  if (d_eta_ <= 0.0) {
    return std::unexpected(CoefficientError("Invalid grid spacing: d_eta must be positive"));
  }
  if (n < 2) {
    return std::unexpected(CoefficientError("Need at least 2 points for derivative"));
  }

  // Cas petits tableaux : schémas simples
  if (n < 5) {
    derivatives[0] = (values[1] - values[0]) / d_eta_;
    for (std::size_t i = 1; i + 1 < n; ++i) {
      derivatives[i] = (values[i + 1] - values[i - 1]) / (2.0 * d_eta_);
    }
    derivatives[n - 1] = (values[n - 1] - values[n - 2]) / d_eta_;
    return derivatives;
  }

  const double dx12 = 12.0 * d_eta_;

  // --- Début : 5 points avant O(h^4) (stencils alignés sur l'ancien code)
  // i = 0
  derivatives[0] =
      (-25.0 * values[0] + 48.0 * values[1] - 36.0 * values[2] + 16.0 * values[3] - 3.0 * values[4]) / dx12;

  // i = 1 : utiliser le stencil de l'ancien code
  // df/dx[1] = (-3*f0 -10*f1 + 18*f2 - 6*f3 + 1*f4)/dx12
  derivatives[1] = (-3.0 * values[0] - 10.0 * values[1] + 18.0 * values[2] - 6.0 * values[3] + 1.0 * values[4]) / dx12;

  // --- Intérieur : 5 points centraux O(h^4)
  for (std::size_t i = 2; i + 2 < n; ++i) {
    derivatives[i] = (values[i - 2] - 8.0 * values[i - 1] + 8.0 * values[i + 1] - values[i + 2]) / dx12;
  }

  // --- Fin : stencils alignés sur l'ancien code
  // i = n-2 : df/dx[n-2] = (3*f[n-6] - 16*f[n-5] + 36*f[n-4] - 48*f[n-3] + 25*f[n-2]) / dx12
  if (n >= 6) {
    derivatives[n - 2] = (3.0 * values[n - 6] - 16.0 * values[n - 5] + 36.0 * values[n - 4]
                           - 48.0 * values[n - 3] + 25.0 * values[n - 2]) / dx12;
  } else {
    // fallback raisonnable (n==5), utiliser le stencil arrière 3 points
    derivatives[n - 2] = (values[n - 2] - values[n - 3]) / d_eta_;
  }

  // i = n-1 : df/dx[n-1] = (3*f[n-5] - 16*f[n-4] + 36*f[n-3] - 48*f[n-2] + 25*f[n-1]) / dx12
  if (n >= 5) {
    derivatives[n - 1] = (3.0 * values[n - 5] - 16.0 * values[n - 4] + 36.0 * values[n - 3]
                           - 48.0 * values[n - 2] + 25.0 * values[n - 1]) / dx12;
  } else {
    derivatives[n - 1] = (values[n - 1] - values[n - 2]) / d_eta_;
  }

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

  const std::size_t n = values.size();
  std::vector<double> second_derivatives(n);

  if (d_eta_ <= 0.0) {
    return std::unexpected(CoefficientError("Invalid grid spacing: d_eta must be positive"));
  }
  if (n < 3) {
    return std::unexpected(CoefficientError("Need at least 3 points for second derivative"));
  }

  // Cas petits tableaux : schéma 3 points classique O(h^2)
  if (n < 5) {
    const double h2 = d_eta_ * d_eta_;
    second_derivatives[0] = (values[0] - 2.0 * values[1] + values[2]) / h2;
    for (std::size_t i = 1; i + 1 < n; ++i) {
      second_derivatives[i] = (values[i - 1] - 2.0 * values[i] + values[i + 1]) / h2;
    }
    second_derivatives[n - 1] = (values[n - 3] - 2.0 * values[n - 2] + values[n - 1]) / h2;
    return second_derivatives;
  }

  const double denom12 = 12.0 * d_eta_ * d_eta_;

  // --- Début : 5 points avant O(h^4)
  // i = 0
  second_derivatives[0] =
      (35.0 * values[0] - 104.0 * values[1] + 114.0 * values[2] - 56.0 * values[3] + 11.0 * values[4]) / denom12;

  // i = 1
  if (n >= 6) {
    second_derivatives[1] =
        (35.0 * values[1] - 104.0 * values[2] + 114.0 * values[3] - 56.0 * values[4] + 11.0 * values[5]) / denom12;
  } else {
    // n == 5 : impossible d'accéder à values[5] → fallback 3 points O(h^2)
    const double h2 = d_eta_ * d_eta_;
    second_derivatives[1] = (values[0] - 2.0 * values[1] + values[2]) / h2;
  }

  // --- Intérieur : 5 points centraux O(h^4)
  for (std::size_t i = 2; i + 2 < n; ++i) {
    second_derivatives[i] =
        (-values[i - 2] + 16.0 * values[i - 1] - 30.0 * values[i] + 16.0 * values[i + 1] - values[i + 2]) / denom12;
  }

  // --- Fin : stencils biaisés O(h^4)
  // i = n-2
  second_derivatives[n - 2] =
      (-1.0 * values[n - 5] + 4.0 * values[n - 4] + 6.0 * values[n - 3] - 20.0 * values[n - 2] +
       11.0 * values[n - 1]) / denom12;

  // i = n-1
  second_derivatives[n - 1] =
      (11.0 * values[n - 5] - 56.0 * values[n - 4] + 114.0 * values[n - 3] - 104.0 * values[n - 2] +
       35.0 * values[n - 1]) / denom12;

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
    // Extraire la ligne j de la matrice c
    std::vector<double> c_row(n_eta);
    for (std::size_t i = 0; i < n_eta; ++i) {
      c_row[i] = solution.c(j, i);
    }

    auto dc_deta_result = eta_derivative_O4(c_row);
    if (!dc_deta_result) {
      return std::unexpected(CoefficientError("Failed to compute dc/deta for species " + std::to_string(j)));
    }
    auto dc_deta = dc_deta_result.value();
    for (std::size_t i = 0; i < n_eta; ++i) {
      derivatives.dc_deta(j, i) = dc_deta[i];
    }

    // Ancien code: d2/deta2 obtenu par dérivation de dc/deta avec le même opérateur
    auto dc_deta2_result = eta_derivative_O4(std::span(dc_deta.data(), dc_deta.size()));
    if (!dc_deta2_result) {
      return std::unexpected(CoefficientError("Failed to compute dc/deta2 (via d/deta of dc/deta) for species " + std::to_string(j)));
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
