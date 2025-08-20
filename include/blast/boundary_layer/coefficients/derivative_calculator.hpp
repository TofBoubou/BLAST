#pragma once
#include "../../core/containers.hpp"
#include "../../core/exceptions.hpp"
#include "../equations/equation_types.hpp"
#include "coefficient_types.hpp"
#include <expected>
#include <span>

namespace blast::boundary_layer::coefficients {

struct UnifiedDerivativeState {
  std::vector<double> dF_deta; // Momentum derivatives
  std::vector<double> dg_deta; // Energy derivatives
  std::vector<double> dV_deta; // Continuity derivatives
  std::vector<double> dT_deta; // Temperature derivatives

  core::Matrix<double> dc_deta;  // First derivatives [n_species x n_eta]
  core::Matrix<double> dc_deta2; // Second derivatives [n_species x n_eta]

  UnifiedDerivativeState(std::size_t n_eta, std::size_t n_species)
      : dF_deta(n_eta), dg_deta(n_eta), dV_deta(n_eta), dT_deta(n_eta), dc_deta(n_species, n_eta),
        dc_deta2(n_species, n_eta) {}
};

class DerivativeCalculator {
private:
  const double d_eta_;

  [[nodiscard]] auto
  eta_derivative_O4(std::span<const double> values) const -> std::expected<std::vector<double>, CoefficientError>;

  [[nodiscard]] auto eta_second_derivative_O4(std::span<const double> values) const
      -> std::expected<std::vector<double>, CoefficientError>;

public:
  explicit DerivativeCalculator(double d_eta) : d_eta_(d_eta) {}

  [[nodiscard]] auto compute_all_derivatives(const equations::SolutionState& solution) const
      -> std::expected<UnifiedDerivativeState, CoefficientError>;

  [[nodiscard]] auto compute_single_derivative(std::span<const double> field) const
      -> std::expected<std::vector<double>, CoefficientError>;

  [[nodiscard]] auto compute_single_second_derivative(std::span<const double> field) const
      -> std::expected<std::vector<double>, CoefficientError>;
};

} // namespace blast::boundary_layer::coefficients