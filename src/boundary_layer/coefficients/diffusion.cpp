#include "blast/boundary_layer/coefficients/diffusion.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

namespace blast::boundary_layer::coefficients::diffusion {

namespace {

// Compute derivative transformation factor
constexpr auto compute_derivative_factor(int station, double xi, const conditions::BoundaryConditions& bc,
                                         const io::SimulationConfig& sim_config) noexcept -> double {

  if (station != 0) {
    return bc.ue() * bc.r_body() / std::sqrt(2.0 * xi);
  }

  // Stagnation point limiting solution
  switch (sim_config.body_type) {
  case io::SimulationConfig::BodyType::Axisymmetric:
    return std::sqrt(2.0 * bc.d_ue_dx() / (bc.rho_e() * bc.mu_e()));
  case io::SimulationConfig::BodyType::TwoD:
    return std::sqrt(bc.d_ue_dx() / (bc.rho_e() * bc.mu_e()));
  case io::SimulationConfig::BodyType::Cone:
  case io::SimulationConfig::BodyType::FlatPlate:
    return 1.0 / std::sqrt(bc.rho_e() * bc.mu_e());
  }
  return 1.0;
}

// Compute average diffusion coefficient D_im
inline auto compute_average_diffusion_coefficient(std::span<const double> x, const Eigen::MatrixXd& D_bin,
                                                  std::size_t species_idx) noexcept -> double {
  const auto n = static_cast<Eigen::Index>(x.size());
  const Eigen::Map<const Eigen::RowVectorXd> x_row(x.data(), n);
  const auto D_row = D_bin.row(static_cast<Eigen::Index>(species_idx));
  const double sum = (x_row.array() / D_row.array()).sum();
  return 1.0 / sum;
}

// Stefan-Maxwell flux calculation at a single eta point using reusable workspace
auto calculate_stefan_maxwell_at_point(std::span<const double> c, std::span<const double> dc_deta,
                                       std::span<const double> x, double rho, double P,
                                       const core::Matrix<double>& D_bin_local, double MW, double dMW_deta,
                                       std::span<const double> tdr_term, double der_fact,
                                       std::span<const double> MW_species, std::span<const double> charges,
                                       Eigen::VectorXd& J_vec,
                                       Eigen::VectorXd& driving_force, Eigen::VectorXd& diffusion_flux,
                                       Eigen::VectorXd& mass_weight, Eigen::VectorXd& electric_mobility,
                                       Eigen::MatrixXd& A, Eigen::PartialPivLU<Eigen::MatrixXd>& lu) -> void {

  const auto n_species = static_cast<std::size_t>(c.size());
  const double full_der_fact = der_fact * rho;
  const double sum_c = std::accumulate(c.begin(), c.end(), 0.0);

  // Use Eigen through the wrapper
  auto& D_bin = D_bin_local.eigen();

  // Ensure correct sizes for workspace (no heap allocation if already sized)
  if (static_cast<std::size_t>(driving_force.size()) != n_species) {
    driving_force.resize(static_cast<Eigen::Index>(n_species));
    diffusion_flux.resize(static_cast<Eigen::Index>(n_species));
    mass_weight.resize(static_cast<Eigen::Index>(n_species));
    electric_mobility.resize(static_cast<Eigen::Index>(n_species));
    J_vec.resize(static_cast<Eigen::Index>(n_species));
  }
  if (A.rows() != static_cast<Eigen::Index>(n_species) || A.cols() != static_cast<Eigen::Index>(n_species)) {
    A.resize(static_cast<Eigen::Index>(n_species), static_cast<Eigen::Index>(n_species));
  }
  A.setIdentity();

  // Build Stefan-Maxwell system
  for (std::size_t i = 0; i < n_species; ++i) {
    const double MW_i = MW_species[i];

    // Driving force: concentration gradient + molecular weight gradient + thermal diffusion
    driving_force[static_cast<Eigen::Index>(i)] = MW / MW_i * dc_deta[i] + c[i] / MW_i * dMW_deta;
    if (!tdr_term.empty()) {
      driving_force[static_cast<Eigen::Index>(i)] += tdr_term[i];
    }

    // Average diffusion coefficient
    const double D_im = compute_average_diffusion_coefficient(x, D_bin, i);

    // Flux and mass weight terms
    diffusion_flux[static_cast<Eigen::Index>(i)] = -rho * MW_i / MW * D_im * driving_force[static_cast<Eigen::Index>(i)] * full_der_fact;
    mass_weight[static_cast<Eigen::Index>(i)] = c[i] * D_im * MW;

    // Build matrix row
    for (std::size_t j = 0; j < n_species; ++j) {
      A(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) -=
          mass_weight[static_cast<Eigen::Index>(i)] /
          (MW_species[j] * D_bin(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)));
    }
  }

  // Apply mass conservation constraint (last row)
  A.row(static_cast<Eigen::Index>(n_species - 1)).setConstant(c[n_species - 1] / sum_c);
  diffusion_flux[static_cast<Eigen::Index>(n_species - 1)] = 0.0;

  // Solve linear system using reusable LU object
  lu.compute(A);
  J_vec = lu.solve(diffusion_flux);

  // Apply ambipolar electric field correction for ionized mixtures
  // Apply ambipolar electric field correction if charges provided (ionized mixtures)
  if (!charges.empty()) {
    // Reuse electric_mobility; MW_species is provided by caller
    for (std::size_t i = 0; i < n_species; ++i) {
      const double D_im = compute_average_diffusion_coefficient(x, D_bin, i);
      electric_mobility[static_cast<Eigen::Index>(i)] = rho / P * D_im * MW_species[i] / MW *
                                                       charges[i] * rho * c[i];
    }

    const Eigen::Map<const Eigen::VectorXd> charges_vec(charges.data(), static_cast<Eigen::Index>(charges.size())) ;
    const double charge_flux = charges_vec.dot(diffusion_flux);
    const double current_flux = charges_vec.dot(J_vec);
    const double charge_mobility = charges_vec.dot(electric_mobility);

    if (std::abs(charge_mobility) > 1e-30) {
      const double E_field = -(charge_flux + current_flux) / charge_mobility;
      J_vec.noalias() += E_field * electric_mobility;
    }
  }

  // Apply mass conservation correction in-place (vectorized)
  const double sum_J = J_vec.sum();
  const Eigen::Map<const Eigen::VectorXd> c_vec(c.data(), static_cast<Eigen::Index>(n_species));
  J_vec.noalias() -= (sum_J / sum_c) * c_vec;
}

} // anonymous namespace

auto compute_stefan_maxwell_fluxes(const CoefficientInputs& inputs,
                                   CoefficientSet& coeffs, // Changed from const to allow modification
                                   const conditions::BoundaryConditions& bc, const XiDerivatives& xi_der,
                                   const io::SimulationConfig& sim_config,
                                   const thermophysics::MixtureInterface& mixture,
                                   double d_eta) -> std::expected<void, CoefficientError> {

  const auto n_species = inputs.c.rows();
  const auto n_eta = inputs.T.size();

  // Access J and dJ_deta directly (no const_cast needed)
  auto& J = coeffs.diffusion.J;
  auto& dJ_deta = coeffs.diffusion.dJ_deta;

  // Ensure matrices are properly sized before use
  J = core::Matrix<double>(n_species, n_eta);
  dJ_deta = core::Matrix<double>(n_species, n_eta);
  J.setZero();
  dJ_deta.setZero();

  // Single species case - no diffusion
  if (n_species == 1) {
    J.setZero();
    dJ_deta.setZero();
    return {};
  }

  // Verify we're using Stefan-Maxwell
  if (sim_config.diffusion_type != io::SimulationConfig::DiffusionType::StefanMaxwell) {
    return std::unexpected(CoefficientError("Only Stefan-Maxwell diffusion is implemented"));
  }

  // Get station and compute derivative factor
  const int station = xi_der.station();
  const double der_fact = compute_derivative_factor(station, inputs.xi, bc, sim_config) * coeffs.transport.K_bl;

  // Pre-allocate work arrays (std::vector + Eigen workspace) to avoid per-η allocations
  std::vector<double> c_local(n_species);
  std::vector<double> dc_deta_local(n_species);
  std::vector<double> x_local(n_species);
  core::Matrix<double> D_bin_local(n_species, n_species);
  // Cache species molecular weights once
  std::vector<double> MW_species(n_species);
  for (std::size_t j = 0; j < n_species; ++j) {
    MW_species[j] = mixture.species_molecular_weight(j);
  }
  // Cache charges once if ionized mixture
  std::span<const double> charges_cached;
  if (mixture.has_electrons()) {
    charges_cached = mixture.species_charges();
  }
  // Reusable Eigen objects
  Eigen::VectorXd driving_force(static_cast<Eigen::Index>(n_species));
  Eigen::VectorXd diffusion_flux(static_cast<Eigen::Index>(n_species));
  Eigen::VectorXd mass_weight(static_cast<Eigen::Index>(n_species));
  Eigen::VectorXd electric_mobility(static_cast<Eigen::Index>(n_species));
  Eigen::VectorXd J_vec(static_cast<Eigen::Index>(n_species));
  Eigen::MatrixXd A(static_cast<Eigen::Index>(n_species), static_cast<Eigen::Index>(n_species));
  Eigen::PartialPivLU<Eigen::MatrixXd> lu;

  // Process each eta point
  for (std::size_t i = 0; i < n_eta; ++i) {
    // Extract local values
    for (std::size_t j = 0; j < n_species; ++j) {
      c_local[j] = inputs.c(j, i);
      dc_deta_local[j] = inputs.dc_deta(j, i);
      x_local[j] = c_local[j] * coeffs.thermodynamic.MW[i] / MW_species[j];
    }

    // Extract binary diffusion coefficients for this eta point
    // Vectorized block copy of Dij for this eta
    D_bin_local.eigen() = coeffs.diffusion.Dij_bin.eigen().block(static_cast<Eigen::Index>(i * n_species), 0,
                                                                 static_cast<Eigen::Index>(n_species),
                                                                 static_cast<Eigen::Index>(n_species));

    // Get thermal diffusion terms if enabled
    std::span<const double> tdr_span;
    if (sim_config.consider_thermal_diffusion && i < coeffs.thermal_diffusion.tdr_term.rows()) {
      auto tdr_row = coeffs.thermal_diffusion.tdr_term.eigen().row(i);
      tdr_span = std::span(tdr_row.data(), n_species);
    }
    // Charges span (cached at function scope)
    std::span<const double> charges_span = charges_cached;

    // Calculate fluxes
    try {
      calculate_stefan_maxwell_at_point(c_local, dc_deta_local, x_local, coeffs.thermodynamic.rho[i], bc.P_e(),
                                        D_bin_local, coeffs.thermodynamic.MW[i],
                                        coeffs.thermodynamic.d_MW_deta[i], tdr_span, der_fact,
                                        std::span<const double>(MW_species.data(), MW_species.size()), charges_span,
                                        J_vec,
                                        driving_force, diffusion_flux, mass_weight, electric_mobility, A, lu);

      // Store results from reusable J_vec
      for (std::size_t j = 0; j < n_species; ++j) {
        J(j, i) = J_vec[static_cast<Eigen::Index>(j)];
      }
    } catch (const std::exception& e) {
      return std::unexpected(
          CoefficientError(std::format("Stefan-Maxwell calculation failed at eta[{}]: {}", i, e.what())));
    }
  }

  // Compute eta derivatives of fluxes: map each row directly (RowMajor) and avoid input copies
  for (std::size_t j = 0; j < n_species; ++j) {
    auto J_row = J.eigen().row(static_cast<Eigen::Index>(j));
    auto dJ_result = derivatives::compute_eta_derivative(std::span<const double>(J_row.data(), n_eta), d_eta);
    if (!dJ_result) {
      return std::unexpected(CoefficientError("Failed to compute diffusion flux derivative"));
    }
    const auto& dJ = dJ_result.value();
    Eigen::Map<const Eigen::RowVectorXd> dJ_map(dJ.data(), static_cast<Eigen::Index>(n_eta));
    dJ_deta.eigen().row(static_cast<Eigen::Index>(j)) = dJ_map;
  }

  return {};
}

} // namespace blast::boundary_layer::coefficients::diffusion
