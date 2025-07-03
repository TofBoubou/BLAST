#include "blast/boundary_layer/coefficients/diffusion.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

namespace blast::boundary_layer::coefficients::diffusion {

namespace {

// Compute derivative transformation factor
constexpr auto compute_derivative_factor(
    int station,
    double xi,
    const conditions::BoundaryConditions& bc,
    const io::SimulationConfig& sim_config
) noexcept -> double {
    
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

// Stefan-Maxwell flux calculation for single eta point
auto calculate_stefan_maxwell_at_point(
    std::span<const double> c,
    std::span<const double> dc_deta,
    std::span<const double> x,
    double rho,
    double P,
    const core::Matrix<double>& D_bin_local,
    double MW,
    double dMW_deta,
    std::span<const double> tdr_term,
    double der_fact,
    double K_bl,
    const thermophysics::MixtureInterface& mixture
) -> std::vector<double> {
    
    const auto n_species = c.size();
    const double full_der_fact = der_fact * rho * K_bl;
    
    double sum_c = std::accumulate(c.begin(), c.end(), 0.0);
    
    // Use Eigen through the wrapper
    auto& D_bin_eigen = D_bin_local.eigen();
    Eigen::VectorXd d(n_species), p1(n_species), p2(n_species);
    Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n_species, n_species);
    
    // Build system
    for (std::size_t i = 0; i < n_species; ++i) {
        // Driving force
        d[i] = MW / mixture.species_molecular_weight(i) * dc_deta[i] + 
               c[i] / mixture.species_molecular_weight(i) * dMW_deta;
        
        if (!tdr_term.empty()) {
            d[i] += tdr_term[i];
        }
        
        // Average diffusion coefficient
        double sum = 0.0;
        for (std::size_t j = 0; j < n_species; ++j) {
            sum += x[j] / D_bin_eigen(i, j);
        }
        double D_im = 1.0 / sum;
        
        p1[i] = -rho * mixture.species_molecular_weight(i) / MW * D_im * d[i] * full_der_fact;
        p2[i] = c[i] * D_im * MW;
        
        // Build matrix
        for (std::size_t j = 0; j < n_species; ++j) {
            A(i, j) -= p2[i] / (mixture.species_molecular_weight(j) * D_bin_eigen(i, j));
        }
    }
    
    // Mass conservation constraint
    A.row(n_species - 1).setConstant(c[n_species - 1] / sum_c);
    p1[n_species - 1] = 0.0;
    
    // Solve
    Eigen::VectorXd J_vec = A.lu().solve(p1);
    
    // Ambipolar electric field
    if (mixture.has_electrons()) {
        auto charges = mixture.species_charges();
        Eigen::VectorXd p3(n_species);
        
        for (std::size_t i = 0; i < n_species; ++i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < n_species; ++j) {
                sum += x[j] / D_bin_eigen(i, j);
            }
            double D_im = 1.0 / sum;
            p3[i] = rho / P * D_im * mixture.species_molecular_weight(i) / MW * 
                    charges[i] * rho * c[i];
        }
        
        double t1 = 0.0, t2 = 0.0, t3 = 0.0;
        for (std::size_t i = 0; i < n_species; ++i) {
            t1 += charges[i] * p1[i];
            t2 += charges[i] * J_vec[i];
            t3 += charges[i] * p3[i];
        }
        
        if (std::abs(t3) > 1e-30) {
            double E = -(t1 + t2) / t3;
            J_vec += E * p3;
        }
    }
    
    // Mass conservation correction
    double sum_J = J_vec.sum();
    std::vector<double> J_result(n_species);
    for (std::size_t i = 0; i < n_species; ++i) {
        J_result[i] = J_vec[i] - c[i] / sum_c * sum_J;
    }
    
    return J_result;
}

} // anonymous namespace

auto compute_stefan_maxwell_fluxes(
    const CoefficientInputs& inputs,
    const CoefficientSet& coeffs,
    const conditions::BoundaryConditions& bc,
    const XiDerivatives& xi_der,
    const io::SimulationConfig& sim_config,
    const thermophysics::MixtureInterface& mixture,
    double d_eta
) -> std::expected<void, CoefficientError> {
    
    const auto n_species = inputs.c.rows();
    const auto n_eta = inputs.T.size();
    
    // Access J and dJ_deta from DiffusionCoefficients (they're mutable through const_cast)
    auto& J = const_cast<core::Matrix<double>&>(coeffs.diffusion.J);
    auto& dJ_deta = const_cast<core::Matrix<double>&>(coeffs.diffusion.dJ_deta);
    
    // Single species case
    if (n_species == 1) {
        J.setZero();
        dJ_deta.setZero();
        return {};
    }
    
    // Verify we're using Stefan-Maxwell
    if (sim_config.diffusion_type != io::SimulationConfig::DiffusionType::StefanMaxwell) {
        return std::unexpected(CoefficientError("Only Stefan-Maxwell diffusion is implemented"));
    }
    
    // Get station from xi_der
    const int station = xi_der.station();
    
    // Compute derivative factor
    const double der_fact = compute_derivative_factor(station, inputs.xi, bc, sim_config);
    const double K_bl = 1.0; // TODO: Get from proper source
    
    // Process each eta point
    for (std::size_t i = 0; i < n_eta; ++i) {
        // Extract local values
        std::vector<double> c_local(n_species);
        std::vector<double> dc_deta_local(n_species);
        std::vector<double> x_local(n_species);
        
        for (std::size_t j = 0; j < n_species; ++j) {
            c_local[j] = inputs.c(j, i);
            dc_deta_local[j] = inputs.dc_deta(j, i);
            x_local[j] = c_local[j] * coeffs.thermodynamic.MW[i] / mixture.species_molecular_weight(j);
        }
        
        // Extract D_bin for this station
        core::Matrix<double> D_bin_local(n_species, n_species);
        for (std::size_t j = 0; j < n_species; ++j) {
            for (std::size_t k = 0; k < n_species; ++k) {
                D_bin_local(j, k) = coeffs.diffusion.Dij_bin(i * n_species + j, k);
            }
        }
        
        // TDR term if thermal diffusion is considered
        std::vector<double> tdr_term_local;
        if (sim_config.consider_thermal_diffusion) {
            tdr_term_local.resize(n_species);
            for (std::size_t j = 0; j < n_species; ++j) {
                tdr_term_local[j] = coeffs.thermal_diffusion.tdr_term(i, j);
            }
        }
        
        // Calculate fluxes
        try {
            auto J_local = calculate_stefan_maxwell_at_point(
                c_local, dc_deta_local, x_local,
                coeffs.thermodynamic.rho[i], bc.P_e(),
                D_bin_local, coeffs.thermodynamic.MW[i],
                coeffs.thermodynamic.d_MW_deta[i],
                tdr_term_local, der_fact, K_bl, mixture
            );
            
            // Store
            for (std::size_t j = 0; j < n_species; ++j) {
                J(j, i) = J_local[j];
            }
        } catch (const std::exception& e) {
            return std::unexpected(CoefficientError(
                std::format("Stefan-Maxwell calculation failed at eta={}: {}", i, e.what())
            ));
        }
    }
    
    // Compute derivatives
    for (std::size_t j = 0; j < n_species; ++j) {
        std::vector<double> J_row(n_eta);
        for (std::size_t i = 0; i < n_eta; ++i) {
            J_row[i] = J(j, i);
        }
        
        auto dJ = derivatives::compute_eta_derivative(J_row, d_eta);
        
        for (std::size_t i = 0; i < n_eta; ++i) {
            dJ_deta(j, i) = dJ[i];
        }
    }
    
    return {};
}

} // namespace blast::boundary_layer::coefficients::diffusion