#include "blast/boundary_layer/coefficients/diffusion.hpp"
#include "blast/boundary_layer/coefficients/coefficient_calculator.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>

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

// Compute average diffusion coefficient D_im
inline auto compute_average_diffusion_coefficient(
    std::span<const double> x,
    const Eigen::MatrixXd& D_bin,
    std::size_t species_idx
) noexcept -> double {
    double sum = 0.0;
    for (std::size_t j = 0; j < x.size(); ++j) {
        sum += x[j] / D_bin(species_idx, j);
    }
    return 1.0 / sum;
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
    const thermophysics::MixtureInterface& mixture
) -> std::vector<double> {
    
    const auto n_species = c.size();
    const double full_der_fact = der_fact * rho;
    const double sum_c = std::accumulate(c.begin(), c.end(), 0.0);

    
/*     std::cout << "[DEBUG] Coefficients à eta = 0" << std::endl;

    for (std::size_t i = 0; i < n_species; ++i) {
        const double MW_i = mixture.species_molecular_weight(i);

        std::cout << "  espèce[" << i << "]" << std::endl;
        std::cout << "    c[i]         = " << std::scientific << c[i] << std::endl;
        std::cout << "    dc_deta[i]   = " << std::scientific << dc_deta[i] << std::endl;
        std::cout << "    x[i]         = " << std::scientific << x[i] << std::endl;
        std::cout << "    MW_i         = " << std::scientific << MW_i << std::endl;

        if (mixture.has_electrons()) {
            auto charges = mixture.species_charges();
            std::cout << "    charge[i]    = " << std::scientific << charges[i] << std::endl;
        }
    }

    std::cout << "  MW (mélange)     = " << std::scientific << MW << std::endl;
    std::cout << "  dMW_deta         = " << std::scientific << dMW_deta << std::endl;
    std::cout << "  rho              = " << std::scientific << rho << std::endl;
    std::cout << "  P                = " << std::scientific << P << std::endl;
    std::cout << "  der_fact         = " << std::scientific << der_fact << std::endl;
    std::cout << "  full_der_fact    = " << std::scientific << der_fact * rho << std::endl;
    std::cout << "  sum_c            = " << std::scientific << sum_c << std::endl; */
    
    // Use Eigen through the wrapper
    auto& D_bin = D_bin_local.eigen();
    
    // Initialize vectors
    Eigen::VectorXd driving_force(n_species);
    Eigen::VectorXd diffusion_flux(n_species);
    Eigen::VectorXd mass_weight(n_species);
    Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n_species, n_species);
    
    // Build Stefan-Maxwell system
    for (std::size_t i = 0; i < n_species; ++i) {
        const double MW_i = mixture.species_molecular_weight(i);
        
        // Driving force: concentration gradient + molecular weight gradient + thermal diffusion
        driving_force[i] = MW / MW_i * dc_deta[i] + c[i] / MW_i * dMW_deta;
        if (!tdr_term.empty()) {
            driving_force[i] += tdr_term[i];
        }
        
        // Average diffusion coefficient
        const double D_im = compute_average_diffusion_coefficient(x, D_bin, i);
        
        // Flux and mass weight terms
        diffusion_flux[i] = -rho * MW_i / MW * D_im * driving_force[i] * full_der_fact;
        mass_weight[i] = c[i] * D_im * MW;
        
        // Build matrix row
        for (std::size_t j = 0; j < n_species; ++j) {
            A(i, j) -= mass_weight[i] / (mixture.species_molecular_weight(j) * D_bin(i, j));
        }
    }
    
    // Apply mass conservation constraint (last row)
    A.row(n_species - 1).setConstant(c[n_species - 1] / sum_c);
    diffusion_flux[n_species - 1] = 0.0;
    
    // Solve linear system
    Eigen::VectorXd J_vec = A.lu().solve(diffusion_flux);
    
    // Apply ambipolar electric field correction for ionized mixtures
    if (mixture.has_electrons()) {
        auto charges = mixture.species_charges();
        Eigen::VectorXd electric_mobility(n_species);
        
        for (std::size_t i = 0; i < n_species; ++i) {
            const double D_im = compute_average_diffusion_coefficient(x, D_bin, i);
            electric_mobility[i] = rho / P * D_im * mixture.species_molecular_weight(i) / MW * 
                                 charges[i] * rho * c[i];
        }
        
        // Calculate electric field strength - Convert span to Eigen vector for dot product
        Eigen::VectorXd charges_vec = Eigen::Map<const Eigen::VectorXd>(charges.data(), charges.size());
        const double charge_flux = charges_vec.dot(diffusion_flux);
        const double current_flux = charges_vec.dot(J_vec);
        const double charge_mobility = charges_vec.dot(electric_mobility);
        
        if (std::abs(charge_mobility) > 1e-30) {
            const double E_field = -(charge_flux + current_flux) / charge_mobility;
            J_vec += E_field * electric_mobility;
        }
    }
    
    // Apply mass conservation correction
    const double sum_J = J_vec.sum();
    std::vector<double> J_result(n_species);
    for (std::size_t i = 0; i < n_species; ++i) {
        J_result[i] = J_vec[i] - c[i] / sum_c * sum_J;
    }

/*     std::cout << "[DEBUG] Stefan-Maxwell flux J_result:" << std::endl;
        for (std::size_t i = 0; i < J_result.size(); ++i) {
            std::cout << "  J[" << i << "] = " << std::scientific << J_result[i] << std::endl;
        } */
    
    return J_result;
}

} // anonymous namespace


auto compute_stefan_maxwell_fluxes(
    const CoefficientInputs& inputs,
    CoefficientSet& coeffs,  // Changed from const to allow modification
    const conditions::BoundaryConditions& bc,
    const XiDerivatives& xi_der,
    const io::SimulationConfig& sim_config,
    const thermophysics::MixtureInterface& mixture,
    double d_eta
) -> std::expected<void, CoefficientError> {
    
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
    const double der_fact = compute_derivative_factor(station, inputs.xi, bc, sim_config);
    
    // Pre-allocate work arrays
    std::vector<double> c_local(n_species);
    std::vector<double> dc_deta_local(n_species);
    std::vector<double> x_local(n_species);
    core::Matrix<double> D_bin_local(n_species, n_species);
    
    // Process each eta point
    for (std::size_t i = 0; i < n_eta; ++i) {
        // Extract local values
        for (std::size_t j = 0; j < n_species; ++j) {
            c_local[j] = inputs.c(j, i);
            dc_deta_local[j] = inputs.dc_deta(j, i);
            x_local[j] = c_local[j] * coeffs.thermodynamic.MW[i] / mixture.species_molecular_weight(j);
        }
        
        // Extract binary diffusion coefficients for this eta point
        for (std::size_t j = 0; j < n_species; ++j) {
            for (std::size_t k = 0; k < n_species; ++k) {
                D_bin_local(j, k) = coeffs.diffusion.Dij_bin(i * n_species + j, k);
            }
        }
        
        // Get thermal diffusion terms if enabled
        std::span<const double> tdr_span;
        if (sim_config.consider_thermal_diffusion && i < coeffs.thermal_diffusion.tdr_term.rows()) {
            auto tdr_row = coeffs.thermal_diffusion.tdr_term.eigen().row(i);
            tdr_span = std::span(tdr_row.data(), n_species);
        }
        
        // Calculate fluxes
        try {
            auto J_local = calculate_stefan_maxwell_at_point(
                c_local, dc_deta_local, x_local,
                coeffs.thermodynamic.rho[i], bc.P_e(),
                D_bin_local, coeffs.thermodynamic.MW[i],
                coeffs.thermodynamic.d_MW_deta[i],
                tdr_span, der_fact, mixture
            );
            
            // Store results
            for (std::size_t j = 0; j < n_species; ++j) {
                J(j, i) = J_local[j];
            }
        } catch (const std::exception& e) {
            return std::unexpected(CoefficientError(
                std::format("Stefan-Maxwell calculation failed at eta[{}]: {}", i, e.what())
            ));
        }
    }
    
    // Compute eta derivatives of fluxes
    for (std::size_t j = 0; j < n_species; ++j) {
        
        // Copy row data to contiguous vector
        std::vector<double> J_species_data(n_eta);
        for (std::size_t i = 0; i < n_eta; ++i) {
            J_species_data[i] = J(j, i);  // Direct access to matrix
        }
        
        // DEBUG: Verify the copied data
/*         std::cout << "[DEBUG] Copied J data for species " << j << ": ";
        for (std::size_t i = 0; i < std::min(n_eta, std::size_t(5)); ++i) {
            std::cout << std::scientific << J_species_data[i] << " ";
        }
        std::cout << std::endl; */
        
        auto dJ_result = derivatives::compute_eta_derivative(J_species_data, d_eta);
        if (!dJ_result) {
            return std::unexpected(CoefficientError("Failed to compute diffusion flux derivative"));
        }
        auto dJ = dJ_result.value();
        
        // DEBUG: Print first few derivatives
/*         std::cout << "[DEBUG] Computed derivatives for species " << j << ": ";
        for (std::size_t i = 0; i < std::min(n_eta, std::size_t(5)); ++i) {
            std::cout << std::scientific << dJ[i] << " ";
        }
        std::cout << std::endl; */
        
        for (std::size_t i = 0; i < n_eta; ++i) {
            dJ_deta(j, i) = dJ[i];
        }
    }
    
    return {};
}

} // namespace blast::boundary_layer::coefficients::diffusion