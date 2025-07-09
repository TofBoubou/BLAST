#include "blast/thermophysics/mutation_mixture.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <format>
#include <iostream>

namespace blast::thermophysics {

namespace {
    // Helper to convert config enums to Mutation++ strings
    [[nodiscard]] auto database_to_string(io::MixtureConfig::Database db) -> std::expected<std::string, ThermophysicsError> {
        switch (db) {
            case io::MixtureConfig::Database::RRHO: return "RRHO";
            case io::MixtureConfig::Database::NASA7: return "NASA-7";
            case io::MixtureConfig::Database::NASA9: return "NASA-9";
        }
        return std::unexpected(ThermophysicsError("Unknown thermodynamic database type"));
    }
    
    [[nodiscard]] auto viscosity_algo_to_string(io::MixtureConfig::ViscosityAlgorithm algo) -> std::expected<std::string, ThermophysicsError> {
        switch (algo) {
            case io::MixtureConfig::ViscosityAlgorithm::chapmanEnskog_CG: 
                return "Chapmann-Enskog_CG";
            case io::MixtureConfig::ViscosityAlgorithm::GuptaYos: 
                return "Gupta-Yos";
            case io::MixtureConfig::ViscosityAlgorithm::chapmanEnskog_LDLT: 
                return "Chapmann-Enskog_LDLT";
            case io::MixtureConfig::ViscosityAlgorithm::Wilke: 
                return "Wilke";
        }
        return std::unexpected(ThermophysicsError("Unknown viscosity algorithm type"));
    }
    
    [[nodiscard]] constexpr auto thermal_cond_algo_to_string(io::MixtureConfig::ViscosityAlgorithm algo) -> std::expected<std::string_view, ThermophysicsError> {
        // Simplified mapping for thermal conductivity
        switch (algo) {
            case io::MixtureConfig::ViscosityAlgorithm::Wilke: 
                return "Wilke";
            case io::MixtureConfig::ViscosityAlgorithm::chapmanEnskog_CG:
            case io::MixtureConfig::ViscosityAlgorithm::chapmanEnskog_LDLT:
            case io::MixtureConfig::ViscosityAlgorithm::GuptaYos:
                return "Chapmann-Enskog_CG";
        }
        return std::unexpected(ThermophysicsError("Unknown thermal conductivity algorithm type"));
    }
}

MutationMixture::MutationMixture(const io::MixtureConfig& config) 
    : n_species_(0), has_electrons_(false) { // n_species_ is a const and needs to be defined at the construction
    
    try {
        // Configure Mutation++ options
        Mutation::MixtureOptions opts(config.name);
        
        // Set state model (always ChemNonEq1T for now)
        opts.setStateModel("ChemNonEq1T");
        
        // Set thermodynamic database
        auto db_result = database_to_string(config.thermodynamic_database);
        if (!db_result) {
            throw db_result.error();
        }
        opts.setThermodynamicDatabase(db_result.value());
        
        // Set viscosity algorithm
        auto visc_result = viscosity_algo_to_string(config.viscosity_algorithm);
        if (!visc_result) {
            throw visc_result.error();
        }
        opts.setViscosityAlgorithm(visc_result.value());
        
        // Set thermal conductivity algorithm
        auto thermal_result = thermal_cond_algo_to_string(config.viscosity_algorithm);
        if (!thermal_result) {
            throw thermal_result.error();
        }
        opts.setThermalConductivityAlgorithm(std::string(thermal_result.value()));
        
        // Create mixture
        mixture_ = std::make_unique<Mutation::Mixture>(opts);
        
        // Initialize const members via const_cast in constructor only
        const_cast<std::size_t&>(n_species_) = mixture_->nSpecies();
        const_cast<bool&>(has_electrons_) = mixture_->hasElectrons();
        
        // Cache species properties
        species_mw_.reserve(n_species_);
        species_charges_.reserve(n_species_);
        species_names_.reserve(n_species_);
        
        for (std::size_t i = 0; i < n_species_; ++i) {
            species_mw_.push_back(mixture_->speciesMw(i));
            
            // Charge in C/kg: q_i = z_i * N_A / M_i
            const double charge = mixture_->speciesCharge(i) * constants::N_Avogadro / species_mw_[i];
            species_charges_.push_back(charge);
            
            species_names_.emplace_back(mixture_->speciesName(i));
        }
        
    } catch (const std::exception& e) {
        throw ThermophysicsError(
            std::format("Failed to create Mutation++ mixture: {}", e.what())
        );
    }
}

auto MutationMixture::validate_composition(std::span<const double> fractions) const 
    -> std::expected<void, ThermophysicsError> {
    
    if (fractions.size() != n_species_) {
        return std::unexpected(ThermophysicsError(
            std::format("Invalid composition size: {} (expected {})", 
                       fractions.size(), n_species_)
        ));
    }
    
    const double sum = std::accumulate(fractions.begin(), fractions.end(), 0.0);
    constexpr double tolerance = 1e-6;
    
    if (std::abs(sum - 1.0) > tolerance) {
        return std::unexpected(ThermophysicsError(
            std::format("Mass fractions sum to {} (should be 1.0)", sum)
        ));
    }
    
    if (std::any_of(fractions.begin(), fractions.end(), 
                    [](double x) { return x < 0.0 || !std::isfinite(x); })) {
        return std::unexpected(ThermophysicsError(
            "Invalid mass fraction values (negative or non-finite)"
        ));
    }
    
    return {};
}

auto MutationMixture::mixture_molecular_weight(
    std::span<const double> mass_fractions
) const -> std::expected<double, ThermophysicsError> {
    
    if (auto validation = validate_composition(mass_fractions); !validation) { // expected si ! means unexpected
        return std::unexpected(validation.error());
    }
    
    // MW = 1 / Σ(Y_i / M_i)
    double sum = 0.0;
    for (std::size_t i = 0; i < n_species_; ++i) {
        sum += mass_fractions[i] / species_mw_[i];
        // std::cout << sum << std::endl;
    }
    
    if (sum <= 0.0 || !std::isfinite(sum)) {
        return std::unexpected(ThermophysicsError("Invalid molecular weight calculation"));
    }
    
    return 1.0 / sum;
}

auto MutationMixture::species_molecular_weight(
    std::size_t species_index
) const noexcept -> double {
    return species_mw_[species_index];
}

auto MutationMixture::species_charges() const noexcept -> std::span<const double> {
    return species_charges_;
}

auto MutationMixture::set_state(
    std::span<const double> mass_fractions,
    double temperature,
    double pressure
) const -> std::expected<void, ThermophysicsError> {
    
    if (auto validation = validate_composition(mass_fractions); !validation) {
        return std::unexpected(validation.error());
    }
    
    if (temperature <= 0.0 || pressure <= 0.0) {
        return std::unexpected(ThermophysicsError(
            std::format("Invalid state: T={}, P={}", temperature, pressure)
        ));
    }
    
    try {
        // Convert span to vector for Mutation++ (C API requirement)
        std::vector<double> c_vec(mass_fractions.begin(), mass_fractions.end());
        double vars[2] = {pressure, temperature};
        
        mixture_->setState(c_vec.data(), vars, 2);
        return {};
        
    } catch (const std::exception& e) {
        return std::unexpected(ThermophysicsError(
            std::format("Failed to set state: {}", e.what())
        ));
    }
}

auto MutationMixture::mixture_enthalpy(
    std::span<const double> mass_fractions,
    double temperature,
    double pressure
) const -> std::expected<double, ThermophysicsError> {
    
    if (auto state_result = set_state(mass_fractions, temperature, pressure); !state_result) {
        return std::unexpected(state_result.error());
    }
    
    try {
        return mixture_->mixtureHMass();
    } catch (const std::exception& e) {
        return std::unexpected(ThermophysicsError(
            std::format("Failed to compute enthalpy: {}", e.what())
        ));
    }
}

auto MutationMixture::species_enthalpies(double temperature) const 
    -> std::expected<std::vector<double>, ThermophysicsError> {
    
    if (temperature <= 0.0) {
        return std::unexpected(ThermophysicsError("Invalid temperature"));
    }
    
    try {
        std::vector<double> h_over_RT(n_species_);
        mixture_->speciesHOverRT(temperature, h_over_RT.data());
        
        // Convert from h/RT to J/kg
        const double RT = constants::R_universal * temperature;
        for (std::size_t i = 0; i < n_species_; ++i) {
            h_over_RT[i] = h_over_RT[i] * RT / species_mw_[i];
        }
        
        return h_over_RT;
        
    } catch (const std::exception& e) {
        return std::unexpected(ThermophysicsError(
            std::format("Failed to compute species enthalpies: {}", e.what())
        ));
    }
}

auto MutationMixture::frozen_cp(
    std::span<const double> mass_fractions,
    double temperature,
    double pressure
) const -> std::expected<double, ThermophysicsError> {
    
    if (auto state_result = set_state(mass_fractions, temperature, pressure); !state_result) {
        return std::unexpected(state_result.error());
    }
    
    try {
        return mixture_->mixtureFrozenCpMass();
    } catch (const std::exception& e) {
        return std::unexpected(ThermophysicsError(
            std::format("Failed to compute frozen Cp: {}", e.what())
        ));
    }
}

auto MutationMixture::viscosity(
    std::span<const double> mass_fractions,
    double temperature,
    double pressure
) const -> std::expected<double, ThermophysicsError> {
    
    if (auto state_result = set_state(mass_fractions, temperature, pressure); !state_result) {
        return std::unexpected(state_result.error());
    }
    
    try {
        return mixture_->viscosity();
    } catch (const std::exception& e) {
        return std::unexpected(ThermophysicsError(
            std::format("Failed to compute viscosity: {}", e.what())
        ));
    }
}

auto MutationMixture::frozen_thermal_conductivity(
    std::span<const double> mass_fractions,
    double temperature,
    double pressure
) const -> std::expected<double, ThermophysicsError> {
    
    if (auto state_result = set_state(mass_fractions, temperature, pressure); !state_result) {
        return std::unexpected(state_result.error());
    }
    
    try {
        return mixture_->frozenThermalConductivity();
    } catch (const std::exception& e) {
        return std::unexpected(ThermophysicsError(
            std::format("Failed to compute thermal conductivity: {}", e.what())
        ));
    }
}

auto MutationMixture::binary_diffusion_coefficients(
    double temperature,
    double pressure
) const -> std::expected<core::Matrix<double>, ThermophysicsError> {
    
    if (temperature <= 0.0 || pressure <= 0.0) {
        return std::unexpected(ThermophysicsError("Invalid T,P for diffusion coefficients"));
    }
    
    try {
/*         // Set state with dummy composition (not needed for binary coefficients)
        double vars[2] = {pressure, temperature};
        mixture_->setState(nullptr, vars, 2);
        
        core::Matrix<double> dij(n_species_, n_species_);
        dij.setZero(); */
        
        // Create uniform composition for setState (required by Mutation++ even if not used)
        std::vector<double> uniform_composition(n_species_, 1.0 / n_species_);
        double vars[2] = {pressure, temperature};
        mixture_->setState(uniform_composition.data(), vars, 2);
        
        core::Matrix<double> dij(n_species_, n_species_);
        dij.setZero();
        
        // Get heavy species diffusion coefficients
        const auto dij_heavy = mixture_->collisionDB().nDij();
        
        // Fill symmetric matrix for heavy species
        std::size_t start_idx = has_electrons_ ? 1 : 0;
        std::size_t counter = 0;
        
        for (std::size_t i = start_idx; i < n_species_; ++i) {
            for (std::size_t j = i; j < n_species_; ++j) {
                dij(i, j) = dij_heavy[counter];
                dij(j, i) = dij_heavy[counter];
                ++counter;
            }
        }
        
        // Add electron diffusion if present
        if (has_electrons_) {
            const auto dei = mixture_->collisionDB().nDei();
            for (std::size_t i = 0; i < n_species_; ++i) {
                dij(0, i) = dei[i];
                dij(i, 0) = dei[i];
            }
        }
        
        // Convert from n*D to D by dividing by number density
        const double n = pressure / (constants::k_Boltzmann * temperature);
        for (std::size_t i = 0; i < n_species_; ++i) {
            for (std::size_t j = 0; j < n_species_; ++j) {
                dij(i, j) /= n;
            }
        }
        
        return dij;
        
    } catch (const std::exception& e) {
        return std::unexpected(ThermophysicsError(
            std::format("Failed to compute binary diffusion: {}", e.what())
        ));
    }
}

auto MutationMixture::production_rates(
    std::span<const double> partial_densities,
    double temperature
) const -> std::expected<std::vector<double>, ThermophysicsError> {
    
    if (partial_densities.size() != n_species_) {
        return std::unexpected(ThermophysicsError("Invalid partial densities size"));
    }
    
    if (temperature <= 0.0) {
        return std::unexpected(ThermophysicsError("Invalid temperature"));
    }
    
    try {
        // Convert span to vector for Mutation++
        std::vector<double> rho_vec(partial_densities.begin(), partial_densities.end());
        
        mixture_->setState(rho_vec.data(), &temperature, 1);
        
        std::vector<double> wi(n_species_);
        mixture_->netProductionRates(wi.data());
        
        return wi;
        
    } catch (const std::exception& e) {
        return std::unexpected(ThermophysicsError(
            std::format("Failed to compute production rates: {}", e.what())
        ));
    }
}

auto MutationMixture::production_rate_jacobian(
    std::span<const double> partial_densities,
    double temperature
) const -> std::expected<core::Matrix<double>, ThermophysicsError> {
    
    if (partial_densities.size() != n_species_) {
        return std::unexpected(ThermophysicsError("Invalid partial densities size"));
    }
    
    try {
        // Set state
        std::vector<double> rho_vec(partial_densities.begin(), partial_densities.end());
        mixture_->setState(rho_vec.data(), &temperature, 1);
        
        // Get linearized Jacobian from Mutation++
        std::vector<double> jac_linear(n_species_ * n_species_);
        mixture_->jacobianRho(jac_linear.data());
        
        // Convert to matrix form
        core::Matrix<double> jacobian(n_species_, n_species_);
        for (std::size_t i = 0; i < n_species_; ++i) {
            for (std::size_t j = 0; j < n_species_; ++j) {
                jacobian(i, j) = jac_linear[i * n_species_ + j];
            }
        }
        
        return jacobian;
        
    } catch (const std::exception& e) {
        return std::unexpected(ThermophysicsError(
            std::format("Failed to compute Jacobian: {}", e.what())
        ));
    }
}

auto MutationMixture::thermal_diffusion_ratios(
    std::span<const double> mass_fractions,
    double temperature,
    double pressure
) const -> std::expected<std::vector<double>, ThermophysicsError> {
    
    if (auto state_result = set_state(mass_fractions, temperature, pressure); !state_result) {
        return std::unexpected(state_result.error());
    }
    
    try {
        std::vector<double> tdr(n_species_);
        mixture_->heavyThermalDiffusionRatios(tdr.data());
        return tdr;
        
    } catch (const std::exception& e) {
        return std::unexpected(ThermophysicsError(
            std::format("Failed to compute TDR: {}", e.what())
        ));
    }
}

auto MutationMixture::equilibrium_composition(
    double temperature,
    double pressure
) const -> std::expected<std::vector<double>, ThermophysicsError> {
    
    if (temperature <= 0.0 || pressure <= 0.0) {
        return std::unexpected(ThermophysicsError("Invalid T,P for equilibrium"));
    }
    
    try {
        std::vector<double> x_mole(n_species_);
        mixture_->equilibriumComposition(temperature, pressure, x_mole.data());
        
        // Convert mole fractions to mass fractions
        std::vector<double> y_mass(n_species_);
        mixture_->convert<Mutation::Thermodynamics::X_TO_Y>(x_mole.data(), y_mass.data());
        
        return y_mass;
        
    } catch (const std::exception& e) {
        return std::unexpected(ThermophysicsError(
            std::format("Failed to compute equilibrium: {}", e.what())
        ));
    }
}

// Factory function implementation
auto create_mixture(const io::MixtureConfig& config) 
    -> std::expected<std::unique_ptr<MixtureInterface>, ThermophysicsError> {
    
    try {
        return std::make_unique<MutationMixture>(config);
    } catch (const ThermophysicsError& e) {
        return std::unexpected(e);
    } catch (const std::exception& e) {
        return std::unexpected(ThermophysicsError(
            std::format("Failed to create mixture: {}", e.what())
        ));
    }
}

} // namespace blast::thermophysics