#pragma once
#include "../../core/containers.hpp"
#include <vector>
#include <span>

namespace blast::boundary_layer::coefficients {

// Structure to manage xi derivatives for finite difference approximations
class XiDerivatives {
private:
    int station_ = 0;
    // Xi values at current and previous stations
    double xi_current_ = 0.0;
    double xi_minus_1_ = 0.0;
    double xi_minus_2_ = 0.0;
    
    // Lambda coefficients for derivative approximation
    double lambda0_ = 0.0;
    double lambda1_ = 0.0;
    double lambda2_ = 0.0;
    
    // Solution at previous stations
    std::vector<double> F_minus_1_;
    std::vector<double> F_minus_2_;
    std::vector<double> g_minus_1_;
    std::vector<double> g_minus_2_;
    core::Matrix<double> c_minus_1_;
    core::Matrix<double> c_minus_2_;
    
    // Computed derivatives
    std::vector<double> F_xi_der_;
    std::vector<double> g_xi_der_;
    core::Matrix<double> c_xi_der_;
    
    std::size_t n_eta_;
    std::size_t n_species_;
    
public:
    explicit XiDerivatives(std::size_t n_eta, std::size_t n_species) 
        : n_eta_(n_eta), n_species_(n_species) {
        initialize();
    }
    
    void initialize() {
        F_minus_1_.resize(n_eta_, 0.0);
        F_minus_2_.resize(n_eta_, 0.0);
        g_minus_1_.resize(n_eta_, 0.0);
        g_minus_2_.resize(n_eta_, 0.0);
        c_minus_1_ = core::Matrix<double>(n_species_, n_eta_);
        c_minus_2_ = core::Matrix<double>(n_species_, n_eta_);
        
        F_xi_der_.resize(n_eta_, 0.0);
        g_xi_der_.resize(n_eta_, 0.0);
        c_xi_der_ = core::Matrix<double>(n_species_, n_eta_);
        
        c_minus_1_.setZero();
        c_minus_2_.setZero();
        c_xi_der_.setZero();
    }
    
    // Update with new station data
    void update_station(int station, double xi, 
                       std::span<const double> F,
                       std::span<const double> g,
                       const core::Matrix<double>& c);
    
    // Handle convergence failure
    void update_failed_convergence(int station, double xi,
                                  std::span<double> F,
                                  std::span<double> g,
                                  core::Matrix<double>& c);
    
    // Accessors
    [[nodiscard]] constexpr auto station() const noexcept { return station_; }
    [[nodiscard]] constexpr auto lambda0() const noexcept { return lambda0_; }
    [[nodiscard]] constexpr auto lambda1() const noexcept { return lambda1_; }
    [[nodiscard]] constexpr auto lambda2() const noexcept { return lambda2_; }
    
    [[nodiscard]] auto F_derivative() const noexcept -> std::span<const double> { 
        return F_xi_der_; 
    }
    [[nodiscard]] auto g_derivative() const noexcept -> std::span<const double> { 
        return g_xi_der_; 
    }
    [[nodiscard]] auto c_derivative() const noexcept -> const core::Matrix<double>& { 
        return c_xi_der_; 
    }

private:
    void compute_derivatives();
};

} // namespace blast::boundary_layer::coefficients