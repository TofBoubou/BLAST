#include "blast/boundary_layer/coefficients/xi_derivatives.hpp"
#include <algorithm>
#include <iostream>

namespace blast::boundary_layer::coefficients {

void XiDerivatives::update_station(int station, double xi,
                                  std::span<const double> F,
                                  std::span<const double> g,
                                  const core::Matrix<double>& c) {
    station_ = station;
    
    if (station == 0) {
    // Stagnation point
    lambda0_ = lambda1_ = lambda2_ = 0.0;
    
    // Initialize xi coordinates for stagnation point
    xi_current_ = xi;  // This should be 0.0 at stagnation
    xi_minus_1_ = xi;  // Initialize to same value to avoid uninitialized values
    xi_minus_2_ = xi;  // Initialize to same value to avoid uninitialized values

    std::ranges::copy(F, F_minus_1_.begin());
    std::ranges::copy(g, g_minus_1_.begin());
    c_minus_1_ = c;

    std::fill(F_xi_der_.begin(), F_xi_der_.end(), 0.0);
    std::fill(g_xi_der_.begin(), g_xi_der_.end(), 0.0);
    c_xi_der_.setZero();

    } else if (station == 1) {
        // First downstream station: use 2-point formula
        // Only update xi values on first call for this station
        if (xi_current_ != xi) {
            xi_minus_1_ = xi_current_;  // Previous xi (should be 0.0 from station 0)  
            xi_current_ = xi;           // Current xi
        }
        
        std::cout << std::scientific << "DEBUG: XiDerivatives station " << station 
                  << " (2-point formula) | xi_current = " << xi_current_ 
                  << " | xi_minus_1 = " << xi_minus_1_ << std::endl;
        
        // Use 2-point backward difference formula
        double d_xi = xi_current_ - xi_minus_1_;
        if (std::abs(d_xi) < 1e-14) {
            std::cout << "ERROR: xi difference too small: " << d_xi << std::endl;
            lambda0_ = lambda1_ = lambda2_ = 0.0;
        } else {
            lambda0_ = 1.0 / d_xi;
            lambda1_ = -1.0 / d_xi;
            lambda2_ = 0.0;  // No second-order term for 2-point formula
        }
        
        std::cout << std::scientific << "DEBUG: Lambda coefficients (2-point): lambda0 = " << lambda0_ 
                  << " | lambda1 = " << lambda1_ 
                  << " | lambda2 = " << lambda2_ << std::endl;
        
        // Shift old solutions
        F_minus_2_ = std::move(F_minus_1_);
        g_minus_2_ = std::move(g_minus_1_);
        c_minus_2_ = std::move(c_minus_1_);
        
        // Save current solution
        F_minus_1_.resize(n_eta_);
        g_minus_1_.resize(n_eta_);
        std::ranges::copy(F, F_minus_1_.begin());
        std::ranges::copy(g, g_minus_1_.begin());
        c_minus_1_ = c;
        
        // Compute derivatives (lambda2 = 0 for 2-point formula)
        compute_derivatives();
        
    } else {
        // Downstream station (>= 2): use 3-point formula
        // Only update xi values on first call for this station
        if (xi_current_ != xi) {
            xi_minus_2_ = xi_minus_1_;
            xi_minus_1_ = xi_current_;
            xi_current_ = xi;
        }
        
        // Update lambda coefficients for three-point formula
        std::cout << std::scientific << "DEBUG: XiDerivatives station " << station 
                  << " (3-point formula) | xi_current = " << xi_current_ 
                  << " | xi_minus_1 = " << xi_minus_1_
                  << " | xi_minus_2 = " << xi_minus_2_ << std::endl;
        std::cout << std::scientific << "DEBUG: Xi differences: (xi_current - xi_minus_1) = " 
                  << (xi_current_ - xi_minus_1_) 
                  << " | (xi_current - xi_minus_2) = " << (xi_current_ - xi_minus_2_) << std::endl;
        
        lambda0_ = 1.0 / (xi_current_ - xi_minus_1_) + 1.0 / (xi_current_ - xi_minus_2_);
        lambda1_ = (xi_current_ - xi_minus_2_) / 
                   ((xi_minus_1_ - xi_current_) * (xi_minus_1_ - xi_minus_2_));
        lambda2_ = (xi_current_ - xi_minus_1_) / 
                   ((xi_minus_2_ - xi_current_) * (xi_minus_2_ - xi_minus_1_));
        
        std::cout << std::scientific << "DEBUG: Lambda coefficients (3-point): lambda0 = " << lambda0_ 
                  << " | lambda1 = " << lambda1_ 
                  << " | lambda2 = " << lambda2_ << std::endl;
        
        // Shift old solutions
        F_minus_2_ = std::move(F_minus_1_);
        g_minus_2_ = std::move(g_minus_1_);
        c_minus_2_ = std::move(c_minus_1_);
        
        // Save current solution
        F_minus_1_.resize(n_eta_);
        g_minus_1_.resize(n_eta_);
        std::ranges::copy(F, F_minus_1_.begin());
        std::ranges::copy(g, g_minus_1_.begin());
        c_minus_1_ = c;
        
        // Compute derivatives
        compute_derivatives();
    }
}

void XiDerivatives::update_failed_convergence(int station, double xi,
                                             std::span<double> F,
                                             std::span<double> g,
                                             core::Matrix<double>& c) {
    station_ = station;
    
    // Restore previous solution
    std::ranges::copy(F_minus_1_, F.begin());
    std::ranges::copy(g_minus_1_, g.begin());
    c = c_minus_1_;
    
    // Update xi
    xi_current_ = xi;
    
    if (station == 1) {
        // First station after stagnation
        const double d_xi = xi;
        lambda0_ = 1.0 / d_xi;
        lambda1_ = -1.0 / d_xi;
        lambda2_ = 0.0;
        
        // Recompute derivatives with backward difference
        for (std::size_t i = 0; i < n_eta_; ++i) {
            F_xi_der_[i] = lambda1_ * F_minus_1_[i];
            g_xi_der_[i] = lambda1_ * g_minus_1_[i];
            
            for (std::size_t j = 0; j < n_species_; ++j) {
                c_xi_der_(j, i) = lambda1_ * c_minus_1_(j, i);
            }
        }
    } else {
        // Recompute lambda coefficients
        lambda0_ = 1.0 / (xi_current_ - xi_minus_1_) + 1.0 / (xi_current_ - xi_minus_2_);
        lambda1_ = (xi_current_ - xi_minus_2_) / 
                   ((xi_minus_1_ - xi_current_) * (xi_minus_1_ - xi_minus_2_));
        lambda2_ = (xi_current_ - xi_minus_1_) / 
                   ((xi_minus_2_ - xi_current_) * (xi_minus_2_ - xi_minus_1_));
        
        // Recompute derivatives
        compute_derivatives();
    }
}

void XiDerivatives::compute_derivatives() {
    std::cout << std::scientific << "DEBUG: compute_derivatives with lambda1 = " << lambda1_ 
              << " | lambda2 = " << lambda2_ << std::endl;
    
    for (std::size_t i = 0; i < n_eta_; ++i) {
        F_xi_der_[i] = lambda1_ * F_minus_1_[i] + lambda2_ * F_minus_2_[i];
        
        if (i < 3) {  // Only print first few values to avoid spam
            std::cout << std::scientific << "DEBUG: F_xi_der[" << i << "] = " << F_xi_der_[i] 
                      << " | F_minus_1[" << i << "] = " << F_minus_1_[i] 
                      << " | F_minus_2[" << i << "] = " << F_minus_2_[i] << std::endl;
        }
        g_xi_der_[i] = lambda1_ * g_minus_1_[i] + lambda2_ * g_minus_2_[i];
        
        for (std::size_t j = 0; j < n_species_; ++j) {
            c_xi_der_(j, i) = lambda1_ * c_minus_1_(j, i) + lambda2_ * c_minus_2_(j, i);
        }
    }
}

} // namespace blast::boundary_layer::coefficients