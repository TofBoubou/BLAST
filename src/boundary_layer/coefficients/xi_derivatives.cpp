#include "blast/boundary_layer/coefficients/xi_derivatives.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace blast::boundary_layer::coefficients {

void XiDerivatives::update_station(int station, double xi,
                                  std::span<const double> F,
                                  std::span<const double> g,
                                  const core::Matrix<double>& c) {
    station_ = station;
    
    if (station == 1) {
    // Stagnation point
    xi_minus_1_ = xi_current_; 
    xi_current_ = xi;    
    double d_xi = xi_current_ - xi_minus_1_;

    std::cout << std::scientific << xi << std::endl;

    lambda0_ = 1.0 / d_xi;
    lambda1_ = -1.0 / d_xi;
    lambda2_ = 0.0;  // No second-order term for 2-point formula

    F_minus_2_ = std::move(F_minus_1_);
    g_minus_2_ = std::move(g_minus_1_);
    c_minus_2_ = std::move(c_minus_1_);
        
    F_minus_1_.resize(n_eta_);
    g_minus_1_.resize(n_eta_);
    std::ranges::copy(F, F_minus_1_.begin());
    std::ranges::copy(g, g_minus_1_.begin());
    c_minus_1_ = c;
        
    compute_derivatives();
        
    } else {
        xi_minus_2_ = xi_minus_1_;
        xi_minus_1_ = xi_current_;
        xi_current_ = xi;
   
        lambda0_ = 1.0 / (xi_current_ - xi_minus_1_) + 1.0 / (xi_current_ - xi_minus_2_);
        lambda1_ = (xi_current_ - xi_minus_2_) / 
                   ((xi_minus_1_ - xi_current_) * (xi_minus_1_ - xi_minus_2_));
        lambda2_ = (xi_current_ - xi_minus_1_) / 
                   ((xi_minus_2_ - xi_current_) * (xi_minus_2_ - xi_minus_1_));
        
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

void XiDerivatives::compute_derivatives() {
    
    for (std::size_t i = 0; i < n_eta_; ++i) {

        F_xi_der_[i] = lambda1_ * F_minus_1_[i] + lambda2_ * F_minus_2_[i];
        g_xi_der_[i] = lambda1_ * g_minus_1_[i] + lambda2_ * g_minus_2_[i];
        
        for (std::size_t j = 0; j < n_species_; ++j) {
            c_xi_der_(j, i) = lambda1_ * c_minus_1_(j, i) + lambda2_ * c_minus_2_(j, i);
        }
    }
}

} // namespace blast::boundary_layer::coefficients