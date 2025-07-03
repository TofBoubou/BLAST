#include "blast/boundary_layer/coefficients/xi_derivatives.hpp"
#include <algorithm>

namespace blast::boundary_layer::coefficients {

void XiDerivatives::update_station(int station, double xi,
                                  std::span<const double> F,
                                  std::span<const double> g,
                                  const core::Matrix<double>& c) {
    station_ = station;
    
    if (station == 0) {
        // Stagnation point
        xi_current_ = xi;
        const double d_xi = xi;
        lambda0_ = 1.0 / d_xi;
        lambda1_ = -1.0 / d_xi;
        lambda2_ = 0.0;
        
        // Save current solution
        std::ranges::copy(F, F_minus_1_.begin());
        std::ranges::copy(g, g_minus_1_.begin());
        c_minus_1_ = c;
        
        // Compute derivatives using backward difference
        for (std::size_t i = 0; i < n_eta_; ++i) {
            F_xi_der_[i] = lambda1_ * F_minus_1_[i];
            g_xi_der_[i] = lambda1_ * g_minus_1_[i];
            
            for (std::size_t j = 0; j < n_species_; ++j) {
                c_xi_der_(j, i) = lambda1_ * c_minus_1_(j, i);
            }
        }
    } else {
        // Downstream station
        xi_minus_2_ = xi_minus_1_;
        xi_minus_1_ = xi_current_;
        xi_current_ = xi;
        
        // Update lambda coefficients for three-point formula
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
    for (std::size_t i = 0; i < n_eta_; ++i) {
        F_xi_der_[i] = lambda1_ * F_minus_1_[i] + lambda2_ * F_minus_2_[i];
        g_xi_der_[i] = lambda1_ * g_minus_1_[i] + lambda2_ * g_minus_2_[i];
        
        for (std::size_t j = 0; j < n_species_; ++j) {
            c_xi_der_(j, i) = lambda1_ * c_minus_1_(j, i) + lambda2_ * c_minus_2_(j, i);
        }
    }
}

} // namespace blast::boundary_layer::coefficients