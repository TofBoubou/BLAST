#include "blast/boundary_layer/solver/adaptive_relaxation_controller.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace blast::boundary_layer::solver {

auto AdaptiveRelaxationController::adapt_relaxation_factor(
    const ConvergenceInfo& conv_info, 
    int iteration
) -> double {
    
    const double current_residual = conv_info.max_residual();
    
    if (first_iteration_) {
        first_iteration_ = false;
        previous_residual_ = current_residual;
        residual_history_.push_back(current_residual);
        return current_factor_;
    }
    
    // Compute residual change ratio
    const double residual_ratio = current_residual / previous_residual_;
    
    // Update history
    residual_history_.push_back(current_residual);
    if (residual_history_.size() > history_size_) {
        residual_history_.erase(residual_history_.begin());
    }
    
    // Detect oscillations
    const bool oscillating = detect_oscillations();
    
    // Adaptation logic — BASED ON CONFIG THRESHOLDS AND FACTORS
    if (oscillating) {
        // If oscillating, reduce aggressively
        current_factor_ *= config_.oscillation_penalty;
        consecutive_improvements_ = 0;
        consecutive_deteriorations_ = 0;
        
        std::cout << "[RELAXATION] Oscillations detected, reducing factor to " 
                  << std::scientific << std::setprecision(3) << current_factor_ << std::endl;
                  
    } else if (residual_ratio > config_.divergence_threshold) {
        // Residual increases → REDUCE
        consecutive_deteriorations_++;
        consecutive_improvements_ = 0;
        
        const double old_factor = current_factor_;
        current_factor_ *= config_.decrease_factor;
        
        std::cout << "[RELAXATION] Residual increased (ratio=" << residual_ratio 
                  << "), reducing factor from " << old_factor << " to " << current_factor_ << std::endl;
                  
    } else if (residual_ratio >= config_.excellent_threshold && residual_ratio <= config_.divergence_threshold) {
        // Moderate convergence → SLIGHT INCREASE
        consecutive_improvements_++;
        consecutive_deteriorations_ = 0;
        
        const double old_factor = current_factor_;
        current_factor_ = std::min(current_factor_ * config_.moderate_increase, config_.max_factor);
        
        std::cout << "[RELAXATION] Moderate convergence (ratio=" << residual_ratio 
                  << "), slight increase from " << old_factor << " to " << current_factor_ << std::endl;
        
    } else if (residual_ratio < config_.excellent_threshold) {
        // Excellent convergence → STRONG INCREASE
        consecutive_improvements_++;
        consecutive_deteriorations_ = 0;
        
        const double old_factor = current_factor_;
        current_factor_ = std::min(current_factor_ * config_.strong_increase, config_.max_factor);
        
        std::cout << "[RELAXATION] Excellent convergence (ratio=" << residual_ratio 
                  << "), strong increase from " << old_factor << " to " << current_factor_ << std::endl;
    }
    
    // Clamp within bounds
    current_factor_ = std::clamp(current_factor_, config_.min_factor, config_.max_factor);
    
    // Extra conservative for initial iterations
    if (iteration < 5) {
        current_factor_ = std::min(current_factor_, 0.3);
    }
    
    // Safety constraint for extremely poor solutions
    if (current_residual > 1e2) {
        current_factor_ = std::min(current_factor_, 0.1);
    }
    
    previous_residual_ = current_residual;
    return current_factor_;
}

auto AdaptiveRelaxationController::detect_oscillations() const -> bool {
    if (residual_history_.size() < 4) return false;
    
    // Check if residuals oscillate around a mean value
    const std::size_t n = residual_history_.size();
    const std::size_t window_size = std::min(n, std::size_t{4});
    
    double mean = 0.0;
    for (std::size_t i = n - window_size; i < n; ++i) {
        mean += residual_history_[i];
    }
    mean /= window_size;
    
    int sign_changes = 0;
    for (std::size_t i = n - window_size + 1; i < n; ++i) {
        const double dev_current = residual_history_[i] - mean;
        const double dev_previous = residual_history_[i - 1] - mean;
        
        if ((dev_current > 0) != (dev_previous > 0)) {
            sign_changes++;
        }
    }
    
    // Consider oscillation if several sign changes
    return sign_changes >= 2;
}

auto AdaptiveRelaxationController::compute_residual_trend() const -> double {
    if (residual_history_.size() < 3) return 0.0;
    
    // Simple linear regression over the last few points
    const std::size_t n = std::min(residual_history_.size(), std::size_t{5});
    const std::size_t start_idx = residual_history_.size() - n;
    
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0;
    
    for (std::size_t i = 0; i < n; ++i) {
        const double x = static_cast<double>(i);
        const double y = std::log(residual_history_[start_idx + i] + 1e-15); // Log for exponential-like trend
        
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_x2 += x * x;
    }
    
    const double slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
    return slope; // Negative slope = improving, positive = worsening
}

// Factory functions
namespace relaxation_factory {

auto create_for_station(int station_number) -> AdaptiveRelaxationController {
    if (station_number == 0) {
        return AdaptiveRelaxationController(AdaptiveRelaxationController::Config::for_stagnation_point());
    } else {
        return AdaptiveRelaxationController(AdaptiveRelaxationController::Config::for_downstream_station());
    }
}

} // namespace relaxation_factory

} // namespace blast::boundary_layer::solver
