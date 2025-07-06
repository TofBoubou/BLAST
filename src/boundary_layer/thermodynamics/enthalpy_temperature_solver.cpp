#include "blast/boundary_layer/thermodynamics/enthalpy_temperature_solver.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace blast::boundary_layer::thermodynamics {

auto EnthalpyTemperatureSolver::solve_single_point(
    std::span<const double> composition,
    double target_enthalpy,
    double pressure,
    double initial_temperature
) const -> std::expected<double, ThermodynamicSolverError> {
    
    const double temp_min = config_.min_temperature;
    const double temp_max = config_.max_temperature;
    
    // Check if initial temperature gives exact solution
    auto h_initial = mixture_.mixture_enthalpy(composition, initial_temperature, pressure);
    if (!h_initial) {
        return std::unexpected(ThermodynamicSolverError("Failed to compute initial enthalpy"));
    }
    
    if (std::abs(*h_initial - target_enthalpy) < config_.tolerance) {
        return initial_temperature;
    }
    
    // Use Brent's method for robust root finding
    return brent_method(composition, target_enthalpy, pressure, temp_min, temp_max);
}

auto EnthalpyTemperatureSolver::brent_method(
    std::span<const double> composition,
    double target_enthalpy,
    double pressure,
    double temp_min,
    double temp_max
) const -> std::expected<double, ThermodynamicSolverError> {
    
    auto enthalpy_residual = [&](double temperature) -> std::expected<double, ThermodynamicSolverError> {
        auto h = mixture_.mixture_enthalpy(composition, temperature, pressure);
        if (!h) {
            return std::unexpected(ThermodynamicSolverError("Failed to compute mixture enthalpy at T={}", std::source_location::current(), temperature));
        }
        return *h - target_enthalpy;
    };

    // DEBUG: std::cout << target_enthalpy << std::endl;
    
    double a = temp_min;
    double b = temp_max;
    
    auto fa_result = enthalpy_residual(a);
    if (!fa_result) return std::unexpected(fa_result.error());
    double fa = *fa_result;
    
    auto fb_result = enthalpy_residual(b);
    if (!fb_result) return std::unexpected(fb_result.error());
    double fb = *fb_result;
    
    if (fa * fb > 0) {
        return std::unexpected(ThermodynamicSolverError(
            "Root not bracketed: f({})={}, f({})={}", std::source_location::current(), a, fa, b, fb));
    }
    
    bool mflag = true;
    double c = a;
    double fc = fa;
    double d = 0.0;
    double s = 0.0;
    
    const double step_tolerance = std::pow(config_.tolerance, 2.0/3.0);
    
    for (int iteration = 0; iteration < config_.max_iterations; ++iteration) {
        
        if (std::abs(fb) < config_.tolerance || std::abs(b - a) < step_tolerance) {
            return b;
        }
        
        // Calculate s using inverse quadratic interpolation or secant method
        if (std::abs(fa - fc) > 1e-10 && std::abs(fb - fc) > 1e-10) {
            // Inverse quadratic interpolation
            double R = fb / fc;
            double S = fb / fa;
            double T = fa / fc;
            double P = S * (T * (R - T) * (c - b) - (1 - R) * (b - a));
            double Q = (T - 1) * (R - 1) * (S - 1);
            s = b + P / Q;
        } else {
            // Secant method
            s = b - fb * (b - a) / (fb - fa);
        }
        
        // Check conditions for using bisection
        bool condition1 = (s < (3 * a + b) / 4) || (s > b);
        bool condition2 = mflag && (std::abs(s - b) >= std::abs(b - c) / 2);
        bool condition3 = !mflag && (std::abs(s - b) >= std::abs(c - d) / 2);
        bool condition4 = mflag && (std::abs(b - c) < step_tolerance);
        bool condition5 = !mflag && (std::abs(c - d) < step_tolerance);
        
        if (condition1 || condition2 || condition3 || condition4 || condition5) {
            s = (a + b) / 2;
            mflag = true;
        } else {
            mflag = false;
        }
        
        auto fs_result = enthalpy_residual(s);
        if (!fs_result) return std::unexpected(fs_result.error());
        double fs = *fs_result;
        
        d = c;
        c = b;
        fc = fb;
        
        if (fa * fs < 0) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        
        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }
    
    return std::unexpected(ThermodynamicSolverError(
        "Maximum iterations ({}) reached without convergence", std::source_location::current(), config_.max_iterations));
}

auto EnthalpyTemperatureSolver::solve(
    std::span<const double> enthalpy_field,
    const core::Matrix<double>& composition,
    const conditions::BoundaryConditions& bc,
    std::span<const double> initial_temperatures
) const -> std::expected<TemperatureField, ThermodynamicSolverError> {
    
    const auto n_points = enthalpy_field.size();

    if (composition.rows() != mixture_.n_species()) {
        return std::unexpected(ThermodynamicSolverError(
            "Composition matrix has {} species, expected {}", 
            std::source_location::current(), composition.rows(), mixture_.n_species()));
    }

    if (composition.cols() != n_points) {
        return std::unexpected(ThermodynamicSolverError(
            "Composition matrix has {} points, expected {}", 
            std::source_location::current(), composition.cols(), n_points));
    }

    if (initial_temperatures.size() != n_points) {
        return std::unexpected(ThermodynamicSolverError(
            "Initial temperature array has {} points, expected {}", 
            std::source_location::current(), initial_temperatures.size(), n_points));
    }

    TemperatureField result;
    result.temperatures.resize(n_points);

    // Fix wall temperature at index 0
    result.temperatures[0] = bc.Tw();

    // Solve for interior points (starting from index 1)
    for (std::size_t i = 1; i < n_points; ++i) {
        // Extract composition for this point
        std::vector<double> point_composition(mixture_.n_species());
        for (std::size_t j = 0; j < mixture_.n_species(); ++j) {
            point_composition[j] = composition(j, i);
        }

        // Solve for temperature
        auto temp_result = solve_single_point(
            point_composition,
            enthalpy_field[i],
            bc.P_e(),
            initial_temperatures[i]
        );

        if (!temp_result) {
            return std::unexpected(ThermodynamicSolverError(
                "Failed to solve temperature at point {}: {}", 
                std::source_location::current(), i, temp_result.error().message()));
        }

        result.temperatures[i] = temp_result.value();
    }

    return result;
}



} // namespace blast::boundary_layer::thermodynamics