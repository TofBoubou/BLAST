/* #include "blast/boundary_layer/thermodynamics/enthalpy_temperature_solver.hpp"
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
    
    double a = temp_min;
    double b = temp_max;

    auto fa_result = enthalpy_residual(a);
    if (!fa_result) return std::unexpected(fa_result.error());
    double fa = *fa_result;

    auto fb_result = enthalpy_residual(b);
    if (!fb_result) return std::unexpected(fb_result.error());
    double fb = *fb_result;

    if (fa * fb > 0) {
        // std::cout << "[DEBUG] Root not initially bracketed, attempting expansion..." << std::endl;
        
        // Attempt to expand the bracketing interval
        for (int i = 0; i < config_.max_bracket_expansions && fa * fb > 0; ++i) {
            // std::cout << "[DEBUG] Expansion attempt " << (i+1) << "/" << config_.max_bracket_expansions << std::endl;
            
            if (fa > 0 && fb > 0) {
                // Root lies below the current lower bound
                double old_a = a;
                a = std::max(1.0, a * 0.5);
                auto fa_res = enthalpy_residual(a);
                if (!fa_res) return std::unexpected(fa_res.error());
                fa = *fa_res;
                // std::cout << "[DEBUG] Expanded lower bound: " << old_a << " -> " << a 
                          // << ", fa=" << fa << std::endl;
            } else if (fa < 0 && fb < 0) {
                // Root lies above the current upper bound
                double old_b = b;
                b *= 2.0;
                auto fb_res = enthalpy_residual(b);
                if (!fb_res) return std::unexpected(fb_res.error());
                fb = *fb_res;
            }
            
            // std::cout << "[DEBUG] After expansion: fa*fb=" << (fa * fb) << std::endl;
        }

        if (fa * fb > 0) {
            // std::cout << "[DEBUG] Failed to bracket root after all expansions!" << std::endl;
            return std::unexpected(ThermodynamicSolverError(
                "Root not bracketed: f({})={}, f({})={}",
                std::source_location::current(), a, fa, b, fb));
        } else {
            // std::cout << "[DEBUG] Successfully bracketed root: a=" << a << ", fa=" << fa 
                      // << ", b=" << b << ", fb=" << fb << std::endl;
        }
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



} // namespace blast::boundary_layer::thermodynamics */


/* #include "blast/boundary_layer/thermodynamics/enthalpy_temperature_solver.hpp"
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

    auto h_initial = mixture_.mixture_enthalpy(composition, initial_temperature, pressure);
    if (!h_initial) {
        // Fallback silencieux
        // std::cout << "BONJOURRRRRRRRRR" << std::endl;
        return temp_max;
    }

    if (std::abs(*h_initial - target_enthalpy) < config_.tolerance) {
        return initial_temperature;
    }

    // Résolution par méthode de Brent avec bornes
    return brent_method(composition, target_enthalpy, pressure, temp_min, temp_max);
}

auto EnthalpyTemperatureSolver::brent_method(
    std::span<const double> composition,
    double target_enthalpy,
    double pressure,
    double temp_min,
    double temp_max
) const -> std::expected<double, ThermodynamicSolverError> {

    auto enthalpy_residual = [&](double temperature) -> std::optional<double> {
        auto h = mixture_.mixture_enthalpy(composition, temperature, pressure);
        if (!h) return std::nullopt;
        // std::cout << "UUUUUUUUUUUUUU : " << *h << std::endl;
        return *h - target_enthalpy;
    };

    double a = temp_min;
    double b = temp_max;

    auto fa_result = enthalpy_residual(a);
    auto fb_result = enthalpy_residual(b);
    if (!fa_result || !fb_result) return temp_max;

    double fa = *fa_result;
    double fb = *fb_result;

    if (fa * fb > 0.0) {
        // Racine non encadrée : fallback silencieux vers la meilleure borne
        // std::cout << "VVVVVVVVVVV : " << std::abs(fa) << "   " << std::abs(fb) << std::endl;
        return (std::abs(fa) < std::abs(fb)) ? a : b;
    }

    bool mflag = true;
    double c = a;
    double fc = fa;
    double d = 0.0;
    double s = 0.0;

    const double step_tolerance = std::pow(config_.tolerance, 2.0 / 3.0);

    // std::cout << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT " << config_.max_iterations << std::endl;

    for (int iteration = 0; iteration < config_.max_iterations; ++iteration) {

        if (std::abs(fb) < config_.tolerance || std::abs(b - a) < step_tolerance) {
            std::cout << "GVYTVV" << std::endl;
            return b;
        }

        if (std::abs(fa - fc) > 1e-10 && std::abs(fb - fc) > 1e-10) {
            // Interpolation quadratique inverse
            double R = fb / fc;
            double S = fb / fa;
            double T = fa / fc;
            double P = S * (T * (R - T) * (c - b) - (1 - R) * (b - a));
            double Q = (T - 1) * (R - 1) * (S - 1);
            s = b + P / Q;
        } else {
            // Méthode de la sécante
            s = b - fb * (b - a) / (fb - fa);
        }

        bool condition1 = (s < (3 * a + b) / 4) || (s > b);
        bool condition2 = mflag && (std::abs(s - b) >= std::abs(b - c) / 2);
        bool condition3 = !mflag && (std::abs(s - b) >= std::abs(c - d) / 2);
        bool condition4 = mflag && (std::abs(b - c) < step_tolerance);
        bool condition5 = !mflag && (std::abs(c - d) < step_tolerance);

        if (condition1 || condition2 || condition3 || condition4 || condition5) {
            s = 0.5 * (a + b);
            mflag = true;
        } else {
            mflag = false;
        }

        auto fs_result = enthalpy_residual(s);
        if (!fs_result) return temp_max;

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

    // Si on atteint le nombre max d'itérations : fallback silencieux
    return b;
}

auto EnthalpyTemperatureSolver::solve(
    std::span<const double> enthalpy_field,
    const core::Matrix<double>& composition,
    const conditions::BoundaryConditions& bc,
    std::span<const double> initial_temperatures
) const -> std::expected<TemperatureField, ThermodynamicSolverError> {

    const auto n_points = enthalpy_field.size();

    TemperatureField result;
    result.temperatures.resize(n_points);
    result.temperatures[0] = bc.Tw();

    std::cout << "Debug dans calcul température : Tw = " << bc.Tw() << " P_e : " << bc.P_e() << " h_e : " << bc.he() << std::endl;

    for (std::size_t i = 1; i < n_points; ++i) {
        std::vector<double> point_composition(mixture_.n_species());
        for (std::size_t j = 0; j < mixture_.n_species(); ++j) {
            point_composition[j] = composition(j, i);
        }

        auto temp_result = solve_single_point(
            point_composition,
            enthalpy_field[i],
            bc.P_e(),
            initial_temperatures[i]
        );

        result.temperatures[i] = temp_result.value_or(config_.max_temperature);
        std::cout << "Température calculée : " << result.temperatures[i] << std::endl;
    }

    return result;
}

} // namespace blast::boundary_layer::thermodynamics */


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


    // std::cout << "Fractions dans température c0: " << composition[0] << "  " << composition[1] << "  " << composition[2] << "  " << composition[3]  << "  " << composition[4] << std::endl;
    // std::cout << "Enthalpie que l'on cherche à avoir : " << target_enthalpy << std::endl;

    // Réplication exacte de get_temp() - bit par bit
    bool mflag = true;
    double a = 1;
    double b = 10000;

    double d = 0;

    double c_brent = a;
    
    // Pas de vérification d'erreur - on assume succès comme dans le code C
    auto fa_result = mixture_.mixture_enthalpy(composition, a, pressure);
    auto fb_result = mixture_.mixture_enthalpy(composition, b, pressure);

    // std::cout << "fa_result : " << *fa_result << std::endl;
    // std::cout << "fb_result : " << *fb_result << std::endl;
    
    double fa = *fa_result - target_enthalpy;
    double fb = *fb_result - target_enthalpy;
    double fc = fb;
    double fs;

    double R, S, T, U, Q;  // Variables nommées exactement comme dans le code C
    double s;

    int count = 0;
    int max_ite = config_.max_iterations;  // Assumé équivalent à num_params.ite_h2t * 1e3
    double conv = std::pow(config_.tolerance, 1.0 / 3.0);  // pow(num_params.conv_h2t, 1.0 / 3.0)
    double step_tol = std::pow(config_.tolerance, 2.0 / 3.0);  // pow(num_params.conv_h2t, 2.0 / 3.0)
    
    // Boucle while EXACTEMENT comme dans le code C
    while(count < max_ite && std::abs(fb) > conv){

        count++;

        if (std::abs(fb) < conv || std::abs(b - a) < step_tol) {
            break;
        }

        if(std::abs(fa - fc) > 1e-10 && std::abs(fb - fc) > 1e-10){

            R = fb / fc;
            S = fb / fa;
            T = fa / fc;

            U = S * (T * (R - T) * (c_brent - b) - (1 - R) * (b - a));
            Q = (T - 1) * (R - 1) * (S - 1);
            s = b + U / Q;
        }
        else{

            s = b - fb * (b - a) / (fb - fa);
        }

        bool condition1 = (s < (3 * a + b) / 4) || (s > b);
        bool condition2 = mflag && (std::abs(s - b) >= std::abs(b - c_brent) / 2);
        bool condition3 = !mflag && (std::abs(s - b) >= std::abs(c_brent - d) / 2);
        bool condition4 = mflag && (std::abs(b - c_brent) < std::abs(step_tol));
        bool condition5 = !mflag && (std::abs(c_brent - d) < std::abs(step_tol));
        
        if (condition1 || condition2 || condition3 || condition4 || condition5) {
            s = (a + b) / 2;
            mflag = true;
        } else {
            mflag = false;
        }

        auto fs_result = mixture_.mixture_enthalpy(composition, s, pressure);
        fs = *fs_result - target_enthalpy;

        d = c_brent;
        c_brent = b;
        fc = fb;

        if(fa * fs < 0){
            b = s;
            fb = fs;
        }
        else{
            a = s;
            fa = fs;
        }

        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }

/*         std::cout << "Enthalpie que l'on cherche à avoir : " << target_enthalpy << std::endl;
        std::cout << "Encadrement enthalpie : " << fa << " " << fb << std::endl;
        std::cout << "count : " << count << std::endl; */
    }

    return b;  // Retour direct, pas de fallback
}

auto EnthalpyTemperatureSolver::solve(
    std::span<const double> enthalpy_field,
    const core::Matrix<double>& composition,
    const conditions::BoundaryConditions& bc,
    std::span<const double> initial_temperatures
) const -> std::expected<TemperatureField, ThermodynamicSolverError> {

    const auto n_points = enthalpy_field.size();

    // Variables exactement comme dans le code C
    double P = bc.P_e(); 
    double h_e = bc.he(); 
    double T_w = bc.Tw(); 

    std::cout << "Debug dans calcul température : Tw = " << T_w << " P_e : " << P << " h_e : " << h_e << std::endl;

    TemperatureField result;
    result.temperatures.resize(n_points);

    // Application inconditionnelle de la condition aux limites (hypothèse posée)
    // Équivalent à thermal_bc == 0 dans le code C
    result.temperatures[0] = T_w;
    int start_index = 1;

    // Boucle for exactement comme dans le code C
    for(int i = start_index; i < static_cast<int>(n_points); i++){

        // Création du vecteur de composition pour ce point
        std::vector<double> c_temp(mixture_.n_species());
        for(std::size_t j = 0; j < mixture_.n_species(); j++){
            c_temp[j] = composition(j, i);  // c_temp[j] = c[j][i];
        }

        // double h = enthalpy_field[i] * h_e;  // double h = g[i] * h_e;

        double h = enthalpy_field[i];
        
        auto temp_result = solve_single_point(
            c_temp,
            h,
            P,
            initial_temperatures[i]  // T[i] comme guess initial
        );

        result.temperatures[i] = *temp_result;  // T[i] = get_temp(c_temp, h, P, T[i]);

        // std::cout << "T[" << i << "] : " << result.temperatures[i] << std::endl;
        // std::cout << "POUR VOIR" << std::endl;
    }

    // std::cout << "Terminer T[19] : " << result.temperatures[19] << std::endl;

    // result.temperatures[19] = 5000.5;

    // Pas de modification de Tw car on applique thermal_bc == 0 équivalent

    return result;
}

} // namespace blast::boundary_layer::thermodynamics
