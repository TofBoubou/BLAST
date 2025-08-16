#include "blast/boundary_layer/thermodynamics/enthalpy_temperature_solver.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace blast::boundary_layer::thermodynamics {

auto EnthalpyTemperatureSolver::solve_single_point(std::span<const double> composition, double target_enthalpy,
                                                   double pressure, double initial_temperature) const
    -> std::expected<double, ThermodynamicSolverError> {

  bool mflag = true;
  double a = 1;
  double b = 10000;

  double d = 0;

  double c_brent = a;

  auto fa_result = mixture_.mixture_enthalpy(composition, a, pressure);
  auto fb_result = mixture_.mixture_enthalpy(composition, b, pressure);

  double fa = *fa_result - target_enthalpy;
  double fb = *fb_result - target_enthalpy;
  double fc = fb;
  double fs;

  double R, S, T, U, Q;
  double s;

  int count = 0;
  int max_ite = config_.max_iterations;
  double conv = std::pow(config_.tolerance, 1.0 / 3.0);
  double step_tol = std::pow(config_.tolerance, 2.0 / 3.0);

  while (count < max_ite && std::abs(fb) > conv) {

    count++;

    if (std::abs(fb) < conv || std::abs(b - a) < step_tol) {
      break;
    }

    if (std::abs(fa - fc) > 1e-10 && std::abs(fb - fc) > 1e-10) {

      R = fb / fc;
      S = fb / fa;
      T = fa / fc;

      U = S * (T * (R - T) * (c_brent - b) - (1 - R) * (b - a));
      Q = (T - 1) * (R - 1) * (S - 1);
      s = b + U / Q;
    } else {

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

  return b;
}

auto EnthalpyTemperatureSolver::solve(std::span<const double> enthalpy_field, const core::Matrix<double>& composition,
                                      const conditions::BoundaryConditions& bc,
                                      std::span<const double> initial_temperatures) const
    -> std::expected<TemperatureField, ThermodynamicSolverError> {

  const auto n_points = enthalpy_field.size();
  
  double P = bc.P_e();
  double h_e = bc.he();
  
  TemperatureField result;
  result.temperatures.resize(n_points);
  
  // Déterminer si on est en mode adiabatique
  bool is_adiabatic = (bc.simulation_config().thermal_bc == io::SimulationConfig::ThermalBC::Adiabatic);
  
  int start_index;
  if (is_adiabatic) {
    // Mode adiabatique : on doit calculer T[0]
    start_index = 0;
  } else {
    // Mode température imposée : T[0] = Tw fixe
    result.temperatures[0] = bc.Tw();
    start_index = 1;
  }
  
  // Calcul des températures
  for (int i = start_index; i < static_cast<int>(n_points); i++) {
    std::vector<double> c_temp(mixture_.n_species());
    for (std::size_t j = 0; j < mixture_.n_species(); j++) {
      c_temp[j] = composition(j, i);
    }
    
    double h = enthalpy_field[i] * h_e;  // Dénormaliser l'enthalpie
    
    auto temp_result = solve_single_point(c_temp, h, P, initial_temperatures[i]);
    if (!temp_result) {
      return std::unexpected(ThermodynamicSolverError(
          std::format("Failed to solve temperature at point {}: {}", i, temp_result.error().message())));
    }
    
    result.temperatures[i] = *temp_result;
  }
  
  // Si adiabatique, mettre à jour la température du mur
  if (is_adiabatic) {
    bc.wall.temperature = result.temperatures[0];  // LA LIGNE CLÉ !
    result.adiabatic_wall_updated = true;
    result.updated_wall_temperature = result.temperatures[0];
  }
  
  return result;
}

} // namespace blast::boundary_layer::thermodynamics
