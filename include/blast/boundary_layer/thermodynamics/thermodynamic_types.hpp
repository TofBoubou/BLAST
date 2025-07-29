#pragma once
#include "../../core/containers.hpp"
#include "../../core/exceptions.hpp"
#include <expected>
#include <vector>

namespace blast::boundary_layer::thermodynamics {

struct EnthalpyTemperatureSolverConfig {
  double tolerance = 1e-8;
  int max_iterations = 100;
  double relaxation_factor = 1.0;
  double min_temperature = 1.0;
  double max_temperature = 10000.0;
  int max_bracket_expansions = 10;
};

struct TemperatureField {
  std::vector<double> temperatures;
  bool adiabatic_wall_updated = false;
  double updated_wall_temperature = 0.0;
};

class ThermodynamicSolverError : public core::BlastException {
public:
  explicit ThermodynamicSolverError(std::string_view message,
                                    std::source_location location = std::source_location::current())
      : BlastException(std::format("Thermodynamic Solver Error: {}", message), location) {}

  template <typename... Args>
  explicit ThermodynamicSolverError(std::string_view format_str, std::source_location location, Args&&... args)
      : BlastException(format_error_message(format_str, args...), location) {}

private:
  template <typename... Args> static std::string format_error_message(std::string_view format_str, Args&&... args) {
    try {
      return std::format("Thermodynamic Solver Error: {}", std::vformat(format_str, std::make_format_args(args...)));
    } catch (...) {
      return "Thermodynamic Solver Error: Exception formatting failed";
    }
  }
};

} // namespace blast::boundary_layer::thermodynamics