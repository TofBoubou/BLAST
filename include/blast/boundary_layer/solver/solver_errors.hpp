#pragma once
#include "../../core/exceptions.hpp"
#include <source_location>
#include <string_view>

namespace blast::boundary_layer::solver {

class SolverError : public core::BlastException {
public:
  explicit SolverError(std::string_view message, std::source_location location = std::source_location::current())
      : BlastException(std::format("Solver Error: {}", message), location) {}
};

class NumericError : public SolverError {
public:
  explicit NumericError(std::string_view message, std::source_location location = std::source_location::current())
      : SolverError(std::format("Numeric Error: {}", message), location) {}
};

class GeometryError : public SolverError {
public:
  explicit GeometryError(std::string_view message, std::source_location location = std::source_location::current())
      : SolverError(std::format("Geometry Error: {}", message), location) {}
};

class RelaxationError : public SolverError {
public:
  explicit RelaxationError(std::string_view message, std::source_location location = std::source_location::current())
      : SolverError(std::format("Relaxation Error: {}", message), location) {}
  RelaxationError(std::string_view message, const SolverError& cause,
                  std::source_location location = std::source_location::current())
      : SolverError(std::format("Relaxation Error: {}: {}", message, cause.message()), location) {}
};

class StepExecutionError : public SolverError {
public:
  explicit StepExecutionError(std::string_view step_name, std::string_view message,
                              std::source_location location = std::source_location::current())
      : SolverError(std::format("Step '{}' failed: {}", step_name, message), location) {}
};

class ConvergenceError : public SolverError {
public:
  explicit ConvergenceError(std::string_view message, std::source_location location = std::source_location::current())
      : SolverError(std::format("Convergence Error: {}", message), location) {}
};

class BoundaryConditionError : public SolverError {
public:
  explicit BoundaryConditionError(std::string_view message,
                                  std::source_location location = std::source_location::current())
      : SolverError(std::format("Boundary Condition Error: {}", message), location) {}
};

class GridError : public SolverError {
public:
  explicit GridError(std::string_view message, std::source_location location = std::source_location::current())
      : SolverError(std::format("Grid Error: {}", message), location) {}
};

class InitializationError : public SolverError {
public:
  explicit InitializationError(std::string_view message,
                               std::source_location location = std::source_location::current())
      : SolverError(std::format("Initialization Error: {}", message), location) {}
};

} // namespace blast::boundary_layer::solver
