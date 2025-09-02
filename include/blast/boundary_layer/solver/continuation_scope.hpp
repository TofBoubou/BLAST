#pragma once

#include "boundary_layer_solver.hpp"

namespace blast::boundary_layer::solver {

class ContinuationScope {
public:
  explicit ContinuationScope(BoundaryLayerSolver& solver) noexcept : solver_(solver), active_(true) {
    solver_.enter_continuation();
  }

  ContinuationScope(const ContinuationScope&) = delete;
  ContinuationScope& operator=(const ContinuationScope&) = delete;

  ContinuationScope(ContinuationScope&& other) noexcept : solver_(other.solver_), active_(other.active_) {
    other.active_ = false;
  }
  ContinuationScope& operator=(ContinuationScope&& other) = delete;

  ~ContinuationScope() {
    if (active_) solver_.exit_continuation();
  }

private:
  BoundaryLayerSolver& solver_;
  bool active_;
};

} // namespace blast::boundary_layer::solver

