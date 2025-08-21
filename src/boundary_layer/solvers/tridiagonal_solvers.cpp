#include "blast/boundary_layer/solvers/tridiagonal_solvers.hpp"
#include "blast/boundary_layer/solver/solver_errors.hpp"
#include "blast/core/constants.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>

namespace blast::boundary_layer::solvers {

using namespace blast::constants::numerical_methods::tridiagonal_solver;

namespace {
// Collocation coefficients structure
struct CollocationCoeffs {
  double p1, z, m1; // Values at collocation points
  double a, b;      // Alpha and beta coefficients
};

// Helper template for monadic error handling
template<typename T, typename Error>
[[nodiscard]] constexpr auto try_unwrap(const std::expected<T, Error>& exp) -> std::expected<T, Error> {
  return exp;
}

// Compute collocation coefficients for a given eta point
[[nodiscard]] auto compute_collocation_coeffs(double a, double b, double c, double t)
    -> std::expected<CollocationCoeffs, SolverError> {
  CollocationCoeffs coeffs;

  if (t == left_point) {
    coeffs.p1 = a - half_coeff * b;
    coeffs.z = minus_two_coeff * a + two_coeff * b;
    coeffs.m1 = a - one_half_coeff * b + c;
    coeffs.a = six_coeff * a - two_coeff * b;
    coeffs.b = minus_ten_coeff * a + two_coeff * b;
  } else if (t == center_point) {
    coeffs.p1 = a + half_coeff * b;
    coeffs.z = minus_two_coeff * a + c;
    coeffs.m1 = a - half_coeff * b;
    coeffs.a = b;
    coeffs.b = two_coeff * a;
  } else if (t == right_point) {
    coeffs.p1 = a + one_half_coeff * b + c;
    coeffs.z = minus_two_coeff * a - two_coeff * b;
    coeffs.m1 = a + half_coeff * b;
    coeffs.a = minus_six_coeff * a - two_coeff * b;
    coeffs.b = minus_ten_coeff * a - two_coeff * b;
  } else {
    return std::unexpected(SolverError(std::format("Invalid collocation point t={}, expected {}, {}, or {}", 
                                                  t, left_point, center_point, right_point)));
  }

  return coeffs;
}

// Helper to compute three collocation coefficients with error propagation
[[nodiscard]] auto compute_three_collocation_coeffs(
    double a_m1, double b_m1, double c_m1,
    double a_0, double b_0, double c_0, 
    double a_p1, double b_p1, double c_p1)
    -> std::expected<std::tuple<CollocationCoeffs, CollocationCoeffs, CollocationCoeffs>, SolverError> {
  
  auto col_m1_result = compute_collocation_coeffs(a_m1, b_m1, c_m1, left_point);
  if (!col_m1_result) return std::unexpected(col_m1_result.error());
  
  auto col_0_result = compute_collocation_coeffs(a_0, b_0, c_0, center_point);
  if (!col_0_result) return std::unexpected(col_0_result.error());
  
  auto col_p1_result = compute_collocation_coeffs(a_p1, b_p1, c_p1, right_point);
  if (!col_p1_result) return std::unexpected(col_p1_result.error());
  
  return std::make_tuple(col_m1_result.value(), col_0_result.value(), col_p1_result.value());
}

// Matrix allocation helper for efficient memory management
template<typename MatrixType>
class MatrixBlockManager {
public:
  explicit MatrixBlockManager(std::size_t n_blocks, std::size_t rows, std::size_t cols) {
    blocks_.reserve(n_blocks);
    for (std::size_t i = 0; i < n_blocks; ++i) {
      blocks_.emplace_back(rows, cols);
    }
  }
  
  [[nodiscard]] auto& operator[](std::size_t index) { return blocks_[index]; }
  [[nodiscard]] const auto& operator[](std::size_t index) const { return blocks_[index]; }
  [[nodiscard]] auto size() const -> std::size_t { return blocks_.size(); }
  
private:
  std::vector<MatrixType> blocks_;
};

using MatrixManager = MatrixBlockManager<core::Matrix<double>>;

} // namespace

// Helper function to compute interior point coefficients
[[nodiscard]] auto compute_interior_coefficients(
    std::size_t i, 
    std::span<const double> a_coeffs, 
    std::span<const double> b_coeffs,
    std::span<const double> c_coeffs, 
    std::span<const double> d_rhs,
    std::span<double> A, 
    std::span<double> B, 
    std::span<double> C, 
    std::span<double> D) -> std::expected<void, SolverError> {
  
  auto coeffs_result = compute_three_collocation_coeffs(
      a_coeffs[i - 1], b_coeffs[i - 1], c_coeffs[i - 1],
      a_coeffs[i], b_coeffs[i], c_coeffs[i],
      a_coeffs[i + 1], b_coeffs[i + 1], c_coeffs[i + 1]);
  
  if (!coeffs_result) {
    return std::unexpected(coeffs_result.error());
  }
  
  auto [col_m1, col_0, col_p1] = coeffs_result.value();

  // Sub-determinants to eliminate alpha and beta
  const double det_m1 = col_0.a * col_p1.b - col_0.b * col_p1.a;
  const double det_0 = col_m1.b * col_p1.a - col_m1.a * col_p1.b;
  const double det_p1 = col_m1.a * col_0.b - col_m1.b * col_0.a;

  // System coefficients
  C[i] = col_m1.m1 * det_m1 + col_0.m1 * det_0 + col_p1.m1 * det_p1;
  B[i] = col_m1.z * det_m1 + col_0.z * det_0 + col_p1.z * det_p1;
  A[i] = col_m1.p1 * det_m1 + col_0.p1 * det_0 + col_p1.p1 * det_p1;
  D[i] = d_rhs[i - 1] * det_m1 + d_rhs[i] * det_0 + d_rhs[i + 1] * det_p1;
  
  return {};
}

// Helper function to setup boundary conditions
[[nodiscard]] auto setup_boundary_conditions(
    double f_bc, double g_bc, double h_bc,
    std::span<const double> a_coeffs, 
    std::span<const double> b_coeffs,
    std::span<const double> c_coeffs, 
    std::span<const double> d_rhs,
    std::span<double> A, 
    std::span<double> B, 
    std::span<double> C, 
    std::span<double> D) -> std::expected<void, SolverError> {
  
  // Wall boundary conditions - BC collocation at t = left_point
  CollocationCoeffs bc_col;
  bc_col.m1 = minus_one_half_coeff * f_bc + g_bc;
  bc_col.z = two_coeff * f_bc;
  bc_col.p1 = minus_half_coeff * f_bc;
  bc_col.a = minus_two_coeff * f_bc;
  bc_col.b = two_coeff * f_bc;

  // Equation collocation at wall
  auto col_1_result = compute_collocation_coeffs(a_coeffs[0], b_coeffs[0], c_coeffs[0], left_point);
  if (!col_1_result) {
    return std::unexpected(col_1_result.error());
  }
  auto col_1 = col_1_result.value();

  auto col_2_result = compute_collocation_coeffs(a_coeffs[1], b_coeffs[1], c_coeffs[1], center_point);
  if (!col_2_result) {
    return std::unexpected(col_2_result.error());
  }
  auto col_2 = col_2_result.value();

  // Sub-determinants for wall
  const double det_m1 = col_1.a * col_2.b - col_1.b * col_2.a;
  const double det_0 = bc_col.b * col_2.a - bc_col.a * col_2.b;
  const double det_p1 = bc_col.a * col_1.b - bc_col.b * col_1.a;

  // Wall coefficients
  const double R0 = bc_col.m1 * det_m1 + col_1.m1 * det_0 + col_2.m1 * det_p1;
  const double R1 = bc_col.z * det_m1 + col_1.z * det_0 + col_2.z * det_p1;
  const double R2 = bc_col.p1 * det_m1 + col_1.p1 * det_0 + col_2.p1 * det_p1;
  const double R = h_bc * det_m1 + d_rhs[0] * det_0 + d_rhs[1] * det_p1;

  // Set boundary coefficient values
  B[0] = R0 * A[1] - R2 * C[1];
  A[0] = R1 * A[1] - R2 * B[1];
  C[0] = 0.0;
  D[0] = R * A[1] - R2 * D[1];
  
  return {};
}

auto solve_momentum_energy_tridiagonal(std::span<const double> prev_solution, double f_bc, double g_bc, double h_bc,
                                       std::span<const double> a_coeffs, std::span<const double> b_coeffs,
                                       std::span<const double> c_coeffs, std::span<const double> d_rhs)
    -> std::expected<std::vector<double>, SolverError> {

  const auto n_eta = prev_solution.size();
  if (n_eta < 3) {
    return std::unexpected(SolverError("Need at least 3 eta points"));
  }

  // Validate coefficient array sizes
  if (a_coeffs.size() != n_eta || b_coeffs.size() != n_eta || c_coeffs.size() != n_eta || d_rhs.size() != n_eta) {
    return std::unexpected(SolverError("Coefficient arrays size mismatch"));
  }

  // Allocate working arrays
  std::vector<double> A(n_eta), B(n_eta), C(n_eta), D(n_eta);

  // Fill coefficients for interior points
  for (std::size_t i = 1; i < n_eta - 1; ++i) {
    if (auto result = compute_interior_coefficients(i, a_coeffs, b_coeffs, c_coeffs, d_rhs, A, B, C, D); !result) {
      return std::unexpected(result.error());
    }
  }

  // Setup boundary conditions
  if (auto result = setup_boundary_conditions(f_bc, g_bc, h_bc, a_coeffs, b_coeffs, c_coeffs, d_rhs, A, B, C, D); !result) {
    return std::unexpected(result.error());
  }

  // Make matrix truly tridiagonal
  D[n_eta - 2] -= A[n_eta - 2] * prev_solution[n_eta - 1];
  A[n_eta - 2] = 0.0;

  // Solve tridiagonal system
  std::vector<double> solution(n_eta);
  if (auto result = detail::thomas_algorithm(C, B, A, D, solution); !result) {
    return std::unexpected(result.error());
  }

  return solution;
}

// Helper to validate species block input parameters
[[nodiscard]] auto validate_species_block_inputs(
    const core::Matrix<double>& prev_solution,
    std::span<const double> f_bc, std::span<const double> g_bc, std::span<const double> h_bc,
    const core::Matrix<double>& a_coeffs, const core::Matrix<double>& b_coeffs,
    const core::Matrix<double>& c_coeffs, const core::Matrix<double>& d_rhs) -> std::expected<void, SolverError> {
  
  const auto n_eta = prev_solution.cols();
  const auto n_species = prev_solution.rows();

  if (n_eta < 3) {
    return std::unexpected(SolverError("Need at least 3 eta points"));
  }

  if (a_coeffs.rows() != n_eta || a_coeffs.cols() != n_species || b_coeffs.rows() != n_eta ||
      b_coeffs.cols() != n_species || c_coeffs.rows() != n_eta || c_coeffs.cols() != n_species ||
      d_rhs.rows() != n_eta || d_rhs.cols() != n_species) {
    return std::unexpected(SolverError("Coefficient matrix dimensions mismatch"));
  }

  if (f_bc.size() != n_species || g_bc.size() != n_species || h_bc.size() != n_species) {
    return std::unexpected(SolverError("Boundary condition arrays size mismatch"));
  }

  return {};
}

// Helper to setup species matrix coefficients for interior points
[[nodiscard]] auto setup_species_interior_matrices(
    std::size_t n_eta, std::size_t n_species, std::size_t start_idx, std::size_t n_heavy,
    const core::Matrix<double>& a_coeffs, const core::Matrix<double>& b_coeffs,
    const core::Matrix<double>& c_coeffs, const core::Matrix<double>& d_rhs,
    MatrixManager& A, MatrixManager& B, MatrixManager& C, core::Matrix<double>& D) -> std::expected<void, SolverError> {

  for (std::size_t i = 1; i < n_eta - 1; ++i) {
    for (std::size_t j = start_idx; j < n_species; ++j) {
      const int jh = j - start_idx;

      auto coeffs_result = compute_three_collocation_coeffs(
          a_coeffs(i - 1, j), b_coeffs(i - 1, j), c_coeffs(i - 1, j),
          a_coeffs(i, j), b_coeffs(i, j), c_coeffs(i, j),
          a_coeffs(i + 1, j), b_coeffs(i + 1, j), c_coeffs(i + 1, j));
      
      if (!coeffs_result) {
        return std::unexpected(coeffs_result.error());
      }
      
      auto [col_m1, col_0, col_p1] = coeffs_result.value();

      const double det_m1 = col_0.a * col_p1.b - col_0.b * col_p1.a;
      const double det_0 = col_m1.b * col_p1.a - col_m1.a * col_p1.b;
      const double det_p1 = col_m1.a * col_0.b - col_m1.b * col_0.a;

      C[i](jh, jh) = col_m1.m1 * det_m1 + col_0.m1 * det_0 + col_p1.m1 * det_p1;
      B[i](jh, jh) = col_m1.z * det_m1 + col_0.z * det_0 + col_p1.z * det_p1;
      A[i](jh, jh) = col_m1.p1 * det_m1 + col_0.p1 * det_0 + col_p1.p1 * det_p1;

      D(i, j) = d_rhs(i - 1, j) * det_m1 + d_rhs(i, j) * det_0 + d_rhs(i + 1, j) * det_p1;
    }
  }
  
  return {};
}

// Helper to setup species boundary conditions
[[nodiscard]] auto setup_species_boundary_matrices(
    std::size_t n_species, std::size_t start_idx, std::size_t n_heavy,
    std::span<const double> f_bc, std::span<const double> g_bc, std::span<const double> h_bc,
    const core::Matrix<double>& a_coeffs, const core::Matrix<double>& b_coeffs,
    const core::Matrix<double>& c_coeffs, const core::Matrix<double>& d_rhs,
    MatrixManager& A, MatrixManager& B, MatrixManager& C, core::Matrix<double>& D) -> std::expected<void, SolverError> {

  core::Matrix<double> R0(n_heavy, n_heavy), R1(n_heavy, n_heavy), R2(n_heavy, n_heavy);
  std::vector<double> R(n_heavy);
  R0.setZero();
  R1.setZero();
  R2.setZero();

  for (std::size_t j = start_idx; j < n_species; ++j) {
    const int jh = j - start_idx;

    // Boundary condition coefficients using named constants
    CollocationCoeffs col_bc;
    col_bc.m1 = minus_one_half_coeff * f_bc[j] + g_bc[j];
    col_bc.z = two_coeff * f_bc[j];
    col_bc.p1 = minus_half_coeff * f_bc[j];
    col_bc.a = minus_two_coeff * f_bc[j];
    col_bc.b = two_coeff * f_bc[j];

    auto col_1_result = compute_collocation_coeffs(a_coeffs(0, j), b_coeffs(0, j), c_coeffs(0, j), left_point);
    if (!col_1_result) {
      return std::unexpected(col_1_result.error());
    }
    auto col_1 = col_1_result.value();

    auto col_2_result = compute_collocation_coeffs(a_coeffs(1, j), b_coeffs(1, j), c_coeffs(1, j), center_point);
    if (!col_2_result) {
      return std::unexpected(col_2_result.error());
    }
    auto col_2 = col_2_result.value();

    const double det_m1 = col_1.a * col_2.b - col_1.b * col_2.a;
    const double det_0 = col_bc.b * col_2.a - col_bc.a * col_2.b;
    const double det_p1 = col_bc.a * col_1.b - col_bc.b * col_1.a;

    R0(jh, jh) = col_bc.m1 * det_m1 + col_1.m1 * det_0 + col_2.m1 * det_p1;
    R1(jh, jh) = col_bc.z * det_m1 + col_1.z * det_0 + col_2.z * det_p1;
    R2(jh, jh) = col_bc.p1 * det_m1 + col_1.p1 * det_0 + col_2.p1 * det_p1;
    R[jh] = h_bc[j] * det_m1 + d_rhs(0, j) * det_0 + d_rhs(1, j) * det_p1;
  }

  try {
    Eigen::PartialPivLU<Eigen::MatrixXd> lu_A1(A[1].eigen());
    if (lu_A1.info() != Eigen::Success) {
      return std::unexpected(SolverError("Singular matrix A[1] in boundary condition setup"));
    }

    auto temp_solve_C = lu_A1.solve(C[1].eigen());
    auto temp_solve_B = lu_A1.solve(B[1].eigen());

    B[0].eigen() = (R0.eigen() - R2.eigen() * temp_solve_C).eval();
    A[0].eigen() = (R1.eigen() - R2.eigen() * temp_solve_B).eval();

    Eigen::VectorXd D1(n_heavy);
    for (std::size_t j = 0; j < n_heavy; ++j) {
      D1[j] = D(1, j + start_idx);
    }
    Eigen::VectorXd temp = R2.eigen() * lu_A1.solve(D1);

    for (std::size_t j = 0; j < n_heavy; ++j) {
      D(0, j + start_idx) = R[j] - temp[j];
    }
  } catch (const std::exception& e) {
    return std::unexpected(SolverError(std::format("Matrix operation failed in boundary setup: {}", e.what())));
  }

  return {};
}

// Helper to normalize solution (ensure non-negative and sum to 1)
void normalize_species_solution(core::Matrix<double>& solution) {
  for (std::size_t j = 0; j < solution.cols(); ++j) {
    // Ensure non-negative values
    for (std::size_t i = 0; i < solution.rows(); ++i) {
      if (solution(i, j) < 0.0) {
        solution(i, j) = 0.0;
      }
    }

    // Normalize to sum = 1
    double sum = 0.0;
    for (std::size_t i = 0; i < solution.rows(); ++i) {
      sum += solution(i, j);
    }

    if (sum > solution_normalization_tolerance) {
      for (std::size_t i = 0; i < solution.rows(); ++i) {
        solution(i, j) /= sum;
      }
    }
  }
}

auto solve_species_block_tridiagonal(const core::Matrix<double>& prev_solution, std::span<const double> f_bc,
                                     std::span<const double> g_bc, std::span<const double> h_bc,
                                     const core::Matrix<double>& a_coeffs, const core::Matrix<double>& b_coeffs,
                                     const core::Matrix<double>& c_coeffs, const core::Matrix<double>& d_rhs,
                                     bool has_electrons) -> std::expected<core::Matrix<double>, SolverError> {

  const auto n_eta = prev_solution.cols();
  const auto n_species = prev_solution.rows();

  // Handle trivial single species case
  if (n_species == 1) {
    core::Matrix<double> result(1, n_eta);
    result.eigen().setOnes();
    return result;
  }

  // Validate inputs
  if (auto validation_result = validate_species_block_inputs(prev_solution, f_bc, g_bc, h_bc, a_coeffs, b_coeffs, c_coeffs, d_rhs); 
      !validation_result) {
    return std::unexpected(validation_result.error());
  }

  const std::size_t start_idx = has_electrons ? 1 : 0;
  const std::size_t n_heavy = n_species - start_idx;

  if (n_heavy == 0) {
    return std::unexpected(SolverError("No heavy species to solve"));
  }

  // Efficiently allocate matrices using the manager
  MatrixManager A(n_eta, n_heavy, n_heavy);
  MatrixManager B(n_eta, n_heavy, n_heavy);
  MatrixManager C(n_eta, n_heavy, n_heavy);

  for (std::size_t i = 0; i < n_eta; ++i) {
    A[i].setZero();
    B[i].setZero();
    C[i].setZero();
  }

  core::Matrix<double> D(n_eta, n_species);
  D.setZero();

  // Setup interior point matrices
  if (auto result = setup_species_interior_matrices(n_eta, n_species, start_idx, n_heavy, a_coeffs, b_coeffs, c_coeffs, d_rhs, A, B, C, D); 
      !result) {
    return std::unexpected(result.error());
  }

  // Setup boundary condition matrices
  if (auto result = setup_species_boundary_matrices(n_species, start_idx, n_heavy, f_bc, g_bc, h_bc, a_coeffs, b_coeffs, c_coeffs, d_rhs, A, B, C, D); 
      !result) {
    return std::unexpected(result.error());
  }

  // Apply far-field boundary condition
  for (std::size_t j = start_idx; j < n_species; ++j) {
    D(n_eta - 2, j) -= A[n_eta - 2](j - start_idx, j - start_idx) * prev_solution(j, n_eta - 1);
  }
  A[n_eta - 2].setZero();

  // Solve the block tridiagonal system
  core::Matrix<double> solution(n_species, n_eta);
  solution.setZero();
  
  std::vector<core::Matrix<double>> A_vec, B_vec, C_vec;
  A_vec.reserve(A.size());
  B_vec.reserve(B.size());
  C_vec.reserve(C.size());
  
  for (std::size_t i = 0; i < A.size(); ++i) {
    A_vec.emplace_back(std::move(A[i]));
    B_vec.emplace_back(std::move(B[i]));
    C_vec.emplace_back(std::move(C[i]));
  }
  
  if (auto result = detail::block_thomas_algorithm(C_vec, B_vec, A_vec, D, solution, start_idx); !result) {
    return std::unexpected(result.error());
  }

  // Normalize the solution
  normalize_species_solution(solution);

  return solution;
}

namespace detail {

// Helper for Thomas algorithm forward sweep
[[nodiscard]] auto thomas_forward_sweep(
    std::span<double> lower_diag, std::span<double> main_diag, 
    std::span<double> upper_diag, std::span<double> rhs) -> std::expected<void, SolverError> {
  
  const auto n = main_diag.size();
  
  // Initial step
  if (std::abs(main_diag[0]) < diagonal_tolerance) {
    return std::unexpected(SolverError("Zero diagonal element at position 0"));
  }

  main_diag[0] = 1.0 / main_diag[0];
  lower_diag[0] = rhs[0] * main_diag[0];

  // Forward elimination
  for (std::size_t i = 1; i < n - 1; ++i) {
    upper_diag[i - 1] = upper_diag[i - 1] * main_diag[i - 1];
    main_diag[i] = main_diag[i] - lower_diag[i] * upper_diag[i - 1];

    if (std::abs(main_diag[i]) < diagonal_tolerance) {
      return std::unexpected(SolverError(std::format("Zero diagonal element at position {}", i)));
    }

    main_diag[i] = 1.0 / main_diag[i];
    lower_diag[i] = (rhs[i] - lower_diag[i] * lower_diag[i - 1]) * main_diag[i];
  }
  
  return {};
}

// Helper for Thomas algorithm backward substitution
void thomas_backward_substitution(
    std::span<const double> lower_diag, std::span<const double> upper_diag, 
    std::span<double> solution) {
  
  const auto n = solution.size();
  
  solution[n - 2] = lower_diag[n - 2];
  for (int i = n - 3; i >= 0; --i) {
    solution[i] = lower_diag[i] - upper_diag[i] * solution[i + 1];
  }
}

auto thomas_algorithm(std::span<double> lower_diag, std::span<double> main_diag, std::span<double> upper_diag,
                      std::span<double> rhs, std::span<double> solution) -> std::expected<void, SolverError> {
  const auto n = main_diag.size();

  if (lower_diag.size() != n || upper_diag.size() != n || rhs.size() != n || solution.size() != n) {
    return std::unexpected(SolverError("Incompatible array sizes in Thomas algorithm"));
  }

  if (n < 2) {
    return std::unexpected(SolverError("System too small for Thomas algorithm"));
  }

  // Forward sweep
  if (auto result = thomas_forward_sweep(lower_diag, main_diag, upper_diag, rhs); !result) {
    return std::unexpected(result.error());
  }

  // Backward substitution
  thomas_backward_substitution(lower_diag, upper_diag, solution);

  return {};
}

// Helper for block Thomas initial step
[[nodiscard]] auto block_thomas_initial_step(
    std::vector<core::Matrix<double>>& main_blocks,
    std::vector<core::Matrix<double>>& upper_blocks,
    core::Matrix<double>& rhs, core::Matrix<double>& solution,
    std::size_t n_heavy, int start_species) -> std::expected<void, SolverError> {

  try {
    Eigen::PartialPivLU<Eigen::MatrixXd> lu_B0(main_blocks[0].eigen());
    if (lu_B0.info() != Eigen::Success) {
      return std::unexpected(SolverError("Singular matrix B[0]"));
    }

    Eigen::MatrixXd beta_0 = lu_B0.solve(upper_blocks[0].eigen());
    upper_blocks[0].eigen() = beta_0;

    Eigen::VectorXd rhs_0(n_heavy);
    for (std::size_t j = 0; j < n_heavy; ++j) {
      rhs_0[j] = rhs(0, j + start_species);
    }

    Eigen::VectorXd cy_0 = lu_B0.solve(rhs_0);

    for (std::size_t j = 0; j < n_heavy; ++j) {
      solution(j + start_species, 0) = cy_0[j];
    }
  } catch (const std::exception& e) {
    return std::unexpected(SolverError(std::format("Matrix operation failed in initial step: {}", e.what())));
  }

  return {};
}

// Helper for block Thomas forward elimination step
[[nodiscard]] auto block_thomas_forward_step(
    std::size_t i, std::vector<core::Matrix<double>>& lower_blocks,
    std::vector<core::Matrix<double>>& main_blocks,
    std::vector<core::Matrix<double>>& upper_blocks,
    core::Matrix<double>& rhs, core::Matrix<double>& solution,
    std::size_t n_heavy, int start_species) -> std::expected<void, SolverError> {

  try {
    auto product = lower_blocks[i].eigen() * upper_blocks[i - 1].eigen();
    Eigen::MatrixXd new_main = main_blocks[i].eigen() - product;
    main_blocks[i].eigen() = new_main;

    for (std::size_t j = 0; j < n_heavy; ++j) {
      double temp = 0.0;
      for (std::size_t k = 0; k < n_heavy; ++k) {
        temp += lower_blocks[i](j, k) * solution(k + start_species, i - 1);
      }
      rhs(i, j + start_species) -= temp;
    }

    Eigen::PartialPivLU<Eigen::MatrixXd> lu_B(main_blocks[i].eigen());
    if (lu_B.info() != Eigen::Success) {
      return std::unexpected(SolverError(std::format("Singular matrix at step {}", i)));
    }

    Eigen::MatrixXd beta_i = lu_B.solve(upper_blocks[i].eigen());
    upper_blocks[i].eigen() = beta_i;

    Eigen::VectorXd rhs_vec(n_heavy);
    for (std::size_t k = 0; k < n_heavy; ++k) {
      rhs_vec[k] = rhs(i, k + start_species);
    }

    Eigen::VectorXd cy = lu_B.solve(rhs_vec);

    for (std::size_t j = 0; j < n_heavy; ++j) {
      solution(j + start_species, i) = cy[j];
    }
  } catch (const std::exception& e) {
    return std::unexpected(SolverError(std::format("Matrix operation failed at forward step {}: {}", i, e.what())));
  }

  return {};
}

// Helper for block Thomas back substitution
void block_thomas_back_substitution(
    const std::vector<core::Matrix<double>>& upper_blocks,
    core::Matrix<double>& solution, std::size_t n_eta,
    std::size_t n_heavy, int start_species) {

  for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n_eta) - 2; i >= 0; --i) {
    for (std::size_t j = 0; j < n_heavy; ++j) {
      double correction = 0.0;
      for (std::size_t k = 0; k < n_heavy; ++k) {
        correction += upper_blocks[i](j, k) * solution(k + start_species, i + 1);
      }

      double cy_val = solution(j + start_species, i);
      solution(j + start_species, i) = cy_val - correction;
    }
  }
}

auto block_thomas_algorithm(std::vector<core::Matrix<double>>& lower_blocks,
                            std::vector<core::Matrix<double>>& main_blocks,
                            std::vector<core::Matrix<double>>& upper_blocks, core::Matrix<double>& rhs,
                            core::Matrix<double>& solution, int start_species) -> std::expected<void, SolverError> {

  const auto n_eta = main_blocks.size();
  const auto n_heavy = main_blocks[0].rows();

  if (n_eta < 2) {
    return std::unexpected(SolverError("System too small"));
  }

  // Initial step (i=0)
  if (auto result = block_thomas_initial_step(main_blocks, upper_blocks, rhs, solution, n_heavy, start_species); !result) {
    return std::unexpected(result.error());
  }

  // Forward elimination (i=1 to n_eta-2)
  for (std::size_t i = 1; i < n_eta - 1; ++i) {
    if (auto result = block_thomas_forward_step(i, lower_blocks, main_blocks, upper_blocks, rhs, solution, n_heavy, start_species); !result) {
      return std::unexpected(result.error());
    }
  }

  // Final step (i = n_eta - 1)
  const std::size_t final_i = n_eta - 1;
  if (auto result = block_thomas_forward_step(final_i, lower_blocks, main_blocks, upper_blocks, rhs, solution, n_heavy, start_species); !result) {
    return std::unexpected(result.error());
  }

  // Back substitution
  block_thomas_back_substitution(upper_blocks, solution, n_eta, n_heavy, start_species);

  return {};
}

} // namespace detail

} // namespace blast::boundary_layer::solvers
