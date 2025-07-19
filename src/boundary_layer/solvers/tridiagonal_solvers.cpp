#include "blast/boundary_layer/solvers/tridiagonal_solvers.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace blast::boundary_layer::solvers {

namespace {
    // Collocation coefficients structure
    struct CollocationCoeffs {
        double p1, z, m1;  // Values at collocation points
        double a, b;       // Alpha and beta coefficients
    };
    
    // Compute collocation coefficients for a given eta point
    [[nodiscard]] auto compute_collocation_coeffs(
        double a, double b, double c, double t
    ) -> std::expected<CollocationCoeffs, SolverError> {
        CollocationCoeffs coeffs;
        
        if (t == -1.0) {
            coeffs.p1 = a - 0.5 * b;
            coeffs.z = -2.0 * a + 2.0 * b;
            coeffs.m1 = a - 1.5 * b + c;
            coeffs.a = 6.0 * a - 2.0 * b;
            coeffs.b = -10.0 * a + 2.0 * b;
        } else if (t == 0.0) {
            coeffs.p1 = a + 0.5 * b;
            coeffs.z = -2.0 * a + c;
            coeffs.m1 = a - 0.5 * b;
            coeffs.a = b;
            coeffs.b = 2.0 * a;
        } else if (t == 1.0) {
            coeffs.p1 = a + 1.5 * b + c;
            coeffs.z = -2.0 * a - 2.0 * b;
            coeffs.m1 = a + 0.5 * b;
            coeffs.a = -6.0 * a - 2.0 * b;
            coeffs.b = -10.0 * a - 2.0 * b;
        } else {
            return std::unexpected(SolverError(
                std::format("Invalid collocation point t={}, expected -1, 0, or 1", t)
            ));
        }
        
        return coeffs;
    }
}

auto solve_momentum_energy_tridiagonal(
    std::span<const double> prev_solution,
    double f_bc,
    double g_bc,
    double h_bc,
    std::span<const double> a_coeffs,
    std::span<const double> b_coeffs,
    std::span<const double> c_coeffs,
    std::span<const double> d_rhs
) -> std::expected<std::vector<double>, SolverError> {
    
    const auto n_eta = prev_solution.size();
    if (n_eta < 3) {
        return std::unexpected(SolverError("Need at least 3 eta points"));
    }
    
    // Validate coefficient array sizes
    if (a_coeffs.size() != n_eta || b_coeffs.size() != n_eta || 
        c_coeffs.size() != n_eta || d_rhs.size() != n_eta) {
        return std::unexpected(SolverError("Coefficient arrays size mismatch"));
    }
    
    // Allocate working arrays
    std::vector<double> A(n_eta), B(n_eta), C(n_eta), D(n_eta);
    
    // Fill coefficients for interior points
    for (std::size_t i = 1; i < n_eta - 1; ++i) {
        // Collocation at t = -1, 0, 1
        auto col_m1_result = compute_collocation_coeffs(a_coeffs[i-1], b_coeffs[i-1], c_coeffs[i-1], -1.0);
        if (!col_m1_result) {
            return std::unexpected(col_m1_result.error());
        }
        auto col_m1 = col_m1_result.value();
        
        auto col_0_result = compute_collocation_coeffs(a_coeffs[i], b_coeffs[i], c_coeffs[i], 0.0);
        if (!col_0_result) {
            return std::unexpected(col_0_result.error());
        }
        auto col_0 = col_0_result.value();
        
        auto col_p1_result = compute_collocation_coeffs(a_coeffs[i+1], b_coeffs[i+1], c_coeffs[i+1], 1.0);
        if (!col_p1_result) {
            return std::unexpected(col_p1_result.error());
        }
        auto col_p1 = col_p1_result.value();
        
        // Sub-determinants to eliminate alpha and beta
        const double det_m1 = col_0.a * col_p1.b - col_0.b * col_p1.a;
        const double det_0 = col_m1.b * col_p1.a - col_m1.a * col_p1.b;
        const double det_p1 = col_m1.a * col_0.b - col_m1.b * col_0.a;
        
        // System coefficients
        C[i] = col_m1.m1 * det_m1 + col_0.m1 * det_0 + col_p1.m1 * det_p1;
        B[i] = col_m1.z * det_m1 + col_0.z * det_0 + col_p1.z * det_p1;
        A[i] = col_m1.p1 * det_m1 + col_0.p1 * det_0 + col_p1.p1 * det_p1;
        D[i] = d_rhs[i-1] * det_m1 + d_rhs[i] * det_0 + d_rhs[i+1] * det_p1;
    }
    
    // Wall boundary conditions
    // BC collocation at t = -1
    CollocationCoeffs bc_col;
    bc_col.m1 = -1.5 * f_bc + g_bc;
    bc_col.z = 2.0 * f_bc;
    bc_col.p1 = -0.5 * f_bc;
    bc_col.a = -2.0 * f_bc;
    bc_col.b = 2.0 * f_bc;
    
    // Equation collocation at wall
    auto col_1_result = compute_collocation_coeffs(a_coeffs[0], b_coeffs[0], c_coeffs[0], -1.0);
    if (!col_1_result) {
        return std::unexpected(col_1_result.error());
    }
    auto col_1 = col_1_result.value();
    
    auto col_2_result = compute_collocation_coeffs(a_coeffs[1], b_coeffs[1], c_coeffs[1], 0.0);
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
    
    // Make matrix truly tridiagonal
    D[n_eta-2] -= A[n_eta-2] * prev_solution[n_eta-1];
    B[0] = R0 * A[1] - R2 * C[1];
    A[0] = R1 * A[1] - R2 * B[1];
    C[0] = 0.0;
    D[0] = R * A[1] - R2 * D[1];
    A[n_eta-2] = 0.0;
    
    // Solve tridiagonal system
    std::vector<double> solution(n_eta);
    if (auto result = detail::thomas_algorithm(C, B, A, D, solution); !result) {
        return std::unexpected(result.error());
    }
    
    return solution;
}

auto solve_species_block_tridiagonal(
    const core::Matrix<double>& prev_solution,
    std::span<const double> f_bc,
    std::span<const double> g_bc,
    std::span<const double> h_bc,
    const core::Matrix<double>& a_coeffs,
    const core::Matrix<double>& b_coeffs,
    const core::Matrix<double>& c_coeffs,
    const core::Matrix<double>& d_rhs,
    bool has_electrons
) -> std::expected<core::Matrix<double>, SolverError> {
    
    const auto n_eta = prev_solution.cols();
    const auto n_species = prev_solution.rows();
    
    if (n_eta < 3) {
        return std::unexpected(SolverError("Need at least 3 eta points"));
    }
    
    if (a_coeffs.rows() != n_eta || a_coeffs.cols() != n_species ||
        b_coeffs.rows() != n_eta || b_coeffs.cols() != n_species ||
        c_coeffs.rows() != n_eta || c_coeffs.cols() != n_species ||
        d_rhs.rows() != n_eta || d_rhs.cols() != n_species) {
        return std::unexpected(SolverError("Coefficient matrix dimensions mismatch"));
    }
    
    if (f_bc.size() != n_species || g_bc.size() != n_species || h_bc.size() != n_species) {
        return std::unexpected(SolverError("Boundary condition arrays size mismatch"));
    }
    
    const std::size_t start_idx = has_electrons ? 1 : 0;
    const std::size_t n_heavy = n_species - start_idx;
    
    if (n_heavy == 0) {
        return std::unexpected(SolverError("No heavy species to solve"));
    }
    
    std::vector<core::Matrix<double>> A(n_eta), B(n_eta), C(n_eta);
    for (auto& mat : A) mat = core::Matrix<double>(n_heavy, n_heavy);
    for (auto& mat : B) mat = core::Matrix<double>(n_heavy, n_heavy);
    for (auto& mat : C) mat = core::Matrix<double>(n_heavy, n_heavy);
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        A[i].setZero();
        B[i].setZero();
        C[i].setZero();
    }
    
    core::Matrix<double> D(n_eta, n_species);
    D.setZero();
    
    std::vector<CollocationCoeffs> col_bc(n_heavy), col_1(n_heavy), col_2(n_heavy);
    std::vector<double> det_m1(n_heavy), det_0(n_heavy), det_p1(n_heavy);
    
    for (std::size_t i = 1; i < n_eta - 1; ++i) {
        for (std::size_t j = start_idx; j < n_species; ++j) {
            const int jh = j - start_idx;
            
            auto col_m1_result = compute_collocation_coeffs(
                a_coeffs(i-1, j), b_coeffs(i-1, j), c_coeffs(i-1, j), -1.0
            );
            if (!col_m1_result) {
                return std::unexpected(col_m1_result.error());
            }
            auto col_m1 = col_m1_result.value();
            
            auto col_0_result = compute_collocation_coeffs(
                a_coeffs(i, j), b_coeffs(i, j), c_coeffs(i, j), 0.0
            );
            if (!col_0_result) {
                return std::unexpected(col_0_result.error());
            }
            auto col_0 = col_0_result.value();
            
            auto col_p1_result = compute_collocation_coeffs(
                a_coeffs(i+1, j), b_coeffs(i+1, j), c_coeffs(i+1, j), 1.0
            );
            if (!col_p1_result) {
                return std::unexpected(col_p1_result.error());
            }
            auto col_p1 = col_p1_result.value();
            
            det_m1[jh] = col_0.a * col_p1.b - col_0.b * col_p1.a;
            det_0[jh] = col_m1.b * col_p1.a - col_m1.a * col_p1.b;
            det_p1[jh] = col_m1.a * col_0.b - col_m1.b * col_0.a;
            
            C[i](jh, jh) = col_m1.m1 * det_m1[jh] + col_0.m1 * det_0[jh] + 
                           col_p1.m1 * det_p1[jh];
            B[i](jh, jh) = col_m1.z * det_m1[jh] + col_0.z * det_0[jh] + 
                           col_p1.z * det_p1[jh];
            A[i](jh, jh) = col_m1.p1 * det_m1[jh] + col_0.p1 * det_0[jh] + 
                           col_p1.p1 * det_p1[jh];
            
            D(i, j) = d_rhs(i-1, j) * det_m1[jh] + d_rhs(i, j) * det_0[jh] + 
                      d_rhs(i+1, j) * det_p1[jh];
        }
    }
    
    core::Matrix<double> R0(n_heavy, n_heavy), R1(n_heavy, n_heavy), R2(n_heavy, n_heavy);
    std::vector<double> R(n_heavy);
    R0.setZero();
    R1.setZero();
    R2.setZero();
    
    for (std::size_t j = start_idx; j < n_species; ++j) {
        const int jh = j - start_idx;
        
        col_bc[jh].m1 = -1.5 * f_bc[j] + g_bc[j];
        col_bc[jh].z = 2.0 * f_bc[j];
        col_bc[jh].p1 = -0.5 * f_bc[j];
        col_bc[jh].a = -2.0 * f_bc[j];
        col_bc[jh].b = 2.0 * f_bc[j];
        
        auto col_1_result = compute_collocation_coeffs(
            a_coeffs(0, j), b_coeffs(0, j), c_coeffs(0, j), -1.0
        );
        if (!col_1_result) {
            return std::unexpected(col_1_result.error());
        }
        col_1[jh] = col_1_result.value();
        
        auto col_2_result = compute_collocation_coeffs(
            a_coeffs(1, j), b_coeffs(1, j), c_coeffs(1, j), 0.0
        );
        if (!col_2_result) {
            return std::unexpected(col_2_result.error());
        }
        col_2[jh] = col_2_result.value();
        
        det_m1[jh] = col_1[jh].a * col_2[jh].b - col_1[jh].b * col_2[jh].a;
        det_0[jh] = col_bc[jh].b * col_2[jh].a - col_bc[jh].a * col_2[jh].b;
        det_p1[jh] = col_bc[jh].a * col_1[jh].b - col_bc[jh].b * col_1[jh].a;
        
        R0(jh, jh) = col_bc[jh].m1 * det_m1[jh] + col_1[jh].m1 * det_0[jh] + 
                     col_2[jh].m1 * det_p1[jh];
        R1(jh, jh) = col_bc[jh].z * det_m1[jh] + col_1[jh].z * det_0[jh] + 
                     col_2[jh].z * det_p1[jh];
        R2(jh, jh) = col_bc[jh].p1 * det_m1[jh] + col_1[jh].p1 * det_0[jh] + 
                     col_2[jh].p1 * det_p1[jh];
        R[jh] = h_bc[j] * det_m1[jh] + d_rhs(0, j) * det_0[jh] + 
                d_rhs(1, j) * det_p1[jh];
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
        return std::unexpected(SolverError(
            std::format("Matrix operation failed in boundary setup: {}", e.what())
        ));
    }
    
    for (std::size_t j = start_idx; j < n_species; ++j) {
        D(n_eta-2, j) -= A[n_eta-2](j-start_idx, j-start_idx) * prev_solution(j, n_eta-1);
    }
    A[n_eta-2].setZero();
    
    core::Matrix<double> solution(n_species, n_eta);
    if (auto result = detail::block_thomas_algorithm(C, B, A, D, solution, start_idx); !result) {
        return std::unexpected(result.error());
    }
    
    for (std::size_t j = 0; j < solution.cols(); ++j) {
        for (std::size_t i = 0; i < solution.rows(); ++i) {
            if (solution(i, j) < 0.0) {
                solution(i, j) = 0.0;
            }
        }
        
        double sum = 0.0;
        for (std::size_t i = 0; i < solution.rows(); ++i) {
            sum += solution(i, j);
        }
        
        if (sum > 1e-9) {
            for (std::size_t i = 0; i < solution.rows(); ++i) {
                solution(i, j) /= sum;
            }
        } else {
            std::cout << "PROBLEME FALLBACK, SUM equal to zero" << std::endl;
            abort();
        }
    }
    
    return solution;
}


namespace detail {

auto thomas_algorithm(
    std::span<double> lower_diag,
    std::span<double> main_diag,
    std::span<double> upper_diag,
    std::span<double> rhs,
    std::span<double> solution
) -> std::expected<void, SolverError> {
    const auto n = main_diag.size();
    
    // Validation des tailles
    if (lower_diag.size() != n || upper_diag.size() != n || 
        rhs.size() != n || solution.size() != n) {
        return std::unexpected(SolverError("Incompatible array sizes in Thomas algorithm"));
    }
    
    if (n < 2) {
        return std::unexpected(SolverError("System too small for Thomas algorithm"));
    }
    
    // Forward sweep
    if (std::abs(main_diag[0]) < 1e-14) {
        return std::unexpected(SolverError("Zero diagonal element at position 0"));
    }
    
    main_diag[0] = 1.0 / main_diag[0];
    lower_diag[0] = rhs[0] * main_diag[0];
    
    for (std::size_t i = 1; i < n - 1; ++i) {
        upper_diag[i-1] = upper_diag[i-1] * main_diag[i-1];
        main_diag[i] = main_diag[i] - lower_diag[i] * upper_diag[i-1];
        
        if (std::abs(main_diag[i]) < 1e-14) {
            return std::unexpected(SolverError(
                std::format("Zero diagonal element at position {}", i)
            ));
        }
        
        main_diag[i] = 1.0 / main_diag[i];
        lower_diag[i] = (rhs[i] - lower_diag[i] * lower_diag[i-1]) * main_diag[i];
    }
    
    // Backward sweep
    solution[n-2] = lower_diag[n-2];
    for (int i = n - 3; i >= 0; --i) {
        solution[i] = lower_diag[i] - upper_diag[i] * solution[i+1];
    }
    
    return {};
}

auto block_thomas_algorithm(
    std::vector<core::Matrix<double>>& lower_blocks,
    std::vector<core::Matrix<double>>& main_blocks,
    std::vector<core::Matrix<double>>& upper_blocks,
    core::Matrix<double>& rhs,
    core::Matrix<double>& solution,
    int start_species
) -> std::expected<void, SolverError> {
    
    const auto n_eta = main_blocks.size();
    const auto n_heavy = main_blocks[0].rows();
    
    if (n_eta < 2) {
        return std::unexpected(SolverError("System too small for block Thomas algorithm"));
    }
    
    if (lower_blocks.size() != n_eta || upper_blocks.size() != n_eta) {
        return std::unexpected(SolverError("Incompatible block matrix sizes"));
    }
    
    if (n_heavy == 0) {
        return std::unexpected(SolverError("Empty block matrices"));
    }
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        if (main_blocks[i].rows() != n_heavy || main_blocks[i].cols() != n_heavy ||
            lower_blocks[i].rows() != n_heavy || lower_blocks[i].cols() != n_heavy ||
            upper_blocks[i].rows() != n_heavy || upper_blocks[i].cols() != n_heavy) {
            return std::unexpected(SolverError(
                std::format("Inconsistent block matrix dimensions at position {}", i)
            ));
        }
    }
    
    try {
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_B0(main_blocks[0].eigen());
        
        if (lu_B0.info() != Eigen::Success) {
            return std::unexpected(SolverError("Singular matrix B[0] in initial conditions"));
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
        return std::unexpected(SolverError("Failed to solve boundary conditions at i=0"));
    }
    
    for (std::size_t i = 1; i < n_eta; ++i) {
        try {
            auto product = lower_blocks[i].eigen() * upper_blocks[i-1].eigen();
            Eigen::MatrixXd new_main = main_blocks[i].eigen() - product;
            main_blocks[i].eigen() = new_main;
            
            for (std::size_t j = 0; j < n_heavy; ++j) {
                double temp = 0.0;
                for (std::size_t k = 0; k < n_heavy; ++k) {
                    temp += lower_blocks[i](j, k) * solution(k + start_species, i-1);
                }
                rhs(i, j + start_species) -= temp;
            }
            
            if (i < n_eta - 1) {
                Eigen::PartialPivLU<Eigen::MatrixXd> lu_B(main_blocks[i].eigen());
                
                if (lu_B.info() != Eigen::Success) {
                    return std::unexpected(SolverError(
                        std::format("Singular matrix at block position {}", i)
                    ));
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
                
            } else {
                Eigen::PartialPivLU<Eigen::MatrixXd> lu_B(main_blocks[i].eigen());
                
                if (lu_B.info() != Eigen::Success) {
                    return std::unexpected(SolverError(
                        std::format("Singular matrix at final block position {}", i)
                    ));
                }
                
                Eigen::VectorXd rhs_vec(n_heavy);
                for (std::size_t k = 0; k < n_heavy; ++k) {
                    rhs_vec[k] = rhs(i, k + start_species);
                }
                
                Eigen::VectorXd cy = lu_B.solve(rhs_vec);
                
                for (std::size_t j = 0; j < n_heavy; ++j) {
                    solution(j + start_species, i) = cy[j];
                }
            }
            
        } catch (const std::exception& e) {
            return std::unexpected(SolverError(
                std::format("Matrix operation failed at block {}: {}", i, e.what())
            ));
        }
    }
    
    try {
        for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n_eta) - 2; i >= 0; --i) {
            for (std::size_t j = 0; j < n_heavy; ++j) {
                double correction = 0.0;
                for (std::size_t k = 0; k < n_heavy; ++k) {
                    correction += upper_blocks[i](j, k) * solution(k + start_species, i+1);
                }
                
                double cy_val = solution(j + start_species, i);
                double new_val = cy_val - correction;
                solution(j + start_species, i) = new_val;
            }
        }
        
    } catch (const std::exception& e) {
        return std::unexpected(SolverError(
            std::format("Back substitution failed: {}", e.what())
        ));
    }
    
    return {};
}

} // namespace detail

} // namespace blast::boundary_layer::solvers


