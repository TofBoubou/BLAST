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
    
    // Validate matrix dimensions
    if (a_coeffs.rows() != n_eta || a_coeffs.cols() != n_species ||
        b_coeffs.rows() != n_eta || b_coeffs.cols() != n_species ||
        c_coeffs.rows() != n_eta || c_coeffs.cols() != n_species ||
        d_rhs.rows() != n_eta || d_rhs.cols() != n_species) {
        return std::unexpected(SolverError("Coefficient matrix dimensions mismatch"));
    }
    
    // Validate boundary condition array sizes
    if (f_bc.size() != n_species || g_bc.size() != n_species || h_bc.size() != n_species) {
        return std::unexpected(SolverError("Boundary condition arrays size mismatch"));
    }
    
    const std::size_t start_idx = has_electrons ? 1 : 0;
    const std::size_t n_heavy = n_species - start_idx;
    
    if (n_heavy == 0) {
        return std::unexpected(SolverError("No heavy species to solve"));
    }
    
    // Allocate block matrices
    std::vector<core::Matrix<double>> A(n_eta), B(n_eta), C(n_eta);
    for (auto& mat : A) mat = core::Matrix<double>(n_heavy, n_heavy);
    for (auto& mat : B) mat = core::Matrix<double>(n_heavy, n_heavy);
    for (auto& mat : C) mat = core::Matrix<double>(n_heavy, n_heavy);
    
    core::Matrix<double> D(n_eta, n_species);
    
    // Working arrays for wall BC
    std::vector<CollocationCoeffs> col_bc(n_heavy), col_1(n_heavy), col_2(n_heavy);
    std::vector<double> det_m1(n_heavy), det_0(n_heavy), det_p1(n_heavy);
    
    // Fill coefficients for interior points
    for (std::size_t i = 1; i < n_eta - 1; ++i) {
        A[i].setZero();
        B[i].setZero();
        C[i].setZero();
        
        // Process each species equation
        for (std::size_t j = start_idx; j < n_species; ++j) {
            const int jh = j - start_idx;  // Heavy species index
            
            // Collocation coefficients
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
            
            // Sub-determinants
            det_m1[jh] = col_0.a * col_p1.b - col_0.b * col_p1.a;
            det_0[jh] = col_m1.b * col_p1.a - col_m1.a * col_p1.b;
            det_p1[jh] = col_m1.a * col_0.b - col_m1.b * col_0.a;
            
            // Fill diagonal blocks
            C[i](jh, jh) = col_m1.m1 * det_m1[jh] + col_0.m1 * det_0[jh] + 
                           col_p1.m1 * det_p1[jh];
            B[i](jh, jh) = col_m1.z * det_m1[jh] + col_0.z * det_0[jh] + 
                           col_p1.z * det_p1[jh];
            A[i](jh, jh) = col_m1.p1 * det_m1[jh] + col_0.p1 * det_0[jh] + 
                           col_p1.p1 * det_p1[jh];
            
            // RHS
            D(i, j) = d_rhs(i-1, j) * det_m1[jh] + d_rhs(i, j) * det_0[jh] + 
                      d_rhs(i+1, j) * det_p1[jh];
        }
    }
    
    // Wall boundary conditions
    core::Matrix<double> R0(n_heavy, n_heavy), R1(n_heavy, n_heavy), R2(n_heavy, n_heavy);
    std::vector<double> R(n_heavy);
    R0.setZero();
    R1.setZero();
    R2.setZero();
    
    for (std::size_t j = start_idx; j < n_species; ++j) {
        const int jh = j - start_idx;
        
        // BC collocation
        col_bc[jh].m1 = -1.5 * f_bc[j] + g_bc[j];
        col_bc[jh].z = 2.0 * f_bc[j];
        col_bc[jh].p1 = -0.5 * f_bc[j];
        col_bc[jh].a = -2.0 * f_bc[j];
        col_bc[jh].b = 2.0 * f_bc[j];
        
        // Equation collocations
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
        
        // Sub-determinants
        det_m1[jh] = col_1[jh].a * col_2[jh].b - col_1[jh].b * col_2[jh].a;
        det_0[jh] = col_bc[jh].b * col_2[jh].a - col_bc[jh].a * col_2[jh].b;
        det_p1[jh] = col_bc[jh].a * col_1[jh].b - col_bc[jh].b * col_1[jh].a;
        
        // Fill R matrices
        R0(jh, jh) = col_bc[jh].m1 * det_m1[jh] + col_1[jh].m1 * det_0[jh] + 
                     col_2[jh].m1 * det_p1[jh];
        R1(jh, jh) = col_bc[jh].z * det_m1[jh] + col_1[jh].z * det_0[jh] + 
                     col_2[jh].z * det_p1[jh];
        R2(jh, jh) = col_bc[jh].p1 * det_m1[jh] + col_1[jh].p1 * det_0[jh] + 
                     col_2[jh].p1 * det_p1[jh];
        R[jh] = h_bc[j] * det_m1[jh] + d_rhs(0, j) * det_0[jh] + 
                d_rhs(1, j) * det_p1[jh];
    }
    
    // Make block matrix truly tridiagonal
    try {
        // Handle top boundary condition
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
    
    // Handle bottom boundary
    for (std::size_t j = start_idx; j < n_species; ++j) {
        D(n_eta-2, j) -= A[n_eta-2](j-start_idx, j-start_idx) * prev_solution(j, n_eta-1);
    }
    A[n_eta-2].setZero();
    
    // Solve block tridiagonal system
    core::Matrix<double> solution(n_species, n_eta);
    if (auto result = detail::block_thomas_algorithm(C, B, A, D, solution, start_idx); !result) {
        return std::unexpected(result.error());
    }
    
    return solution;
}


/* auto solve_species_block_tridiagonal(
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
    
    std::cout << "=== DEBUG: Entering solve_species_block_tridiagonal ===" << std::endl;
    
    const auto n_eta = prev_solution.cols();
    const auto n_species = prev_solution.rows();
    
    if (n_eta < 3) {
        std::cout << "DEBUG: ERROR - Need at least 3 eta points" << std::endl;
        return std::unexpected(SolverError("Need at least 3 eta points"));
    }
    
    // Validate matrix dimensions
    if (a_coeffs.rows() != n_eta || a_coeffs.cols() != n_species ||
        b_coeffs.rows() != n_eta || b_coeffs.cols() != n_species ||
        c_coeffs.rows() != n_eta || c_coeffs.cols() != n_species ||
        d_rhs.rows() != n_eta || d_rhs.cols() != n_species) {
        std::cout << "DEBUG: ERROR - Coefficient matrix dimensions mismatch" << std::endl;
        return std::unexpected(SolverError("Coefficient matrix dimensions mismatch"));
    }
    
    // Validate boundary condition array sizes
    if (f_bc.size() != n_species || g_bc.size() != n_species || h_bc.size() != n_species) {
        std::cout << "DEBUG: ERROR - Boundary condition arrays size mismatch" << std::endl;
        return std::unexpected(SolverError("Boundary condition arrays size mismatch"));
    }
    
    const std::size_t start_idx = has_electrons ? 1 : 0;
    const std::size_t n_heavy = n_species - start_idx;
    
    std::cout << "DEBUG: start_idx=" << start_idx << ", n_heavy=" << n_heavy 
              << ", n_species=" << n_species << ", n_eta=" << n_eta << std::endl;
    
    if (n_heavy == 0) {
        std::cout << "DEBUG: ERROR - No heavy species to solve" << std::endl;
        return std::unexpected(SolverError("No heavy species to solve"));
    }
    
    std::cout << "DEBUG: Allocating block matrices..." << std::endl;
    
    // Allocate block matrices
    std::vector<core::Matrix<double>> A(n_eta), B(n_eta), C(n_eta);
    for (auto& mat : A) mat = core::Matrix<double>(n_heavy, n_heavy);
    for (auto& mat : B) mat = core::Matrix<double>(n_heavy, n_heavy);
    for (auto& mat : C) mat = core::Matrix<double>(n_heavy, n_heavy);
    
    core::Matrix<double> D(n_eta, n_species);
    
    // Working arrays for wall BC
    std::vector<CollocationCoeffs> col_bc(n_heavy), col_1(n_heavy), col_2(n_heavy);
    std::vector<double> det_m1(n_heavy), det_0(n_heavy), det_p1(n_heavy);
    
    std::cout << "DEBUG: Filling coefficients for interior points (i=1 to " << n_eta-2 << ")..." << std::endl;
    
    // Fill coefficients for interior points
    for (std::size_t i = 1; i < n_eta - 1; ++i) {
        
        A[i].setZero();
        B[i].setZero();
        C[i].setZero();
        
        // Process each species equation
        for (std::size_t j = start_idx; j < n_species; ++j) {

                    // Dans la boucle principale, avant compute_collocation_coeffs
        std::cout << "DEBUG: Input coeffs at i=" << i << ", j=" << j 
          << ": a[i-1]=" << a_coeffs(i-1, j) 
          << ", b[i-1]=" << b_coeffs(i-1, j) 
          << ", c[i-1]=" << c_coeffs(i-1, j) << std::endl;
        std::cout << "DEBUG: Processing interior point i=" << i << std::endl;
            const int jh = j - start_idx;  // Heavy species index
            
            std::cout << "DEBUG:   Species j=" << j << " (heavy index=" << jh << ")" << std::endl;
            
            // Collocation coefficients
            auto col_m1_result = compute_collocation_coeffs(
                a_coeffs(i-1, j), b_coeffs(i-1, j), c_coeffs(i-1, j), -1.0
            );
            if (!col_m1_result) {
                std::cout << "DEBUG: ERROR - Failed to compute collocation coeffs for col_m1" << std::endl;
                return std::unexpected(col_m1_result.error());
            }
            auto col_m1 = col_m1_result.value();
            
            auto col_0_result = compute_collocation_coeffs(
                a_coeffs(i, j), b_coeffs(i, j), c_coeffs(i, j), 0.0
            );
            if (!col_0_result) {
                std::cout << "DEBUG: ERROR - Failed to compute collocation coeffs for col_0" << std::endl;
                return std::unexpected(col_0_result.error());
            }
            auto col_0 = col_0_result.value();
            
            auto col_p1_result = compute_collocation_coeffs(
                a_coeffs(i+1, j), b_coeffs(i+1, j), c_coeffs(i+1, j), 1.0
            );
            if (!col_p1_result) {
                std::cout << "DEBUG: ERROR - Failed to compute collocation coeffs for col_p1" << std::endl;
                return std::unexpected(col_p1_result.error());
            }
            auto col_p1 = col_p1_result.value();
            
            // Sub-determinants
            det_m1[jh] = col_0.a * col_p1.b - col_0.b * col_p1.a;
            det_0[jh] = col_m1.b * col_p1.a - col_m1.a * col_p1.b;
            det_p1[jh] = col_m1.a * col_0.b - col_m1.b * col_0.a;
            
            std::cout << "DEBUG:     Determinants: det_m1=" << std::scientific << std::setprecision(6) 
                      << det_m1[jh] << ", det_0=" << det_0[jh] << ", det_p1=" << det_p1[jh] << std::endl;
            
            // Fill diagonal blocks
            C[i](jh, jh) = col_m1.m1 * det_m1[jh] + col_0.m1 * det_0[jh] + 
                           col_p1.m1 * det_p1[jh];
            B[i](jh, jh) = col_m1.z * det_m1[jh] + col_0.z * det_0[jh] + 
                           col_p1.z * det_p1[jh];
            A[i](jh, jh) = col_m1.p1 * det_m1[jh] + col_0.p1 * det_0[jh] + 
                           col_p1.p1 * det_p1[jh];
            
            // RHS
            D(i, j) = d_rhs(i-1, j) * det_m1[jh] + d_rhs(i, j) * det_0[jh] + 
                      d_rhs(i+1, j) * det_p1[jh];
            
            std::cout << "DEBUG:     Diagonal elements: A=" << A[i](jh, jh) 
                      << ", B=" << B[i](jh, jh) << ", C=" << C[i](jh, jh) 
                      << ", D=" << D(i, j) << std::endl;
        }
        
        std::cout << "DEBUG:   Adding df coupling to A,B,C..." << std::endl;
        // NOTE: df coupling would be added here if df matrices were available
        // This is where the C++ version differs - it doesn't have df coupling in this snippet
        
        std::cout << "DEBUG:   Adding df coupling to D..." << std::endl;
        // NOTE: df coupling to D would be added here if df matrices were available
    }
    
    std::cout << "DEBUG: Computing wall boundary conditions..." << std::endl;
    
    // Wall boundary conditions
    core::Matrix<double> R0(n_heavy, n_heavy), R1(n_heavy, n_heavy), R2(n_heavy, n_heavy);
    std::vector<double> R(n_heavy);
    R0.setZero();
    R1.setZero();
    R2.setZero();
    
    for (std::size_t j = start_idx; j < n_species; ++j) {
        const int jh = j - start_idx;
        
        std::cout << "DEBUG:   Wall BC for species j=" << j << std::endl;
        
        // BC collocation
        col_bc[jh].m1 = -1.5 * f_bc[j] + g_bc[j];
        col_bc[jh].z = 2.0 * f_bc[j];
        col_bc[jh].p1 = -0.5 * f_bc[j];
        col_bc[jh].a = -2.0 * f_bc[j];
        col_bc[jh].b = 2.0 * f_bc[j];
        
        // Equation collocations
        auto col_1_result = compute_collocation_coeffs(
            a_coeffs(0, j), b_coeffs(0, j), c_coeffs(0, j), -1.0
        );
        if (!col_1_result) {
            std::cout << "DEBUG: ERROR - Failed to compute wall collocation coeffs for col_1" << std::endl;
            return std::unexpected(col_1_result.error());
        }
        col_1[jh] = col_1_result.value();
        
        auto col_2_result = compute_collocation_coeffs(
            a_coeffs(1, j), b_coeffs(1, j), c_coeffs(1, j), 0.0
        );
        if (!col_2_result) {
            std::cout << "DEBUG: ERROR - Failed to compute wall collocation coeffs for col_2" << std::endl;
            return std::unexpected(col_2_result.error());
        }
        col_2[jh] = col_2_result.value();
        
        // Sub-determinants
        det_m1[jh] = col_1[jh].a * col_2[jh].b - col_1[jh].b * col_2[jh].a;
        det_0[jh] = col_bc[jh].b * col_2[jh].a - col_bc[jh].a * col_2[jh].b;
        det_p1[jh] = col_bc[jh].a * col_1[jh].b - col_bc[jh].b * col_1[jh].a;
        
        std::cout << "DEBUG:     Wall determinants: det_m1=" << det_m1[jh] 
                  << ", det_0=" << det_0[jh] << ", det_p1=" << det_p1[jh] << std::endl;
        
        // Fill R matrices
        R0(jh, jh) = col_bc[jh].m1 * det_m1[jh] + col_1[jh].m1 * det_0[jh] + 
                     col_2[jh].m1 * det_p1[jh];
        R1(jh, jh) = col_bc[jh].z * det_m1[jh] + col_1[jh].z * det_0[jh] + 
                     col_2[jh].z * det_p1[jh];
        R2(jh, jh) = col_bc[jh].p1 * det_m1[jh] + col_1[jh].p1 * det_0[jh] + 
                     col_2[jh].p1 * det_p1[jh];
        R[jh] = h_bc[j] * det_m1[jh] + d_rhs(0, j) * det_0[jh] + 
                d_rhs(1, j) * det_p1[jh];
        
        std::cout << "DEBUG:     R matrices diagonal: R0=" << R0(jh, jh) 
                  << ", R1=" << R1(jh, jh) << ", R2=" << R2(jh, jh) 
                  << ", R=" << R[jh] << std::endl;
    }
    
    std::cout << "DEBUG: Adding df coupling to R0, R1..." << std::endl;
    // NOTE: df coupling would be added here if df matrices were available
    
    std::cout << "DEBUG: Adding dh_bc coupling to R0..." << std::endl;
    // NOTE: dh_bc coupling would be added here if dh_bc_dc matrix was available
    
    std::cout << "DEBUG: Adding df coupling to R..." << std::endl;
    // NOTE: df coupling to R would be added here if df matrices were available
    
    std::cout << "DEBUG: Adding dh_bc coupling to R..." << std::endl;
    // NOTE: dh_bc coupling to R would be added here if dh_bc_dc matrix was available
    
    std::cout << "DEBUG: Making matrix truly tridiagonal..." << std::endl;
    
    // Make block matrix truly tridiagonal
    try {
        // Handle top boundary condition
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_A1(A[1].eigen());
        if (lu_A1.info() != Eigen::Success) {
            std::cout << "DEBUG: ERROR - Singular matrix A[1] in boundary condition setup" << std::endl;
            return std::unexpected(SolverError("Singular matrix A[1] in boundary condition setup"));
        }
        
        std::cout << "DEBUG: Computing LU decomposition and matrix operations..." << std::endl;
        
        auto temp_solve_C = lu_A1.solve(C[1].eigen());
        auto temp_solve_B = lu_A1.solve(B[1].eigen());

        B[0].eigen() = (R0.eigen() - R2.eigen() * temp_solve_C).eval();
        A[0].eigen() = (R1.eigen() - R2.eigen() * temp_solve_B).eval();
 
        Eigen::VectorXd D1(n_heavy);
        for (std::size_t j = 0; j < n_heavy; ++j) {
            D1[j] = D(1, j + start_idx);
        }
        Eigen::VectorXd temp = R2.eigen() * lu_A1.solve(D1);
        
        std::cout << "DEBUG: Finalizing boundary conditions..." << std::endl;
        
        for (std::size_t j = 0; j < n_heavy; ++j) {
            D(0, j + start_idx) = R[j] - temp[j];
        }
    } catch (const std::exception& e) {
        std::cout << "DEBUG: ERROR - Matrix operation failed in boundary setup: " << e.what() << std::endl;
        return std::unexpected(SolverError(
            std::format("Matrix operation failed in boundary setup: {}", e.what())
        ));
    }
    
    std::cout << "DEBUG: Handling bottom boundary..." << std::endl;
    
    // Handle bottom boundary
    for (std::size_t j = start_idx; j < n_species; ++j) {
        double correction = A[n_eta-2](j-start_idx, j-start_idx) * prev_solution(j, n_eta-1);
        D(n_eta-2, j) -= correction;
        std::cout << "DEBUG:   Bottom correction for species " << j << ": " << correction << std::endl;
    }
    A[n_eta-2].setZero();
    
    std::cout << "DEBUG: Solving block tridiagonal system..." << std::endl;
    
    // Solve block tridiagonal system
    core::Matrix<double> solution(n_species, n_eta);
    if (auto result = detail::block_thomas_algorithm(C, B, A, D, solution, start_idx); !result) {
        std::cout << "DEBUG: ERROR - Block Thomas algorithm failed" << std::endl;
        return std::unexpected(result.error());
    }
    
    std::cout << "=== DEBUG: Exiting solve_species_block_tridiagonal ===" << std::endl;
    
    return solution;
} */

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

/* auto block_thomas_algorithm(
    std::vector<core::Matrix<double>>& lower_blocks,
    std::vector<core::Matrix<double>>& main_blocks,
    std::vector<core::Matrix<double>>& upper_blocks,
    core::Matrix<double>& rhs,
    core::Matrix<double>& solution,
    int start_species
) -> std::expected<void, SolverError> {
    const auto n_eta = main_blocks.size();
    
    if (n_eta < 2) {
        return std::unexpected(SolverError("System too small for block Thomas algorithm"));
    }
    
    if (lower_blocks.size() != n_eta || upper_blocks.size() != n_eta) {
        return std::unexpected(SolverError("Incompatible block matrix sizes"));
    }
    
    const auto n_heavy = main_blocks[0].rows();
    if (n_heavy == 0) {
        return std::unexpected(SolverError("Empty block matrices"));
    }
    
    // Validation des tailles de matrices
    for (std::size_t i = 0; i < n_eta; ++i) {
        if (main_blocks[i].rows() != n_heavy || main_blocks[i].cols() != n_heavy ||
            lower_blocks[i].rows() != n_heavy || lower_blocks[i].cols() != n_heavy ||
            upper_blocks[i].rows() != n_heavy || upper_blocks[i].cols() != n_heavy) {
            return std::unexpected(SolverError(
                std::format("Inconsistent block matrix dimensions at position {}", i)
            ));
        }
    }
    
    // Forward elimination
    for (std::size_t i = 1; i < n_eta - 1; ++i) {
        try {
            Eigen::PartialPivLU<Eigen::MatrixXd> lu_B(main_blocks[i-1].eigen());
            
            if (lu_B.info() != Eigen::Success) {
                return std::unexpected(SolverError(
                    std::format("Singular matrix at block position {}", i-1)
                ));
            }
            
            upper_blocks[i-1].eigen() = (lu_B.solve(upper_blocks[i-1].eigen())).eval();
            main_blocks[i].eigen() = (main_blocks[i].eigen() - 
                                     lower_blocks[i].eigen() * upper_blocks[i-1].eigen()).eval();
            
            // Update RHS
            for (std::size_t j = 0; j < n_heavy; ++j) {
                double temp = 0.0;
                for (std::size_t k = 0; k < n_heavy; ++k) {
                    temp += lower_blocks[i](j, k) * rhs(i-1, k + start_species);
                }
                rhs(i, j + start_species) -= temp;
            }
        } catch (const std::exception& e) {
            return std::unexpected(SolverError(
                std::format("Matrix operation failed at block {}: {}", i, e.what())
            ));
        }
    }
    
    // Back substitution
    try {
        // Last interior point
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_B(main_blocks[n_eta - 2].eigen());
        
        if (lu_B.info() != Eigen::Success) {
            return std::unexpected(SolverError(
                std::format("Singular matrix at final block position {}", n_eta - 2)
            ));
        }
        
        Eigen::VectorXd rhs_vec(n_heavy);
        for (std::size_t k = 0; k < n_heavy; ++k) {
            rhs_vec[k] = rhs(n_eta - 2, k + start_species);
        }
        
        Eigen::VectorXd sol = lu_B.solve(rhs_vec);
        for (std::size_t j = 0; j < n_heavy; ++j) {
            solution(j + start_species, n_eta - 2) = sol[j];
        }
        
        // Remaining points
        for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n_eta) - 3; i >= 0; --i) {
            for (std::size_t j = 0; j < n_heavy; ++j) {
                solution(j + start_species, i) = rhs(i, j + start_species);
                for (std::size_t k = 0; k < n_heavy; ++k) {
                    solution(j + start_species, i) -= 
                        upper_blocks[i](j, k) * solution(k + start_species, i+1);
                }
            }
        }
    } catch (const std::exception& e) {
        return std::unexpected(SolverError(
            std::format("Back substitution failed: {}", e.what())
        ));
    }
    
    return {};
} */


/* auto block_thomas_algorithm(
    std::vector<core::Matrix<double>>& lower_blocks,
    std::vector<core::Matrix<double>>& main_blocks,
    std::vector<core::Matrix<double>>& upper_blocks,
    core::Matrix<double>& rhs,
    core::Matrix<double>& solution,
    int start_species
) -> std::expected<void, SolverError> {
    
    std::cout << "=== DEBUG BLOCK_THOMAS: Starting block Thomas algorithm ===" << std::endl;
    
    const auto n_eta = main_blocks.size();
    const auto n_heavy = main_blocks[0].rows();
    
    std::cout << "DEBUG: n_eta=" << n_eta << ", n_heavy=" << n_heavy 
              << ", start_species=" << start_species << std::endl;
    
    if (n_eta < 2) {
        std::cout << "DEBUG: ERROR - System too small" << std::endl;
        return std::unexpected(SolverError("System too small for block Thomas algorithm"));
    }
    
    if (lower_blocks.size() != n_eta || upper_blocks.size() != n_eta) {
        std::cout << "DEBUG: ERROR - Incompatible block matrix sizes" << std::endl;
        return std::unexpected(SolverError("Incompatible block matrix sizes"));
    }
    
    if (n_heavy == 0) {
        std::cout << "DEBUG: ERROR - Empty block matrices" << std::endl;
        return std::unexpected(SolverError("Empty block matrices"));
    }
    
    // Validation des tailles de matrices
    for (std::size_t i = 0; i < n_eta; ++i) {
        if (main_blocks[i].rows() != n_heavy || main_blocks[i].cols() != n_heavy ||
            lower_blocks[i].rows() != n_heavy || lower_blocks[i].cols() != n_heavy ||
            upper_blocks[i].rows() != n_heavy || upper_blocks[i].cols() != n_heavy) {
            std::cout << "DEBUG: ERROR - Inconsistent block matrix dimensions at position " << i << std::endl;
            return std::unexpected(SolverError(
                std::format("Inconsistent block matrix dimensions at position {}", i)
            ));
        }
    }
    
    std::cout << "DEBUG: Initial conditions at i=0" << std::endl;
    std::cout << "DEBUG: B[0] matrix (main_blocks[0]):" << std::endl;
    for (std::size_t j = 0; j < n_heavy; ++j) {
        for (std::size_t k = 0; k < n_heavy; ++k) {
            std::cout << "  B[0][" << j << "][" << k << "] = " 
                      << std::scientific << std::setprecision(12) 
                      << main_blocks[0](j, k) << std::endl;
        }
    }
    std::cout << "DEBUG: C[0] matrix (upper_blocks[0]):" << std::endl;
    for (std::size_t j = 0; j < n_heavy; ++j) {
        for (std::size_t k = 0; k < n_heavy; ++k) {
            std::cout << "  C[0][" << j << "][" << k << "] = " 
                      << std::scientific << std::setprecision(12) 
                      << upper_blocks[0](j, k) << std::endl;
        }
    }
    std::cout << "DEBUG: RHS[0] before processing:" << std::endl;
    for (std::size_t j = 0; j < n_heavy; ++j) {
        std::cout << "  RHS[0][" << j << "] = " 
                  << std::scientific << std::setprecision(12) 
                  << rhs(0, j + start_species) << std::endl;
    }
    
    // *** Traitement initial à i=0 : calculer beta[0] et cy[0] ***
    std::cout << "DEBUG: Computing initial beta[0] and cy[0] at i=0" << std::endl;
    try {
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_B0(main_blocks[0].eigen());
        
        if (lu_B0.info() != Eigen::Success) {
            std::cout << "DEBUG: ERROR - Singular matrix B[0]" << std::endl;
            return std::unexpected(SolverError("Singular matrix B[0] in initial conditions"));
        }
        
        // Calculer beta[0] = inv(B[0]) * C[0]
        upper_blocks[0].eigen() = (lu_B0.solve(upper_blocks[0].eigen())).eval();
        std::cout << "DEBUG: upper_blocks[0] after solve (beta[0]):" << std::endl;
        for (std::size_t j = 0; j < n_heavy; ++j) {
            for (std::size_t k = 0; k < n_heavy; ++k) {
                std::cout << "  upper_blocks[0][" << j << "][" << k << "] = " 
                          << std::scientific << std::setprecision(12) 
                          << upper_blocks[0](j, k) << std::endl;
            }
        }
        
        // Calculer cy[0] = inv(B[0]) * RHS[0]
        Eigen::VectorXd rhs_0(n_heavy);
        for (std::size_t j = 0; j < n_heavy; ++j) {
            rhs_0[j] = rhs(0, j + start_species);
        }
        
        Eigen::VectorXd cy_0 = lu_B0.solve(rhs_0);
        std::cout << "DEBUG: cy[0] (= inv(B[0]) * RHS[0]):" << std::endl;
        for (std::size_t j = 0; j < n_heavy; ++j) {
            solution(j + start_species, 0) = cy_0[j];
            std::cout << "  cy[0][" << j << "] = " 
                      << std::scientific << std::setprecision(12) 
                      << cy_0[j] << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cout << "DEBUG: ERROR - Failed to solve initial conditions: " << e.what() << std::endl;
        return std::unexpected(SolverError("Failed to solve boundary conditions at i=0"));
    }
    
    // Forward elimination
    std::cout << "\n=== DEBUG: Forward elimination phase ===" << std::endl;
    for (std::size_t i = 1; i < n_eta; ++i) {
        // Only debug first few iterations to avoid terminal overflow
        bool debug_this_step = (i <= 2);
        
        if (debug_this_step) {
            std::cout << "\n=== DEBUG: Forward elimination step i=" << i << " ===" << std::endl;
        }
        
        if (debug_this_step) {
            std::cout << "DEBUG: A[" << i << "] matrix (lower_blocks[" << i << "]):" << std::endl;
            for (std::size_t j = 0; j < n_heavy; ++j) {
                for (std::size_t k = 0; k < n_heavy; ++k) {
                    std::cout << "  A[" << i << "][" << j << "][" << k << "] = " 
                              << std::scientific << std::setprecision(12) 
                              << lower_blocks[i](j, k) << std::endl;
                }
            }
            
            std::cout << "DEBUG: B[" << i << "] matrix before update (main_blocks[" << i << "]):" << std::endl;
            for (std::size_t j = 0; j < n_heavy; ++j) {
                for (std::size_t k = 0; k < n_heavy; ++k) {
                    std::cout << "  B[" << i << "][" << j << "][" << k << "] = " 
                              << std::scientific << std::setprecision(12) 
                              << main_blocks[i](j, k) << std::endl;
                }
            }
        }
        
        try {
            // Store original values for comparison
            auto original_main = main_blocks[i];
            
            if (debug_this_step) {
                std::cout << "DEBUG: Computing main_blocks[" << i << "] = B[" << i << "] - A[" << i << "] * upper_blocks[" << i-1 << "]" << std::endl;
            }
            auto product = lower_blocks[i].eigen() * upper_blocks[i-1].eigen();
            main_blocks[i].eigen() = (main_blocks[i].eigen() - product).eval();
            
            if (debug_this_step) {
                std::cout << "DEBUG: Product A[" << i << "] * upper_blocks[" << i-1 << "]:" << std::endl;
                for (std::size_t j = 0; j < n_heavy; ++j) {
                    for (std::size_t k = 0; k < n_heavy; ++k) {
                        std::cout << "  product[" << j << "][" << k << "] = " 
                                  << std::scientific << std::setprecision(12) 
                                  << product(j, k) << std::endl;
                    }
                }
                
                std::cout << "DEBUG: main_blocks[" << i << "] after update:" << std::endl;
                for (std::size_t j = 0; j < n_heavy; ++j) {
                    for (std::size_t k = 0; k < n_heavy; ++k) {
                        std::cout << "  B[" << i << "][" << j << "][" << k << "]: " 
                                  << std::scientific << std::setprecision(12) 
                                  << original_main(j, k) << " - " << product(j, k) 
                                  << " = " << main_blocks[i](j, k) << std::endl;
                    }
                }
            }
            
            // Update RHS
            if (debug_this_step) {
                std::cout << "DEBUG: RHS[" << i << "] before update:" << std::endl;
                for (std::size_t j = 0; j < n_heavy; ++j) {
                    std::cout << "  RHS[" << i << "][" << j << "] = " 
                              << std::scientific << std::setprecision(12) 
                              << rhs(i, j + start_species) << std::endl;
                }
            }
            
            for (std::size_t j = 0; j < n_heavy; ++j) {
                double temp = 0.0;
                for (std::size_t k = 0; k < n_heavy; ++k) {
                    temp += lower_blocks[i](j, k) * solution(k + start_species, i-1);
                }
                double old_val = rhs(i, j + start_species);
                rhs(i, j + start_species) -= temp;
                if (debug_this_step) {
                    std::cout << "  RHS[" << i << "][" << j << "]: " 
                              << std::scientific << std::setprecision(12) 
                              << old_val << " - " << temp << " = " << rhs(i, j + start_species) << std::endl;
                }
            }
            
            // Pour les points intérieurs (pas le dernier), calculer beta[i] et cy[i]
            if (i < n_eta - 1) {
                if (debug_this_step) {
                    std::cout << "DEBUG: Computing inv_prod for i=" << i << " (interior point)" << std::endl;
                }
                
                Eigen::PartialPivLU<Eigen::MatrixXd> lu_B(main_blocks[i].eigen());
                
                if (lu_B.info() != Eigen::Success) {
                    std::cout << "DEBUG: ERROR - Singular matrix at block position " << i << std::endl;
                    return std::unexpected(SolverError(
                        std::format("Singular matrix at block position {}", i)
                    ));
                }
                
                // Calculer beta[i] = inv(B[i]) * C[i]
                upper_blocks[i].eigen() = (lu_B.solve(upper_blocks[i].eigen())).eval();
                
                if (debug_this_step) {
                    std::cout << "DEBUG: upper_blocks[" << i << "] after solve (beta[" << i << "]):" << std::endl;
                    for (std::size_t j = 0; j < n_heavy; ++j) {
                        for (std::size_t k = 0; k < n_heavy; ++k) {
                            std::cout << "  upper_blocks[" << i << "][" << j << "][" << k << "] = " 
                                      << std::scientific << std::setprecision(12) 
                                      << upper_blocks[i](j, k) << std::endl;
                        }
                    }
                }
                
                // Calculer cy[i] = inv(B[i]) * RHS[i]
                Eigen::VectorXd rhs_vec(n_heavy);
                for (std::size_t k = 0; k < n_heavy; ++k) {
                    rhs_vec[k] = rhs(i, k + start_species);
                }
                
                Eigen::VectorXd cy = lu_B.solve(rhs_vec);
                for (std::size_t j = 0; j < n_heavy; ++j) {
                    solution(j + start_species, i) = cy[j];
                }
                
                if (debug_this_step) {
                    std::cout << "DEBUG: cy[" << i << "] after solve:" << std::endl;
                    for (std::size_t j = 0; j < n_heavy; ++j) {
                        std::cout << "  cy[" << i << "][" << j << "] = " 
                                  << std::scientific << std::setprecision(12) 
                                  << cy[j] << std::endl;
                    }
                }
            } else {
                // Dernier point : seulement calculer cy[i] = inv(B[i]) * RHS[i]
                std::cout << "DEBUG: Computing inv_vec for i=" << i << " (last point)" << std::endl;
                
                Eigen::PartialPivLU<Eigen::MatrixXd> lu_B(main_blocks[i].eigen());
                
                if (lu_B.info() != Eigen::Success) {
                    std::cout << "DEBUG: ERROR - Singular matrix at final block position " << i << std::endl;
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
                
                std::cout << "DEBUG: cy[" << i << "] after inv_vec (final solution at last point):" << std::endl;
                for (std::size_t j = 0; j < n_heavy; ++j) {
                    std::cout << "  cy[" << i << "][" << j << "] = " 
                              << std::scientific << std::setprecision(12) 
                              << cy[j] << std::endl;
                }
            }
            
        } catch (const std::exception& e) {
            std::cout << "DEBUG: ERROR - Matrix operation failed at block " << i << ": " << e.what() << std::endl;
            return std::unexpected(SolverError(
                std::format("Matrix operation failed at block {}: {}", i, e.what())
            ));
        }
    }
    
    // Back substitution
    std::cout << "\n=== DEBUG: Back substitution phase ===" << std::endl;
    try {
        std::cout << "DEBUG: Final solution at i=" << (n_eta-1) << ":" << std::endl;
        for (std::size_t j = 0; j < n_heavy; ++j) {
            std::cout << "  solution[" << j + start_species << "][" << (n_eta-1) << "] = " 
                      << std::scientific << std::setprecision(12) 
                      << solution(j + start_species, n_eta-1) << std::endl;
        }
        
        // Back substitution pour les points restants
        std::cout << "DEBUG: Back substitution for remaining points" << std::endl;
        for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n_eta) - 2; i >= 0; --i) {
            std::cout << "DEBUG: Back substitution step i=" << i << std::endl;
            for (std::size_t j = 0; j < n_heavy; ++j) {
                double correction = 0.0;
                for (std::size_t k = 0; k < n_heavy; ++k) {
                    correction += upper_blocks[i](j, k) * solution(k + start_species, i+1);
                }
                double cy_val = solution(j + start_species, i);
                solution(j + start_species, i) = cy_val - correction;
                std::cout << "  solution[" << j + start_species << "][" << i << "] = " 
                          << std::scientific << std::setprecision(12) 
                          << cy_val << " - " << correction << " = " 
                          << solution(j + start_species, i) << std::endl;
            }
        }
        
    } catch (const std::exception& e) {
        std::cout << "DEBUG: ERROR - Back substitution failed: " << e.what() << std::endl;
        return std::unexpected(SolverError(
            std::format("Back substitution failed: {}", e.what())
        ));
    }
    
    std::cout << "\n=== DEBUG: Final solution vector ===" << std::endl;
    for (std::size_t i = 0; i < n_eta; ++i) {
        for (std::size_t j = 0; j < n_heavy; ++j) {
            std::cout << "SOLUTION[" << i << "][" << j + start_species << "] = " 
                      << std::scientific << std::setprecision(12) 
                      << solution(j + start_species, i) << std::endl;
        }
    }
    
    std::cout << "=== DEBUG BLOCK_THOMAS: Completed ===" << std::endl;
    return {};
} */

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
    
    // Validation des tailles de matrices
    for (std::size_t i = 0; i < n_eta; ++i) {
        if (main_blocks[i].rows() != n_heavy || main_blocks[i].cols() != n_heavy ||
            lower_blocks[i].rows() != n_heavy || lower_blocks[i].cols() != n_heavy ||
            upper_blocks[i].rows() != n_heavy || upper_blocks[i].cols() != n_heavy) {
            return std::unexpected(SolverError(
                std::format("Inconsistent block matrix dimensions at position {}", i)
            ));
        }
    }
    
    // *** Traitement initial à i=0 : calculer beta[0] et cy[0] ***
    try {
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_B0(main_blocks[0].eigen());
        
        if (lu_B0.info() != Eigen::Success) {
            return std::unexpected(SolverError("Singular matrix B[0] in initial conditions"));
        }
        
        // Calculer beta[0] = inv(B[0]) * C[0]
        upper_blocks[0].eigen() = (lu_B0.solve(upper_blocks[0].eigen())).eval();
        
        // Calculer cy[0] = inv(B[0]) * RHS[0]
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
    
    // Forward elimination
    for (std::size_t i = 1; i < n_eta; ++i) {
        try {
            // Update main diagonal block
            auto product = lower_blocks[i].eigen() * upper_blocks[i-1].eigen();
            main_blocks[i].eigen() = (main_blocks[i].eigen() - product).eval();
            
            // Update RHS
            for (std::size_t j = 0; j < n_heavy; ++j) {
                double temp = 0.0;
                for (std::size_t k = 0; k < n_heavy; ++k) {
                    temp += lower_blocks[i](j, k) * solution(k + start_species, i-1);
                }
                rhs(i, j + start_species) -= temp;
            }
            
            // Pour les points intérieurs (pas le dernier), calculer beta[i] et cy[i]
            if (i < n_eta - 1) {
                Eigen::PartialPivLU<Eigen::MatrixXd> lu_B(main_blocks[i].eigen());
                
                if (lu_B.info() != Eigen::Success) {
                    return std::unexpected(SolverError(
                        std::format("Singular matrix at block position {}", i)
                    ));
                }
                
                // Calculer beta[i] = inv(B[i]) * C[i]
                upper_blocks[i].eigen() = (lu_B.solve(upper_blocks[i].eigen())).eval();
                
                // Calculer cy[i] = inv(B[i]) * RHS[i]
                Eigen::VectorXd rhs_vec(n_heavy);
                for (std::size_t k = 0; k < n_heavy; ++k) {
                    rhs_vec[k] = rhs(i, k + start_species);
                }
                
                Eigen::VectorXd cy = lu_B.solve(rhs_vec);
                for (std::size_t j = 0; j < n_heavy; ++j) {
                    solution(j + start_species, i) = cy[j];
                }
            } else {
                // Dernier point : seulement calculer cy[i] = inv(B[i]) * RHS[i]
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
    
    // Back substitution
    try {
        // Back substitution pour les points restants
        for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n_eta) - 2; i >= 0; --i) {
            for (std::size_t j = 0; j < n_heavy; ++j) {
                double correction = 0.0;
                for (std::size_t k = 0; k < n_heavy; ++k) {
                    correction += upper_blocks[i](j, k) * solution(k + start_species, i+1);
                }
                double cy_val = solution(j + start_species, i);
                solution(j + start_species, i) = cy_val - correction;
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