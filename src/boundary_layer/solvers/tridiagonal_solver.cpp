#include "blast/boundary_layer/solvers/tridiagonal_solvers.hpp"
#include <algorithm>
#include <numeric>

namespace blast::boundary_layer::solvers {

namespace {
    // Collocation coefficients structure
    struct CollocationCoeffs {
        double p1, z, m1;  // Values at collocation points
        double a, b;       // Alpha and beta coefficients
    };
    
    // Compute collocation coefficients for a given eta point
    // This function now exactly matches the original C code formulas
    [[nodiscard]] auto compute_collocation_coeffs(
        double a, double b, double c, double t
    ) noexcept -> CollocationCoeffs {
        CollocationCoeffs coeffs;
        
        if (t == -1.0) {
            // Formulas for t = -1 (from original C code)
            coeffs.p1 = a - 0.5 * b;
            coeffs.z = -2.0 * a + 2.0 * b;
            coeffs.m1 = a - 1.5 * b + c;
            coeffs.a = 6.0 * a - 2.0 * b;
            coeffs.b = -10.0 * a + 2.0 * b;
        } else if (t == 0.0) {
            // Formulas for t = 0 (from original C code)
            coeffs.p1 = a + 0.5 * b;
            coeffs.z = -2.0 * a + c;
            coeffs.m1 = a - 0.5 * b;
            coeffs.a = b;
            coeffs.b = 2.0 * a;
        } else if (t == 1.0) {
            // Formulas for t = 1 (from original C code)
            coeffs.p1 = a + 1.5 * b + c;
            coeffs.z = -2.0 * a - 2.0 * b;
            coeffs.m1 = a + 0.5 * b;
            coeffs.a = -6.0 * a - 2.0 * b;
            coeffs.b = -10.0 * a - 2.0 * b;
        } else {
            // This should not happen in the tridiagonal solver context
            // where only t = -1, 0, 1 are used
            std::terminate(); // Or handle error appropriately
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
    
    // Allocate working arrays
    std::vector<double> A(n_eta), B(n_eta), C(n_eta), D(n_eta);
    
    // Fill coefficients for interior points
    for (std::size_t i = 1; i < n_eta - 1; ++i) {
        // Collocation at t = -1, 0, 1
        auto col_m1 = compute_collocation_coeffs(a_coeffs[i-1], b_coeffs[i-1], c_coeffs[i-1], -1.0);
        auto col_0 = compute_collocation_coeffs(a_coeffs[i], b_coeffs[i], c_coeffs[i], 0.0);
        auto col_p1 = compute_collocation_coeffs(a_coeffs[i+1], b_coeffs[i+1], c_coeffs[i+1], 1.0);
        
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
    auto col_1 = compute_collocation_coeffs(a_coeffs[0], b_coeffs[0], c_coeffs[0], -1.0);
    auto col_2 = compute_collocation_coeffs(a_coeffs[1], b_coeffs[1], c_coeffs[1], 0.0);
    
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
    detail::thomas_algorithm(C, B, A, D, solution);
    
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
    
    const int start_idx = has_electrons ? 1 : 0;
    const int n_heavy = n_species - start_idx;
    
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
        for (int j = start_idx; j < n_species; ++j) {
            const int jh = j - start_idx;  // Heavy species index
            
            // Collocation coefficients
            auto col_m1 = compute_collocation_coeffs(
                a_coeffs(i-1, j), b_coeffs(i-1, j), c_coeffs(i-1, j), -1.0
            );
            auto col_0 = compute_collocation_coeffs(
                a_coeffs(i, j), b_coeffs(i, j), c_coeffs(i, j), 0.0
            );
            auto col_p1 = compute_collocation_coeffs(
                a_coeffs(i+1, j), b_coeffs(i+1, j), c_coeffs(i+1, j), 1.0
            );
            
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
    
    for (int j = start_idx; j < n_species; ++j) {
        const int jh = j - start_idx;
        
        // BC collocation
        col_bc[jh].m1 = -1.5 * f_bc[j] + g_bc[j];
        col_bc[jh].z = 2.0 * f_bc[j];
        col_bc[jh].p1 = -0.5 * f_bc[j];
        col_bc[jh].a = -2.0 * f_bc[j];
        col_bc[jh].b = 2.0 * f_bc[j];
        
        // Equation collocations
        col_1[jh] = compute_collocation_coeffs(
            a_coeffs(0, j), b_coeffs(0, j), c_coeffs(0, j), -1.0
        );
        col_2[jh] = compute_collocation_coeffs(
            a_coeffs(1, j), b_coeffs(1, j), c_coeffs(1, j), 0.0
        );
        
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
    // Handle top boundary condition

    Eigen::PartialPivLU<Eigen::MatrixXd> lu_A1(A[1].eigen());
    B[0] = (R0 - R2 * lu_A1.solve(C[1].eigen())).eval();
    A[0] = (R1 - R2 * lu_A1.solve(B[1].eigen())).eval();
        
    Eigen::VectorXd D1(n_heavy);
    for (int j = 0; j < n_heavy; ++j) {
        D1[j] = D(1, j + start_idx);
    }
    Eigen::VectorXd temp = R2.eigen() * lu_A1.solve(D1);
        
    for (int j = 0; j < n_heavy; ++j) {
        D(0, j + start_idx) = R[j] - temp[j];
    }
    
    // Handle bottom boundary
    for (int j = start_idx; j < n_species; ++j) {
        D(n_eta-2, j) -= A[n_eta-2](j-start_idx, j-start_idx) * prev_solution(j, n_eta-1);
    }
    A[n_eta-2].setZero();
    
    // Solve block tridiagonal system
    core::Matrix<double> solution(n_species, n_eta);
    detail::block_thomas_algorithm(C, B, A, D, solution, start_idx);
    
    return solution;
}

namespace detail {

void thomas_algorithm(
    std::span<double> lower_diag,
    std::span<double> main_diag,
    std::span<double> upper_diag,
    std::span<double> rhs,
    std::span<double> solution
) {
    const auto n = main_diag.size();
    
    // Forward sweep
    main_diag[0] = 1.0 / main_diag[0];
    lower_diag[0] = rhs[0] * main_diag[0];
    
    for (std::size_t i = 1; i < n - 1; ++i) {
        upper_diag[i-1] = upper_diag[i-1] * main_diag[i-1];
        main_diag[i] = main_diag[i] - lower_diag[i] * upper_diag[i-1];
        main_diag[i] = 1.0 / main_diag[i];
        lower_diag[i] = (rhs[i] - lower_diag[i] * lower_diag[i-1]) * main_diag[i];
    }
    
    // Backward sweep
    solution[n-2] = lower_diag[n-2];
    for (int i = n - 3; i >= 0; --i) {
        solution[i] = lower_diag[i] - upper_diag[i] * solution[i+1];
    }
}

void block_thomas_algorithm(
    std::vector<core::Matrix<double>>& lower_blocks,
    std::vector<core::Matrix<double>>& main_blocks,
    std::vector<core::Matrix<double>>& upper_blocks,
    core::Matrix<double>& rhs,
    core::Matrix<double>& solution,
    int start_species
) {
    const auto n_eta = main_blocks.size();
    const auto n_heavy = main_blocks[0].rows();
    
    // Forward elimination
    for (std::size_t i = 1; i < n_eta - 1; ++i) {
            Eigen::PartialPivLU<Eigen::MatrixXd> lu_B(main_blocks[i-1].eigen());
            upper_blocks[i-1] = (lu_B.solve(upper_blocks[i-1].eigen())).eval();
            main_blocks[i] = (main_blocks[i].eigen() -
                            lower_blocks[i].eigen() * upper_blocks[i-1].eigen()).eval();
            
            // Update RHS
            for (int j = 0; j < n_heavy; ++j) {
                double temp = 0.0;
                for (int k = 0; k < n_heavy; ++k) {
                    temp += lower_blocks[i](j, k) * rhs(i-1, k + start_species);
                }
                rhs(i, j + start_species) -= temp;
            }
    }
    
    // Back substitution
    // Last interior point
    Eigen::PartialPivLU<Eigen::MatrixXd> lu_B(main_blocks[n_eta - 2].eigen());
    Eigen::VectorXd rhs_vec(n_heavy);
    for (int k = 0; k < n_heavy; ++k) {
        rhs_vec[k] = rhs(n_eta - 2, k + start_species);
    }
    Eigen::VectorXd sol = lu_B.solve(rhs_vec);
    for (int j = 0; j < n_heavy; ++j) {
        solution(j + start_species, n_eta - 2) = sol[j];
    }
    
    // Remaining points
    for (int i = n_eta - 3; i >= 0; --i) {
        for (int j = 0; j < n_heavy; ++j) {
            solution(j + start_species, i) = rhs(i, j + start_species);
            for (int k = 0; k < n_heavy; ++k) {
                solution(j + start_species, i) -= 
                    upper_blocks[i](j, k) * solution(k + start_species, i+1);
            }
        }
    }
}

} // namespace detail

} // namespace blast::boundary_layer::solvers