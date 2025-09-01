#include "gasp2/catalysis/finite_rate/steady_surface_coverage.hpp"
#include "gasp2/catalysis/finite_rate/stoichiometry.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace gasp2::catalysis::finite_rate {
namespace details {
namespace {

//=========================== DEBUG UTILITIES ===============================//

/**
 * @brief Print a vector with a label when debug mode is enabled.
 *
 * @param label Text label preceding the vector values.
 * @param v     Vector to be printed.
 */
static void debug_print_vector(const std::string &label,
                               const std::vector<double> &v) {
  if (!ctx.debug)
    return;
  std::cout << label;
  for (double value : v)
    std::cout << ' ' << value;
  std::cout << '\n';
}

/**
 * @brief Print a matrix with a label when debug mode is enabled.
 *
 * Rows are printed one per line following the supplied label.
 *
 * @param label Text label preceding the matrix.
 * @param M     Matrix to be printed.
 */
static void debug_print_matrix(const std::string &label,
                               const std::vector<std::vector<double>> &M) {
  if (!ctx.debug)
    return;
  std::cout << label << '\n';
  for (const auto &row : M) {
    for (double value : row)
      std::cout << value << ' ';
    std::cout << '\n';
  }
}

/**
 * @brief Print an integer matrix with a label when debug mode is enabled.
 *
 * Rows are printed one per line following the supplied label.
 *
 * @param label Text label preceding the matrix.
 * @param M     Integer matrix to be printed.
 */
static void debug_print_int_matrix(const std::string &label,
                                   const std::vector<std::vector<int>> &M) {
  if (!ctx.debug)
    return;
  std::cout << label << '\n';
  for (const auto &row : M) {
    for (int value : row)
      std::cout << value << ' ';
    std::cout << '\n';
  }
}

//============================== RESIDUAL ===================================//

/**
 * @brief Evaluate the residual vector for the current coverages.
 *
 * The reaction rate products (forward and backward) are assumed to be
 * precomputed for the supplied coverage state. This avoids recomputing the
 * expensive mass-action products when both the residual and the Jacobian are
 * required during the Newton iteration.
 *
 * @param props           Global reaction properties (provides site density).
 * @param nu_reactants    Stoichiometric coefficients of reactants
 *                        (indexed by species then reaction).
 * @param nu_products     Stoichiometric coefficients of products
 *                        (indexed by species then reaction).
 * @param surface_index   Indices in the global species array that correspond to
 *                        surface species (last entry represents vacant sites).
 * @param coverage        Current surface coverages followed by the vacant-site
 *                        density.
 * @param rate_f          Precomputed forward reaction rates for each reaction.
 * @param rate_b          Precomputed backward reaction rates for each
 *                        reaction.
 * @return Residual vector where each entry enforces zero net production for a
 *         surface species and the final entry enforces the site-balance
 *         constraint.
 */
static std::vector<double>
residual(const ReactionProperties &props,
         const std::vector<std::vector<int>> &nu_reactants,
         const std::vector<std::vector<int>> &nu_products,
         const std::vector<std::size_t> &surface_index,
         const std::vector<double> &coverage, const std::vector<double> &rate_f,
         const std::vector<double> &rate_b) {
  const std::size_t n_reactions = props.reactions.size();
  const std::size_t n_ads = coverage.size() - 1; // last entry: empty sites

  // Compute residuals for adsorbed species.
  std::vector<double> res(coverage.size(), 0.0);
  for (std::size_t j = 0; j < n_ads; ++j) {
    double ws = 0.0;
    std::size_t idx = surface_index[j];
    for (std::size_t r = 0; r < n_reactions; ++r) {
      int nu_p = nu_products[idx][r];
      int nu_r = nu_reactants[idx][r];
      ws += static_cast<double>(nu_p - nu_r) * (rate_f[r] - rate_b[r]);
    }
    res[j] = ws;
  }

  // Site balance: sum of all surface species equals total site density.
  res.back() = -props.site_density;
  for (double theta : coverage)
    res.back() += theta;

  return res;
}

//============================= ANALYTIC JACOBIAN ===========================//

/**
 * @brief Construct the analytic Jacobian for the steady-state residual.
 *
 * For each reaction \(i\) and coverage variable \(X_j\) the derivative is given
 * by
 *
 * \f[
 * \frac{\partial R_i}{\partial X_j} = k_{f,i}\nu'_{ji}X_j^{\nu'_{ji}-1}
 * \prod_{m\ne j}X_m^{\nu'_{mi}} -
 * k_{b,i}\nu''_{ji}X_j^{\nu''_{ji}-1}\prod_{m\ne j}X_m^{\nu''_{mi}}.
 * \f]
 *
 * Special cases are handled to avoid singularities:
 * - When the reaction order with respect to \(X_j\) is zero, the derivative
 *   is zero.
 * - If \(X_j = 0\) and the exponent is greater than one, the derivative is
 *   zero.
 * - If \(X_j = 0\) and the exponent equals one, the derivative reduces to
 *   the rate constant multiplied by the product of the remaining species
 *   terms.
 *
 * @param props           Global reaction properties containing rate constants.
 * @param nu_reactants    Reactant stoichiometric matrix.
 * @param nu_products     Product stoichiometric matrix.
 * @param species_order   Ordering of gas and surface species (gas first).
 * @param surface_index   Indices of surface species in the global ordering.
 * @param coverage        Current surface coverages followed by the vacant-site
 *                        density.
 * @param conc            Concentration vector (gas and surface) consistent with
 *                        @p species_order.
 * @param rate_f          Precomputed forward reaction rates.
 * @param rate_b          Precomputed backward reaction rates.
 * @return Jacobian matrix of size `coverage.size()` squared.
 */
static std::vector<std::vector<double>> analytic_jacobian(
    const ReactionProperties &props,
    const std::vector<std::vector<int>> &nu_reactants,
    const std::vector<std::vector<int>> &nu_products,
    const std::vector<std::string> &species_order,
    const std::vector<std::size_t> &surface_index,
    const std::vector<double> &coverage, const std::vector<double> &conc,
    const std::vector<double> &rate_f, const std::vector<double> &rate_b) {
  const std::size_t n_species = species_order.size();
  const std::size_t n_reactions = props.reactions.size();
  const std::size_t n_unknowns = coverage.size();
  const std::size_t n_ads = n_unknowns - 1;

  std::vector<std::vector<double>> J(n_unknowns,
                                     std::vector<double>(n_unknowns, 0.0));

  // Columns correspond to derivatives with respect to coverage[j].
  for (std::size_t j = 0; j < n_unknowns; ++j) {
    std::size_t idx_j = surface_index[j];
    double conc_j = conc[idx_j];

    // Derivative of reaction rates with respect to X_j.
    std::vector<double> drate_f(n_reactions, 0.0);
    std::vector<double> drate_b(n_reactions, 0.0);
    for (std::size_t r = 0; r < n_reactions; ++r) {
      int nu_rj = nu_reactants[idx_j][r];
      if (nu_rj > 0) {
        if (conc_j == 0.0) {
          if (nu_rj == 1) {
            double prod = props.reactions[r].kf;
            for (std::size_t m = 0; m < n_species; ++m) {
              if (m == idx_j)
                continue;
              int nu_r = nu_reactants[m][r];
              if (nu_r > 0)
                prod *= std::pow(conc[m], nu_r);
            }
            drate_f[r] = prod;
          }
        } else {
          drate_f[r] = rate_f[r] * static_cast<double>(nu_rj) / conc_j;
        }
      }

      int nu_pj = nu_products[idx_j][r];
      if (nu_pj > 0) {
        if (conc_j == 0.0) {
          if (nu_pj == 1) {
            double prod = props.reactions[r].kb;
            for (std::size_t m = 0; m < n_species; ++m) {
              if (m == idx_j)
                continue;
              int nu_p = nu_products[m][r];
              if (nu_p > 0)
                prod *= std::pow(conc[m], nu_p);
            }
            drate_b[r] = prod;
          }
        } else {
          drate_b[r] = rate_b[r] * static_cast<double>(nu_pj) / conc_j;
        }
      }
    }

    // Assemble column j of the Jacobian.
    for (std::size_t i = 0; i < n_ads; ++i) {
      std::size_t idx_i = surface_index[i];
      for (std::size_t r = 0; r < n_reactions; ++r) {
        int nu_p = nu_products[idx_i][r];
        int nu_r = nu_reactants[idx_i][r];
        J[i][j] += static_cast<double>(nu_p - nu_r) * (drate_f[r] - drate_b[r]);
      }
    }

    // Site balance row: derivative is 1 for all coverages.
    J[n_ads][j] = 1.0;
  }

  return J;
}

//============================ LINEAR SOLVER ================================//

/**
 * @brief Solve the linear system \(Ax = b\) using Eigen's LU decomposition.
 *
 * The std::vector-based matrix and right-hand side are converted to Eigen
 * objects, the system is solved via `fullPivLu`, and the resulting solution is
 * returned. Debug output of the system and solution is emitted when enabled.
 *
 * @param A Coefficient matrix stored as a row-major `std::vector` of rows.
 * @param b Right-hand side vector.
 * @return Solution vector \(x\).
 */
static std::vector<double>
solve_linear(const std::vector<std::vector<double>> &A,
             const std::vector<double> &b) {
  const std::size_t n = b.size();
  if (ctx.debug) {
    debug_print_matrix("Linear system matrix A", A);
    debug_print_vector("Linear system vector b", b);
  }
  Eigen::MatrixXd mat(n, n);
  Eigen::VectorXd rhs(n);
  for (std::size_t i = 0; i < n; ++i) {
    rhs(i) = b[i];
    for (std::size_t j = 0; j < n; ++j)
      mat(i, j) = A[i][j];
  }
  Eigen::VectorXd sol = mat.fullPivLu().solve(rhs);
  std::vector<double> x(n);
  for (std::size_t i = 0; i < n; ++i)
    x[i] = sol(i);
  if (ctx.debug)
    debug_print_vector("Linear solution x", x);
  return x;
}

} // namespace

//======================== NEWTON SOLVER INTERFACE ==========================//

/**
 * @brief Solve for steady-state surface coverages using Newton's method.
 *
 * This function computes the steady state surface coverages for a catalytic
 * system by solving a nonlinear system of equations. The system enforces:
 * 1. Zero net production for each adsorbed species (steady state)
 * 2. Conservation of total surface sites (site balance)
 *
 * The Newton-Raphson iterative method is used with an analytic Jacobian
 * to achieve quadratic convergence near the solution.
 *
 * @param props Reaction properties containing:
 *              - reactions: Vector of reaction data (rate constants, etc.)
 *              - concentration: Gas-phase species concentrations
 *              - site_density: Total density of surface sites
 * @return Vector of surface coverages followed by vacant site density.
 *         Returns the last iterate if convergence is not achieved.
 */
std::vector<double> solve_surface_coverage(const ReactionProperties &props) {
  // Precompute stoichiometric data shared by the solver and flux evaluation.
  // This extracts reaction stoichiometry matrices and species ordering.
  auto stoich = build_stoichiometry(props);
  const auto &nu_r = stoich.nu_reactants;    // Reactant stoichiometric matrix
  const auto &nu_p = stoich.nu_products;     // Product stoichiometric matrix
  const auto &species_order = stoich.species_order;  // Global species ordering
  const auto &surface_index = stoich.surface_index;  // Surface species indices
  // ===== DEBUG: Print stoichiometric matrices and species information =====
  if(ctx.debug) {
    std::cout << "=== Stoichiometric Matrices ===\n";
    debug_print_int_matrix("Nu Reactants", nu_r);
    debug_print_int_matrix("Nu Products", nu_p);
    std::cout << "Species Order:\n";
    for (const auto &sp : species_order)
      std::cout << sp << ' ';
    std::cout << '\n';
    std::cout << "Surface Species Indices:\n";
    for (std::size_t idx : surface_index)
      std::cout << idx << ' ';
    std::cout << '\n';
  }

  // Initial guess: all adsorbed species empty, all sites available.
  // theta[0..n-2] = surface coverages, theta[n-1] = vacant site density
  std::vector<double> theta(ctx.surface_species.size() + 1, 0.0);
  theta.back() = props.site_density;
  // ===== DEBUG: Print initial guess =====
  if (ctx.debug) {
    std::cout << "=== Initial Guess ===\n";
    debug_print_vector("theta", theta);
  }

  // Newton solver parameters
  const std::size_t n_unknowns = theta.size();
  constexpr std::size_t max_iter = 50;      // Maximum Newton iterations
  constexpr double tol = 1e-12;             // Convergence tolerance (Euclidean norm)

  if (ctx.debug)
    std::cout << "=== Beginning surface coverage solve ===\n";

  // Main Newton iteration loop
  for (std::size_t iter = 0; iter < max_iter; ++iter) {
    if (ctx.debug) {
      std::cout << "-- Iteration " << iter << " --\n";
      debug_print_vector("theta", theta);
    }

    // ========== CONCENTRATION & RATES ==================================
    // Step 1: Assemble full concentration vector (gas + surface species)
    // This combines gas-phase concentrations with current surface coverages
    std::vector<double> conc(species_order.size(), 0.0);
    
    // Fill gas-phase concentrations from input properties
    for (std::size_t i = 0; i < ctx.species_order.size(); ++i) {
      const auto &name = ctx.species_order[i];
      conc[i] = props.concentration.at(name);
    }
    
    // Fill surface concentrations from current Newton iterate
    for (std::size_t j = 0; j < n_unknowns; ++j)
      conc[surface_index[j]] = theta[j];
    // DEBUG
    if (ctx.debug) {
      std::cout << "Concentration vector (iteration " << iter << "):\n";
      debug_print_vector("conc", conc);
    }
    // Step 2: Precompute forward/backward reaction rates
    // These are the mass-action rate expressions for each elementary reaction
    std::vector<double> rate_f(props.reactions.size(), 0.0);  // Forward rates
    std::vector<double> rate_b(props.reactions.size(), 0.0);  // Backward rates
    
    for (std::size_t r = 0; r < props.reactions.size(); ++r) {
      const auto &rx = props.reactions[r];
      rate_f[r] = rx.kf;  // Start with forward rate constant
      rate_b[r] = rx.kb;  // Start with backward rate constant
      
      // Apply mass-action kinetics: rate = k * ∏(concentrations^stoichiometry)
      for (std::size_t m = 0; m < species_order.size(); ++m) {
        // Forward rate contribution
        int nu_r_coeff = nu_r[m][r];
        if (nu_r_coeff > 0)
          rate_f[r] *= std::pow(conc[m], nu_r_coeff);
          
        // Backward rate contribution
        int nu_p_coeff = nu_p[m][r];
        if (nu_p_coeff > 0)
          rate_b[r] *= std::pow(conc[m], nu_p_coeff);
      }
      
      if (ctx.debug)
        std::cout << "  Reaction " << r << ": kf=" << rx.kf << ", kb=" << rx.kb
                  << ", Kc=" << rx.Kc << ", Ka=" << rx.Ka
                  << ", rate_f=" << rate_f[r] << ", rate_b=" << rate_b[r]
                  << '\n';
    }

    // ========== RESIDUAL ==============================================
    // Step 3: Evaluate the nonlinear residual vector
    // Each equation enforces zero net production for surface species
    // plus a site balance constraint
    auto res =
        residual(props, nu_r, nu_p, surface_index, theta, rate_f, rate_b);
    debug_print_vector("residual", res);

    // Step 4: Check convergence using Euclidean (L2) norm
    // If ||residual||_2 < tolerance, we have found the solution
    double norm = 0.0;
    for (double v : res)
      norm += v * v;  // Sum of squares
    norm = std::sqrt(norm);  // Take square root for Euclidean norm
    if (ctx.debug)
      std::cout << "Residual norm = " << norm << '\n';
    if (norm < tol) {
      if (ctx.debug)
        std::cout << "Converged after " << iter << " iterations\n";
      return theta;  // Solution found!
    }

    // ========== JACOBIAN ==============================================
    // Step 5: Compute the analytic Jacobian matrix
    // J[i][j] = ∂(residual_i)/∂(theta_j)
    auto J = analytic_jacobian(props, nu_r, nu_p, species_order, surface_index,
                               theta, conc, rate_f, rate_b);
    debug_print_matrix("Jacobian", J);

    // Step 6: Solve the Newton linear system J * delta = -residual
    // This gives the Newton step direction
    for (double &v : res)
      v = -v;  // Negate residual to get right-hand side
    debug_print_vector("Right-hand side b", res);
    auto delta = solve_linear(J, res);
    debug_print_vector("Newton step", delta);

    // Step 7: Update the solution estimate
    // theta_new = theta_old + delta
    for (std::size_t i = 0; i < n_unknowns; ++i)
      theta[i] += delta[i];
    debug_print_vector("Updated theta", theta);
  }
  
  // If we reach here, Newton's method did not converge within max_iter
  if (ctx.debug)
    std::cout << "Newton solver did not converge after " << max_iter
              << " iterations. "<<std::endl;
  
  // Throw exception for non-convergence
  throw std::runtime_error("Newton solver failed to converge for surface coverage "
                          "calculation after " + std::to_string(max_iter) + 
                          " iterations. ");
}

} // namespace details
} // namespace gasp2::catalysis::finite_rate
