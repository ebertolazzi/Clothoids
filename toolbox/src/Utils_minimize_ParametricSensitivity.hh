/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2025                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: Utils_minimize_ParametricSensitivity.hh
//

/*--------------------------------------------------------------------------*\
 |  Parametric Sensitivity Analysis for Optimization Problems               |
 |                                                                          |
 |  Given: f(x, p) where p are fixed parameters                             |
 |  At optimum x*(p): ∇_x f(x*, p) = 0                                      |
 |                                                                          |
 |  Computes: dx/dp using the Implicit Function Theorem                     |
 |                                                                          |
 |  === UNCONSTRAINED CASE ===                                              |
 |  Formula: dx/dp = -[∇²_xx f]^(-1) · ∇²_xp f                              |
 |                                                                          |
 |  === BOX CONSTRAINED CASE: l ≤ x ≤ u ===                                 |
 |  At optimum, KKT conditions:                                             |
 |    ∇_x f + λ = 0                                                         |
 |    λᵢ ≥ 0 if xᵢ = lᵢ (lower active)                                      |
 |    λᵢ ≤ 0 if xᵢ = uᵢ (upper active)                                      |
 |    λᵢ = 0 if lᵢ < xᵢ < uᵢ (free)                                         |
 |                                                                          |
 |  Differentiating w.r.t. p (assuming bounds don't depend on p):           |
 |    For FREE variables i ∈ F:                                             |
 |      ∑_j H_ij (dx_j/dp_k) + G_ik = 0                                     |
 |      where H = ∇²_xx f, G = ∇²_xp f                                      |
 |                                                                          |
 |    For ACTIVE variables i ∈ A:                                           |
 |      dx_i/dp_k = 0  (remains at bound)                                   |
 |                                                                          |
 |  Reduced system (solve only for free variables):                         |
 |    H_FF · (dx/dp)_F = -G_F                                               |
 |                                                                          |
 |  === REGULARIZED CASE: min f(x,p) + ε‖x‖² ===                            |
 |  Effective objective: f_ε(x,p) = f(x,p) + ε‖x‖²                          |
 |  At optimum: ∇_x f + 2εx* = 0                                            |
 |                                                                          |
 |  Differentiating:                                                        |
 |    ∇²_xx f · (dx/dp) + ∇²_xp f + 2ε(dx/dp) = 0                           |
 |    [∇²_xx f + 2εI] · (dx/dp) = -∇²_xp f                                  |
 |                                                                          |
 |  Formula: dx/dp = -[H + 2εI]^(-1) · G                                    |
 |                                                                          |
 |  Effect of regularization:                                               |
 |  - Makes Hessian more positive definite (better conditioned)             |
 |  - Reduces sensitivity magnitude: ‖dx/dp‖ decreases with ε               |
 |  - Solution x*(p) becomes "smoother" function of p                       |
 |                                                                          |
 |  === IMPORTANT NOTES ===                                                 |
 |  1. Active set discontinuities:                                          |
 |     - dx/dp is NOT continuous when active set changes!                   |
 |     - Need active set identification at each p                           |
 |     - Sensitivity valid only for small Δp where active set stable        |
 |                                                                          |
 |  2. If bounds depend on p: l(p), u(p):                                   |
 |     - For active at lower: x_i* = l_i(p)                                 |
 |       ⟹ dx_i/dp = dl_i/dp                                               |
 |     - For active at upper: x_i* = u_i(p)                                 |
 |       ⟹ dx_i/dp = du_i/dp                                               |
 |                                                                          |
 |  References:                                                             |
 |  - Nocedal & Wright, "Numerical Optimization", Section 19.3              |
 |  - Bonnans & Shapiro, "Perturbation Analysis of Optimization Problems"   |
 |  - Fiacco, "Introduction to Sensitivity and Stability Analysis"          |
\*--------------------------------------------------------------------------*/

#pragma once

#ifndef UTILS_MINIMIZE_PARAMETRIC_SENSITIVITY_dot_HH
#define UTILS_MINIMIZE_PARAMETRIC_SENSITIVITY_dot_HH

#include "Utils_minimize.hh"
#include "Utils_minimize_Newton.hh"
#include "Utils_fmt.hh"

namespace Utils
{

  // ===========================================================================
  // Parametric Sensitivity Analyzer
  // ===========================================================================

  /**
   * @brief Computes sensitivity of optimal solution to parameter changes
   *
   * For optimization problem:
   *   min_x f(x, p)
   *
   * At optimal x*(p), computes dx/dp using implicit function theorem.
   *
   * @tparam Scalar Floating point type (double/float)
   */
  template <typename Scalar = double> class ParametricSensitivity
  {
  public:
    using Vector       = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Matrix       = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;

    /**
     * @brief Callback type for parametric objective function
     *
     * @param x Decision variables
     * @param p Parameters
     * @param grad_x [out] Gradient w.r.t. x (∇_x f)
     * @param hess_xx [out] Hessian w.r.t. x (∇²_xx f)
     * @param grad_xp [out] Mixed derivatives (∇²_xp f), size n_x × n_p
     * @return Function value f(x, p)
     */
    using ParametricCallback = std::function<Scalar(
      Vector const & x,
      Vector const & p,
      Vector *       grad_x,
      SparseMatrix * hess_xx,
      Matrix *       grad_xp  // ∂²f/∂x_i∂p_j
      )>;

    /**
     * @brief Active set information at optimum
     */
    struct ActiveSet
    {
      std::vector<size_t> lower_active;  // Indices where x_i = l_i
      std::vector<size_t> upper_active;  // Indices where x_i = u_i
      std::vector<size_t> free;          // Indices where l_i < x_i < u_i

      bool is_active( size_t i ) const
      {
        return std::find( lower_active.begin(), lower_active.end(), i ) != lower_active.end() ||
               std::find( upper_active.begin(), upper_active.end(), i ) != upper_active.end();
      }

      size_t n_free() const { return free.size(); }
      size_t n_active() const { return lower_active.size() + upper_active.size(); }
    };

    struct Options
    {
      // Numerical differentiation for ∇²_xp if not provided analytically
      bool   use_finite_differences{ false };
      Scalar fd_epsilon{ 1e-7 };

      // Solver options for sensitivity system
      Scalar solver_tolerance{ 1e-12 };

      // Bound constraints
      bool   has_bounds{ false };
      Scalar active_set_tolerance{ 1e-8 };  // Tolerance for detecting active constraints

      // Regularization
      bool   account_for_regularization{ false };
      Scalar regularization_epsilon{ 0 };  // ε in f(x,p) + ε‖x‖²

      // Verbosity
      size_t verbosity_level{ 1 };
    };

  private:
    Options m_opts;
    Scalar  m_epsi{ std::numeric_limits<Scalar>::epsilon() };

    // Results storage
    Matrix      m_sensitivity;  // dx/dp, size n_x × n_p
    ActiveSet   m_active_set;   // Active constraints at optimum
    Scalar      m_condition_number{ 0 };
    bool        m_success{ false };
    bool        m_active_set_changed{ false };  // Flag for potential discontinuity
    std::string m_error_message{ "" };

    /**
     * @brief Identify active constraints at optimum
     *
     * For box constraints l ≤ x ≤ u:
     * - Lower active: x_i ≈ l_i AND ∇f_i ≥ 0 (pushing against lower bound)
     * - Upper active: x_i ≈ u_i AND ∇f_i ≤ 0 (pushing against upper bound)
     * - Free: l_i < x_i < u_i (interior)
     */
    ActiveSet identify_active_set( Vector const & x, Vector const & grad, Vector const & lower, Vector const & upper )
      const
    {
      ActiveSet    aset;
      Eigen::Index n = x.size();

      for ( Eigen::Index i = 0; i < n; ++i )
      {
        bool at_lower = std::abs( x( i ) - lower( i ) ) < m_opts.active_set_tolerance;
        bool at_upper = std::abs( x( i ) - upper( i ) ) < m_opts.active_set_tolerance;

        // Check complementarity: gradient must push against active bound
        bool lower_active = at_lower && grad( i ) > -m_opts.active_set_tolerance;
        bool upper_active = at_upper && grad( i ) < m_opts.active_set_tolerance;

        if ( lower_active ) { aset.lower_active.push_back( static_cast<size_t>( i ) ); }
        else if ( upper_active ) { aset.upper_active.push_back( static_cast<size_t>( i ) ); }
        else
        {
          aset.free.push_back( static_cast<size_t>( i ) );
        }
      }

      return aset;
    }

    /**
     * @brief Extract submatrix for free variables
     */
    Matrix extract_free_submatrix( SparseMatrix const & A, std::vector<size_t> const & free_indices ) const
    {
      Eigen::Index n_free = static_cast<Eigen::Index>( free_indices.size() );
      Matrix       A_free( n_free, n_free );

      for ( Eigen::Index i = 0; i < n_free; ++i )
      {
        for ( Eigen::Index j = 0; j < n_free; ++j )
        {
          A_free( i, j ) = A.coeff(
            static_cast<Eigen::Index>( free_indices[i] ),
            static_cast<Eigen::Index>( free_indices[j] ) );
        }
      }

      return A_free;
    }

    Vector extract_free_subvector( Vector const & v, std::vector<size_t> const & free_indices ) const
    {
      Eigen::Index n_free = static_cast<Eigen::Index>( free_indices.size() );
      Vector       v_free( n_free );

      for ( Eigen::Index i = 0; i < n_free; ++i ) { v_free( i ) = v( static_cast<Eigen::Index>( free_indices[i] ) ); }

      return v_free;
    }

    Matrix extract_free_rows( Matrix const & M, std::vector<size_t> const & free_indices ) const
    {
      Eigen::Index n_free = static_cast<Eigen::Index>( free_indices.size() );
      Eigen::Index n_cols = M.cols();
      Matrix       M_free( n_free, n_cols );

      for ( Eigen::Index i = 0; i < n_free; ++i )
      {
        M_free.row( i ) = M.row( static_cast<Eigen::Index>( free_indices[i] ) );
      }

      return M_free;
    }

    /**
     * @brief Compute ∇²_xp f using finite differences
     *
     * For each parameter p_j:
     *   ∇²_xp[:, j] ≈ [∇_x f(x, p + ε e_j) - ∇_x f(x, p)] / ε
     */
    Matrix compute_mixed_derivatives_fd( Vector const & x, Vector const & p, ParametricCallback const & callback ) const
    {
      Eigen::Index n_x = x.size();
      Eigen::Index n_p = p.size();

      Matrix grad_xp( n_x, n_p );

      // Base gradient
      Vector grad_base( n_x );
      callback( x, p, &grad_base, nullptr, nullptr );

      // Perturb each parameter
      Vector p_pert = p;
      Vector grad_pert( n_x );

      for ( Eigen::Index j = 0; j < n_p; ++j )
      {
        Scalar h = m_opts.fd_epsilon * std::max( Scalar( 1 ), std::abs( p( j ) ) );

        p_pert( j ) = p( j ) + h;
        callback( x, p_pert, &grad_pert, nullptr, nullptr );

        grad_xp.col( j ) = ( grad_pert - grad_base ) / h;

        p_pert( j ) = p( j );  // Restore
      }

      return grad_xp;
    }

    /**
     * @brief Solve sensitivity system: H · S = -G
     * where H = ∇²_xx f, S = dx/dp, G = ∇²_xp f
     *
     * For constrained problems, only solve for free variables.
     * Active variables remain at bounds: dx_active/dp = 0
     */
    bool solve_sensitivity_system( SparseMatrix const & hess_xx, Matrix const & grad_xp, ActiveSet const & active_set )
    {
      m_active_set = active_set;
      m_success    = false;

      Eigen::Index n_x = hess_xx.rows();
      Eigen::Index n_p = grad_xp.cols();

      m_sensitivity.resize( n_x, n_p );
      m_sensitivity.setZero();

      // If no free variables, all sensitivity is zero (stuck at bounds)
      if ( active_set.n_free() == 0 )
      {
        m_success       = true;
        m_error_message = "All variables at bounds - zero sensitivity";

        if ( m_opts.verbosity_level >= 1 )
        {
          fmt::print( fmt::fg( fmt::color::yellow ), "[Sensitivity] WARNING: No free variables - all at bounds\n" );
        }

        return true;
      }

      // Extract reduced system for free variables only
      Matrix H_free = extract_free_submatrix( hess_xx, active_set.free );
      Matrix G_free = extract_free_rows( grad_xp, active_set.free );

      Eigen::Index n_free = static_cast<Eigen::Index>( active_set.n_free() );

      if ( m_opts.verbosity_level >= 2 )
      {
        fmt::print( "[Sensitivity] Solving reduced system: {} free / {} total\n", n_free, n_x );
        fmt::print(
          "  Lower active: {}, Upper active: {}\n",
          active_set.lower_active.size(),
          active_set.upper_active.size() );
      }

      // Solve reduced system: H_free · S_free = -G_free
      Eigen::LDLT<Matrix> ldlt( H_free );

      if ( ldlt.info() != Eigen::Success )
      {
        m_success       = false;
        m_error_message = "Reduced Hessian factorization failed";

        if ( m_opts.verbosity_level >= 1 )
        {
          fmt::print( fmt::fg( fmt::color::red ), "[Sensitivity] ERROR: Reduced Hessian not positive definite\n" );
        }

        return false;
      }

      // Solve for each parameter
      Matrix S_free( n_free, n_p );
      for ( Eigen::Index j = 0; j < n_p; ++j )
      {
        Vector rhs      = -G_free.col( j );
        S_free.col( j ) = ldlt.solve( rhs );
      }

      // Insert back into full sensitivity matrix
      for ( Eigen::Index i = 0; i < n_free; ++i )
      {
        size_t idx                                            = active_set.free[static_cast<size_t>( i )];
        m_sensitivity.row( static_cast<Eigen::Index>( idx ) ) = S_free.row( i );
      }

      // Active variables have zero sensitivity (remain at bounds)
      // Already set to zero by initialization

      m_success = true;

      // Estimate condition number of reduced Hessian
      Eigen::SelfAdjointEigenSolver<Matrix> eigensolver( H_free );
      if ( eigensolver.info() == Eigen::Success )
      {
        Vector eigenvalues = eigensolver.eigenvalues();
        Scalar lambda_min  = eigenvalues.minCoeff();
        Scalar lambda_max  = eigenvalues.maxCoeff();

        if ( lambda_min > m_epsi ) { m_condition_number = lambda_max / lambda_min; }
      }

      return true;
    }

  public:
    // =======================================================================
    // Getter methods for results
    // =======================================================================

    /**
     * @brief Get sensitivity matrix dx/dp
     */
    Matrix const & sensitivity() const { return m_sensitivity; }

    /**
     * @brief Get active set information
     */
    ActiveSet const & active_set() const { return m_active_set; }

    /**
     * @brief Get condition number of reduced Hessian
     */
    Scalar condition_number() const { return m_condition_number; }

    /**
     * @brief Check if sensitivity computation was successful
     */
    bool success() const { return m_success; }

    /**
     * @brief Check if active set changed during sensitivity computation
     */
    bool active_set_changed() const { return m_active_set_changed; }

    /**
     * @brief Get error message if computation failed
     */
    std::string const & error_message() const { return m_error_message; }

    /**
     * @brief Get the norm of the sensitivity matrix
     */
    Scalar sensitivity_norm() const { return m_sensitivity.norm(); }

    /**
     * @brief Get maximum absolute sensitivity value
     */
    Scalar max_sensitivity() const { return m_sensitivity.cwiseAbs().maxCoeff(); }

  public:
    ParametricSensitivity( Options opts = Options() ) : m_opts( opts ) {}

    /**
     * @brief Compute dx/dp at optimal point x*(p)
     *
     * @param x_opt Optimal solution x* for given p
     * @param p Current parameter values
     * @param callback Function providing f, ∇_x f, ∇²_xx f, and optionally ∇²_xp f
     * @param lower Lower bounds (empty if unconstrained)
     * @param upper Upper bounds (empty if unconstrained)
     *
     * For box-constrained problems l ≤ x ≤ u:
     * - At optimum, some variables may be at bounds (active)
     * - KKT conditions: ∇f + λ = 0, where λᵢ ≠ 0 only if xᵢ at bound
     * - Sensitivity only for free variables (interior points)
     * - Active variables: dx_i/dp = 0 (remain at bound, assuming bound doesn't depend on p)
     *
     * For regularized problems f(x,p) + ε‖x‖²:
     * - Effective Hessian: H_eff = H + 2εI
     * - Mixed derivatives unchanged: ∇²_xp f
     * - Automatically makes problem more stable (ε > 0)
     */
    void compute_sensitivity(
      Vector const &             x_opt,
      Vector const &             p,
      ParametricCallback const & callback,
      Vector const &             lower = Vector(),
      Vector const &             upper = Vector() )
    {
      Eigen::Index n_x = x_opt.size();
      Eigen::Index n_p = p.size();

      bool has_bounds = ( lower.size() == n_x && upper.size() == n_x );

      if ( m_opts.verbosity_level >= 1 )
      {
        fmt::print(
          fmt::fg( fmt::color::cyan ),
          "\n"
          "╔══════════════════════════════════════════════════════════╗\n"
          "║           Parametric Sensitivity Analysis                ║\n"
          "╠══════════════════════════════════════════════════════════╣\n"
          "║ Variables (x): {:<41} ║\n"
          "║ Parameters (p): {:<40} ║\n"
          "║ Constraints: {:<43} ║\n"
          "║ Regularization: ε = {:<36.2e} ║\n"
          "║ Method: {:<48} ║\n"
          "╚══════════════════════════════════════════════════════════╝\n",
          static_cast<int>( n_x ),
          static_cast<int>( n_p ),
          has_bounds ? "Box constraints" : "Unconstrained",
          m_opts.regularization_epsilon,
          m_opts.use_finite_differences ? "Finite Differences" : "Analytical" );
      }

      // Reset results
      m_success = false;
      m_error_message.clear();
      m_condition_number   = 0;
      m_active_set_changed = false;

      // Evaluate at optimal point
      Vector       grad_x( n_x );
      SparseMatrix hess_xx( n_x, n_x );
      Matrix       grad_xp( n_x, n_p );

      Scalar f_opt = callback( x_opt, p, &grad_x, &hess_xx, &grad_xp );
      (void) f_opt;  // Suppress unused variable warning

      // Add regularization to Hessian if specified
      if ( m_opts.account_for_regularization && m_opts.regularization_epsilon > 0 )
      {
        std::vector<Eigen::Triplet<Scalar>> reg_triplets;

        // Extract existing triplets
        for ( Eigen::Index k = 0; k < hess_xx.outerSize(); ++k )
        {
          for ( typename SparseMatrix::InnerIterator it( hess_xx, k ); it; ++it )
          {
            reg_triplets.emplace_back( it.row(), it.col(), it.value() );
          }
        }

        // Add 2ε to diagonal
        for ( Eigen::Index i = 0; i < n_x; ++i )
        {
          reg_triplets.emplace_back( i, i, 2 * m_opts.regularization_epsilon );
        }

        hess_xx.setFromTriplets( reg_triplets.begin(), reg_triplets.end() );

        if ( m_opts.verbosity_level >= 2 )
        {
          fmt::print(
            "[Sensitivity] Added regularization 2ε = {:.2e} to Hessian diagonal\n",
            2 * m_opts.regularization_epsilon );
        }
      }

      // Identify active set if constrained
      ActiveSet active_set;
      if ( has_bounds )
      {
        active_set = identify_active_set( x_opt, grad_x, lower, upper );

        if ( m_opts.verbosity_level >= 1 )
        {
          fmt::print( "[Sensitivity] Active set analysis:\n" );
          fmt::print( "  Free variables:  {} / {}\n", active_set.n_free(), static_cast<size_t>( n_x ) );
          fmt::print( "  Lower active:    {}\n", active_set.lower_active.size() );
          fmt::print( "  Upper active:    {}\n", active_set.upper_active.size() );

          if ( m_opts.verbosity_level >= 3 && active_set.n_active() > 0 )
          {
            fmt::print( "  Lower active indices: " );
            for ( auto i : active_set.lower_active ) fmt::print( "{} ", i );
            fmt::print( "\n" );
            fmt::print( "  Upper active indices: " );
            for ( auto i : active_set.upper_active ) fmt::print( "{} ", i );
            fmt::print( "\n" );
          }
        }
      }
      else
      {
        // Unconstrained: all variables are free
        active_set.free.resize( static_cast<size_t>( n_x ) );
        std::iota( active_set.free.begin(), active_set.free.end(), 0 );
      }

      // Check optimality for free variables
      Vector grad_free = extract_free_subvector( grad_x, active_set.free );
      Scalar grad_norm = grad_free.norm();

      if ( grad_norm > 1e-4 )
      {
        if ( m_opts.verbosity_level >= 1 )
        {
          fmt::print(
            fmt::fg( fmt::color::yellow ),
            "[Sensitivity] WARNING: ‖∇_x f (free)‖ = {:.2e} > 1e-4\n"
            "              Point may not be optimal for free variables!\n",
            grad_norm );
        }
      }

      // Compute mixed derivatives if not provided
      bool grad_xp_provided = ( grad_xp.rows() == n_x && grad_xp.cols() == n_p );

      if ( !grad_xp_provided || m_opts.use_finite_differences )
      {
        if ( m_opts.verbosity_level >= 2 )
        {
          fmt::print( "[Sensitivity] Computing ∇²_xp using finite differences...\n" );
        }
        grad_xp = compute_mixed_derivatives_fd( x_opt, p, callback );

        // Regularization doesn't affect mixed derivatives (no x in ε‖x‖²'s p-dependence)
      }

      // Solve sensitivity system (handles active set internally)
      bool solved = solve_sensitivity_system( hess_xx, grad_xp, active_set );

      // Print results
      if ( m_opts.verbosity_level >= 1 )
      {
        if ( solved )
        {
          fmt::print( fmt::fg( fmt::color::green ), "\n[Sensitivity] ✓ Successfully computed dx/dp\n" );

          if ( has_bounds )
          {
            fmt::print( "  Sensitivity computed for {} free variables\n", active_set.n_free() );
            fmt::print( "  {} active variables have zero sensitivity\n", active_set.n_active() );
          }

          fmt::print( "  Condition number (free vars): {:.2e}\n", m_condition_number );

          Scalar max_sens = m_sensitivity.cwiseAbs().maxCoeff();
          fmt::print( "  Max |∂x_i*/∂p_j|: {:.2e}\n", max_sens );

          if ( m_opts.verbosity_level >= 2 )
          {
            fmt::print( "\n  Sensitivity matrix (dx/dp) - showing free variables:\n" );
            size_t n_show = std::min( size_t( 5 ), active_set.n_free() );

            for ( size_t k = 0; k < n_show; ++k )
            {
              size_t i = active_set.free[k];
              fmt::print( "    x[{:2d}]: ", static_cast<int>( i ) );
              for ( Eigen::Index j = 0; j < n_p; ++j )
              {
                fmt::print( "{:12.4e} ", m_sensitivity( static_cast<Eigen::Index>( i ), j ) );
              }
              fmt::print( "\n" );
            }
            if ( active_set.n_free() > 5 ) { fmt::print( "    ... ({} free variables total)\n", active_set.n_free() ); }

            if ( active_set.n_active() > 0 && m_opts.verbosity_level >= 3 )
            {
              fmt::print( "\n  Active variables (zero sensitivity):\n" );
              for ( auto i : active_set.lower_active )
              {
                fmt::print( "    x[{:2d}] at lower bound\n", static_cast<int>( i ) );
              }
              for ( auto i : active_set.upper_active )
              {
                fmt::print( "    x[{:2d}] at upper bound\n", static_cast<int>( i ) );
              }
            }
          }
        }
        else
        {
          fmt::print( fmt::fg( fmt::color::red ), "\n[Sensitivity] ✗ Failed: {}\n", m_error_message );
        }
      }
    }

    /**
     * @brief Predict optimal solution at perturbed parameters
     *
     * Linear approximation: x*(p + δp) ≈ x*(p) + (dx/dp) · δp
     *
     * @param x_opt Current optimal x*(p)
     * @param delta_p Parameter perturbation δp
     * @return Predicted new optimum
     */
    Vector predict_new_optimum( Vector const & x_opt, Vector const & delta_p ) const
    {
      return x_opt + m_sensitivity * delta_p;
    }
  };

  // ===========================================================================
  // Parametric Newton Minimizer
  // ===========================================================================

  /**
   * @brief Newton minimizer with built-in parametric sensitivity analysis
   */
  template <typename Scalar = double> class ParametricNewtonMinimizer
  {
  public:
    using Vector             = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Matrix             = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMatrix       = Eigen::SparseMatrix<Scalar>;
    using ParametricCallback = typename ParametricSensitivity<Scalar>::ParametricCallback;
    using Status             = typename Newton_minimizer<Scalar>::Status;
    using integer            = Eigen::Index;

    struct Options
    {
      typename Newton_minimizer<Scalar>::Options      optimizer_opts;
      typename ParametricSensitivity<Scalar>::Options sensitivity_opts;
      bool                                            compute_sensitivity{ true };
    };

  private:
    Options                       m_opts;
    Vector                        m_current_p;
    ParametricSensitivity<Scalar> m_sens_analyzer;
    Newton_minimizer<Scalar>      m_optimizer;

    // Storage for optimization results
    Vector  m_x_opt;
    Scalar  m_f_opt{ 0 };
    Status  m_status;
    integer m_iterations{ 0 };
    integer m_function_evals{ 0 };
    integer m_hessian_evals{ 0 };
    Scalar  m_final_grad_norm{ 0 };
    bool    m_optimization_success{ false };

  public:
    ParametricNewtonMinimizer( Options opts = Options() )
      : m_opts( opts ), m_sens_analyzer( m_opts.sensitivity_opts ), m_optimizer( m_opts.optimizer_opts )
    {
    }

    /**
     * @brief Minimize f(x, p) and compute dx/dp
     */
    void minimize( Vector const & x0, Vector const & p, ParametricCallback const & parametric_callback )
    {
      m_current_p = p;

      // Reset results
      m_optimization_success = false;

      // Create standard callback for optimizer (fixing p)
      auto fixed_p_callback = [&]( Vector const & x, Vector * grad, SparseMatrix * hess ) -> Scalar
      { return parametric_callback( x, m_current_p, grad, hess, nullptr ); };

      // Optimize
      m_optimizer.minimize( x0, fixed_p_callback );

      // Store optimization results using getters
      m_x_opt           = m_optimizer.solution();
      m_f_opt           = m_optimizer.final_f();
      m_status          = m_optimizer.status();
      m_iterations      = m_optimizer.iterations();
      m_function_evals  = m_optimizer.function_evals();
      m_hessian_evals   = m_optimizer.hessian_evals();
      m_final_grad_norm = m_optimizer.final_grad_norm();

      m_optimization_success = m_status == Status::CONVERGED;

      // Compute sensitivity if requested
      if ( m_opts.compute_sensitivity && m_optimization_success )
      {
        m_sens_analyzer.compute_sensitivity( m_x_opt, m_current_p, parametric_callback );
      }
    }

    // =======================================================================
    // Getter methods for optimization results
    // =======================================================================

    /**
     * @brief Get optimal solution
     */
    Vector const & solution() const { return m_x_opt; }

    /**
     * @brief Get final function value
     */
    Scalar final_f() const { return m_f_opt; }

    /**
     * @brief Get optimization status
     */
    Status status() const { return m_status; }

    /**
     * @brief Get number of iterations
     */
    integer iterations() const { return m_iterations; }

    /**
     * @brief Get number of function evaluations
     */
    integer function_evals() const { return m_function_evals; }

    /**
     * @brief Get number of Hessian evaluations
     */
    integer hessian_evals() const { return m_hessian_evals; }

    /**
     * @brief Get final gradient norm
     */
    Scalar final_grad_norm() const { return m_final_grad_norm; }

    /**
     * @brief Check if optimization was successful
     */
    bool optimization_success() const { return m_optimization_success; }

    // =======================================================================
    // Getter methods for sensitivity results
    // =======================================================================

    /**
     * @brief Get sensitivity matrix dx/dp
     */
    Matrix const & sensitivity() const { return m_sens_analyzer.sensitivity(); }

    /**
     * @brief Get active set information
     */
    typename ParametricSensitivity<Scalar>::ActiveSet const & active_set() const
    {
      return m_sens_analyzer.active_set();
    }

    /**
     * @brief Get condition number of reduced Hessian
     */
    Scalar condition_number() const { return m_sens_analyzer.condition_number(); }

    /**
     * @brief Check if sensitivity computation was successful
     */
    bool sensitivity_success() const { return m_sens_analyzer.success(); }

    /**
     * @brief Check if active set changed during sensitivity computation
     */
    bool active_set_changed() const { return m_sens_analyzer.active_set_changed(); }

    /**
     * @brief Get error message if sensitivity computation failed
     */
    std::string const & sensitivity_error_message() const { return m_sens_analyzer.error_message(); }

    // =======================================================================
    // Utility methods
    // =======================================================================

    /**
     * @brief Get reference to the sensitivity analyzer for direct access
     */
    ParametricSensitivity<Scalar> const & sensitivity_analyzer() const { return m_sens_analyzer; }

    /**
     * @brief Get reference to the sensitivity analyzer for direct access (non-const)
     */
    ParametricSensitivity<Scalar> & sensitivity_analyzer() { return m_sens_analyzer; }

    /**
     * @brief Get reference to the optimizer for direct access
     */
    Newton_minimizer<Scalar> const & optimizer() const { return m_optimizer; }

    /**
     * @brief Get reference to the optimizer for direct access (non-const)
     */
    Newton_minimizer<Scalar> & optimizer() { return m_optimizer; }

    /**
     * @brief Predict optimal solution at perturbed parameters
     *
     * Linear approximation: x*(p + δp) ≈ x*(p) + (dx/dp) · δp
     *
     * @param delta_p Parameter perturbation δp
     * @return Predicted new optimum
     */
    Vector predict_new_optimum( Vector const & delta_p ) const { return m_x_opt + sensitivity() * delta_p; }
  };

}  // namespace Utils

#endif
