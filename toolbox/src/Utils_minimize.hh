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
// file: Utils_minimize_Newton.hh
//
/*--------------------------------------------------------------------------*\
 |  Improved Newton Minimizer with Enhanced Recovery Mechanisms             |
 |  Based on: Nocedal & Wright, "Numerical Optimization", 2nd Ed.           |
\*--------------------------------------------------------------------------*/

#pragma once

#ifndef UTILS_MINIMIZE_dot_HH
#define UTILS_MINIMIZE_dot_HH

#include <set>
#include <optional>
#include <limits>
#include "Utils_fmt.hh"
#include "Utils_Linesearch.hh"

namespace Utils
{

  /**
   * @class BoxConstraintHandler
   * @brief Efficient handler for box constraints with QP subproblem solving capabilities
   *
   * This class provides comprehensive handling of box constraints for optimization problems:
   * \f[
   *   \ell \leq x \leq u
   * \f]
   *
   * Features include:
   * - Projection operations
   * - Active set identification
   * - KKT condition checking
   * - QP subproblem solving with both sparse and dense Hessians
   * - Efficient vectorized operations using Eigen3
   *
   * @tparam Scalar Floating-point type (default: double)
   */
  template <typename Scalar = double> class BoxConstraintHandler
  {
  public:
    /// @brief Integer type for indexing
    using integer = Eigen::Index;

    /// @brief Vector type
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    /// @brief Integer vector type
    using VectorI = Eigen::Matrix<integer, Eigen::Dynamic, 1>;

    /// @brief Dense matrix type
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    /// @brief Sparse matrix type
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;

  private:
    Vector m_lower;            ///< Lower bounds ℓ
    Vector m_upper;            ///< Upper bounds u
    Vector m_tol_lower;        ///< Adaptive tolerances for lower bounds
    Vector m_tol_upper;        ///< Adaptive tolerances for upper bounds
    Vector m_lambda_lower;     ///< Multipliers for lower bounds
    Vector m_lambda_upper;     ///< Multipliers for upper bounds
    bool   m_active{ false };  ///< Whether constraints are active
    Scalar m_epsi;             ///< Machine epsilon or user-defined tolerance

  public:
    // =========================================================================
    //  Constructors and Destructor
    // =========================================================================

    /**
     * @brief Default constructor (no constraints)
     */
    BoxConstraintHandler() : m_epsi( std::numeric_limits<Scalar>::epsilon() ) {}

    /**
     * @brief Constructor with explicit bounds
     * @param lower Lower bounds vector
     * @param upper Upper bounds vector
     */
    BoxConstraintHandler( Vector const & lower, Vector const & upper )
      : m_lower( lower ), m_upper( upper ), m_active( true ), m_epsi( std::numeric_limits<Scalar>::epsilon() )
    {
      UTILS_ASSERT(
        lower.size() == upper.size(),
        "Lower and upper bounds must have same size (#lower={},#upper={})",
        lower.size(),
        upper.size() );
      UTILS_ASSERT( ( lower.array() <= upper.array() ).all(), "Lower bounds must be <= upper bounds" );
      update_tolerances();
      m_lambda_lower.resize( m_lower.size() );
      m_lambda_upper.resize( m_upper.size() );
    }

    /**
     * @brief Destructor
     */
    ~BoxConstraintHandler() = default;

    // =========================================================================
    //  Constraint Management
    // =========================================================================

    /**
     * @brief Set box constraints
     *
     * Resets all existing constraints and applies new bounds
     *
     * @param lower Lower bounds ℓ
     * @param upper Upper bounds u
     */
    void set_bounds( Vector const & lower, Vector const & upper )
    {
      UTILS_ASSERT(
        lower.size() == upper.size(),
        "Lower and upper bounds must have same size (#lower={},#upper={})",
        lower.size(),
        upper.size() );
      UTILS_ASSERT( ( lower.array() <= upper.array() ).all(), "Lower bounds must be <= upper bounds" );

      m_lower  = lower;
      m_upper  = upper;
      m_active = true;
      update_tolerances();
    }

    /**
     * @brief Update adaptive tolerances based on current epsilon
     *
     * Computes: tol = ε * (1 + |bound|) for numerical stability
     */
    void update_tolerances()
    {
      if ( !m_active ) return;

      // Vectorized tolerance computation
      m_tol_lower = m_epsi * ( Vector::Ones( m_lower.size() ).array() + m_lower.array().abs() );
      m_tol_upper = m_epsi * ( Vector::Ones( m_upper.size() ).array() + m_upper.array().abs() );
      m_lambda_lower.resize( m_tol_lower.size() );
      m_lambda_upper.resize( m_tol_upper.size() );
    }

    /**
     * @brief Set epsilon value for adaptive tolerances
     * @param eps New epsilon value
     */
    void set_epsilon( Scalar eps )
    {
      m_epsi = eps;
      update_tolerances();
    }

    /**
     * @brief Clear all constraints
     */
    void clear()
    {
      m_active = false;
      m_lower.resize( 0 );
      m_upper.resize( 0 );
      m_tol_lower.resize( 0 );
      m_tol_upper.resize( 0 );
      m_lambda_lower.resize( 0 );
      m_lambda_upper.resize( 0 );
    }

    // =========================================================================
    //  Accessors
    // =========================================================================

    /// @brief Check if constraints are active
    bool is_active() const { return m_active; }

    /// @brief Get problem dimension (0 if inactive)
    integer dimension() const { return m_active ? m_lower.size() : 0; }

    /// @brief Get current epsilon value
    Scalar epsilon() const { return m_epsi; }

    /// @brief Get lower bounds (read-only)
    Vector const & lower() const { return m_lower; }

    /// @brief Get upper bounds (read-only)
    Vector const & upper() const { return m_upper; }

    /// @brief Get lower bound tolerances (read-only)
    Vector const & tolerance_lower() const { return m_tol_lower; }

    /// @brief Get upper bound tolerances (read-only)
    Vector const & tolerance_upper() const { return m_tol_upper; }

    /**
     * @brief Get specific lower bound
     * @param i Index
     * @return Lower bound at index i
     */
    Scalar lower( integer i ) const { return m_lower( i ); }

    /**
     * @brief Get specific upper bound
     * @param i Index
     * @return Upper bound at index i
     */
    Scalar upper( integer i ) const { return m_upper( i ); }

    /// @brief Get lower bound multipliers (read-only)
    Vector const & lambda_lower() const { return m_lambda_lower; }

    /// @brief Get upper bound multipliers (read-only)
    Vector const & lambda_upper() const { return m_lambda_upper; }

    // =========================================================================
    //  Feasibility Checking
    // =========================================================================

    /**
     * @brief Check if point is within bounds
     *
     * Vectorized implementation using Eigen array operations
     *
     * @param x Point to check
     * @param extra_tol Additional tolerance margin
     * @return True if ℓ - tol ≤ x ≤ u + tol
     */
    bool is_feasible( Vector const & x, Scalar extra_tol = 0 ) const
    {
      if ( !m_active ) return true;

      // Vectorized feasibility check
      auto lower_ok = ( x.array() >= ( m_lower.array() - m_tol_lower.array() - extra_tol ) ).all();
      auto upper_ok = ( x.array() <= ( m_upper.array() + m_tol_upper.array() + extra_tol ) ).all();

      return lower_ok && upper_ok;
    }

    /**
     * @brief Count constraint violations
     *
     * Vectorized implementation using coefficient-wise operations
     *
     * @param x Point to check
     * @return Number of violated constraints
     */
    integer count_violations( Vector const & x ) const
    {
      if ( !m_active ) return 0;

      // Vectorized violation counting
      auto lower_viol = ( x.array() < ( m_lower.array() - m_tol_lower.array() ) ).template cast<integer>();
      auto upper_viol = ( x.array() > ( m_upper.array() + m_tol_upper.array() ) ).template cast<integer>();

      return lower_viol.sum() + upper_viol.sum();
    }

    /**
     * @brief Minimum distance to boundary
     *
     * Computes minimum perpendicular distance to any constraint boundary
     *
     * @param x Reference point
     * @return Minimum distance to boundary (infinity if no constraints)
     */
    Scalar min_distance_to_boundary( Vector const & x ) const
    {
      if ( !m_active ) return std::numeric_limits<Scalar>::infinity();

      // Vectorized distance computation
      auto dist_lower = ( x - m_lower ).array().abs();
      auto dist_upper = ( m_upper - x ).array().abs();

      return std::min( dist_lower.minCoeff(), dist_upper.minCoeff() );
    }

    // =========================================================================
    //  Projection Operations
    // =========================================================================

    /**
     * @brief Project point onto box (in-place)
     *
     * Vectorized projection: x = max(ℓ, min(x, u))
     *
     * @param x Point to project (modified in-place)
     */
    void project( Vector & x ) const
    {
      if ( m_active )
      {
        // Vectorized projection using coefficient-wise operations
        x = x.cwiseMax( m_lower ).cwiseMin( m_upper );
      }
    }

    /**
     * @brief Project point onto box (non-destructive)
     * @param x Original point
     * @return Projected point
     */
    Vector projected( Vector const & x ) const
    {
      Vector y = x;
      project( y );
      return y;
    }

    /**
     * @brief Project gradient according to optimality conditions
     *
     * For active constraints, zeroes gradient components pointing outward
     * Vectorized implementation using masks
     *
     * @param x Current point
     * @param g Gradient to project (modified in-place)
     */
    void project_gradient( Vector const & x, Vector & g ) const
    {
      if ( !m_active ) return;

      // Vectorized gradient projection using masks
      auto lower_mask = ( x.array() <= ( m_lower.array() + m_tol_lower.array() ) ) && ( g.array() > 0 );
      auto upper_mask = ( x.array() >= ( m_upper.array() - m_tol_upper.array() ) ) && ( g.array() < 0 );

      // Apply masks using select()
      g.array() = lower_mask.select( 0, g.array() );
      g.array() = upper_mask.select( 0, g.array() );
    }

    /**
     * @brief Project direction to maintain feasibility
     *
     * Zeroes direction components that would violate active constraints
     * Vectorized implementation
     *
     * @param x Current point
     * @param d Direction to project (modified in-place)
     */
    void project_direction( Vector const & x, Vector & d ) const
    {
      if ( !m_active ) return;

      // Vectorized direction projection using masks
      auto lower_mask = ( x.array() <= ( m_lower.array() + m_tol_lower.array() ) ) && ( d.array() < 0 );
      auto upper_mask = ( x.array() >= ( m_upper.array() - m_tol_upper.array() ) ) && ( d.array() > 0 );

      // Apply masks using select()
      d.array() = lower_mask.select( 0, d.array() );
      d.array() = upper_mask.select( 0, d.array() );
    }

    /**
     * @brief Projected gradient (non-destructive)
     * @param x Current point
     * @param g Original gradient
     * @return Projected gradient
     */
    Vector projected_gradient( Vector const & x, Vector const & g ) const
    {
      Vector pg = g;
      project_gradient( x, pg );
      return pg;
    }

    // =========================================================================
    //  Active Set and Multipliers
    // =========================================================================

    /**
     * @brief Identify active constraints using Bertsekas' method
     *
     * Vectorized identification of active constraints:
     * - -1: lower bound active
     * -  1: upper bound active
     * -  0: free variable
     *
     * @param x Current point
     * @param g Gradient ∇f(x)
     * @param active Output vector indicating active constraints
     * @return Number of active constraints
     */
    integer identify_active_set( Vector const & x, Vector const & g, VectorI & active ) const
    {
      active.setZero();
      if ( !m_active ) return 0;

      // Vectorized active set identification
      auto dist_lower = ( x - m_lower ).array();
      auto dist_upper = ( m_upper - x ).array();

      auto lower_active = ( dist_lower <= m_tol_lower.array() ) && ( g.array() > 0 );
      auto upper_active = ( dist_upper <= m_tol_upper.array() ) && ( g.array() < 0 );

      // Combine masks using select()
      active.array() = lower_active.select( -1, upper_active.select( 1, active.array() ) );

      return lower_active.template cast<integer>().sum() + upper_active.template cast<integer>().sum();
    }

    /**
     * @brief Counts the number of active constraints at a given point
     * @param x Point to check
     * @return Number of active constraints (variables at lower or upper bounds)
     */
    integer num_active( Vector const & x ) const
    {
      if ( !m_active ) return 0;
      return ( ( x.array() <= ( m_lower + m_tol_lower ).array() ) ||
               ( x.array() >= ( m_upper - m_tol_upper ).array() ) )
        .count();
    }

    /**
     * @brief Check if a specific variable is at a bound
     * @param i Index of the variable
     * @param x Current point
     * @return 1 if at upper bound, -1 if at lower bound, 0 if free
     */
    integer active_status( integer i, Vector const & x ) const
    {
      if ( !m_active ) return 0;

      if ( x( i ) <= ( m_lower( i ) + m_tol_lower( i ) ) ) return -1;
      if ( x( i ) >= ( m_upper( i ) - m_tol_upper( i ) ) ) return 1;
      return 0;
    }

    /**
     * @brief Compute optimal Lagrange multipliers
     *
     * Computes KKT multipliers for box constraints:
     * λ_lower = max(0, -g) when x at lower bound
     * λ_upper = max(0, g) when x at upper bound
     *
     * @param x Current point
     * @param g Gradient ∇f(x)
     * @param lambda_lower Output multipliers for lower bounds
     * @param lambda_upper Output multipliers for upper bounds
     */
    void update_multipliers( Vector const & x, Vector const & g )
    {
      m_lambda_lower.setZero();
      m_lambda_upper.setZero();

      if ( !m_active ) return;

      // Vectorized multiplier computation - CORRETTO
      auto at_lower = ( x.array() <= ( m_lower.array() + m_tol_lower.array() ) );
      auto at_upper = ( x.array() >= ( m_upper.array() - m_tol_upper.array() ) );

      // Per limiti inferiori: λ = max(0, ∇f(x))
      m_lambda_lower.array() = at_lower.select( g.array().cwiseMax( 0 ), 0 );

      // Per limiti superiori: λ = max(0, -∇f(x))
      m_lambda_upper.array() = at_upper.select( ( -g.array() ).cwiseMax( 0 ), 0 );
    }

    /**
     * @brief Compute complementarity error
     *
     * Computes maximum violation of complementarity conditions:
     * λ_lower·(x - ℓ) = 0 and λ_upper·(u - x) = 0
     *
     * @param x Current point
     * @return Maximum complementarity error
     */
    Scalar complementarity_error( Vector const & x ) const
    {
      if ( !m_active ) return 0;

      // Vectorized complementarity error computation
      auto comp_lower = m_lambda_lower.array() * ( x.array() - m_lower.array() );
      auto comp_upper = m_lambda_upper.array() * ( m_upper.array() - x.array() );

      return ( comp_lower.abs() + comp_upper.abs() ).maxCoeff();
    }

    // =========================================================================
    //  KKT Analysis
    // =========================================================================

    /**
     * @brief Projected gradient norm (∞-norm)
     *
     * Computes ‖P(x - ∇f(x)) - x‖_∞ for optimality measure
     *
     * @param x Current point
     * @param g Gradient ∇f(x)
     * @return Infinity norm of projected gradient
     */
    Scalar projected_gradient_norm_inf( Vector const & x, Vector const & g ) const
    {
      Vector pg = projected_gradient( x, g );
      return pg.template lpNorm<Eigen::Infinity>();
    }

    /**
     * @brief Projected gradient norm (2-norm)
     * @param x Current point
     * @param g Gradient ∇f(x)
     * @return 2-norm of projected gradient
     */
    Scalar projected_gradient_norm2( Vector const & x, Vector const & g ) const
    {
      Vector pg = projected_gradient( x, g );
      return pg.norm();
    }

    /**
     * @brief KKT stationarity error (∞-norm)
     *
     * Computes ‖∇f(x) - λ_lower + λ_upper‖_∞
     *
     * @param g Gradient ∇f(x)
     * @return Stationarity error
     */
    Scalar kkt_stationarity_error_inf( Vector const & g ) const
    {
      if ( !m_active ) return g.template lpNorm<Eigen::Infinity>();
      Vector stationarity = g - m_lambda_lower + m_lambda_upper;
      return stationarity.template lpNorm<Eigen::Infinity>();
    }

    /**
     * @brief KKT violation norm (2-norm)
     *
     * Computes weighted violation of KKT conditions
     *
     * @param x Current point
     * @param v Violation vector
     * @return 2-norm of KKT violation
     */
    Scalar kkt_violation_norm2( Vector const & x, Vector const & v ) const
    {
      if ( !m_active ) return 0;

      // Vectorized KKT violation computation
      auto xa = x.array();
      auto va = v.array();

      auto at_lower = ( xa <= m_lower.array() + m_tol_lower.array() );
      auto at_upper = ( xa >= m_upper.array() - m_tol_upper.array() );
      auto free_var = !( at_lower || at_upper );

      // Compute violations for different constraint types
      auto free_viol  = free_var.select( va.abs2(), 0 );
      auto lower_viol = at_lower.select( va.min( 0 ).abs2(), 0 );
      auto upper_viol = at_upper.select( va.max( 0 ).abs2(), 0 );

      return ( free_viol + lower_viol + upper_viol ).sum();
    }

    /**
     * @brief Total KKT error
     *
     * Computes maximum of stationarity, complementarity, and feasibility errors
     *
     * @param x Current point
     * @param g Gradient ∇f(x)
     * @return Maximum KKT error component
     */
    Scalar total_kkt_error( Vector const & x, Vector const & g ) const
    {
      Scalar stationarity    = kkt_stationarity_error_inf( g );
      Scalar complementarity = complementarity_error( x );

      // Vectorized feasibility error computation
      auto   lower_viol  = ( m_lower - x ).cwiseMax( 0 );
      auto   upper_viol  = ( x - m_upper ).cwiseMax( 0 );
      Scalar feasibility = std::max(
        lower_viol.template lpNorm<Eigen::Infinity>(),
        upper_viol.template lpNorm<Eigen::Infinity>() );

      return std::max( { stationarity, complementarity, feasibility } );
    }

    // =========================================================================
    //  Line Search Optimization
    // =========================================================================

    /**
     * @brief Maximum feasible step along direction
     *
     * Computes α_max = max{α : ℓ ≤ x + αd ≤ u}
     * Vectorized computation of minimum ratio
     *
     * @param x Current point
     * @param d Search direction
     * @param safety_margin Safety factor (default: 0.99)
     * @return Maximum feasible step
     */
    Scalar max_feasible_step( Vector const & x, Vector const & d, Scalar safety_margin = 0.99 ) const
    {
      if ( !m_active ) return std::numeric_limits<Scalar>::infinity();

      // Vectorized step computation using masks
      auto pos_mask = ( d.array() > m_epsi );
      auto neg_mask = ( d.array() < -m_epsi );

      Vector ratios = Vector::Constant( x.size(), std::numeric_limits<Scalar>::infinity() );

      // Compute upper bound ratios for positive directions
      ratios.array() = pos_mask.select( ( m_upper - x ).array() / d.array(), ratios.array() );

      // Compute lower bound ratios for negative directions
      ratios.array() = neg_mask.select( ( m_lower - x ).array() / d.array(), ratios.array() );

      // Find minimum positive ratio
      Scalar alpha_max = ratios.minCoeff();

      return safety_margin * ( alpha_max < 0 ? 0 : alpha_max );
    }

    /**
     * @brief Optimal projection step for steepest descent
     *
     * Computes step for projected gradient method
     *
     * @param x Current point
     * @param g Gradient ∇f(x)
     * @param max_step Maximum step length
     * @return Optimal projection step
     */
    Scalar optimal_projection_step( Vector const & x, Vector const & g, Scalar max_step = 1.0 ) const
    {
      if ( !m_active ) return max_step;

      Vector d         = -projected_gradient( x, g );
      Scalar alpha_max = max_feasible_step( x, d, 1.0 );

      return std::min( max_step, alpha_max );
    }

    // =========================================================================
    //  QP Subproblem Solving
    // =========================================================================

    /**
     * @brief Solve QP subproblem with sparse Hessian
     *
     * Solves the box-constrained QP:
     * \f[
     * \min_{p} \frac{1}{2} p^T (H + \lambda I) p + g^T p
     * \quad \text{s.t.} \quad \ell - x \leq p \leq u - x
     * \f]
     *
     * Uses:
     * - Sparse LDLT decomposition for large sparse problems
     * - Dense LDLT for smaller or denser problems
     * - Efficient projection with incremental Hessian-vector product updates
     *
     * @param x Current point
     * @param g Gradient ∇f(x)
     * @param H Hessian ∇²f(x) (sparse)
     * @param lambda Regularization parameter
     * @param verbosity Verbosity level (0: silent, ≥3: detailed)
     * @return QPSolution structure
     */
    bool solve_qp_subproblem_sparse(
      Vector const &       x,
      Vector const &       g,
      SparseMatrix const & H,
      Scalar               lambda,
      Vector &             p,            ///< Search direction (step)
      Vector &             multipliers,  ///< Lagrange multipliers (size 2*n: lower then upper)
      integer              verbosity = 0 ) const
    {
      integer n = x.size();
      p.resize( n );
      p.setZero();
      multipliers.resize( 2 * n );
      multipliers.setZero();

      try
      {
        // Choose solver based on problem characteristics
        integer nnz               = H.nonZeros();
        Scalar  sparsity_ratio    = Scalar( nnz ) / ( n * n );
        bool    use_sparse_solver = ( n > 100 && sparsity_ratio < 0.1 );

        // Regularize Hessian matrix
        SparseMatrix A = H;
        A.diagonal().array() += lambda;

        if ( use_sparse_solver )
        {
          // Sparse LDLT for large sparse problems
          Eigen::SimplicialLDLT<SparseMatrix> solver;
          solver.compute( A );

          if ( solver.info() != Eigen::Success )
          {
            if ( verbosity >= 3 ) fmt::print( "  QP: Sparse decomposition failed\n" );
            return false;
          }
          p = -solver.solve( g );
        }
        else
        {
          // Dense LDLT for smaller/denser problems
          Matrix              A_dense = A;
          Eigen::LDLT<Matrix> solver;
          solver.compute( A_dense );

          if ( solver.info() != Eigen::Success )
          {
            if ( verbosity >= 3 ) fmt::print( "  QP: Dense decomposition failed\n" );
            return false;
          }
          p = -solver.solve( g );
        }

        // Check for numerical issues
        if ( !p.allFinite() )
        {
          if ( verbosity >= 3 ) fmt::print( "  QP: Non-finite solution\n" );
          return false;
        }

        // Project step to satisfy box constraints
        if ( m_active )
        {
          // Compute H*p once
          Vector Hp = H * p;

          // Vectorized clamping with efficient Hp updates
          for ( integer i = 0; i < n; ++i )
          {
            Scalar x_new = x( i ) + p( i );

            if ( x_new < m_lower( i ) )
            {
              // Clamp to lower bound
              Scalar delta = m_lower( i ) - x( i ) - p( i );
              p( i )       = m_lower( i ) - x( i );

              if ( sparsity_ratio > 0.01 )
              {
                // Incremental update of Hp using sparse column access
                for ( typename SparseMatrix::InnerIterator it( H, i ); it; ++it )
                {
                  integer row = it.row();
                  Scalar  val = it.value();

                  if ( row == i )
                    Hp( i ) += val * delta;  // Diagonal
                  else
                    Hp( row ) += val * delta;  // Off-diagonal
                }
              }
              else
              {
                // Recompute for very sparse matrices
                Hp = H * p;
              }

              multipliers( i ) = std::max<Scalar>( 0, -g( i ) - Hp( i ) );
            }
            else if ( x_new > m_upper( i ) )
            {
              // Clamp to upper bound
              Scalar delta = m_upper( i ) - x( i ) - p( i );
              p( i )       = m_upper( i ) - x( i );

              if ( sparsity_ratio > 0.01 )
              {
                for ( typename SparseMatrix::InnerIterator it( H, i ); it; ++it )
                {
                  integer row = it.row();
                  Scalar  val = it.value();

                  if ( row == i )
                    Hp( i ) += val * delta;
                  else
                    Hp( row ) += val * delta;
                }
              }
              else
              {
                Hp = H * p;
              }

              multipliers( i + n ) = std::max<Scalar>( 0, g( i ) + Hp( i ) );
            }
          }
        }
        if ( verbosity >= 3 ) { fmt::print( "  QP: ‖p‖={:.2e}, λ={:.2e}\n", p.norm(), lambda ); }
      }
      catch ( const std::exception & e )
      {
        if ( verbosity >= 3 ) fmt::print( "  QP: Exception: {}\n", e.what() );
        return false;
      }

      return true;
    }

    /**
     * @brief Solve QP subproblem with dense Hessian
     *
     * Optimized version for dense Hessian matrices
     * Uses dense LDLT decomposition and efficient matrix-vector operations
     *
     * @param x Current point
     * @param g Gradient ∇f(x)
     * @param H Hessian ∇²f(x) (dense)
     * @param lambda Regularization parameter
     * @param verbosity Verbosity level
     * @return QPSolution structure
     */
    bool solve_qp_subproblem_dense(
      Vector const & x,
      Vector const & g,
      Matrix const & H,
      Scalar         lambda,
      Vector &       p,            ///< Search direction (step)
      Vector &       multipliers,  ///< Lagrange multipliers (size 2*n: lower then upper)
      integer        verbosity = 0 ) const
    {
      integer n = x.size();
      p.resize( n );
      p.setZero();
      multipliers.resize( 2 * n );
      multipliers.setZero();

      try
      {
        // Regularize Hessian
        Matrix A = H;
        A.diagonal().array() += lambda;

        // Dense LDLT decomposition
        Eigen::LDLT<Matrix> solver;
        solver.compute( A );

        if ( solver.info() != Eigen::Success )
        {
          if ( verbosity >= 3 ) fmt::print( "  QP: Dense decomposition failed\n" );
          return false;
        }

        p = -solver.solve( g );

        if ( !p.allFinite() )
        {
          if ( verbosity >= 3 ) fmt::print( "  QP: Non-finite solution\n" );
          return false;
        }

        // Project step with efficient Hp updates
        if ( m_active )
        {
          Vector Hp = H * p;

          // Vectorized constraint checking with manual loop for Hp updates
          for ( integer i = 0; i < n; ++i )
          {
            Scalar x_new = x( i ) + p( i );

            if ( x_new < m_lower( i ) )
            {
              // Efficient Hp update using matrix column
              Scalar delta = m_lower( i ) - x( i ) - p( i );
              p( i )       = m_lower( i ) - x( i );
              Hp.noalias() += H.col( i ) * delta;

              multipliers( i ) = std::max<Scalar>( 0, -g( i ) - Hp( i ) );
            }
            else if ( x_new > m_upper( i ) )
            {
              Scalar delta = m_upper( i ) - x( i ) - p( i );
              p( i )       = m_upper( i ) - x( i );
              Hp.noalias() += H.col( i ) * delta;

              multipliers( i + n ) = std::max<Scalar>( 0, g( i ) + Hp( i ) );
            }
          }
        }

        if ( verbosity >= 3 ) { fmt::print( "  QP: ‖p‖={:.2e}, λ={:.2e}\n", p.norm(), lambda ); }
      }
      catch ( const std::exception & e )
      {
        if ( verbosity >= 3 ) fmt::print( "  QP: Exception: {}\n", e.what() );
        return false;
      }

      return true;
    }

    /**
     * @brief Template QP solver with automatic Hessian type dispatch
     *
     * Compile-time dispatch to appropriate implementation
     * Supports both sparse and dense Hessian matrices
     *
     * @tparam HessianType Type of Hessian matrix (SparseMatrix or Matrix)
     * @param x Current point
     * @param g Gradient ∇f(x)
     * @param H Hessian ∇²f(x)
     * @param lambda Regularization parameter
     * @param verbosity Verbosity level
     * @return QPSolution structure
     */
    template <typename HessianType> bool solve_qp_subproblem(
      Vector const &      x,
      Vector const &      g,
      HessianType const & H,
      Scalar              lambda,
      Vector &            p,            ///< Search direction (step)
      Vector &            multipliers,  ///< Lagrange multipliers (size 2*n: lower then upper)
      integer             verbosity = 0 ) const
    {
      // Compile-time dispatch based on Hessian type
      if constexpr ( std::is_same_v<HessianType, SparseMatrix> )
      {
        return solve_qp_subproblem_sparse( x, g, H, lambda, p, multipliers, verbosity );
      }
      else if constexpr ( std::is_same_v<HessianType, Matrix> )
      {
        return solve_qp_subproblem_dense( x, g, H, lambda, p, multipliers, verbosity );
      }
      else
      {
        // Convert to dense matrix for unknown types
        Matrix H_dense = H;
        return solve_qp_subproblem_dense( x, g, H_dense, lambda, p, multipliers, verbosity );
      }
    }

    // =========================================================================
    //  Debugging and Reporting
    // =========================================================================

    /**
     * @brief Print constraint bounds
     *
     * Displays first 10 bounds for large problems
     */
    void print_bounds() const
    {
      if ( !m_active )
      {
        fmt::print( "No box constraints active\n" );
        return;
      }

      integer n = dimension();
      fmt::print( "Box constraints (dimension={}):\n", n );
      fmt::print( "{:>5} {:>15} {:>15} {:>15}\n", "i", "lower", "upper", "active" );

      integer display_count = std::min( n, integer( 10 ) );
      for ( integer i = 0; i < display_count; ++i )
      {
        fmt::print( "{:5d} {:15.6e} {:15.6e}\n", i, m_lower( i ), m_upper( i ) );
      }

      if ( n > 10 ) fmt::print( "... (showing first 10 of {})\n", n );
    }

    /**
     * @brief Get constraint summary
     * @return String with basic constraint information
     */
    std::string summary() const
    {
      if ( !m_active ) return "No box constraints";

      integer n      = dimension();
      Scalar  volume = ( m_upper - m_lower ).prod();

      return fmt::format( "Box[dim={}, volume={:.3e}, active={}]", n, volume, m_active ? "yes" : "no" );
    }

    /**
     * @brief Get detailed constraint information
     * @return String with comprehensive constraint statistics
     */
    std::string detailed_info() const
    {
      if ( !m_active ) return "No box constraints active";

      integer n         = dimension();
      Scalar  min_lower = m_lower.minCoeff();
      Scalar  max_lower = m_lower.maxCoeff();
      Scalar  min_upper = m_upper.minCoeff();
      Scalar  max_upper = m_upper.maxCoeff();
      Scalar  avg_width = ( m_upper - m_lower ).mean();
      Scalar  max_width = ( m_upper - m_lower ).maxCoeff();
      Scalar  min_width = ( m_upper - m_lower ).minCoeff();

      return fmt::format(
        "Box Constraints Statistics:\n"
        "  Dimension: {}\n"
        "  Lower bounds: [{:.3e}, {:.3e}]\n"
        "  Upper bounds: [{:.3e}, {:.3e}]\n"
        "  Width range: [{:.3e}, {:.3e}]\n"
        "  Average width: {:.3e}",
        n,
        min_lower,
        max_lower,
        min_upper,
        max_upper,
        min_width,
        max_width,
        avg_width );
    }
  };

}  // namespace Utils

#endif
