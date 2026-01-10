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
 |      Università degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: Utils_TrustRegion.hh
//

#pragma once

#ifndef UTILS_THRUST_REGION_dot_HH
#define UTILS_THRUST_REGION_dot_HH

#include "Utils_fmt.hh"
#include "Utils_ssolver.hh"

namespace Utils
{

  // ===========================================================================
  // Trust Region Subproblem Solver
  // ===========================================================================

  template <typename Scalar = double> class TrustRegionSolver
  {
  public:
    using Vector       = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using integer      = Eigen::Index;
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;
    using DenseMatrix  = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    enum class Method
    {
      DOGLEG,         // Powell's dogleg method
      STEIHAUG_CG,    // Steihaug's conjugate gradient
      MORE_SORENSEN,  // More-Sorensen with exact solution
      REGULARIZED     // Simple regularization (fallback)
    };

    struct Options
    {
      Method  method           = Method::MORE_SORENSEN;
      Scalar  cg_tol           = 0.1;  // Tolerance for CG termination
      integer cg_max_iter      = 20;   // Max CG iterations
      Scalar  tol              = 1e-9;
      Scalar  ms_tol           = 1e-8;  // Tolerance for More-Sorensen Newton
      integer ms_max_iter      = 10;    // Max Newton iterations for More-Sorensen
      bool    use_precondition = true;  // Use preconditioning in CG
      Scalar  min_radius_ratio = 0.1;   // Minimum p_norm/Δ for dogleg
      Scalar  max_radius_ratio = 0.8;   // Maximum p_norm/Δ for dogleg
    };

  private:
    Options m_opts;

    /**
     * @brief Compute Cauchy point (steepest descent within trust region)
     *
     * Solves: min α s.t. ‖αg‖ ≤ Δ
     * Solution: α = min(Δ/‖g‖, gᵀg/gᵀHg)
     */
    Vector compute_cauchy_point( Vector const & g, SparseMatrix const & H, Scalar delta ) const
    {
      Scalar g_norm = g.norm();
      if ( g_norm < std::numeric_limits<Scalar>::epsilon() ) return Vector::Zero( g.size() );

      // Maximum step to trust region boundary
      Scalar alpha_max = delta / g_norm;

      // Optimal step for quadratic model along gradient
      Scalar gHg = g.dot( H * g );
      Scalar alpha_opt;

      if ( gHg > 0 ) { alpha_opt = ( g.dot( g ) / gHg ); }
      else
      {
        // Negative curvature: take maximum step
        alpha_opt = alpha_max;
      }

      // Take minimum
      Scalar alpha = std::min( alpha_max, alpha_opt );
      return -alpha * g;
    }

    /**
     * @brief Dogleg method (Powell 1970)
     *
     * Path: p(τ) = { τ·p_u,           0 ≤ τ ≤ 1
     *               p_u + (τ-1)(p_b - p_u), 1 ≤ τ ≤ 2 }
     * where p_u = - (gᵀg/gᵀHg)·g (Cauchy point)
     *       p_b = -H⁻¹g (Newton point)
     */
    bool solve_dogleg( Vector const & g, SparseMatrix const & H, Scalar delta, Vector & p ) const
    {
      try
      {
        // Handle zero gradient
        if ( g.norm() < m_opts.tol )
        {
          p.setZero( g.size() );
          return true;
        }

        // Compute Newton point
        SymmetricSolver<Scalar> solver( H, 0 );
        Vector                  p_b = solver.solve( -g );

        Scalar p_b_norm = p_b.norm();

        // If Newton point inside trust region, accept it
        if ( p_b_norm <= delta )
        {
          p = p_b;
          return true;
        }

        // Compute Cauchy point
        Vector p_u      = compute_cauchy_point( g, H, delta );
        Scalar p_u_norm = p_u.norm();

        // If Cauchy point on boundary, return it
        if ( std::abs( p_u_norm - delta ) < m_opts.tol * delta )
        {
          p = p_u;
          return true;
        }

        // Find τ where ‖p(τ)‖ = Δ
        Vector d = p_b - p_u;
        Scalar a = d.dot( d );
        Scalar b = 2.0 * p_u.dot( d );
        Scalar c = p_u.dot( p_u ) - delta * delta;

        // Solve quadratic equation aτ² + bτ + c = 0
        Scalar disc = b * b - 4.0 * a * c;
        if ( disc < 0 )
        {
          // Fallback to Cauchy point scaled to boundary
          p = p_u * ( delta / p_u_norm );
          return true;
        }

        Scalar sqrt_disc = std::sqrt( disc );
        Scalar tau1      = ( -b + sqrt_disc ) / ( 2.0 * a );
        Scalar tau2      = ( -b - sqrt_disc ) / ( 2.0 * a );

        // Choose τ in [1,2]
        Scalar tau = 1.0;
        if ( tau1 >= 1.0 && tau1 <= 2.0 )
          tau = tau1;
        else if ( tau2 >= 1.0 && tau2 <= 2.0 )
          tau = tau2;
        else
        {
          // Fallback: take point on boundary along dogleg path
          Vector p1 = p_u;
          Vector p2 = p_b;
          Scalar f1 = p1.norm() - delta;
          Scalar f2 = p2.norm() - delta;

          if ( f1 * f2 > 0 )
          {
            // Same sign, use interpolation
            tau = 1.5;
          }
          else
          {
            // Bisection
            for ( int i = 0; i < 10; ++i )
            {
              tau          = ( tau1 + tau2 ) / 2.0;
              Vector p_tau = p_u + ( tau - 1.0 ) * d;
              Scalar f_tau = p_tau.norm() - delta;

              if ( std::abs( f_tau ) < 1e-6 * delta ) break;

              if ( f1 * f_tau < 0 )
              {
                tau2 = tau;
                f2   = f_tau;
              }
              else
              {
                tau1 = tau;
                f1   = f_tau;
              }
            }
          }
        }

        p = p_u + ( tau - 1.0 ) * d;
        return true;
      }
      catch ( ... )
      {
        // Fallback to steepest descent
        p             = compute_cauchy_point( g, H, delta );
        Scalar p_norm = p.norm();
        if ( p_norm > delta ) p *= delta / p_norm;
        return true;
      }
    }

    /**
     * @brief Steihaug's conjugate gradient method
     *
     * Solves trust-region subproblem approximately using CG.
     * Terminates when:
     * 1. ‖p‖ ≥ Δ (hit boundary)
     * 2. Negative curvature detected
     * 3. Residual small enough
     */
    bool solve_steihaug_cg( Vector const & g, SparseMatrix const & H, Scalar delta, Vector & p ) const
    {
      p.setZero( g.size() );

      // Handle zero gradient case
      Scalar g_norm = g.norm();
      if ( g_norm < m_opts.tol )
      {
        return true;  // p = 0 is optimal
      }

      Vector r = -g;  // Initial residual: r = -g - H·p (p=0)
      Vector d = r;   // Initial search direction

      Scalar rTr = r.dot( r );
      Scalar rTr_old;
      Scalar alpha, beta;

      for ( integer k = 0; k < m_opts.cg_max_iter; ++k )
      {
        // Compute Hd
        Vector Hd = H * d;

        // Check for negative curvature
        Scalar dHd = d.dot( Hd );
        if ( dHd <= 0 )
        {
          // Negative curvature: move to boundary along d
          Scalar p_norm  = p.norm();
          Scalar p_dot_d = p.dot( d );
          Scalar d_norm2 = d.dot( d );

          if ( d_norm2 < m_opts.tol ) break;  // Avoid division by zero

          Scalar a = d_norm2;
          Scalar b = 2.0 * p_dot_d;
          Scalar c = p_norm * p_norm - delta * delta;

          Scalar disc = b * b - 4.0 * a * c;
          if ( disc < 0 ) break;

          Scalar sqrt_disc = std::sqrt( disc );
          Scalar tau1      = ( -b + sqrt_disc ) / ( 2.0 * a );
          Scalar tau2      = ( -b - sqrt_disc ) / ( 2.0 * a );

          Scalar tau = std::max( tau1, tau2 );
          if ( tau > 0 && tau < 1e6 )
          {  // Reasonable tau
            p += tau * d;
            return true;
          }
          break;
        }

        // Standard CG step
        alpha = rTr / dHd;

        // Check if step would exceed trust region
        Vector p_new      = p + alpha * d;
        Scalar p_new_norm = p_new.norm();

        if ( p_new_norm >= delta )
        {
          // Move to boundary
          Scalar p_norm  = p.norm();
          Scalar p_dot_d = p.dot( d );
          Scalar d_norm2 = d.dot( d );

          if ( d_norm2 < m_opts.tol ) break;

          Scalar a = d_norm2;
          Scalar b = 2.0 * p_dot_d;
          Scalar c = p_norm * p_norm - delta * delta;

          Scalar disc = b * b - 4.0 * a * c;
          if ( disc < 0 ) break;

          Scalar sqrt_disc = std::sqrt( disc );
          Scalar tau1      = ( -b + sqrt_disc ) / ( 2.0 * a );
          Scalar tau2      = ( -b - sqrt_disc ) / ( 2.0 * a );

          Scalar tau = std::max( tau1, tau2 );
          if ( tau > 0 && tau < 1e6 )
          {
            p += tau * d;
            return true;
          }
          break;
        }

        // Accept step
        p = p_new;
        r -= alpha * Hd;

        // Check convergence
        rTr_old = rTr;
        rTr     = r.dot( r );

        if ( std::sqrt( rTr ) < m_opts.cg_tol * g_norm ) return true;

        // Update search direction
        beta = rTr / rTr_old;
        d    = r + beta * d;
      }

      // If we get here, scale to trust region boundary if necessary
      Scalar p_norm = p.norm();
      if ( p_norm > delta ) { p *= delta / p_norm; }
      return p_norm > 0;
    }

    /**
     * @brief More-Sorensen method for exact trust-region solution
     *
     * Solves (H + λI)p = -g with ‖p‖ = Δ using Newton's method on:
     * φ(λ) = 1/‖p(λ)‖ - 1/Δ = 0
     */
    bool solve_more_sorensen( Vector const & g, SparseMatrix const & H, Scalar delta, Vector & p, Scalar & lambda )
      const
    {
      // First try with lambda = 0 (Newton point)
      try
      {
        SymmetricSolver<Scalar> solver0( H, 0 );
        Vector                  p_newton    = solver0.solve( -g );
        Scalar                  newton_norm = p_newton.norm();

        if ( newton_norm <= delta + m_opts.tol * delta )
        {
          // Newton point is inside trust region
          p      = p_newton;
          lambda = 0;
          return true;
        }
      }
      catch ( ... )
      {
        // Newton point failed, proceed with iterative method
      }

      // If gradient is zero, solution is p = 0
      Scalar g_norm = g.norm();
      if ( g_norm < m_opts.tol )
      {
        p.setZero( g.size() );
        lambda = 0;
        return true;
      }

      // Compute initial bounds for λ
      Scalar lambda_low = 0;
      Scalar lambda_up  = 0;

      // Try to estimate bounds using eigenvalue computation
      try
      {
        // Create a dense copy for eigenvalue computation
        DenseMatrix                                H_dense = H;
        Eigen::SelfAdjointEigenSolver<DenseMatrix> eigensolver( H_dense );
        if ( eigensolver.info() == Eigen::Success )
        {
          auto   eigenvals  = eigensolver.eigenvalues();
          Scalar lambda_min = eigenvals.minCoeff();
          lambda_low        = std::max( 0.0, -lambda_min );

          // Estimate upper bound: lambda such that (H + λI) is sufficiently positive definite
          lambda_up = lambda_low + 100 * std::max( 1.0, std::abs( lambda_min ) );

          // If gradient is large, we might need larger lambda to reduce step size
          if ( g_norm > 0 ) { lambda_up = std::max( lambda_up, g_norm / delta - lambda_min ); }
        }
        else
        {
          lambda_low = 0;
          lambda_up  = 1e8;
        }
      }
      catch ( ... )
      {
        lambda_low = 0;
        lambda_up  = 1e8;
      }

      // Initialize lambda with a reasonable value
      lambda = lambda_low;

      // Use a bracketing approach for robustness
      Scalar left_lambda  = lambda_low;
      Scalar right_lambda = lambda_up;

      // Try to find good initial bracket
      for ( int i = 0; i < 10; ++i )
      {
        try
        {
          SymmetricSolver<Scalar> solver( H, lambda );
          p             = solver.solve( -g );
          Scalar p_norm = p.norm();

          if ( p_norm > delta )
          {
            left_lambda = lambda;
            if ( right_lambda <= lambda ) { right_lambda = lambda * 2.0; }
          }
          else
          {
            right_lambda = lambda;
            if ( left_lambda >= lambda ) { left_lambda = std::max( 0.0, lambda / 2.0 ); }
          }

          if ( left_lambda < right_lambda ) break;
        }
        catch ( ... )
        {
          // Increase lambda and try again
          lambda *= 2.0;
        }
      }

      lambda = ( left_lambda + right_lambda ) / 2.0;

      for ( integer iter = 0; iter < m_opts.ms_max_iter; ++iter )
      {
        try
        {
          // Solve (H + λI)p = -g
          SymmetricSolver<Scalar> solver( H, lambda );
          p = solver.solve( -g );

          Scalar p_norm   = p.norm();
          Scalar diff     = p_norm - delta;
          Scalar rel_diff = std::abs( diff ) / ( delta + m_opts.tol );

          if ( rel_diff < m_opts.ms_tol )
          {
            // Converged
            return true;
          }

          // Update brackets based on whether p_norm is too large or too small
          if ( p_norm > delta )
          {
            // Need to increase lambda to reduce step size
            left_lambda = lambda;
          }
          else
          {
            // Can decrease lambda to get larger step
            right_lambda = lambda;
          }

          // Try Newton step if we have gradient information
          if ( p_norm > 0 )
          {
            try
            {
              Vector q       = solver.solve( -p );
              Scalar p_dot_q = p.dot( q );
              if ( std::abs( p_dot_q ) > m_opts.tol )
              {
                Scalar phi   = 1.0 / p_norm - 1.0 / delta;
                Scalar deriv = -p_dot_q / ( p_norm * p_norm * p_norm );

                if ( std::abs( deriv ) > m_opts.tol )
                {
                  Scalar lambda_new = lambda - phi / deriv;

                  // Safeguard the Newton step
                  if ( lambda_new > left_lambda && lambda_new < right_lambda ) { lambda = lambda_new; }
                  else
                  {
                    lambda = ( left_lambda + right_lambda ) / 2.0;
                  }
                }
                else
                {
                  lambda = ( left_lambda + right_lambda ) / 2.0;
                }
              }
              else
              {
                lambda = ( left_lambda + right_lambda ) / 2.0;
              }
            }
            catch ( ... )
            {
              lambda = ( left_lambda + right_lambda ) / 2.0;
            }
          }
          else
          {
            lambda = ( left_lambda + right_lambda ) / 2.0;
          }

          // Check bracket width
          if ( right_lambda - left_lambda < m_opts.ms_tol * ( 1.0 + lambda ) ) { return true; }
        }
        catch ( ... )
        {
          // Solver failed, increase λ and adjust bracket
          left_lambda = lambda;
          lambda      = std::min( right_lambda, lambda * 2.0 );
        }
      }

      // Fallback: use the best lambda we found and ensure constraint satisfaction
      try
      {
        SymmetricSolver<Scalar> solver( H, lambda );
        p             = solver.solve( -g );
        Scalar p_norm = p.norm();
        if ( p_norm > delta ) p *= delta / p_norm;
        return true;
      }
      catch ( ... )
      {
        // Last resort: gradient descent step
        p = -( delta / g_norm ) * g;
        return true;
      }
    }

  public:
    TrustRegionSolver( Options opts = Options() ) : m_opts( opts ) {}

    /**
     * @brief Solve trust-region subproblem: min ½pᵀHp + gᵀp s.t. ‖p‖ ≤ Δ
     *
     * @param g Gradient
     * @param H Hessian
     * @param delta Trust region radius
     * @param p [out] Solution
     * @param lambda [out] Regularization parameter used
     * @return true if successful
     */
    bool solve( Vector const & g, SparseMatrix const & H, Scalar delta, Vector & p, Scalar & lambda )
    {
      if ( delta <= 0 ) return false;

      switch ( m_opts.method )
      {
        case Method::DOGLEG: return solve_dogleg( g, H, delta, p );
        case Method::STEIHAUG_CG: return solve_steihaug_cg( g, H, delta, p );
        case Method::MORE_SORENSEN: return solve_more_sorensen( g, H, delta, p, lambda );
        case Method::REGULARIZED:
        default:
          try
          {
            SymmetricSolver<Scalar> solver( H, lambda );
            p             = solver.solve( -g );
            Scalar p_norm = p.norm();
            if ( p_norm > delta ) p *= delta / p_norm;
            return true;
          }
          catch ( ... )
          {
            return false;
          }
      }
    }

    void set_method( Method m ) { m_opts.method = m; }
  };

}  // namespace Utils

#endif
