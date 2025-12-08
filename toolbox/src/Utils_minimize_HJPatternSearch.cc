/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2022                                                      |
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
// file: Utils_minimize_HJPatternSearch.cc
//

#include "Utils_minimize_HJPatternSearch.hh"

namespace Utils
{

  // =================================================================
  // Setup
  // =================================================================

  template <typename Real>
  void
  HJPatternSearch<Real>::set_tolerance( Real tol )
  {
    UTILS_ASSERT( tol > 0, "set_tolerance({}) argument must be >0\n", tol );
    m_tolerance = tol;
  }

  // =================================================================

  template <typename Real>
  void
  HJPatternSearch<Real>::set_max_iterations( integer mit )
  {
    UTILS_ASSERT( mit > 0, "set_max_iterations({}) argument must be >0\n", mit );
    m_max_iteration = mit;
  }

  // =================================================================

  template <typename Real>
  void
  HJPatternSearch<Real>::set_max_fun_evaluation( integer mfev )
  {
    UTILS_ASSERT( mfev > 0, "set_max_fun_evaluation({}) argument must be >0\n", mfev );
    m_max_fun_evaluation = mfev;
  }

  // =================================================================

  template <typename Real>
  void
  HJPatternSearch<Real>::set_max_num_stagnation( integer nstg )
  {
    UTILS_ASSERT( nstg > 0, "set_max_num_stagnation({}) argument must be >0\n", nstg );
    m_max_num_stagnation = nstg;
  }

  // =================================================================
  // Allocate
  // =================================================================
  template <typename Real>
  void
  HJPatternSearch<Real>::allocate( integer n )
  {
    m_dim        = n;
    integer ntot = n * 10;
    m_base_value.reallocate( ntot );

    new ( &this->m_dir ) MapVec( m_base_value( n ), n );
    new ( &this->m_new_x ) MapVec( m_base_value( n ), n );
    new ( &this->m_search_sign ) MapVec( m_base_value( n ), n );
    new ( &this->m_x_old ) MapVec( m_base_value( n ), n );
    new ( &this->m_x_best ) MapVec( m_base_value( n ), n );
    new ( &this->m_p ) MapVec( m_base_value( n ), n );
    new ( &this->m_p1 ) MapVec( m_base_value( n ), n );
  }

  // =================================================================
  // Setup
  // =================================================================
  template <typename Real>
  void
  HJPatternSearch<Real>::setup( integer const dim, HJFunc & fun, Console const * console )
  {
    // HJPatternSearch The constructor initialize the solver
    // parameters and check the inputs when the class is instanciated.
    m_fun     = fun;
    m_console = console;
    allocate( dim );
  }

  // =================================================================
  // info
  // =================================================================

  template <typename Real>
  string
  HJPatternSearch<Real>::info() const
  {
    string res = fmt::format( "\nOptimization stats:\n" );
    if ( m_iteration_count >= m_max_iteration )
    {
      res += fmt::format( "iteration number reached  Max Iteration Limit = {}\n", m_max_iteration );
    }
    else if ( m_fun_evaluation >= m_max_fun_evaluation )
    {
      res += fmt::format( "function evaluations {} exceeded the maximum limit [{}]\n", m_fun_evaluation,
                          m_max_fun_evaluation );
    }
    else
    {
      res += fmt::format( "mesh size h = {} less than tolerance = {}\n", m_h, m_tolerance );
    }
    res += fmt::format( "[#iterations/#feval] F(x_best): [{}/{}]  {:6}\n\n", m_iteration_count, m_fun_evaluation,
                        m_f_best );
    return res;
  }

  // =================================================================
  // Search
  // =================================================================

  template <typename Real>
  void
  HJPatternSearch<Real>::search()
  {
    /*
    // SEARCH This method call the explore method on the first
    // iteration and then continue to call explore until a stencil
    // fails. In the case of a stencil failure, it tries once to go
    // back of half a step along the search direction by setting x_center
    // equal to the base point x_best.
    // If the stencil fails again, it exits the while loop and stencil_failure
    // is set to zero in order to signal that a reduction of h is necessary.
    */
    ++m_iteration_count;  // augment counter
    // Print info
    if ( m_verbose > 0 && m_console != nullptr && m_console->get_level() >= 3 )
    {
      string line =
          "--------------------------------------------------------------------"
          "-----";
      string msg = fmt::format( "Iteration={} f(x_best)/#f/|h| = {:.6} / {} / {:.6}\n", m_iteration_count, m_f_best,
                                m_fun_evaluation, m_h );
      if ( m_verbose > 2 )
      {
        for ( integer ii = 0; ii < m_dim; ++ii ) msg += fmt::format( "x[{}] = {:.6}\n", ii, m_x_best( ii ) );
      }
      msg += fmt::format( "{}\n", line );
      m_console->message( msg, 3 );
    }

    best_nearby();

    while ( m_stencil_failure )
    {
      // reduce the scale
      m_h *= m_rho;
      best_nearby();
      if ( m_h <= m_tolerance ) return;
      if ( m_fun_evaluation_count >= m_max_fun_evaluation ) return;
    }

    m_dir.noalias() = m_x_best - m_x_old;  // Compute search direction

    // Continue exploring until stencil failure or exceed of
    Real lambda  = 1;
    Real max_der = 0;
    while ( m_fun_evaluation_count < m_max_fun_evaluation && lambda > 0.1 )
    {
      m_new_x.noalias() = m_x_best + lambda * m_dir;
      Real new_f        = eval_function( m_new_x );
      if ( new_f < m_f_best - ( 0.25 * lambda ) * max_der )
      {
        Real der = ( m_f_best - new_f ) / lambda;
        if ( der > max_der ) max_der = der;
        m_x_best.noalias() = m_new_x;
        m_f_best           = new_f;
        // lambda *= 2;
      }
      lambda /= 2;
    }
  }

  // =================================================================
  // Explore
  // =================================================================

  template <typename Real>
  void
  HJPatternSearch<Real>::best_nearby()
  {
    /*
    // EXPLORE This method explore all points on the stencil center at
    // x_temporary = x_center and updates the current iteration x to the current
    // best point x_current_best. If the current best point x_current_best is
    worse than the
    // base point x_best, the current iteration x will remain constant
    // (x = x_best) and stencil failure flag stencil_failure will be set to
    zero.
    */
    // Initialize
    m_stencil_failure = true;

    // ----------------------------------------------------------------------------------------
    // Cycle on all stencil directions

    for ( integer j = 0; j < m_dim; ++j )
    {
      Real s_dirh   = m_search_sign( j ) * m_h;
      m_p.noalias() = m_x_best;
      m_p( j ) += s_dirh;
      Real fp = eval_function( m_p );
      if ( fp >= m_f_best )
      {
        m_p1.noalias() = m_x_best;
        m_p1( j ) -= s_dirh;  // try the opposite direction
        Real fp1 = eval_function( m_p1 );
        if ( fp1 < fp )
        {
          m_p.noalias() = m_p1;
          fp            = fp1;
          // change priority of search direction to the opposite verse
          m_search_sign( j ) = -m_search_sign( j );
        }
      }
      // Update temporary and current best point before checking
      // the remaining directions j
      if ( fp < m_f_best )
      {
        m_x_best.noalias() = m_p;    // move temporary point
        m_f_best           = fp;     // new current best point
        m_stencil_failure  = false;  // update stencil failure flag
      }
    }
  }

  // =================================================================
  // Run
  // =================================================================

  template <typename Real>
  void
  HJPatternSearch<Real>::run( Real const x[], Real h )
  {
    // RUN This method run the whole Hooke-Jeeves algorithm. Search
    // is repeated until it fails, then the scal h is reduced. When h
    // is less than a threshold, the method returns the solution.

    // initialize current iteration to guess for the first iteration
    std::copy_n( x, m_dim, m_x_best.data() );
    m_f_best          = eval_function( m_x_best );
    m_stencil_failure = false;  // initialize stencil failure flag
    m_search_sign.setOnes();    // Initialize search verse vector to all ones (first verse
                                // will be positive for each direction)
    m_h = h;

    m_iteration_count      = 0;  // initialize explore iteration counter
    m_fun_evaluation_count = 0;  // initialize function evaluation counter

    // RUN This method run the whole Hooke-Jeeves algorithm.
    // Search is repeated until it fails, then the scal h is reduced.
    // When h is less than a threshold, the method returns the solution.

    integer num_stagnation = 0;
    while ( m_h > m_tolerance && m_fun_evaluation_count < m_max_fun_evaluation )
    {
      m_x_old.noalias() = m_x_best;
      m_f_old           = m_f_best;

      search();
      if ( m_stencil_failure ) break;

      // If iteration limit is reached,stop.
      if ( m_iteration_count >= m_max_iteration ) break;

      // reduce the scale
      m_h *= m_rho;

      // check stagnation
      if ( m_f_old <= m_f_best )
      {
        ++num_stagnation;
        if ( num_stagnation > m_max_num_stagnation ) break;
      }
      else
      {
        num_stagnation = 0;
      }
    }
  }

  // =================================================================

  template class HJPatternSearch<double>;
  template class HJPatternSearch<float>;
}  // namespace Utils

//
// eof: Utils_HJPatternSearch.c
//
