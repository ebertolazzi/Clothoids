 /*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2003                                                      |
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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Utils_zeros.hh"

#include <vector>

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

namespace Utils {

  using std::isfinite;

  /*
  //   __________ ____   ___  ____
  //  |__  / ____|  _ \ / _ \/ ___|
  //    / /|  _| | |_) | | | \___ \
  //   / /_| |___|  _ <| |_| |___) |
  //  /____|_____|_| \_\\___/|____/
  */

  // =================================================================
  // set_max_iterations
  // =================================================================

  template <typename Real>
  void
  Zeros<Real>::set_max_iterations( Integer mit ) {
    UTILS_ASSERT(
      mit > 0,
      "Zeros::set_max_iterations({}) argument must be >0\n", mit
    );
    m_max_iteration = mit;
  }

  // =================================================================
  // set_max_fun_evaluation
  // =================================================================

  template <typename Real>
  void
  Zeros<Real>::set_max_fun_evaluation( Integer mfev ) {
    UTILS_ASSERT(
      mfev > 0,
      "Zeros::set_max_fun_evaluation({}) argument must be >0\n", mfev
    );
    m_max_fun_evaluation = mfev;
  }

  template <typename Real>
  void
  Zeros<Real>::set_tolerance( Real tol ) {
    UTILS_ASSERT(
      tol > 0,
      "Zeros::set_tolerance({}) argument must be >0\n", tol
    );
    m_tolerance = tol;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  Zeros<Real>::solve_Newton( Real x_guess, Zeros_base_fun<Real> * fun ) {
    m_iteration_count      = 0;
    m_fun_evaluation_count = 0;

    Real x{x_guess};
    while ( true ) {
      Real f = fun->eval( x ); ++m_fun_evaluation_count;
      m_converged = abs(f) < m_tolerance;
      if ( m_converged ) return x;
      ++m_iteration_count;
      if ( m_iteration_count > m_max_iteration || m_fun_evaluation_count > m_max_fun_evaluation ) break;

      Real df = fun->eval_D( x ); ++m_fun_evaluation_count;
      x -= f/df;
      if ( !Utils::is_finite(x) ) break;
    }
    return 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  Zeros<Real>::solve_Chebyshev( Real x_guess, Zeros_base_fun<Real> * fun ) {
    m_iteration_count      = 0;
    m_fun_evaluation_count = 0;

    Real x{x_guess};
    while ( true ) {
      Real f = fun->eval( x ); ++m_fun_evaluation_count;
      m_converged = abs(f) < m_tolerance;
      if ( m_converged ) return x;
      ++m_iteration_count;
      if ( m_iteration_count > m_max_iteration || m_fun_evaluation_count > m_max_fun_evaluation ) break;

      Real df  = fun->eval_D( x );  ++m_fun_evaluation_count;
      Real ddf = fun->eval_DD( x ); ++m_fun_evaluation_count;
      x -= (f/df)*(1+(f*ddf)/(2*df*df));
      if ( !Utils::is_finite(x) ) {
        fmt::print("found x[{}] = {}\n",m_iteration_count,x);
        break;
      }
    }
    return 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Real>
  Real
  Zeros<Real>::solve_Halley( Real x_guess, Zeros_base_fun<Real> * fun ) {
    m_iteration_count      = 0;
    m_fun_evaluation_count = 0;

    Real x{x_guess};
    while ( true ) {
      Real f = fun->eval( x ); ++m_fun_evaluation_count;
      m_converged = abs(f) < m_tolerance;
      if ( m_converged ) return x;
      ++m_iteration_count;
      if ( m_iteration_count > m_max_iteration || m_fun_evaluation_count > m_max_fun_evaluation ) break;

      Real df  = fun->eval_D( x );  ++m_fun_evaluation_count;
      Real ddf = fun->eval_DD( x ); ++m_fun_evaluation_count;
      x -= (f/df)/(1-(f*ddf)/(2*df*df));
      if ( !Utils::is_finite(x) ) {
        fmt::print("found x[{}] = {}\n",m_iteration_count,x);
        break;
      }
    }
    return 0;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /*

  An Optimal Thirty-Second-Order Iterative Method for
  Solving Nonlinear Equations and a Conjecture
  Juan Luis Varona
  Qualitative Theory of Dynamical Systems (2022)
  https://doi.org/10.1007/s12346-022-00572-3
  */

  template <typename Real>
  Real
  Zeros<Real>::solve_Order4( Real x_guess, Zeros_base_fun<Real> * fun ) {
    m_iteration_count      = 0;
    m_fun_evaluation_count = 0;

    Real x{x_guess};
    while ( true ) {
      Real f = fun->eval( x ); ++m_fun_evaluation_count;
      m_converged = abs(f) < m_tolerance;
      if ( m_converged ) return x;
      ++m_iteration_count;
      if ( m_iteration_count > m_max_iteration || m_fun_evaluation_count > m_max_fun_evaluation ) break;

      Real df = fun->eval_D( x );  ++m_fun_evaluation_count;
      Real y  = x - (f/df);
      Real fy = fun->eval( y ); ++m_fun_evaluation_count;
      Real t  = fy/f;
      x = y-Q(t)*(fy/df);
      if ( !Utils::is_finite(x) ) {
        fmt::print("found x[{}] = {}\n",m_iteration_count,x);
        break;
      }
    }
    return 0;
  }

  template <typename Real>
  Real
  Zeros<Real>::solve_Order8( Real x_guess, Zeros_base_fun<Real> * fun ) {
    m_iteration_count      = 0;
    m_fun_evaluation_count = 0;

    Real x{x_guess};
    while ( true ) {
      Real f = fun->eval( x ); ++m_fun_evaluation_count;
      m_converged = abs(f) < m_tolerance;
      if ( m_converged ) return x;
      ++m_iteration_count;
      if ( m_iteration_count > m_max_iteration || m_fun_evaluation_count > m_max_fun_evaluation ) break;

      Real df = fun->eval_D( x );  ++m_fun_evaluation_count;
      Real y  = x - (f/df);
      Real fy = fun->eval( y ); ++m_fun_evaluation_count;
      Real t  = fy/f;
      Real z  = y - Q(t)*(fy/df);
      Real fz = fun->eval( z );   ++m_fun_evaluation_count;
      Real s  = fz/fy;
      x = z-W(t,s)*(fz/df);
      if ( !Utils::is_finite(x) ) {
        fmt::print("found x[{}] = {}\n",m_iteration_count,x);
        break;
      }
    }
    return 0;
  }

  /*

  An Optimal Thirty-Second-Order Iterative Method for
  Solving Nonlinear Equations and a Conjecture
  Juan Luis Varona
  Qualitative Theory of Dynamical Systems (2022)
  https://doi.org/10.1007/s12346-022-00572-3
  */

  template <typename Real>
  Real
  Zeros<Real>::solve_Order16( Real x_guess, Zeros_base_fun<Real> * fun ) {
    m_iteration_count      = 0;
    m_fun_evaluation_count = 0;

    Real x{x_guess};
    while ( true ) {
      Real f = fun->eval( x ); ++m_fun_evaluation_count;
      m_converged = abs(f) < m_tolerance;
      if ( m_converged ) return x;
      ++m_iteration_count;
      if ( m_iteration_count > m_max_iteration || m_fun_evaluation_count > m_max_fun_evaluation ) break;

      Real df = fun->eval_D( x );  ++m_fun_evaluation_count;
      Real y  = x - (f/df);
      Real fy = fun->eval( y ); ++m_fun_evaluation_count;

      m_converged = abs(fy) < m_tolerance;
      if ( m_converged ) return y;

      Real t  = fy/f;
      Real z  = y - Q(t)*(fy/df);
      Real fz = fun->eval( z ); ++m_fun_evaluation_count;

      m_converged = abs(fz) < m_tolerance;
      if ( m_converged ) return z;

      Real s  = fz/fy;
      Real w  = z - W(t,s)*(fz/df);
      Real fw = fun->eval( w ); ++m_fun_evaluation_count;

      Real u = fw/fz;
      x = w-H(t,s,u)*(fw/df);
      if ( !Utils::is_finite(x) ) {
        fmt::print("found x[{}] = {}\n",m_iteration_count,x);
        break;
      }
    }
    return 0;
  }

  template <typename Real>
  Real
  Zeros<Real>::solve_Order32( Real x_guess, Zeros_base_fun<Real> * fun ) {
    m_iteration_count      = 0;
    m_fun_evaluation_count = 0;

    Real x{x_guess};
    while ( true ) {
      Real f = fun->eval( x ); ++m_fun_evaluation_count;
      m_converged = abs(f) < m_tolerance;
      if ( m_converged ) return x;

      ++m_iteration_count;
      if ( m_iteration_count > m_max_iteration || m_fun_evaluation_count > m_max_fun_evaluation ) break;

      Real df = fun->eval_D( x ); ++m_fun_evaluation_count;
      Real y  = x - (f/df);

      if ( !Utils::is_finite(y) ) {
        fmt::print( "found y[{}] = {}, f'(x={}) = {}\n", m_iteration_count, y, x, df );
        break;
      }

      Real fy = fun->eval( y ); ++m_fun_evaluation_count;
      Real t  = fy/f;
      Real z  = y - Q(t)*(fy/df);
      Real fz = fun->eval( z ); ++m_fun_evaluation_count;

      m_converged = abs(fz) < m_tolerance;
      if ( m_converged ) return z;

      Real s  = fz/fy;
      Real w  = z - W(t,s)*(fz/df);
      Real fw = fun->eval( w ); ++m_fun_evaluation_count;

      m_converged = abs(fw) < m_tolerance;
      if ( m_converged ) return w;

      Real u  = fw/fz;
      Real h  = w - H(t,s,u)*(fz/df);
      Real fh = fun->eval( h ); ++m_fun_evaluation_count;

      Real v  = fh/fw;
      x = h - J(t,s,u,v)*(fw/df);

      if ( !Utils::is_finite(x) ) {
        fmt::print("found x[{}] = {}\n",m_iteration_count,x);
        break;
      }
    }
    return 0;
  }

  template class Zeros<float>;
  template class Zeros<double>;
}
