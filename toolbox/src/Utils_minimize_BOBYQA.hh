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
// file: Utils_minimize_BOBYQA.hh
//

#pragma once

#ifndef UTILS_BOBYQA_dot_HH
#define UTILS_BOBYQA_dot_HH

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wswitch-enum"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma clang diagnostic ignored "-Wweak-vtables"
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wdouble-promotion"
#pragma clang diagnostic ignored "-Wsigned-enum-bitfield"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-vtables"
#pragma clang diagnostic ignored "-Wunused-template"
#pragma clang diagnostic ignored "-Wnon-virtual-dtor"
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored "-Wmissing-noreturn"
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include "Utils.hh"
#include "Utils_fmt.hh"
#include "Utils_eigen.hh"

namespace Utils
{

  using std::abs;
  using std::max;
  using std::min;
  using std::pow;
  using std::sqrt;

  /*
   * bobyqa.h -
   *
   * Definitions for Mike Powell's BOBYQA algorithm for minimizing a function of
   * many variables.  The method is "derivatives free" (only the function values
   * are needed) and accounts for bound constraints on the variables.  The
   * algorithm is described in:
   *
   *   M.J.D. Powell, "The BOBYQA Algorithm for Bound Constrained Optimization
   *   Without Derivatives."  Technical report, Department of Applied Mathematics
   *   and Theoretical Physics, University of Cambridge (2009).
   *
   * The present code is based on the original FORTRAN version written by Mike
   * Powell who kindly provides his code on demand (at mjdp@cam.ac.uk) and has
   * been converted to C by É. Thiébaut.
   *
   * Copyright (c) 2009, Mike Powell (FORTRAN version).
   * Copyright (c) 2015, Éric Thiébaut (C version).
   *
   * Read the accompanying `LICENSE` file for details.
   */

  /* BOBYQA seeks the least value of a function of many variables, by applying a
     trust region method that forms quadratic models by interpolation.  There is
     usually some freedom in the interpolation conditions, which is taken up by
     minimizing the Frobenius norm of the change to the second derivative of the
     model, beginning with the zero matrix.  The values of the variables are
     constrained by upper and lower bounds.  The arguments of the subroutine are
     as follows.

     N must be set to the number of variables and must be at least two.  NPT is
     the number of interpolation conditions.  Its value must be in the interval
     [N+2,(N+1)(N+2)/2].  Choices that exceed 2*N+1 are not recommended.

     OBJFUN is provided by the user to compute the objective function value at
     the values of the variables X(1),X(2),...,X(N), which are generated
     automatically by BOBYQA in a way that satisfies the bounds given in XL and
     XU.  DATA is anything needed by the function and which is passed as is to
     OBJFUN by BOBYQA.

     Initial values of the variables must be set in X(1),X(2),...,X(N).  They
     will be changed to the values that give the least calculated F.  For
     I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper bounds,
     respectively, on X(I).  The construction of quadratic models requires XL(I)
     to be strictly less than XU(I) for each I.  Further, the contribution to a
     model from changes to the I-th variable is damaged severely by rounding
     errors if XU(I)-XL(I) is too small.

     RHOBEG and RHOEND must be set to the initial and final values of a trust
     region radius, so both must be positive with RHOEND no greater than RHOBEG.
     Typically, RHOBEG should be about one tenth of the greatest expected change
     to a variable, while RHOEND should indicate the accuracy that is required in
     the final values of the variables.  An error return occurs if any of the
     differences XU(I)-XL(I), I=1,...,N, is less than 2*RHOBEG.

     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the amount
     of printing.  Specifically, there is no output if IPRINT=0 and there is
     output only at the return if IPRINT=1.  Otherwise, each new value of RHO is
     printed, with the best vector of variables so far and the corresponding
     value of the objective function.  Further, each new value of F with its
     variables are output if IPRINT=3.

     MAXFUN must be set to an upper bound on the number of calls of OBJFUN.

     The array W will be used for working space.  Its length must be at least
     (NPT+5)*(NPT+N)+3*N*(N+5)/2.  Upon successful return, the first element of W
     will be set to the function value at the solution.
  */

  template <typename Scalar>
  class BOBYQA_minimizer
  {
  public:
    using integer = Eigen::Index;
    using Vector  = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using VectorI = Eigen::Matrix<integer, Eigen::Dynamic, 1>;
    using Matrix  = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    /* Prototype of the objective function assumed by the BOBYQA routine.  The
       returned value is the function value at X the current variables, N is the
       number of variables and DATA is anything needed by the function (unused by
       BOBYQA itself). */

    using bobyqa_objfun = std::function<Scalar( Vector const & x )>;

    using Status = enum class Status : int {
      BOBYQA_SUCCESS              = 0,   // algorithm converged
      BOBYQA_BAD_NPT              = -1,  // NPT is not in the required interval
      BOBYQA_TOO_CLOSE            = -2,  // insufficient space between the bounds
      BOBYQA_ROUNDING_ERRORS      = -3,  // too much cancellation in a denominator
      BOBYQA_TOO_MANY_EVALUATIONS = -4,  // maximum number of function evaluations exceeded
      BOBYQA_STEP_FAILED          = -5   // a trust region step has failed to reduce Q
    };

    static string
    to_string( Status status )
    {
      switch ( status )
      {
        case Status::BOBYQA_SUCCESS:
          return "algorithm converged";
        case Status::BOBYQA_BAD_NPT:
          return "NPT is not in the required interval";
        case Status::BOBYQA_TOO_CLOSE:
          return "insufficient space between the bounds";
        case Status::BOBYQA_ROUNDING_ERRORS:
          return "too much cancellation in a denominator";
        case Status::BOBYQA_TOO_MANY_EVALUATIONS:
          return "maximum number of function evaluations exceeded";
        case Status::BOBYQA_STEP_FAILED:
          return "a trust region step has failed to reduce Q";
      }
      return "";
    }

  private:
    bobyqa_objfun m_fun;
    integer       m_num_f_eval   = 0;
    integer       m_max_f_eval   = 10000;
    integer       n_num_f_saved  = 0;
    integer       n_num_f_rescue = 0;
    
    Status        m_status;
    string        m_reason;

    Scalar eval( Vector const & x );

    integer m_nv;
    integer m_npt;
    integer m_dim;
    integer m_nptm;
    integer m_print_level = 3;

    integer m_kopt;
    integer m_knew;

    // ========================================================================
    // COSTANTI NUMERICHE
    // ========================================================================

    Scalar m_tol_convergence = Scalar( 1e-4 );  ///< Tolleranza per convergenza
    Scalar m_tol_step        = Scalar( 1e-2 );  ///< Tolleranza passo minimo
    Scalar m_eps             = Scalar( 1e-20 );

    Scalar m_rho;
    Scalar m_rhobeg = Scalar( 0.1 );
    Scalar m_rhoend = Scalar( 1e-6 );
    Scalar m_crvmin;
    Scalar m_dsq;
    Scalar m_alpha;
    Scalar m_beta;
    Scalar m_delta;
    Scalar m_adelt;
    Scalar m_cauchy;
    Scalar m_denom;
    Scalar m_x_opt_square;
    Scalar m_distsq;

    VectorI m_xbdi;

    Vector m_x_lower;
    Vector m_x_upper;

    Vector m_x_base;
    Vector m_x_new;
    Vector m_g_new;
    Vector m_x_opt;
    Vector m_g_opt;

    Vector m_x_alt;
    Vector m_f_val;
    Vector m_pq;
    Vector m_s_lower;
    Vector m_s_upper;
    Vector m_d;
    Vector m_v_lag;
    Vector m_g_lag;
    Vector m_hcol;
    Vector m_curv;

    Matrix m_xpt;
    Matrix m_B;
    Matrix m_Z;
    Matrix m_ptsaux;
    Vector m_ptsid;
    Matrix m_HQ;

    Vector m_s;
    Vector m_hs;
    Vector m_hred;

    Vector m_WNPT;

    Vector m_VN0;
    Vector m_VN1;

    static Scalar
    power2( Scalar const x )
    {
      return x * x;
    }
    static bool
    is_zero( Scalar const x )
    {
      return std::abs( x ) < 1e-20;
    }

  public:
    void
    set_rho( Scalar beg, Scalar end )
    {
      m_rhobeg = beg;
      m_rhoend = end;
    }
    void
    set_verbosity( integer v )
    {
      m_print_level = v;
    }
    void
    set_maxfun( integer m )
    {
      m_max_f_eval = m;
    }

  private:
    std::string print_vec( Vector const & x, integer max_elem ) const;

    void
    print_error( string const & reason )
    {
      fmt::print( "\n    Return from BOBYQA because {}.\n", reason );
    }

    void compute_hessian_product( Vector const & s, Vector & hs ) const;
    void add_hessian_product( Vector const & s, Vector & hs ) const;

    void altmov();
    void prelim( Vector & x );
    void trsbox();
    void rescue();
    void update();

  public:
    // Debug accessor to validate BMAT symmetry post-update (used by tests)
    bool debug_check_bmat_symmetry( Scalar tol = Scalar( 1e-12 ) ) const;

    // Test and debug helpers (initialize internal sizes and set/get internals)
    void debug_init( integer n, integer npt );
    void
    debug_set_B( Matrix const & B )
    {
      m_B = B;
    }
    void
    debug_set_Z( Matrix const & Z )
    {
      m_Z = Z;
    }
    void
    debug_set_vlag( Vector const & vlag )
    {
      m_v_lag = vlag;
    }
    void
    debug_set_pq( Vector const & pq )
    {
      m_pq = pq;
    }
    void
    debug_set_knew( integer knew )
    {
      m_knew = knew;
    }
    void
    debug_set_denom( Scalar denom )
    {
      m_denom = denom;
    }
    void
    debug_set_beta( Scalar beta )
    {
      m_beta = beta;
    }
    Matrix
    debug_B() const
    {
      return m_B;
    }
    Matrix
    debug_Z() const
    {
      return m_Z;
    }
    Vector
    debug_vlag() const
    {
      return m_v_lag;
    }
    Vector
    debug_pq() const
    {
      return m_pq;
    }
    integer
    debug_knew() const
    {
      return m_knew;
    }
    Scalar
    debug_denom() const
    {
      return m_denom;
    }
    Scalar
    debug_beta() const
    {
      return m_beta;
    }
    // Execute the private update() method for testing
    void
    debug_update()
    {
      update();
    }

    Status minimize( integer const n,
                     integer const npt,
                     bobyqa_objfun const &,
                     Vector &       x,
                     Vector const & xlower,
                     Vector const & xupper );
  };

}  // namespace Utils

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif

//
// eof: Utils_minimize_BOBYQA.hh
//
