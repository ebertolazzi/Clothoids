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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: Utils_NelderMead.hh
///

#pragma once

#ifndef UTILS_NELDER_MEAD_dot_HH
#define UTILS_NELDER_MEAD_dot_HH

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <functional>
#include <cmath>

#include "Utils.hh"
#include "Utils_eigen.hh"

namespace Utils {

  /*
  //
  //  J.A.Nelder and R.Mead
  //  Computer Journal, vol. 7, pp. 308-313, 1965.
  //
  //  L.A.Yarbro and S.N.Deming
  //  Analytica Chimica Acta, vol. 73, pp. 391-398, 1974.
  //
  //  S.L.S.Jacoby, J.S.Kowalik and J.T.Pizzo
  //  Iterative Methods for Nonlinear Optimization Problems
  //  (Englewood Cliffs, NJ: Prentice-Hall), 1972
  */
  template <typename Real>
  class NelderMead {

    typedef Eigen::Matrix<Real,Eigen::Dynamic,1>              Vec_t;
    typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> Mat_t;
    typedef Eigen::Map<Vec_t>                                 MapVec;
    typedef Eigen::Map<Mat_t>                                 MapMat;
    typedef int                                               integer;
    typedef std::function<Real(Real const[])>                 NMFunc;
    typedef enum {
      NM_INIT=0,
      NM_REFLECT,
      NM_EXPAND_FE,
      NM_EXPAND_FR,
      NM_CONTRACT_O,
      NM_CONTRACT_I,
      NM_SHRINK,
      NM_RESTART,
      NM_WORSE
    } NMstype;

  private:

    string const m_name;

    Malloc<Real> m_base_value;

    NMFunc m_fun; // handle to the value function to minimize
    Console const * m_console{nullptr}; //!< pointer to the message stream class

    integer m_dim; // problem dimension (number of variables)
    Real    m_r_dim;
    Real    m_dim_factorial;
    Real    m_dim_regular_simplex_volume;
    Real    m_tolerance{Real(1e-8)};

    integer m_low;
    integer m_0high;
    integer m_high;

    MapVec m_f{nullptr,0};
    MapMat m_p{nullptr,0,0};
    MapVec m_psum{nullptr,0};
    MapMat m_dist{nullptr,0,0};

    MapMat m_p_work{nullptr,0,0};
    MapVec m_f_work{nullptr,0};

    Eigen::PartialPivLU<Mat_t> m_lu;

    MapVec m_pr{nullptr,0};
    MapVec m_pe{nullptr,0};
    MapVec m_pc{nullptr,0};

    MapVec m_grad{nullptr,0};
    Real   m_diameter{0};
    Real   m_simplex_volume{0};

    Real   m_rho{Real(1.0)};
    Real   m_chi{Real(2.0)};
    Real   m_gamma{Real(0.5)};
    Real   m_sigma{Real(0.25)};
    Real   m_prob{Real(0.3)};
    Real   m_volume_tolerance{Real(0.01)};

    integer m_max_fun_evaluation{5000}; // max number of function evaluations
    integer m_max_iteration{1000};      // max number of iterations

    mutable integer m_iteration_count{0};    // explore iteration counter
    mutable integer m_fun_evaluation_count{0};

    NMstype m_which_step{NM_INIT};
    integer m_verbose{1}; // flag to activate info printing

    void allocate( integer n );

    Real extrapolate( Real alpha, integer j, MapVec & pe ) const;
    void spendley( MapVec const & X0, Real len );
    void gradient( MapVec & G ) const;
    void shrink();
    void dist_init();
    void dist_update( integer jpos );
    void grad_update( integer jpos );
    Real reflect( MapVec & pr ) { return this->extrapolate( m_rho,         m_high, pr ); }
    Real expand ( MapVec & pe ) { return this->extrapolate( m_rho*m_chi,   m_high, pe ); }
    Real outside( MapVec & po ) { return this->extrapolate( m_rho*m_gamma, m_high, po ); }
    Real inside ( MapVec & pi ) { return this->extrapolate( -m_gamma,      m_high, pi ); }

    void
    replace_point( Real fj, MapVec const & pj, integer jpos ) {
      m_f(jpos)               = fj;
      m_psum                 += pj - m_p.col(jpos);
      m_p.col(jpos).noalias() = pj;
      this->dist_update( jpos );
    }

    void spendley( Real const X0[], Real delta );
    void diamond( Real const X0[], Real delta );

    Real
    eval_function( Real const x[] ) const {
      ++m_fun_evaluation_count;
      return m_fun( x );
    }

  public:

    explicit
    NelderMead( string const & name )
    : m_name(name)
    , m_base_value("NelderMead_"+name)
    {}

    //!
    //! Get the name of the class.
    //!
    string const & name(void) const { return m_name; }

    void setup( integer dim, NMFunc & fun, Console const * console );

    void change_console( Console const * console ) { m_console = console; }

    void set_verbose( integer level ) { m_verbose = level; }
    void set_tolerance( Real tol );
    void set_max_iterations( integer mit );
    void set_max_fun_evaluation( integer mfev );

    void message( Real tol ) const;

    string info() const;
    void print_info( ostream_type & stream ) const { stream << info(); }
    void explore();
    bool search();
    bool run( Real const x_sol[], Real h );

    void
    get_last_solution( Real x[] ) const {
      std::copy_n( m_p.col(m_low).data(), m_dim, x );
    }
  };
}

#endif

///
/// eof: Utils_NelderMead.hh
///
