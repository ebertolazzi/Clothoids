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
  public:

    using Vec_t   = Eigen::Matrix<Real,Eigen::Dynamic,1>;
    using Mat_t   = Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic>;
    using MapVec  = Eigen::Map<Vec_t>;
    using MapMat  = Eigen::Map<Mat_t>;
    using integer = int;
    using NMFunc  = std::function<Real(Real const[])>;
    using NM_move = enum class NelderMead_move : integer {
      INIT,
      REFLECT,
      EXPAND_FE,
      EXPAND_FR,
      CONTRACT_O,
      CONTRACT_I,
      SHRINK,
      RESTART,
      WORSE
    };

    static string to_string( NM_move n ) {
      string res = "";
      switch ( n ) {
      case NM_move::INIT:       res = "INIT";       break;
      case NM_move::REFLECT:    res = "REFLECT";    break;
      case NM_move::EXPAND_FE:  res = "EXPAND_FE";  break;
      case NM_move::EXPAND_FR:  res = "EXPAND_FR";  break;
      case NM_move::CONTRACT_O: res = "CONTRACT_O"; break;
      case NM_move::CONTRACT_I: res = "CONTRACT_I"; break;
      case NM_move::SHRINK:     res = "SHRINK";     break;
      case NM_move::RESTART:    res = "RESTART";    break;
      case NM_move::WORSE:      res = "WORSE";      break;
      }
      return res;
    }

  private:

    string const m_name;

    Malloc<Real> m_base_value;

    NMFunc m_fun; // handle to the value function to minimize
    Console const * m_console{nullptr}; //!< pointer to the message stream class

    integer m_dim{0}; // problem dimension (number of variables)
    Real    m_r_dim{0};
    Real    m_dim_factorial{0};
    Real    m_dim_regular_simplex_volume{0};
    Real    m_tolerance{Real(1e-8)};

    integer m_low{0};
    integer m_0high{0};
    integer m_high{0};

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

    NM_move m_which_step{NM_move::INIT};
    integer m_verbose{1}; // flag to activate info printing

    void allocate( integer n );

    Real extrapolate( Real alpha, integer j, MapVec & pe ) const;
    //void spendley( MapVec const & X0, Real len );
    //void gradient( MapVec & G ) const;
    void shrink();
    void dist_init();
    void dist_update( integer jpos );
    void grad_update( integer jpos );
    Real reflect( MapVec & pr ) { return this->extrapolate( m_rho,         m_high, pr ); }
    Real expand ( MapVec & pe ) { return this->extrapolate( m_rho*m_chi,   m_high, pe ); }
    Real outside( MapVec & po ) { return this->extrapolate( m_rho*m_gamma, m_high, po ); }
    Real inside ( MapVec & pi ) { return this->extrapolate( -m_gamma,      m_high, pi ); }

    void replace_point( Real fj, MapVec const & pj, integer jpos );
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
    //void explore();
    bool search();
    bool run( Real const x_sol[], Real h );

    Real
    get_last_solution( Real x[] ) const {
      std::copy_n( m_p.col(m_low).data(), m_dim, x );
      return m_f.coeff(m_low);
    }

    Real get_better_value() const { return m_f.coeff(m_low); }
    Real get_worst_value() const { return m_f.coeff(m_high); }
  };
}

#endif

///
/// eof: Utils_NelderMead.hh
///
