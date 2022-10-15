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
/// file: Utils_HJPatternSearch.hh
///

#pragma once

#ifndef UTILS_HJ_PATTERN_SEARCH_dot_HH
#define UTILS_HJ_PATTERN_SEARCH_dot_HH

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
  // HJPatternSearch class decription
  //
  // This class implements the Hooke-Jeeves Pattern Search Algorithm.
  // The references can be found in:
  //
  // R. Hooke, T. A. Jeeves,
  // "Direct Search" Solution of Numerical and Statistical Problems,
  // Westinghouse Research Laboratories,
  // Pittsburg, Pennsylvania
  //
  // Arthur Kaupe,
  // Algorithm 178: Direct Search,
  // Communications of the ACM,
  // Volume 6, Number 6, June 1963, page 313.
  //
  // M Bell, Malcolm Pike,
  // Remark on Algorithm 178: Direct Search,
  // Communications of the ACM,
  // Volume 9, Number 9, September 1966, page 684.
  //
  // FK Tomlin, LB Smith,
  // Remark on Algorithm 178: Direct Search,
  // Communications of the ACM,
  // Volume 12, Number 11, November 1969, page 637-638.
  //
  */
  template <typename Real>
  class HJPatternSearch {
    typedef Eigen::Matrix<Real,Eigen::Dynamic,1> Vec_t;
    typedef Eigen::Map<Vec_t>                    MapVec;
    typedef int                                  integer;
    typedef std::function<Real(Real const[])>    HJFunc;
  private:
    string const m_name;

    Malloc<Real> m_base_value;

    HJFunc m_fun; // handle to the value function to minimize
    Console const * m_console{nullptr}; //!< pointer to the message stream class

    Real m_rho{Real(0.9)}; // stencil step decreasing factor (must be 0 < rho < 1)
    MapVec m_h{nullptr,0}; // scale of the stencil (Remark 1966)

    bool m_stencil_failure; // stencil failure flag - used to shrink h,
                            // stencil_failure = true means failure

    integer m_N{0};           // problem dimension (number of variables)
    MapVec m_x{nullptr,0};    // current iteration
    MapVec m_pars{nullptr,0}; // current iteration

    MapVec m_x_best{nullptr,0};    // base point
    MapVec m_pars_best{nullptr,0}; // base point

    MapVec m_x_center{nullptr,0};    // stencil center
    MapVec m_pars_center{nullptr,0}; // stencil center

    MapVec m_d{nullptr,0};           // search direction
    MapVec m_search_sign{nullptr,0}; // vector to keep in memory the direction of function value descent from the previous iteration in each direction j

    MapVec m_x_current_best{nullptr,0};
    MapVec m_pars_current_best{nullptr,0};

    MapVec m_x_temporary{nullptr,0}; // temporary point representing the center of the stencil
    MapVec m_p{nullptr,0};
    MapVec m_p1{nullptr,0};

    Real m_tolerance{Real(1e-8)};                  // tolerance on the scale h
    Real m_sqrt_tolerance{Real(sqrt(Real(1e-8)))}; // tolerance on the scale h
    integer m_max_iteration{50};                   // max number of iterations
    integer m_max_fun_evaluation{1000};            // max number of function evaluations
    integer m_max_num_stagnation{10};
    integer m_verbose{1};                          // flag to activate info printing

    Real m_fun_best;           // best value function (fun evaluated in x_best_current)
    integer m_iteration_count; // explore iteration counter
    integer m_num_stagnation;

    void allocate( integer n );
    Real h_norm_inf() const { return m_h.template lpNorm<Eigen::Infinity>(); }
    Real eval_function( MapVec const & x, MapVec & pars ) const;

    mutable integer m_fun_evaluation;

  public:

    explicit
    HJPatternSearch( string const & name )
    : m_name(name)
    , m_base_value("HJ_"+name)
    {}

    //!
    //! Get the name of the class.
    //!
    string const & name(void) const { return m_name; }

    void setup( HJFunc & fun, Console const * console );

    void change_console( Console const * console ) { m_console = console; }

    void set_verbose( integer level ) { m_verbose = level; }
    void set_tolerance( Real tol );
    void set_max_iterations( integer mit );
    void set_max_fun_evaluation( integer mfev );
    void set_max_num_stagnation( integer nstg );

    string info() const;
    void print_info( ostream_type & stream ) const { stream << info(); }
    void explore();
    void search();
    void run( Real x_sol[], Real h );
  };
}

#endif

///
/// eof: Utils_HJPatternSearch.hh
///
