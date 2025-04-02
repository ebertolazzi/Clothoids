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
// file: Utils_NelderMead.hh
//

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
  /*!
   * \addtogroup Minimize
   * @{
   */

  //!
  //!  \class NelderMead
  //!  \brief Implements the Nelder-Mead optimization algorithm for nonlinear optimization problems.
  //!
  //!  The Nelder-Mead algorithm is a widely used simplex-based method for minimizing an objective function
  //!  in multiple dimensions. It is especially useful for problems where derivative information is not readily
  //!  available or when the objective function is noisy.
  //!
  //!  \tparam Real The numeric type used for calculations, typically `float` or `double`.
  //!
  //!  \details
  //!  This implementation of the Nelder-Mead algorithm is based on several foundational works:
  //!
  //!  - J.A. Nelder and R. Mead, "A Simplex Method for Function Minimization," _Computer Journal_, vol. 7, pp. 308-313, 1965.
  //!  - L.A. Yarbro and S.N. Deming, "Optimization by the Simplex Method," _Analytica Chimica Acta_, vol. 73, pp. 391-398, 1974.
  //!  - S.L.S. Jacoby, J.S. Kowalik, and J.T. Pizzo, _Iterative Methods for Nonlinear Optimization Problems_ (Englewood Cliffs, NJ: Prentice-Hall), 1972.
  //!
  //!  \note The Nelder-Mead algorithm does not require derivatives, making it suitable for non-smooth or noisy functions.
  //!
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

    //!
    //! \brief Constructor that initializes the NelderMead instance with a specific name.
    //!
    //! \param name The name to be assigned to this instance of the Nelder-Mead optimizer.
    //!
    explicit
    NelderMead( string_view name )
    : m_name(name)
    , m_base_value( string("NelderMead_")+string(name) )
    {}

    //!
    //! \brief Get the name of the class.
    //!
    //! \return A string representing the name of the Nelder-Mead optimizer instance.
    //!
    string_view name(void) const { return m_name; }

    //!
    //! \brief Setup the optimizer with the function and console interface.
    //!
    //! \param dim     The dimensionality of the optimization problem.
    //! \param fun     A reference to the function to be minimized.
    //! \param console A pointer to the console for output logging during optimization.
    //!
    void setup( integer dim, NMFunc & fun, Console const * console );

    //!
    //! \brief Change the console output interface.
    //!
    //! \param console A pointer to the new console object for logging.
    //!
    void change_console( Console const * console ) { m_console = console; }

    //!
    //! \brief Set the verbosity level for logging during the optimization process.
    //!
    //! \param level The verbosity level (higher values produce more detailed output).
    //!
    void set_verbose( integer level ) { m_verbose = level; }

    //!
    //! \brief Set the tolerance for stopping criteria in the optimization process.
    //!
    //! \param tol The tolerance value, typically a small positive number (e.g., 1e-6).
    //!
    void set_tolerance( Real tol );

    //!
    //! \brief Set the maximum number of iterations allowed in the optimization process.
    //!
    //! \param mit The maximum number of iterations to be performed.
    //!
    void set_max_iterations( integer mit );

    //!
    //! \brief Set the maximum number of function evaluations allowed in the optimization process.
    //!
    //! \param mfev The maximum number of function evaluations.
    //!
    void set_max_fun_evaluation( integer mfev );

    //!
    //! \brief Display a message with information about the current tolerance level.
    //!
    //! \param tol The tolerance level for which a message should be displayed.
    //!
    void message( Real tol ) const;

    //!
    //! \brief Provide information about the current state of the Nelder-Mead optimization.
    //!
    //! \return A string with detailed information about the current state of the optimizer.
    //!
    string info() const;

    //!
    //! \brief Print detailed information about the current state of the Nelder-Mead optimizer.
    //!
    //! \param stream The output stream where the information will be printed.
    //!
    void print_info( ostream_type & stream ) const { stream << info(); }

    //void explore();

    //!
    //! \brief Perform the search step of the optimization process.
    //!
    //! \return True if the search was successful; otherwise, false.
    //!
    bool search();

    //!
    //! \brief Run the optimization process starting from an initial guess and step size.
    //!
    //! \param x_sol An array containing the initial guess for the solution.
    //! \param h     The initial step size for the simplex search.
    //!
    //! \return True if the optimization converged successfully; otherwise, false.
    //!
    bool run( Real const x_sol[], Real h );

    //!
    //! \brief Get the solution found by the Nelder-Mead algorithm after the last run.
    //!
    //! \param x An array to store the solution values (the optimizer writes the solution here).
    //! \return The function value at the solution point.
    //!
    Real
    get_last_solution( Real x[] ) const {
      std::copy_n( m_p.col(m_low).data(), m_dim, x );
      return m_f.coeff(m_low);
    }

    //!
    //! \brief Get the best function value found during the optimization process.
    //!
    //! \return The best (lowest) function value found.
    //!
    Real get_better_value() const { return m_f.coeff(m_low); }

    //!
    //! \brief Get the worst function value encountered during the optimization process.
    //!
    //! \return The worst (highest) function value found.
    //!
    Real get_worst_value() const { return m_f.coeff(m_high); }
  };

  /*! @} */

}

#endif

//
// eof: Utils_NelderMead.hh
//
