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

//
// file: Utils_Trichotomy.hh
//

#pragma once

#ifndef UTILS_TRICHOTOMY_dot_HH
#define UTILS_TRICHOTOMY_dot_HH

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>

#include "Utils.hh"

namespace Utils {

  using std::pow;
  using std::abs;

  /*!
   * \addtogroup Minimize
   * @{
   */

  //!
  //! \class Trichotomy_base_fun
  //! \brief Abstract base class for defining mathematical functions used in the minimization algorithm.
  //!
  //! This class serves as a base interface for user-defined functions that can be evaluated.
  //! It allows for the implementation the numerical method to
  //! find the minimum of the one dimensional equation \f$ f(x) \f$.
  //! Users must inherit from this class and implement the virtual method to define their specific functions.
  //!
  //! **Template Parameter:**
  //! - `Real`: A numeric type representing the data type of the function's input and output,
  //!   such as `float`, `double`, etc.
  //!
  //! **Usage Example:**
  //! To create a custom function, derive from this class and implement the required methods.
  //! Here is an example for the function \f$ f(x) = x^2 - 2 \f$:
  //!
  //! \code{cpp}
  //! class Fun1 : public Trichotomy_base_fun<double> {
  //! public:
  //!     double eval(double x) const override { return x*x - 2; }
  //! };
  //! \endcode
  template <typename Real>
  class Trichotomy_base_fun {
  public:
    //!
    //! Evaluate the function \f$ f(x) \f$
    //!
    //! \param[in] x the point to evaluate \f$ f(x) \f$
    //! \return the value of \f$ f(x) \f$
    //!
    virtual Real eval( Real x ) const = 0;
    //!
    //! Evaluate the function \f$ f(x) \f$
    //!
    //! \param[in] x the point to evaluate \f$ f(x) \f$
    //! \return the value of \f$ f(x) \f$
    //!
    Real operator () ( Real x ) const { return this->eval(x); }
  };

  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename Real, typename PFUN>
  class Trichotomy_fun : public Trichotomy_base_fun<Real> {
    PFUN m_fun;
  public:
    explicit Trichotomy_fun( PFUN pfun ) : m_fun(pfun) {}
    Real eval( Real x ) const override { return m_fun(x); };
  };
  #endif

  //!
  //! \class Trichotomy
  //! \brief Class for minimize the function \f$ f(x) \f$ without the usew of derivative
  //!
  //! \note The used algorithm is described in:
  //!       *A new zero-order 1-D optimization algorithm: trichotomy method*,
  //!       by **Alena Antonova, Olga Ibryaeva**, [link](https://doi.org/10.48550/arXiv.1903.07117)
  //!
  //! **Usage Example**
  //!
  //! To use this class, first wrap your function in a derived class. For instance, for the function \f$ f(x) = x^2 - 2 \f$, you can define:
  //!
  //! \code{cpp}
  //! class Fun1 : public Trichotomy_base_fun<double> {
  //! public:
  //!   double eval(double x) const override { return x*x - 2; }
  //! };
  //! \endcode
  //!
  //! Next, instantiate the function and the solver. Then, call the desired method to find the root:
  //!
  //! \code{cpp}
  //! Trichotomy<real_type> solver;
  //! Fun1 f;
  //! real_type a=-1,b=2;
  //! real_type x_solution = solver.eval2(a,b,f);
  //! \endcode
  //!
  //! If the method converges, `x_solution` will contain the computed solution.
  //!
  template <typename Real>
  class Trichotomy {

    using Integer = int;

    Integer m_num_iter_done{0};
    Integer m_num_fun_eval{0};
    Integer m_max_iteration{100};
    Integer m_max_fun_evaluation{1000};
    Real    m_tolerance{pow(machine_eps<Real>(),Real(2./3.))};

    bool m_converged{false};

    Real m_a{0},  m_fa{0};
    Real m_x1{0}, m_f1{0};
    Real m_x2{0}, m_f2{0};
    Real m_x3{0}, m_f3{0};
    Real m_x4{0}, m_f4{0};
    Real m_x5{0}, m_f5{0};
    Real m_b{0},  m_fb{0};

    Trichotomy_base_fun<Real> * m_function{nullptr};

    bool bracketing();
    Real minimize();
    //void set_tolerance( Real tol );
    Real evaluate( Real x ) { ++m_num_fun_eval; return m_function->eval(x); };
    Real eval( Real a, Real b );
    Real search( Real x, Real delta );

  public:

    Trichotomy() = default;
    ~Trichotomy() = default;

    //!
    //! Find the minimum of a function wrapped in the class `Trichotomy_base_fun<Real>`
    //! starting from guess interval `[a,b]`
    //!
    //! \param a    lower bound search interval
    //! \param b    upper bound search interval
    //! \param fun  the pointer to base class `Trichotomy_base_fun<Real>` wrapping the user function
    //!
    Real
    eval( Real a, Real b, Trichotomy_base_fun<Real> * fun ) {
      m_function = fun;
      return this->eval( a, b );
    }

    //!
    //! Find the minimum of a function `pfun`
    //! starting from guess interval `[a,b]`
    //!
    //! \param a    lower bound search interval
    //! \param b    upper bound search interval
    //! \param pfun the pointer to base class `Trichotomy_base_fun<Real>` wrapping the user function
    //!
    template <typename PFUN>
    Real
    eval2( Real a, Real b, PFUN pfun ) {
      Trichotomy_fun<Real,PFUN> fun( pfun );
      m_function = &fun;
      return this->eval(a,b);
    }

    //!
    //! Find the minimum of a function wrapped in the class `Trichotomy_base_fun<Real>`
    //! starting from guess interval `[x-delta,x+delta]`
    //!
    //! \param x     center of the search interval
    //! \param delta ray of th eseaerch interval
    //! \param fun   the pointer to base class `Trichotomy_base_fun<Real>` wrapping the user function
    //!
    Real
    search( Real x, Real delta, Trichotomy_base_fun<Real> * fun ) {
      m_function = fun;
      return this->search( x, delta );
    }

    //!
    //! Find the minimum of a function object
    //! starting from guess interval `[x-delta,x+delta]`
    //!
    //! \param x     center of the search interval
    //! \param delta ray of th eseaerch interval
    //! \param pfun  the function object
    //!
    template <typename PFUN>
    Real
    search2( Real x, Real delta, PFUN pfun ) {
      Trichotomy_fun<Real,PFUN> fun( pfun );
      m_function = &fun;
      return this->search( x, delta );
    }

    //!
    //! \return the number of iterations used in the last computation
    //!
    Integer used_iter() const { return m_num_iter_done; }
    //!
    //! \return the number of evaluation used in the last computation
    //!
    Integer num_fun_eval() const { return m_num_fun_eval; }
    //!
    //! \return the tolerance set for computation
    //!
    Real tolerance() const { return m_tolerance; }
    //!
    //! \return true if the last computation was successfull
    //!
    bool converged() const { return m_converged; }

    //!
    //! Fix the maximum number of iteration.
    //!
    //! \param mit the maximum number of iteration
    //!
    void set_max_iterations( Integer mit );
    //!
    //! Fix the maximum number of evaluation.
    //!
    //! \param mfev the maximum number of evaluation of \f$ f(x) \f$
    //!
    void set_max_fun_evaluation( Integer mfev );

  };

  #ifndef UTILS_OS_WINDOWS
  extern template class Trichotomy<float>;
  extern template class Trichotomy<double>;
  #endif

  /*! @} */

}

#endif

//
// eof: Utils_Trichotomy.hh
//
