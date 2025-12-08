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
// file: Utils_HNewton.hh
//

#pragma once

#ifndef UTILS_ALGO_HNEWTON_dot_HH
#define UTILS_ALGO_HNEWTON_dot_HH

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "Utils.hh"

namespace Utils
{

  using std::abs;
  using std::pow;

  /*!
   * \addtogroup Zeros
   * @{
   */

  /*
  //   _   _ _   _               _
  //  | | | | \ | | _____      _| |_ ___  _ __
  //  | |_| |  \| |/ _ \ \ /\ / / __/ _ \| '_ \
  //  |  _  | |\  |  __/\ V  V /| || (_) | | | |
  //  |_| |_|_| \_|\___| \_/\_/  \__\___/|_| |_|
  */

  //!
  //! \class AlgoHNewton_base_fun
  //! \brief Abstract base class for defining mathematical functions used in the
  //! zero search algorithm.
  //!
  //! This class serves as a base interface for user-defined functions that can
  //! be evaluated. It allows for the implementation the numerical method to
  //! find the solution of the one dimensional equation \f$ f(x) = 0 \f$.
  //! Users must inherit from this class and implement the virtual method to
  //! define their specific functions.
  //!
  //! **Template Parameter:**
  //! - `Real`: A numeric type representing the data type of the function's
  //! input and output,
  //!   such as `float`, `double`, etc.
  //!
  //! **Usage Example:**
  //! To create a custom function, derive from this class and implement the
  //! required methods. Here is an example for the function \f$ f(x) = x^2 - 2
  //! \f$:
  //!
  //! \code{cpp}
  //! class Fun1 : public AlgoHNewton_base_fun<double> {
  //! public:
  //!   double eval(double x) const override { return x*x - 2; }
  //!   double D   (double x) const override { return 2*x; }
  //! };
  //! \endcode
  //!
  template <typename Real>
  class AlgoHNewton_base_fun
  {
  public:
    //!
    //! Evaluate the function \f$ f(x) \f$
    //!
    //! \param[in] x the point to evaluate \f$ f(x) \f$
    //! \return the value of \f$ f(x) \f$
    //!
    virtual Real eval( Real x ) const = 0;
    //!
    //! Evaluate the function \f$ f'(x) \f$
    //!
    //! \param[in] x the point to evaluate \f$ f(x) \f$
    //! \return the value of \f$ f'(x) \f$
    //!
    virtual Real D( Real x ) const = 0;
  };

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  template <typename Real, typename PFUN, typename PFUN_D>
  class AlgoHNewton_fun : public AlgoHNewton_base_fun<Real>
  {
    PFUN   m_fun;
    PFUN_D m_fun_D;

  public:
    explicit AlgoHNewton_fun( PFUN f, PFUN_D Df ) : m_fun( f ), m_fun_D( Df ) {}
    Real
    eval( Real x ) const override
    {
      return m_fun( x );
    };
    Real
    D( Real x ) const override
    {
      return m_fun_D( x );
    };
  };
#endif

  //!
  //!  \class AlgoHNewton
  //!  \brief Class for solving \f$ f(x)=0 \f$ without the usew of derivative
  //!
  //!  ## Usage simple example:
  //!
  //!  To use this class, first wrap your function in a derived class. For
  //!  instance, for the function \f$ f(x) = x^2 - 2 \f$, you can define:
  //!
  //!  \code{cpp}
  //!  class Fun1 : public AlgoHNewton_base_fun<double> {
  //!  public:
  //!    double eval(double x) const override { return x*x - 2; }
  //!    double D   (double x) const override { return 2*x; }
  //!  };
  //!  \endcode
  //!
  //!  Next, instantiate the function and the solver. Then, call the desired
  //!  method to find the root:
  //!
  //!  \code{cpp}
  //!  HNewton<real_type> solver;
  //!  Fun1      f;
  //!  real_type a=-1, b=2;
  //!  real_type x_solution = solver.eval2(a,b,f);
  //!  \endcode
  //!
  //!  If the method converges, `x_solution` will contain the computed solution.
  //!
  //!  ## Usage Detailed Example:
  //!
  //!  To create a custom function, derive from this class and implement the
  //!  required methods. Here is an example for the function \f$ f(x) = x^2 - 2
  //!  \f$:
  //!
  //!  ### Step 1: including headers
  //!
  //!  To use the `HNewton` class, include the necessary headers:
  //!
  //!  \code{cpp}
  //!  #include "Utils_HNewton.hh"
  //!  #include "Utils_fmt.hh" // For formatted output
  //!  \endcode
  //!
  //!  ### defining your function
  //!
  //!  Define the function for which you want to find the root. The function
  //!  must take a single argument (the variable for which you're solving) and
  //!  return a value.
  //!
  //!  For example, to find the root of \f$ \sin(x) - \frac{x}{2} = 0 \f$:
  //!
  //!  \code{cpp}
  //!  class Fun1 : public AlgoHNewton_base_fun<double> {
  //!  public:
  //!    double eval(double x) const override { return sin(x) - x / 2; }
  //!    double D   (double x) const override { return cos(x) - 0.5; }
  //!  };
  //!  static Fun1 myFunction;
  //!  \endcode
  //!
  //!  ## Calling the solver
  //!
  //!  To find a root, call the `eval` method with the initial guesses for the
  //!  root (interval `[a, b]`).
  //!  Here's an example of how to set up a program to find a root:
  //!
  //!  \code{cpp}
  //!  int
  //!  main() {
  //!    // Create an instance of the solver
  //!    HNewton<real_type> solver;
  //!
  //!    // Define the interval [a, b]
  //!    real_type a = 0.0; // lower bound
  //!    real_type b = 2.0; // upper bound
  //!
  //!    // Solve for the root
  //!    real_type root = solver.eval(a, b, myFunction);
  //!
  //!    // Print the results
  //!    cout << "Root found: " << root << endl;
  //!    cout << "f(root) = " << myFunction.eval(root) << endl;
  //!
  //!    return 0;
  //!  }
  //!  \endcode
  //!
  //!  ### Analyzing results
  //!
  //!  After calling the `eval` method, you can check the number of iterations,
  //!  number of function evaluations, and whether the algorithm converged:
  //!
  //!  \code{cpp}
  //!  cout << "Iterations: " << solver.used_iter() << endl;
  //!  cout << "Function Evaluations: " << solver.num_fun_eval() << endl;
  //!  cout << "Function Derivative Evaluations: " << solver.num_fun_D_eval() <<
  //!  endl; cout << "Converged: " << (solver.converged() ? "Yes" : "No") <<
  //!  endl;
  //!  \endcode
  //!
  template <typename Real>
  class AlgoHNewton
  {
    using Integer = int;

    Real m_tolerance{ pow( machine_eps<Real>(), Real( 2. / 3. ) ) };
    bool m_converged{ false };

    Real m_a{ 0 }, m_fa{ 0 };
    Real m_b{ 0 }, m_fb{ 0 };
    Real m_c{ 0 }, m_fc{ 0 };
    Real m_d{ 0 }, m_fd{ 0 };
    Real m_ba{ 0 };
    Real m_kappa{ 0.05 };

    AlgoHNewton_base_fun<Real> const * m_function{ nullptr };

    Integer m_max_iteration{ 200 };  // max number of iterations

    mutable Integer m_iteration_count{ 0 };  // explore iteration counter
    mutable Integer m_fun_evaluation_count{ 0 };
    mutable Integer m_fun_D_evaluation_count{ 0 };

    Real
    evaluate( Real x ) const
    {
      ++m_fun_evaluation_count;
      return m_function->eval( x );
    };
    Real
    evaluate_D( Real x ) const
    {
      ++m_fun_D_evaluation_count;
      return m_function->D( x );
    };

    Real eval();
    Real eval( Real a, Real b );
    Real eval( Real a, Real b, Real amin, Real bmax );

    void set_tolerance( Real tol );

    Real p_zero2() const;
    Real invp_zero2() const;
#if 0
    Real invp_zero3() const;
#endif

  public:
    AlgoHNewton()  = default;
    ~AlgoHNewton() = default;

    //!
    //! Find the solution for a function wrapped in the class
    //! `AlgoHNewton_base_fun<Real>` starting from guess interval `[a,b]`
    //!
    //! \param a    lower bound search interval
    //! \param b    upper bound search interval
    //! \param fun  the pointer to base class `AlgoHNewton_base_fun<Real>`
    //! wrapping the user function
    //!
    Real
    eval( Real a, Real b, AlgoHNewton_base_fun<Real> const * fun )
    {
      m_function = fun;
      return this->eval( a, b );
    }

    //!
    //! Find the solution for a function wrapped in the class
    //! `AlgoHNewton_base_fun<Real>` starting from guess interval `[a,b]`
    //!
    //! \param a    guess interval lower bound
    //! \param b    guess interval upper bound
    //! \param amin lower bound search interval
    //! \param bmax upper bound search interval
    //! \param fun  the pointer to base class `AlgoHNewton_base_fun<Real>`
    //! wrapping the user function
    //!
    Real
    eval( Real a, Real b, Real amin, Real bmax, AlgoHNewton_base_fun<Real> const * fun )
    {
      m_function = fun;
      return this->eval( a, b, amin, bmax );
    }

    //!
    //! Fix the maximum number of iteration.
    //!
    //! \param mit the maximum number of iteration
    //!
    void set_max_iterations( Integer mit );

    //!
    //! \return the number of iterations used in the last computation
    //!
    Integer
    used_iter() const
    {
      return m_iteration_count;
    }

    //!
    //! \return the number of evaluation used in the last computation
    //!
    Integer
    num_fun_eval() const
    {
      return m_fun_evaluation_count;
    }

    //!
    //! \return the number of evaluation used in the last computation
    //!
    Integer
    num_fun_D_eval() const
    {
      return m_fun_D_evaluation_count;
    }

    //!
    //! \return the tolerance set for computation
    //!
    Real
    tolerance() const
    {
      return m_tolerance;
    }

    //!
    //! \return true if the last computation was successfull
    //!
    bool
    converged() const
    {
      return m_converged;
    }

    Real
    a() const
    {
      return m_a;
    }
    Real
    b() const
    {
      return m_b;
    }
    Real
    fa() const
    {
      return m_fa;
    }
    Real
    fb() const
    {
      return m_fb;
    }
  };

#ifndef UTILS_OS_WINDOWS
  extern template class AlgoHNewton<float>;
  extern template class AlgoHNewton<double>;
#endif

  /*! @} */

}  // namespace Utils

#endif

//
// EOF: Utils_AlgoHNewton.hh
//
