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
// file: Utils_zeros.hh
//

#pragma once

#ifndef UTILS_ZEROS_dot_HH
#define UTILS_ZEROS_dot_HH

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
   * \addtogroup Zeros
   * @{
   */

  /*
  //   __________ ____   ___  ____
  //  |__  / ____|  _ \ / _ \/ ___|
  //    / /|  _| | |_) | | | \___ \
  //   / /_| |___|  _ <| |_| |___) |
  //  /____|_____|_| \_\\___/|____/
  */
  //!
  //! \ingroup Zeros
  //! \class Zeros_base_fun
  //! \brief Abstract base class for defining mathematical functions used in root-finding algorithms.
  //!
  //! This class serves as a base interface for user-defined functions that can be evaluated
  //! and differentiated. It allows for the implementation of various numerical methods to
  //! find the roots of the equation \f$ f(x) = 0 \f$. Users must inherit from this class and
  //! implement the virtual methods to define their specific functions.
  //!
  //! The class provides methods to evaluate the function itself as well as its first,
  //! second, and third derivatives, which are essential for many root-finding algorithms
  //! such as Newton-Raphson and Halley methods.
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
  //! class Fun1 : public Zeros_base_fun<double> {
  //! public:
  //!     double eval(double x) const override { return x*x - 2; }
  //!     double eval_D(double x) const override { return 2*x; }
  //!     double eval_DD(double x) const override { return 2; }
  //!     double eval_DDD(double x) const override { return 0; }
  //! };
  //! \endcode

  template <typename Real>
  class Zeros_base_fun {
  public:

    //!
    //! Evaluate the function \f$ f(x) \f$
    //!
    //! \param[in] x the point to evaluate \f$ f(x) \f$
    //! \return the value of \f$ f(x) \f$
    //!
    virtual Real eval( Real x ) const = 0;

    //!
    //! Evaluate the first derivative of \f$ f(x) \f$
    //!
    //! \param[in] x the point to evaluate \f$ f'(x) \f$
    //! \return the value of \f$ f'(x) \f$
    //!
    virtual Real eval_D( Real x ) const = 0;

    //!
    //! Evaluate the second derivative of \f$ f(x) \f$
    //!
    //! \param[in] x the point to evaluate \f$ f''(x) \f$
    //! \return the value of \f$ f''(x) \f$
    //!
    virtual Real eval_DD( Real x ) const = 0;

    //!
    //! Evaluate the third derivative of \f$ f(x) \f$
    //!
    //! \param[in] x the point to evaluate \f$ f'''(x) \f$
    //! \return the value of \f$ f'''(x) \f$
    //!
    virtual Real eval_DDD( Real x ) const = 0;

    //!
    //! Evaluate the function \f$ f(x) \f$
    //!
    //! This operator allows for a more intuitive usage of the function object,
    //! enabling the evaluation of the function using the call operator.
    //!
    //! \param[in] x the point to evaluate \f$ f(x) \f$
    //! \return the value of \f$ f(x) \f$
    //!
    Real operator () ( Real x ) const { return this->eval(x); }

    //!
    //! Evaluate the first derivative of \f$ f(x) \f$
    //!
    //! This operator provides a convenient way to evaluate the first derivative
    //! using the call operator.
    //!
    //! \param[in] x the point to evaluate \f$ f'(x) \f$
    //! \return the value of \f$ f'(x) \f$
    //!
    Real D( Real x ) const { return this->eval_D(x); }

    //!
    //! Evaluate the second derivative of \f$ f(x) \f$
    //!
    //! This operator allows for an easy evaluation of the second derivative
    //! using the call operator.
    //!
    //! \param[in] x the point to evaluate \f$ f''(x) \f$
    //! \return the value of \f$ f''(x) \f$
    //!
    Real DD( Real x ) const { return this->eval_DD(x); }

    //!
    //! Evaluate the third derivative of \f$ f(x) \f$
    //!
    //! This operator enables straightforward evaluation of the third derivative
    //! using the call operator.
    //!
    //! \param[in] x the point to evaluate \f$ f'''(x) \f$
    //! \return the value of \f$ f'''(x) \f$
    //!
    Real DDD( Real x ) const { return this->eval_DDD(x); }
  };

  //!
  //! \class Zeros
  //! \brief Class for solving the equation \f$ f(x) = 0 \f$ using various numerical methods.
  //!
  //! This class implements multiple solvers to find the roots of a given function. The available methods include:
  //!
  //! - **Newton-Raphson Method**: A widely used iterative method for finding successively better approximations to the roots of a real-valued function.
  //!   [Learn more](https://en.wikipedia.org/wiki/Newton%27s_method).
  //! - **Chebyshev Method**: A higher-order root-finding method that offers faster convergence compared to the Newton-Raphson method.
  //! - **Halley Method**: An iterative method that is a generalization of the Newton-Raphson method and provides faster convergence.
  //!   [Learn more](https://en.wikipedia.org/wiki/Halley%27s_method).
  //! - **Methods by Juan Luis Varona**: A series of methods developed for enhanced convergence properties:
  //!   - Order 4 method
  //!   - Order 8 method
  //!   - Order 16 method
  //!   - Order 32 method
  //!
  //! For a detailed exploration of these methods, refer to the paper:
  //! - *An Optimal Thirty-Second-Order Iterative Method for Solving Nonlinear Equations and a Conjecture*,
  //!   **Juan Luis Varona**, Qualitative Theory of Dynamical Systems (2022).
  //!   [Link to the paper](https://link.springer.com/article/10.1007/s12346-022-00572-3).
  //!
  //! \note This class is designed to work with user-defined functions that extend from `Utils::Zeros_base_fun`.
  //!
  //! **Usage Example**
  //!
  //! To use this class, first wrap your function in a derived class. For instance, for the function \f$ f(x) = x^2 - 2 \f$, you can define:
  //!
  //! \code{cpp}
  //! class Fun1 : public Utils::Zeros_base_fun<real_type> {
  //! public:
  //!     real_type eval(real_type x) const override { return x*x - 2; }
  //!     real_type eval_D(real_type x) const override { return 2*x; }
  //!     real_type eval_DD(real_type x) const override { return 2; }
  //!     real_type eval_DDD(real_type x) const override { return 0; }
  //! };
  //! \endcode
  //!
  //! Next, instantiate the function and the solver. Then, call the desired method to find the root:
  //!
  //! \code{cpp}
  //! Zeros<real_type> solver;
  //! Fun1 f;
  //! real_type x_guess = 1.0;  // Initial guess
  //! real_type x_solution = solver.solve_Newton(x_guess, f);
  //! \endcode
  //!
  //! If the method converges, `x_solution` will contain the computed solution.

  template <typename Real>
  class Zeros {

    using Integer = int;

    Integer m_max_fun_evaluation{200};  //< max number of function evaluations
    Integer m_max_iteration{100};       //< max number of iterations
    Real    m_tolerance{pow(machine_eps<Real>(),Real(2./3.))};
    bool    m_converged{false};

    mutable Integer m_iteration_count{0};    // explore iteration counter
    mutable Integer m_fun_evaluation_count{0};

    //static Real Q( Real t ) { return 1/(1-2*t); }
    static Real Q( Real t ) { return 1+2*t; }

    static Real W( Real t, Real s ) {
      Real t2{t*t};
      return t2*(1-4*t)+(4*s+2)*t+s+1;
    }

    static Real H( Real t, Real s, Real u ) {
      Real t1  = t*t;
      Real t2  = t1*t1;
      Real t8  = s*s;
      Real t17 = s*t8;
      Real t23 = 2.0*u;
      return ( (8*u+6*t2+4)*s-(6*t8+4*(s+u+1))*t1 + 2*t8 - 4*t17 + t23 + 2 )*t +
               t1*(t8+s+u+1) + (1-3*t2+t23)*s + u - t17 + 1;
    }

    static Real J( Real t, Real s, Real u, Real v ) {
      Real t1  = s*s;
      Real t2  = t1*t1;
      Real t17 = t*t;
      Real t22 = u*u;
      Real t32 = t17*t17;
      Real t34 = t*t32;
      Real t37 = t*t17;
      Real t46 = 1+v;
      Real t65 = u+1+v;
      Real t76 = (-2*t22+u+4.0*v+2)*u;
      return (2*t-1)*(2+5*t)*u*t*t2 +
             (4*t+1)*u*s*t2+
             (u*t22-2*u*v-u-v-1)*(4*t17+3*t+1)*(t-1)-
             8.0*(t22*(t17/2.0-1.0/4.0)+u*(t17*t32-5.0/8.0*t34-3.0/4.0*t32+
             3.0/8.0*t37+3.0/4.0*t17-t/8.0-1.0/4.0)+3.0/4.0*t46*(t+1.0/2.0)*(t-2.0/3.0))*t*t1+
             4.0*(t22*(-3.0/2.0*t-1.0/4.0)+u*(t34-t32-3.0/2.0*t37+t17/4.0-t-1.0/4.0)-
             t46*(t+1.0/4.0))*s*t1+(1.0+v+t65*t17-4.0*t65*t37-3.0*t65*t32+6.0*t65*t34+t76+4.0*(1.0+v+t76)*t)*s;
    }

  public:

    Zeros() = default;
    ~Zeros() = default;

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
    //!
    //! Fix the requested tolerance for ieration stop.
    //! Stop when \f$ |f(x)| < \epsilon \f$
    //!
    //! \param tol the requested tolerance
    //!
    void set_tolerance( Real tol );

    //!
    //! Find the zero of a function wrapped in the class `Zeros_base_fun<Real>`
    //! starting from guess value `x_guess`
    //!
    //! \param x_guess starting value for iterative method
    //! \param fun     the pointer to base class `Zeros_base_fun<Real>` wrapping the user function
    //!
    Real solve_Newton( Real x_guess, Zeros_base_fun<Real> * fun );
    //!
    //! Find the zero of a function wrapped in the class `Zeros_base_fun<Real>`
    //! starting from guess value `x_guess`.
    //!
    //! \param x_guess starting value for iterative method
    //! \param fun     the pointer to base class `Zeros_base_fun<Real>` wrapping the user function
    //!
    Real solve_Chebyshev( Real x_guess, Zeros_base_fun<Real> * fun );
    //!
    //! Find the zero of a function wrapped in the class `Zeros_base_fun<Real>`
    //! starting from guess value `x_guess`
    //!
    //! \param x_guess starting value for iterative method
    //! \param fun     the pointer to base class `Zeros_base_fun<Real>` wrapping the user function
    //!
    Real solve_Halley( Real x_guess, Zeros_base_fun<Real> * fun );
    //!
    //! Find the zero of a function wrapped in the class `Zeros_base_fun<Real>`
    //! starting from guess value `x_guess`
    //!
    //! \param x_guess starting value for iterative method
    //! \param fun     the pointer to base class `Zeros_base_fun<Real>` wrapping the user function
    //!
    Real solve_Order4( Real x_guess, Zeros_base_fun<Real> * fun );
    //!
    //! Find the zero of a function wrapped in the class `Zeros_base_fun<Real>`
    //! starting from guess value `x_guess`
    //!
    //! \param x_guess starting value for iterative method
    //! \param fun     the pointer to base class `Zeros_base_fun<Real>` wrapping the user function
    //!
    Real solve_Order8( Real x_guess, Zeros_base_fun<Real> * fun );
    //!
    //! Find the zero of a function wrapped in the class `Zeros_base_fun<Real>`
    //! starting from guess value `x_guess`
    //!
    //! \param x_guess starting value for iterative method
    //! \param fun     the pointer to base class `Zeros_base_fun<Real>` wrapping the user function
    //!
    Real solve_Order16( Real x_guess, Zeros_base_fun<Real> * fun );
    //!
    //! Find the zero of a function wrapped in the class `Zeros_base_fun<Real>`
    //! starting from guess value `x_guess`
    //!
    //! \param x_guess starting value for iterative method
    //! \param fun     the pointer to base class `Zeros_base_fun<Real>` wrapping the user function
    //!
    Real solve_Order32( Real x_guess, Zeros_base_fun<Real> * fun );

    //!
    //! \return the number of iterations used in the last computation
    //!
    Integer used_iter() const { return m_iteration_count; }
    //!
    //! \return the number of evaluation used in the last computation
    //!
    Integer num_fun_eval() const { return m_fun_evaluation_count; }
    //!
    //! \return the tolerance set for computation
    //!
    Real tolerance() const { return m_tolerance; }
    //!
    //! \return true if the last computation was successfull
    //!
    bool converged() const { return m_converged; }

  };

  #ifndef UTILS_OS_WINDOWS
  extern template class Zeros<float>;
  extern template class Zeros<double>;
  #endif

}

#endif

//
// eof: Utils_zeros.hh
//
