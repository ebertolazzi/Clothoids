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
// file: Utils_Algo748.hh
//

#pragma once

#ifndef UTILS_ALGO748_dot_HH
#define UTILS_ALGO748_dot_HH

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
  //      _    _           _____ _  _    ___
  //     / \  | | __ _  __|___  | || |  ( _ )
  //    / _ \ | |/ _` |/ _ \ / /| || |_ / _ \
  //   / ___ \| | (_| | (_) / / |__   _| (_) |
  //  /_/   \_\_|\__, |\___/_/     |_|  \___/
  //             |___/
  */

  //!
  //! \class Algo748_base_fun
  //! \brief Abstract base class for defining mathematical functions used in the zero search algorithm.
  //!
  //! This class serves as a base interface for user-defined functions that can be evaluated.
  //! It allows for the implementation the numerical method to
  //! find the solution of the one dimensional equation \f$ f(x) = 0 \f$.
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
  //!
  template <typename Real>
  class Algo748_base_fun {
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
  class Algo748_fun : public Algo748_base_fun<Real> {
    PFUN m_fun;
  public:
    explicit Algo748_fun( PFUN pfun ) : m_fun(pfun) {}
    Real eval( Real x ) const override { return m_fun(x); };
  };
  #endif

  //!
  //!  \class Algo748
  //!  \brief Class for solving \f$ f(x)=0 \f$ without the usew of derivative
  //!
  //!  \note The used algorithm is described in:
  //!        **G. E. Alefeld, Florian A Potra, Yixun Shi**,
  //!        *Algorithm 748: enclosing zeros of continuous functions*,
  //!        ACM Transactions on Mathematical Software, vol 21, N.3, 1995
  //!
  //!  ## Usage simple example:
  //!
  //!  To use this class, first wrap your function in a derived class. For instance, for the function \f$ f(x) = x^2 - 2 \f$, you can define:
  //!
  //!  \code{cpp}
  //!  class Fun1 : public Algo748_base_fun<double> {
  //!  public:
  //!    double eval(double x) const override { return x*x - 2; }
  //!  };
  //!  \endcode
  //!
  //!  Next, instantiate the function and the solver. Then, call the desired method to find the root:
  //!
  //!  \code{cpp}
  //!  Algo748<real_type> solver;
  //!  Fun1 f;
  //!  real_type a=-1,b=2;
  //!  real_type x_solution = solver.eval2(a,b,f);
  //!  \endcode
  //!
  //!  If the method converges, `x_solution` will contain the computed solution.
  //!
  //!  ## Usage Detailed Example:
  //!
  //!  To create a custom function, derive from this class and implement the required methods.
  //!  Here is an example for the function \f$ f(x) = x^2 - 2 \f$:
  //!
  //!  ### Step 1: including headers
  //!
  //!  To use the `Algo748` class, include the necessary headers:
  //!
  //!  \code{cpp}
  //!  #include "Utils_Algo748.hh"
  //!  #include "Utils_fmt.hh" // For formatted output
  //!  \endcode
  //!
  //!  ### Step 2: defining your function
  //!
  //!  Define the function for which you want to find the root. The function must
  //!  take a single argument (the variable for which you're solving) and return
  //!  a value.
  //!
  //!  For example, to find the root of \f$ \sin(x) - \frac{x}{2} = 0 \f$:
  //!
  //!  \code{cpp}
  //!  static real_type myFunction(real_type x) {
  //!    return sin(x) - x / 2;
  //!  }
  //!  \endcode
  //!
  //!  ### Step 3: Creating solver instance
  //!
  //!  Create an instance of the `Algo748` class to solve your function. The
  //!  class provides methods to evaluate functions and find roots.
  //!
  //!  ### Step 4: Calling the solver
  //!
  //!  To find a root, call the `eval` method with the initial guesses for the
  //!  root (interval `[a, b]`).
  //!  Here's an example of how to set up a program to find a root:
  //!
  //!  \code{cpp}
  //!  int
  //!  main() {
  //!    // Create an instance of the solver
  //!    Algo748<real_type> solver;
  //!
  //!    // Define the interval [a, b]
  //!    real_type a = 0.0; // lower bound
  //!    real_type b = 2.0; // upper bound
  //!
  //!    // Solve for the root
  //!    real_type root = solver.eval2(a, b, myFunction);
  //!
  //!    // Print the results
  //!    cout << "Root found: " << root << endl;
  //!    cout << "f(root) = " << myFunction(root) << endl;
  //!
  //!    return 0;
  //!  }
  //!  \endcode
  //!
  //!  ### Step 5: Using lambda functions
  //!
  //!  You can also use lambda functions to define your function inline, which
  //!  can simplify your code. Here's how to modify the previous example:
  //!
  //!  \code{cpp}
  //!  int
  //!  main() {
  //!    // Create an instance of the solver
  //!    Algo748<real_type> solver;
  //!
  //!    // Define the interval [a, b]
  //!    real_type a = 0.0; // lower bound
  //!    real_type b = 2.0; // upper bound
  //!
  //!    // Solve for the root using a lambda function
  //!    real_type root = solver.eval2(a, b, [](real_type x) { return sin(x) - x / 2; });
  //!
  //!    // Print the results
  //!    cout << "Root found: " << root << endl;
  //!    cout << "f(root) = " << sin(root) - root / 2 << endl;
  //!
  //!    return 0;
  //!  }
  //!  \endcode
  //!
  //!  ### Step 6: Advanced usage with function parameters
  //!
  //!  If your function requires additional parameters, you can wrap it in a
  //!  lambda or use `std::bind`. Here's an example using `std::bind`:
  //!
  //!  \code{cpp}
  //!  #include <functional>
  //!
  //!  static real_type myParameterizedFunction(real_type x, real_type a) {
  //!    return a * x * exp(-x);
  //!  }
  //!
  //!  int
  //!  main() {
  //!    // Create an instance of the solver
  //!    Algo748<real_type> solver;
  //!
  //!    // Define the interval [a, b]
  //!    real_type a = 0.0;
  //!    real_type b = 5.0;
  //!    real_type parameter = -1.0;
  //!
  //!    // Solve for the root using std::bind
  //!    real_type root = solver.eval2(a, b, std::bind(myParameterizedFunction, std::placeholders::_1, parameter));
  //!
  //!    // Print the results
  //!    cout << "Root found: " << root << endl;
  //!    cout << "f(root) = " << myParameterizedFunction(root, parameter) << endl;
  //!
  //!    return 0;
  //!  }
  //!  \endcode
  //!
  //!  ### Step 7: Analyzing results
  //!
  //!  After calling the `eval2` method, you can check the number of iterations,
  //!  number of function evaluations, and whether the algorithm converged:
  //!
  //!  \code{cpp}
  //!  cout << "Iterations: " << solver.used_iter() << endl;
  //!  cout << "Function Evaluations: " << solver.num_fun_eval() << endl;
  //!  cout << "Converged: " << (solver.converged() ? "Yes" : "No") << endl;
  //!  \endcode
  //!
  template <typename Real>
  class Algo748 {

    using Integer = int;

    Real m_mu{Real(0.5)};
    Real m_tolerance{pow(machine_eps<Real>(),Real(2./3.))};
    Real m_interval_shink{Real(0.025)};

    bool m_converged{false};

    Real m_a{0}, m_fa{0};
    Real m_b{0}, m_fb{0};
    Real m_c{0}, m_fc{0};
    Real m_d{0}, m_fd{0};
    Real m_e{0}, m_fe{0};

    Algo748_base_fun<Real> * m_function{nullptr};

    Integer m_max_fun_evaluation{1000}; // max number of function evaluations
    Integer m_max_iteration{200};       // max number of iterations

    mutable Integer m_iteration_count{0};    // explore iteration counter
    mutable Integer m_fun_evaluation_count{0};

    bool bracketing();
    void set_tolerance( Real tol );
    Real pzero();
    bool newton_quadratic( Integer niter, Real & c );
    Real evaluate( Real x ) { ++m_fun_evaluation_count; return m_function->eval(x); };
    bool all_different( Real a, Real b, Real c, Real d ) const;

    Real eval();
    Real eval( Real a, Real b );
    Real eval( Real a, Real b, Real amin, Real bmax );

  public:

    Algo748() = default;
    ~Algo748() = default;

    //!
    //! Find the solution for a function wrapped in the class `Algo748_base_fun<Real>`
    //! starting from guess interval `[a,b]`
    //!
    //! \param a    lower bound search interval
    //! \param b    upper bound search interval
    //! \param fun  the pointer to base class `Algo748_base_fun<Real>` wrapping the user function
    //!
    Real
    eval( Real a, Real b, Algo748_base_fun<Real> * fun ) {
      m_function = fun;
      return this->eval( a, b );
    }

    //!
    //! Find the solution for a function wrapped in the class `Algo748_base_fun<Real>`
    //! starting from guess interval `[a,b]`
    //!
    //! \param a    guess interval lower bound
    //! \param b    guess interval upper bound
    //! \param amin lower bound search interval
    //! \param bmax upper bound search interval
    //! \param fun  the pointer to base class `Algo748_base_fun<Real>` wrapping the user function
    //!
    Real
    eval( Real a, Real b, Real amin, Real bmax, Algo748_base_fun<Real> * fun ) {
      m_function = fun;
      return this->eval( a, b, amin, bmax );
    }

    //!
    //! Find the solution for a function stored in `pfun`
    //! starting from guess interval `[a,b]`
    //!
    //! \param a    lower bound search interval
    //! \param b    upper bound search interval
    //! \param pfun object storing the function
    //!
    template <typename PFUN>
    Real
    eval2( Real a, Real b, PFUN pfun ) {
      Algo748_fun<Real,PFUN> fun( pfun );
      m_function = &fun;
      return this->eval( a, b );
    }

    //!
    //! Find the solution for a function stored in `pfun`
    //! starting from guess interval `[a,b]`
    //!
    //! \param a    guess interval lower bound
    //! \param b    guess interval upper bound
    //! \param amin lower bound search interval
    //! \param bmax upper bound search interval
    //! \param pfun object storing the function
    //!
    template <typename PFUN>
    Real
    eval2( Real a, Real b, Real amin, Real bmax, PFUN pfun ) {
      Algo748_fun<Real,PFUN> fun( pfun );
      m_function = &fun;
      return this->eval( a, b, amin, bmax );
    }

    //!
    //! Find the solution for a function wrapped into `pfun`
    //! starting from guess interval `[a,b]`
    //!
    //! \param a    lower bound search interval
    //! \param b    upper bound search interval
    //! \param fa   the value \f$ f(a) \f$
    //! \param fb   the value \f$ f(b) \f$
    //! \param pfun object storing the function
    //!
    template <typename PFUN>
    Real
    eval3( Real a, Real b, Real fa, Real fb, PFUN pfun ) {
      Algo748_fun<Real,PFUN> fun( pfun );
      m_function             = &fun;
      m_iteration_count      = 0;
      m_fun_evaluation_count = 0;
      m_a = a; m_fa = fa;
      m_b = b; m_fb = fb;
      return eval();
    }

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
  extern template class Algo748<float>;
  extern template class Algo748<double>;
  #endif

  /*! @} */

}

#endif

//
// EOF: Utils_Algo748.hh
//
