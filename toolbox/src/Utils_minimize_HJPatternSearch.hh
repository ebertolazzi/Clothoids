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
 |      email: enrico.bertolazzi\unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: Utils_HJPatternSearch.hh
//

#pragma once

#ifndef UTILS_HJ_PATTERN_SEARCH_dot_HH
#define UTILS_HJ_PATTERN_SEARCH_dot_HH

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

#include "Utils.hh"
#include "Utils_eigen.hh"

namespace Utils
{

  /*!
   * \addtogroup Minimize
   * @{
   */

  //!
  //!  \brief Class for implementing the Hooke-Jeeves Pattern Search Algorithm.
  //!
  //!  The `HJPatternSearch` class provides an implementation of the
  //!  Hooke-Jeeves Pattern Search Algorithm, which is a direct search method
  //!  for optimization problems. This method does not require the computation
  //!  of gradients and is suitable for minimizing functions in arbitrary
  //!  dimensions.
  //!
  //!  \tparam Real The data type used for the optimization (e.g., float,
  //!  double).
  //!
  //!  ## References
  //!
  //!  The implementation is based on the following references:
  //!
  //!  - R. Hooke, T. A. Jeeves, "Direct Search" Solution of Numerical and
  //!  Statistical Problems,
  //!    Westinghouse Research Laboratories, Pittsburgh, Pennsylvania.
  //!
  //!  - Arthur Kaupe, "Algorithm 178: Direct Search," Communications of the
  //!  ACM,
  //!    Volume 6, Number 6, June 1963, page 313.
  //!
  //!  - M. Bell, Malcolm Pike, "Remark on Algorithm 178: Direct Search,"
  //!    Communications of the ACM, Volume 9, Number 9, September 1966, page
  //!    684.
  //!
  //!  - F.K. Tomlin, L.B. Smith, "Remark on Algorithm 178: Direct Search,"
  //!    Communications of the ACM, Volume 12, Number 11, November 1969, pages
  //!    637-638.
  //!
  //!
  //!  ## Using the **HJPatternSearch** Class for Function Minimization
  //!
  //!  The `HJPatternSearch` class is a powerful tool for finding the minimum of
  //!  a function in an arbitrary number of dimensions using the Hooke-Jeeves
  //!  Pattern Search algorithm. This tutorial will guide you through the
  //!  process of using this class effectively.
  //!
  //!  ### Prerequisites
  //!
  //!  Before you begin, ensure you have the following:
  //!
  //!  1. **C++ Compiler**: A C++ compiler that supports C++11 or later.
  //!  2. **Eigen Library**: Make sure the Eigen library is included in your
  //!  project for matrix operations.
  //!  3. **Your Implementation of `HJPatternSearch`**: Ensure you have the
  //!  implementation of the `HJPatternSearch` class available in your project.
  //!
  //!  ### Step 1: Include Necessary Headers
  //!
  //!  Start by including the headers for the `HJPatternSearch` class and any
  //!  necessary libraries.
  //!
  //!  \code{cpp}
  //!  #include "Utils_HJPatternSearch.hh" // Adjust the path according to your
  //!  project structure #include "Utils_fmt.hh" // For formatted output
  //!  (optional)
  //!
  //!  #include <cmath>
  //!  #include <iostream>
  //!  #include <functional>
  //!
  //!  using namespace std;
  //!  using Utils::HJPatternSearch; // Adjust according to your namespace
  //!  using real_type = double; // Define the real type as double for precision
  //!  \endcode
  //!
  //!  ### Step 2: Define the Objective Function
  //!
  //!  You need to define the function you want to minimize.
  //!  The function should take an array of doubles as input and return a double
  //!  representing the function value.
  //!
  //!  Here’s an example of a simple quadratic function:
  //!
  //!  \code{cpp}
  //!  static real_type fun3(real_type const X[]) {
  //!    real_type x = X[0];
  //!    real_type y = X[1];
  //!    return 100 * pow(y - x * x, 2) + pow(1 - x, 2); // Rosenbrock function
  //!  }
  //!  \endcode
  //!
  //!  ### Step 3: Set Up the Solver
  //!
  //!  Create a function to set up and run the **HJPatternSearch** solver.
  //!
  //!  \code{cpp}
  //!  template <typename FUN>
  //!  void do_solve(FUN f, real_type const X0[], real_type delta) {
  //!    Utils::Console console(&cout, 4); // Create a console for logging
  //!    HJPatternSearch<real_type> solver("HJPatternSearch"); // Create the
  //!    solver solver.setup(2, f, &console); // Setup for 2D optimization
  //!    solver.set_tolerance(1e-20); // Set the desired tolerance for
  //!    convergence solver.run(X0, delta); // Run the optimization from initial
  //!    guess X0 with step size delta
  //!
  //!    real_type X[2]; // Array to store the result
  //!    solver.get_last_solution(X); // Retrieve the best solution found
  //!    fmt::print("X={}, Y={}\n", X[0], X[1]); // Print the results
  //!  }
  //!  \endcode
  //!
  //!  ### Step 4: Main Function to Execute the Solver
  //!
  //!  In your main function, define your initial guess and the step size for
  //!  the search. Then call the do_solve function with your defined objective
  //!  function.
  //!
  //!  \code{cpp}
  //!  int main() {
  //!    real_type X0[2]{-1, 1}; // Initial guess for [x, y]
  //!    real_type delta = 0.1; // Initial step size for search
  //!    std::function<real_type(real_type const[])> F(fun3); // Create a
  //!    function wrapper do_solve(F, X0, delta); // Execute the optimization
  //!
  //!    cout << "\nAll Done Folks!\n"; // Indicate completion
  //!    return 0; // Return success
  //!  }
  //!  \endcode
  //!
  //!  ### Step 5: Compile and Run
  //!
  //!  Compile your program using a suitable C++ compiler.
  //!  For example, using g++ you can compile with:
  //!
  //!  \code{bash}
  //!  g++ -std=c++11 -o HJPatternSearchExample your_file.cpp
  //!  \endcode
  //!
  //!  Run the executable:
  //!
  //!  \code{bash}
  //!  ./HJPatternSearchExample
  //!  \endcode
  //!
  template <typename Real>
  class HJPatternSearch
  {
    using Vec_t   = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
    using MapVec  = Eigen::Map<Vec_t>;
    using integer = int;
    using HJFunc  = std::function<Real( Real const[] )>;

  private:
    string const m_name;  //!< Name of the pattern search instance.

    Malloc<Real> m_base_value;

    HJFunc          m_fun;                 //!< Handle to the value function to minimize.
    Console const * m_console{ nullptr };  //!< Pointer to the message stream class for logging.

    Real m_rho{ Real( 0.9 ) };  //!< Stencil step decreasing factor (must be 0 < rho < 1).
    Real m_h{ Real( 0.1 ) };
    bool m_stencil_failure{ false };  // stencil failure flag - used to shrink
                                      // h, stencil_failure = true means failure

    integer m_dim{ 0 };
    MapVec  m_x_old{ nullptr, 0 };
    MapVec  m_x_best{ nullptr, 0 };
    Real    m_f_old;
    Real    m_f_best;

    MapVec m_dir{ nullptr, 0 };          //!< Current search direction.
    MapVec m_search_sign{ nullptr, 0 };  //!< vector to keep in memory the direction of function value descent
                                         //!< from the previous iteration in each direction j

    MapVec m_p{ nullptr, 0 };
    MapVec m_p1{ nullptr, 0 };
    MapVec m_new_x{ nullptr, 0 };

    Real    m_tolerance{ Real( 1e-8 ) };   //!< Tolerance for convergence.
    integer m_max_iteration{ 500 };        //!< Maximum number of iterations.
    integer m_max_fun_evaluation{ 1000 };  //!< Maximum number of function evaluations.
    integer m_max_num_stagnation{ 10 };    //!< Maximum number of stagnations before stopping.
    integer m_verbose{ 1 };                //!< Verbosity level for output.

    integer m_iteration_count;       //!< Current iteration counter.
    integer m_fun_evaluation_count;  //!< Current function evaluation counter.

    void allocate( integer n );

    //! Evaluates the objective function.
    Real
    eval_function( MapVec const & x ) const
    {
      ++m_fun_evaluation;
      return m_fun( x.data() );
    }

    mutable integer m_fun_evaluation;

  public:
    //!
    //! \brief Constructs a Hooke-Jeeves Pattern Search instance.
    //!
    //! \param name The name of this instance for identification.
    //!
    explicit HJPatternSearch( string_view name ) : m_name( name ), m_base_value( string( "HJ_" ) + string( name ) ) {}

    //!
    //! \brief Retrieves the name of the pattern search instance.
    //!
    //! \return The name of the instance.
    //!
    string_view
    name( void ) const
    {
      return m_name;
    }

    //!
    //! \brief Sets up the pattern search with the specified dimension and
    //! objective function.
    //!
    //! \param dim The dimension of the problem.
    //! \param fun A reference to the objective function to minimize.
    //! \param console Optional pointer to a console for output messages.
    //!
    void setup( integer dim, HJFunc & fun, Console const * console );

    //!
    //! \brief Changes the console used for output messages.
    //!
    //! \param console Pointer to the new console for output.
    //!
    void
    change_console( Console const * console )
    {
      m_console = console;
    }

    //!
    //! \brief Sets the verbosity level for output messages.
    //!
    //! \param level The verbosity level (0: none, 1: basic info, 2: detailed
    //! info).
    //!
    void
    set_verbose( integer level )
    {
      m_verbose = level;
    }

    //!
    //! \brief Sets the convergence tolerance.
    //!
    //! \param tol The new tolerance value.
    //!
    void set_tolerance( Real tol );

    //!
    //! \brief Sets the maximum number of iterations.
    //!
    //! \param mit The maximum iterations allowed.
    //!
    void set_max_iterations( integer mit );

    //!
    //! \brief Sets the maximum number of function evaluations.
    //!
    //! \param mfev The maximum number of function evaluations allowed.
    //!
    void set_max_fun_evaluation( integer mfev );

    //!
    //! \brief Sets the maximum number of allowed stagnation iterations.
    //!
    //! \param nstg The maximum number of allowed stagnation iterations.
    //!
    void set_max_num_stagnation( integer nstg );

    //!
    //! \brief Provides information about the current state of the algorithm.
    //!
    //! \return A string containing the current state information.
    //!
    string info() const;

    //!
    //! \brief Prints the current state information to the specified output
    //! stream.
    //!
    //! \param stream The output stream to which the information will be
    //! printed.
    //!
    void
    print_info( ostream_type & stream ) const
    {
      stream << info();
    }

    //!
    //! \brief Searches for the best nearby solution.
    //!
    //! This method performs a search around the current best solution to find a
    //! nearby minimum.
    //!
    void best_nearby();

    //!
    //! \brief Performs the main search algorithm.
    //!
    //! This method executes the Hooke-Jeeves pattern search algorithm
    //! iteratively.
    //!
    void search();

    //!
    //! \brief Runs the optimization algorithm starting from the provided
    //! solution.
    //!
    //! \param x_sol Initial solution to start the optimization.
    //! \param h Initial step size for the search.
    //!
    void run( Real const x_sol[], Real h );

    //!
    //! \brief Retrieves the last solution found by the algorithm.
    //!
    //! \param x Array to store the last solution.
    //! \return The objective function value at the last solution.
    //!
    Real
    get_last_solution( Real x[] ) const
    {
      std::copy_n( m_x_best.data(), m_dim, x );
      return m_f_best;
    }
  };

  /*! @} */

}  // namespace Utils

#endif

//
// eof: Utils_HJPatternSearch.hh
//
