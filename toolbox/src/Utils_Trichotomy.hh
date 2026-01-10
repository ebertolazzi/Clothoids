/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2022-2026                                                 |
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

/**
 * \file Utils_Trichotomy.hh
 * \brief Implementation of the trichotomy method for 1D minimization
 *
 * This file contains the complete implementation of the trichotomy algorithm,
 * a zero-order optimization method for finding the minimum of unimodal
 * one-dimensional functions.
 *
 * \author Enrico Bertolazzi
 * \date 2022-2026
 */

#pragma once

#ifndef UTILS_TRICHOTOMY_dot_HH
#define UTILS_TRICHOTOMY_dot_HH

#include "Utils.hh"

namespace Utils
{

  using std::abs;
  using std::isfinite;
  using std::pow;

  /**
   * \defgroup Minimization 1D Minimization Algorithms
   * \brief Algorithms for one-dimensional function minimization
   *
   * This group contains implementations of algorithms for finding
   * the minimum of one-dimensional functions. The algorithms are designed
   * to work with unimodal functions (functions with a single minimum
   * in the search interval).
   *
   * @{
   */

  /**
   * \class Trichotomy_base_fun
   * \brief Abstract base class for defining mathematical functions
   *
   * This class provides an abstract interface for one-dimensional functions
   * that can be minimized using the trichotomy algorithm. Users must derive
   * from this class and implement the pure virtual method `eval`.
   *
   * \tparam Real Numeric type for input and output (float, double, long double)
   *
   * \par Design Rationale:
   * Using an abstract base class allows:
   * - Polymorphism: handle different function types through a common interface
   * - Extensibility: easily add new function types
   * - Type safety: compile-time verification of data types
   *
   * \par Usage Example:
   * \code{.cpp}
   * // Define a quadratic function: f(x) = (x-2)² + 3
   * class QuadraticFunction : public Trichotomy_base_fun<double> {
   * public:
   *     double eval(double x) const override {
   *         return (x - 2.0) * (x - 2.0) + 3.0;
   *     }
   * };
   *
   * // Use the function
   * QuadraticFunction func;
   * double result = func(1.5);  // Uses operator()
   * double direct = func.eval(1.5);  // Direct call
   * \endcode
   *
   * \see Trichotomy for the minimization algorithm
   */
  template <typename Real> class Trichotomy_base_fun
  {
  public:
    /**
     * \brief Virtual destructor
     *
     * The virtual destructor is essential to ensure proper destruction
     * of derived objects when managed through base class pointers.
     */
    virtual ~Trichotomy_base_fun() = default;

    /**
     * \brief Evaluate the function at the specified point
     *
     * Pure virtual method that must be implemented by derived classes
     * to define the behavior of the function \f$ f(x) \f$.
     *
     * \param x Point at which to evaluate the function
     * \return The value \f$ f(x) \f$
     *
     * \note Implementations should properly handle:
     *       - Input values that might cause overflow/underflow
     *       - Discontinuity points (if applicable)
     *       - Special values (NaN, infinity)
     */
    virtual Real eval( Real x ) const = 0;

    /**
     * \brief Function call operator
     *
     * Provides a more natural syntax for function evaluation,
     * allowing the object to be used as if it were a function.
     *
     * \param x Point at which to evaluate the function
     * \return The value \f$ f(x) \f$
     *
     * \par Example:
     * \code{.cpp}
     * QuadraticFunction f;
     * double y = f(2.5);  // Equivalent to f.eval(2.5)
     * \endcode
     */
    Real operator()( Real x ) const { return this->eval( x ); }
  };

  /**
   * \brief Internal wrapper class for function objects and function pointers
   *
   * This class allows the use of lambdas, free functions, and function objects
   * with the trichotomy algorithm without requiring the user to explicitly
   * derive from Trichotomy_base_fun.
   *
   * \tparam Real Floating-point type
   * \tparam PFUN Type of the function object or function pointer
   *
   * \internal
   */
  template <typename Real, typename PFUN> class Trichotomy_fun : public Trichotomy_base_fun<Real>
  {
    PFUN m_fun;  ///< Oggetto funzione o puntatore wrappato

  public:
    /**
     * \brief Construct a wrapper around the provided function
     * \param pfun Function object or pointer to wrap
     */
    explicit Trichotomy_fun( PFUN pfun ) : m_fun( pfun ) {}

    /**
     * \brief Evaluate the wrapped function
     * \param x Evaluation point
     * \return Function value at x
     */
    Real eval( Real x ) const override { return m_fun( x ); }
  };

  /**
   * \class Trichotomy
   * \brief Implementation of the trichotomy method for 1D minimization
   *
   * This class implements the trichotomy method, a zero-order optimization
   * algorithm for finding the minimum of unimodal one-dimensional functions.
   *
   * \par Algorithm Description:
   * The trichotomy algorithm is a generalization of the golden section method
   * that divides the search interval into three parts instead of two. At each iteration:
   *
   * 1. **Bracketing**: The current interval [a, b] is divided using internal points
   * 2. **Evaluation**: The function is evaluated at the division points
   * 3. **Reduction**: The interval is reduced based on function values
   * 4. **Convergence**: The process continues until the interval is sufficiently small
   *
   * \par Complexity:
   * - Convergence: Superlinear, order \f$ \phi \approx 1.839 \f$
   * - Evaluations per iteration: 2-4 (adaptive)
   * - Reduction factor: Approximately 1/3 per iteration
   *
   * \par Advantages:
   * - Does not require derivatives (zero-order method)
   * - Guaranteed convergence for unimodal functions
   * - Faster than golden section method
   * - Robust against numerical noise
   *
   * \par Limitations:
   * - Requires unimodal functions (single minimum in the interval)
   * - Not suitable for multimodal functions
   * - Degraded performance with very flat functions
   *
   * \tparam Real Floating-point type (float, double, long double)
   *
   * \par Reference:
   * **A new zero-order 1-D optimization algorithm: trichotomy method**,\n
   * by **Alena Antonova, Olga Ibryaeva**,\n
   * [arXiv:1903.07117](https://doi.org/10.48550/arXiv.1903.07117)
   *
   * \par Complete Example:
   * \code{.cpp}
   * #include "Utils_Trichotomy.hh"
   * #include <iostream>
   * #include <cmath>
   *
   * int main() {
   *     // Minimize f(x) = (x-2)² + 3 on interval [-5, 10]
   *     auto quadratic = [](double x) {
   *         return (x - 2.0) * (x - 2.0) + 3.0;
   *     };
   *
   *     Utils::Trichotomy<double> solver;
   *     solver.set_tolerance(1e-6);
   *     solver.set_max_iterations(50);
   *
   *     double a = -5.0, b = 10.0;
   *     double x_min = solver.eval2(a, b, quadratic);
   *
   *     if (solver.converged()) {
   *         double a_final, b_final;
   *         solver.get_interval(a_final, b_final);
   *
   *         std::cout << "Minimum found at x = " << x_min << "\n";
   *         std::cout << "Minimum value: " << solver.min_value() << "\n";
   *         std::cout << "Iterations: " << solver.used_iter() << "\n";
   *         std::cout << "Evaluations: " << solver.num_fun_eval() << "\n";
   *         std::cout << "Final interval: [" << a_final
   *                   << ", " << b_final << "]\n";
   *     } else {
   *         std::cout << "Convergence not reached\n";
   *     }
   *
   *     return 0;
   * }
   * \endcode
   *
   * \see Trichotomy_base_fun for defining custom functions
   */
  template <typename Real> class Trichotomy
  {
  public:
    using Integer = int;  ///< Tipo intero per conteggi

    /**
     * \brief Default constructor
     *
     * Initializes the solver with conservative default parameters:
     * - Maximum iterations: 100
     * - Maximum evaluations: 1000
     * - Tolerance: \f$ \epsilon^{2/3} \f$ where \f$ \epsilon \f$ is machine epsilon
     *
     * The default tolerance \f$ \epsilon^{2/3} \f$ balances precision and speed,
     * providing approximately 10-11 significant digits for double precision.
     */
    Trichotomy() = default;

    /**
     * \brief Destructor
     */
    ~Trichotomy() = default;
    // =================================================================
    // Public interface methods
    // =================================================================

    /**
     * \brief Find minimum in interval [a, b]
     *
     * This method searches for the minimum of a function in the specified interval
     * using a pointer to an object derived from Trichotomy_base_fun.
     *
     * \param a Lower bound of the search interval
     * \param b Upper bound of the search interval
     * \param fun Pointer to the function to minimize
     * \return Estimated location of the minimum \f$ x^* \f$
     *
     * \pre \f$ a < b \f$ (the interval must be valid)
     * \pre The function must be unimodal on [a, b]
     * \pre `fun` must be a valid non-null pointer
     *
     * \post If converged() == true, the result is in the final interval
     * \post The returned value approximately minimizes the function
     *
     * \note After the call, use converged() to verify success
     *
     * \par Example:
     * \code{.cpp}
     * class MyFunction : public Trichotomy_base_fun<double> {
     * public:
     *     double eval(double x) const override { return x*x - 4*x + 7; }
     * };
     *
     * MyFunction func;
     * Trichotomy<double> solver;
     * double x_min = solver.eval(-10.0, 10.0, &func);
     * \endcode
     *
     * \see eval2() to use lambdas or free functions
     * \see converged() to verify convergence
     * \see get_interval() to obtain the final interval
     */
    Real eval( Real a, Real b, Trichotomy_base_fun<Real> * fun )
    {
      m_function = fun;
      return this->eval_impl( a, b );
    }

    /**
     * \brief Find minimum in interval [a, b] using a function object
     *
     * Template version that accepts any callable object: lambda, functor,
     * function pointer. More flexible than eval() but with the same semantics.
     *
     * \tparam PFUN Type of the function object (automatically deduced)
     * \param a Lower bound of the search interval
     * \param b Upper bound of the search interval
     * \param pfun Callable function object with signature `Real(Real)`
     * \return Estimated location of the minimum
     *
     * \pre \f$ a < b \f$
     * \pre The function must be unimodal on [a, b]
     * \pre `pfun` must be callable as `Real pfun(Real)`
     *
     * \par Examples:
     * \code{.cpp}
     * Trichotomy<double> solver;
     *
     * // With lambda
     * auto f1 = [](double x) { return x*x - 4*x + 7; };
     * double x1 = solver.eval2(-10.0, 10.0, f1);
     *
     * // With free function
     * double my_func(double x) { return std::sin(x); }
     * double x2 = solver.eval2(0.0, 3.14159, my_func);
     *
     * // With functor
     * struct Polynomial {
     *     double operator()(double x) const { return x*x*x - 2*x + 1; }
     * };
     * double x3 = solver.eval2(-2.0, 2.0, Polynomial{});
     * \endcode
     *
     * \see eval() for the base class version
     */
    template <typename PFUN> Real eval2( Real a, Real b, PFUN pfun )
    {
      Trichotomy_fun<Real, PFUN> fun( pfun );
      m_function = &fun;
      return this->eval_impl( a, b );
    }

    /**
     * \brief Search for minimum with automatic interval expansion
     *
     * This method starts with an interval centered at `x` with radius `delta`,
     * then automatically expands the interval if necessary to bracket the minimum.
     *
     * \param x Center of the initial search interval
     * \param delta Radius of the initial interval (must be > 0)
     * \param fun Pointer to the function to minimize
     * \return Estimated location of the minimum
     *
     * \pre \f$ \delta > 0 \f$
     * \pre `fun` must be a valid pointer
     * \pre The function must have a reachable minimum
     *
     * \post The final interval contains the minimum
     *
     * \par Expansion Strategy:
     * 1. Start with interval [x-delta, x+delta]
     * 2. Evaluate the function at three points: x-delta, x, x+delta
     * 3. If the minimum appears to be outside the interval, expand in the appropriate direction
     * 4. Continue expansion until the minimum is bracketed
     * 5. Apply the standard trichotomy algorithm
     *
     * \par When to Use search():
     * - When you have an approximate estimate of the minimum's location
     * - When the optimal interval is not known a priori
     * - For problems where the minimum might be far from the initial point
     *
     * \par Example:
     * \code{.cpp}
     * class RosenbrockSlice : public Trichotomy_base_fun<double> {
     * public:
     *     double eval(double x) const override {
     *         // Rosenbrock: (1-x)² + 100*(x²-x)²
     *         return (1-x)*(1-x) + 100*(x*x-x)*(x*x-x);
     *     }
     * };
     *
     * RosenbrockSlice func;
     * Trichotomy<double> solver;
     *
     * // Start near x=0 with small radius
     * double x_min = solver.search(0.0, 0.1, &func);
     * // The algorithm will automatically expand to find the minimum near x=1
     * \endcode
     *
     * \warning Expansion can require many evaluations if the minimum is very far away
     *
     * \see search2() for the function object version
     * \see set_max_fun_evaluation() to limit evaluations during expansion
     */
    Real search( Real x, Real delta, Trichotomy_base_fun<Real> * fun )
    {
      m_function = fun;
      return this->search_impl( x, delta );
    }

    /**
     * \brief Search with automatic expansion using a function object
     *
     * Template version of search() that accepts lambdas, functors, and function pointers.
     *
     * \tparam PFUN Type of the function object
     * \param x Center of the initial interval
     * \param delta Radius of the initial interval
     * \param pfun Callable function object
     * \return Estimated location of the minimum
     *
     * \see search() for complete details on behavior
     *
     * \par Example:
     * \code{.cpp}
     * Trichotomy<double> solver;
     *
     * // Search for minimum of sin(x) starting from x=1, radius=0.5
     * auto sine = [](double x) { return std::sin(x); };
     * double x_min = solver.search2(1.0, 0.5, sine);
     * // Will find the minimum near 3π/2
     * \endcode
     */
    template <typename PFUN> Real search2( Real x, Real delta, PFUN pfun )
    {
      Trichotomy_fun<Real, PFUN> fun( pfun );
      m_function = &fun;
      return this->search_impl( x, delta );
    }

    // =================================================================
    // Configuration methods
    // =================================================================

    /**
     * \brief Set maximum number of iterations
     *
     * Limits the number of algorithm iterations to prevent infinite loops
     * on problematic functions or to control execution time.
     *
     * \param mit Maximum number of iterations (must be > 0)
     * \throws std::invalid_argument if mit <= 0
     *
     * \par Guidelines:
     * - Smooth and well-behaved functions: 20-50 iterations
     * - Functions with noise: 50-100 iterations
     * - High precision required: 100-200 iterations
     * - Default (100): suitable for most cases
     *
     * \note Each iteration performs 2-4 function evaluations
     *
     * \par Example:
     * \code{.cpp}
     * Trichotomy<double> solver;
     * solver.set_max_iterations(50);  // Limit to 50 iterations
     * \endcode
     *
     * \see set_max_fun_evaluation() to limit evaluations instead of iterations
     */
    void set_max_iterations( Integer mit )
    {
      if ( mit <= 0 ) { throw std::invalid_argument( "Trichotomy::set_max_iterations: argument must be > 0" ); }
      m_max_iteration = mit;
    }

    /**
     * \brief Set maximum number of function evaluations
     *
     * Limits the total number of calls to the objective function. Useful when
     * function evaluation is computationally expensive.
     *
     * \param mfev Maximum number of evaluations (must be > 0)
     * \throws std::invalid_argument if mfev <= 0
     *
     * \par Relationship with Iterations:
     * On average, each iteration requires about 2-3 evaluations, so:
     * - max_evaluations ≈ 2.5 × max_iterations is a good rule
     * - The default value (1000) corresponds to about 400 iterations
     *
     * \par When to Use:
     * - Functions that require expensive simulations
     * - Functions with I/O or database access
     * - Limited computational budget
     * - Execution time control in real-time systems
     *
     * \par Example:
     * \code{.cpp}
     * Trichotomy<double> solver;
     * solver.set_max_fun_evaluation(100);  // Maximum 100 function calls
     *
     * // Expensive function we want to minimize
     * auto expensive = [](double x) {
     *     // Complex simulation...
     *     return compute_expensive_result(x);
     * };
     * \endcode
     *
     * \see num_fun_eval() to get the number of evaluations performed
     */
    void set_max_fun_evaluation( Integer mfev )
    {
      if ( mfev <= 0 ) { throw std::invalid_argument( "Trichotomy::set_max_fun_evaluation: argument must be > 0" ); }
      m_max_fun_evaluation = mfev;
    }

    /**
     * \brief Set convergence tolerance
     *
     * The tolerance defines when the algorithm considers the minimum found.
     * The algorithm terminates when the interval [a, b] satisfies: \f$ b - a < tol \f$
     *
     * \param tol Tolerance (must be > 0)
     * \throws std::invalid_argument if tol <= 0
     *
     * \par Tolerance Selection:
     * - \f$ 10^{-3} \f$: low precision, fast
     * - \f$ 10^{-6} \f$: medium precision, balanced
     * - \f$ 10^{-12} \f$: high precision (near limit for double)
     * - \f$ \epsilon^{2/3} \f$ (default): excellent compromise, about \f$ 10^{-11} \f$ for double
     *
     * \par Numerical Limitations:
     * For double precision (64-bit), the practical minimum tolerance is about \f$ 10^{-14} \f$
     * due to machine epsilon (\f$ \approx 2.22 \times 10^{-16} \f$).
     *
     * \par Considerations:
     * - Very small tolerances drastically increase iterations
     * - With numerical noise, tolerances < \f$ 10^{-12} \f$ may not improve the result
     * - For flat functions, a larger tolerance may be necessary
     *
     * \par Example:
     * \code{.cpp}
     * Trichotomy<double> solver;
     *
     * // High precision for critical applications
     * solver.set_tolerance(1e-12);
     *
     * // Low precision for fast prototypes
     * solver.set_tolerance(1e-3);
     * \endcode
     *
     * \see tolerance() to get the current tolerance
     */
    void set_tolerance( Real tol )
    {
      if ( tol <= Real( 0 ) ) { throw std::invalid_argument( "Trichotomy::set_tolerance: argument must be > 0" ); }
      m_tolerance = tol;
    }

    // =================================================================
    // Status query methods
    // =================================================================

    /**
     * \brief Get the number of iterations performed
     *
     * Returns the number of complete bracketing iterations performed
     * in the last call to eval(), eval2(), search(), or search2().
     *
     * \return Number of iterations performed
     *
     * \note This value is updated with each call to the minimization methods
     *
     * \par Interpretation:
     * - Few iterations (< 10): simple function or already tight initial interval
     * - Medium iterations (10-50): typical behavior
     * - Many iterations (> 100): complex function, tight tolerance, or convergence issues
     *
     * \par Example:
     * \code{.cpp}
     * Trichotomy<double> solver;
     * auto f = [](double x) { return x*x; };
     *
     * solver.eval2(-10.0, 10.0, f);
     * std::cout << "Iterations: " << solver.used_iter() << "\n";
     * std::cout << "Evaluations: " << solver.num_fun_eval() << "\n";
     * std::cout << "Eval/iter ratio: "
     *           << double(solver.num_fun_eval()) / solver.used_iter() << "\n";
     * \endcode
     *
     * \see num_fun_eval() for the total number of function evaluations
     */
    Integer used_iter() const { return m_num_iter_done; }

    /**
     * \brief Get the number of function evaluations performed
     *
     * Returns the total count of calls to the objective function
     * during the last algorithm execution.
     *
     * \return Number of function evaluations
     *
     * \par Typical Usage:
     * - Algorithm performance analysis
     * - Profiling expensive functions
     * - Verifying computational budget was not exceeded
     * - Comparing different algorithm configurations
     *
     * \par Typical Values:
     * For n iterations:
     * - Minimum: 3 + 2n (optimal case)
     * - Typical: 3 + 2.5n
     * - Maximum: 3 + 4n (worst case)
     *
     * \note The search() method may use extra evaluations for interval expansion
     *
     * \see set_max_fun_evaluation() to set a limit
     */
    Integer num_fun_eval() const { return m_num_fun_eval; }

    /**
     * \brief Get the current convergence tolerance
     *
     * \return Value of the tolerance used by the algorithm
     *
     * \note If the tolerance was not explicitly set, it is automatically
     *       initialized to \f$ \epsilon^{2/3} \f$ on the first algorithm execution
     *
     * \par Example:
     * \code{.cpp}
     * Trichotomy<double> solver;
     * std::cout << "Default tolerance: " << solver.tolerance() << "\n";
     *
     * solver.set_tolerance(1e-8);
     * std::cout << "New tolerance: " << solver.tolerance() << "\n";
     * \endcode
     */
    Real tolerance() const { return m_tolerance; }

    /**
     * \brief Check if the last execution converged
     *
     * Returns true if the algorithm terminated successfully by reaching
     * the required tolerance, false if it terminated due to reaching
     * iteration or evaluation limits.
     *
     * \return true if converged, false otherwise
     *
     * \par Convergence Conditions:
     * The algorithm converges when: \f$ |b - a| < \text{tolerance} \f$
     *
     * \par Reasons for Non-Convergence:
     * - Reached max_iterations
     * - Reached max_fun_evaluation
     * - Function is not unimodal on the interval
     * - Tolerance is too tight relative to numerical precision
     *
     * \warning If converged() returns false, the result may be
     *          imprecise. Consider increasing the limits or relaxing the tolerance.
     *
     * \par Example:
     * \code{.cpp}
     * Trichotomy<double> solver;
     * auto f = [](double x) { return std::abs(x); };
     *
     * double x_min = solver.eval2(-1.0, 1.0, f);
     *
     * if (solver.converged()) {
     *     std::cout << "Success! x_min = " << x_min << "\n";
     * } else {
     *     std::cerr << "Warning: convergence not reached\n";
     * }
     * \endcode
     *
     * \see used_iter() to verify how many iterations were used
     */
    bool converged() const { return m_converged; }

    /**
     * \brief Get the final interval containing the minimum
     *
     * After algorithm execution, this method returns the final interval
     * [a, b] that brackets the found minimum. If the algorithm converged,
     * the interval width will be less than the tolerance.
     *
     * \param[out] a Lower bound of the final interval
     * \param[out] b Upper bound of the final interval
     *
     * \post If converged() == true: \f$ b - a < \text{tolerance} \f$
     * \post The estimated minimum lies in [a, b]
     *
     * \par Typical Usage:
     * - Assess uncertainty in the found minimum
     * - Verify the final interval width
     * - Algorithm debugging and analysis
     * - Refinement with other methods
     *
     * \par Example:
     * \code{.cpp}
     * Trichotomy<double> solver;
     * auto f = [](double x) { return (x-3.14)*(x-3.14); };
     *
     * double x_min = solver.eval2(0.0, 10.0, f);
     *
     * double a, b;
     * solver.get_interval(a, b);
     *
     * std::cout << "Estimated minimum: " << x_min << "\n";
     * std::cout << "Interval: [" << a << ", " << b << "]\n";
     * std::cout << "Uncertainty: ±" << (b-a)/2.0 << "\n";
     * std::cout << "Width: " << (b-a) << "\n";
     *
     * // Verification
     * if (x_min >= a && x_min <= b) {
     *     std::cout << "Minimum in interval: OK\n";
     * }
     * \endcode
     *
     * \see converged() to verify if the interval respects the tolerance
     */
    void get_interval( Real & a, Real & b ) const
    {
      a = m_a;
      b = m_b;
    }

    /**
     * \brief Get the function value at the estimated minimum
     *
     * Returns \f$ f(x^*) \f$ where \f$ x^* \f$ is the minimum estimated
     * by the algorithm (the center point of the final interval).
     *
     * \return Function value at the estimated minimum
     *
     * \note This is the value already computed during the algorithm,
     *       no new function evaluation is performed
     *
     * \par Example:
     * \code{.cpp}
     * Trichotomy<double> solver;
     * auto rosenbrock = [](double x) {
     *     return (1-x)*(1-x) + 100*(x*x-x)*(x*x-x);
     * };
     *
     * double x_min = solver.search2(0.0, 0.5, rosenbrock);
     * double f_min = solver.min_value();
     *
     * std::cout << "Minimum found at x = " << x_min << "\n";
     * std::cout << "Minimum value f(x*) = " << f_min << "\n";
     *
     * // Theoretical value for Rosenbrock is 0 at x=1
     * std::cout << "Error: " << std::abs(f_min - 0.0) << "\n";
     * \endcode
     *
     * \see eval(), eval2(), search(), search2() for the methods that compute the minimum
     */
    Real min_value() const { return m_f3; }

  private:
    // =================================================================
    // Internal state variables
    // =================================================================

    Integer m_num_iter_done{ 0 };          ///< Iterations performed in last call
    Integer m_num_fun_eval{ 0 };           ///< Function evaluations performed
    Integer m_max_iteration{ 100 };        ///< Maximum iteration limit
    Integer m_max_fun_evaluation{ 1000 };  ///< Maximum evaluation limit
    Real    m_tolerance{ 0 };              ///< Convergence tolerance

    bool m_converged{ false };  ///< Convergence flag

    // Interval state and function values
    Real m_a{ 0 }, m_fa{ 0 };   ///< Left bound and f(a)
    Real m_b{ 0 }, m_fb{ 0 };   ///< Right bound and f(b)
    Real m_x1{ 0 }, m_f1{ 0 };  ///< First interior point and f(x1)
    Real m_x2{ 0 }, m_f2{ 0 };  ///< Second interior point and f(x2)
    Real m_x3{ 0 }, m_f3{ 0 };  ///< Third interior point (best estimate) and f(x3)
    Real m_x4{ 0 }, m_f4{ 0 };  ///< Fourth interior point and f(x4)
    Real m_x5{ 0 }, m_f5{ 0 };  ///< Fifth interior point and f(x5)

    Trichotomy_base_fun<Real> * m_function{ nullptr };  ///< Pointer to function to minimize

    // =================================================================
    // Private helper methods
    // =================================================================

    /**
     * \brief Evaluate function with automatic counting
     *
     * Internal wrapper that evaluates the function and increments the
     * evaluation counter. All evaluation points in the algorithm
     * must use this method to ensure accurate counting.
     *
     * \param x Point at which to evaluate the function
     * \return Function value at x
     *
     * \note This method does not validate the result.
     *       It is the user's responsibility to provide a function that
     *       returns finite values.
     *
     * \internal
     */
    Real evaluate( Real x )
    {
      ++m_num_fun_eval;
      return m_function->eval( x );
    }

    /**
     * \brief Internal implementation of eval
     *
     * Initializes the algorithm with three points: the endpoints a, b and the
     * midpoint, then calls minimize() to perform the optimization.
     *
     * \param a Lower bound
     * \param b Upper bound
     * \return Estimated minimum location
     *
     * \note IMPORTANT: This version uses m_function->eval() directly
     *       for the initial evaluations instead of evaluate() to maintain
     *       correct counting (3 initial evaluations).
     *
     * \internal
     */
    Real eval_impl( Real a, Real b )
    {
      // Inizializza con tre punti: a, midpoint, b
      m_a  = a;
      m_fa = m_function->eval( m_a );
      m_b  = b;
      m_fb = m_function->eval( m_b );
      m_x3 = ( a + b ) / Real( 2 );
      m_f3 = m_function->eval( m_x3 );

      m_num_iter_done = 0;
      m_num_fun_eval  = 3;

      // Inizializza tolleranza se non impostata
      if ( m_tolerance == Real( 0 ) ) { m_tolerance = pow( machine_eps<Real>(), Real( 2.0 / 3.0 ) ); }

      return minimize();
    }

    /**
     * \brief Perform one bracketing step of the trichotomy algorithm
     *
     * This is the heart of the algorithm. At each call:
     * 1. Divide the current interval [a, b] into subintervals
     * 2. Evaluate the function at 2-4 new points (adaptive)
     * 3. Identify the subinterval containing the minimum
     * 4. Update [a, b] to the new narrower subinterval
     *
     * \return true if |b-a| < tolerance (convergence reached), false otherwise
     *
     * \par Adaptive Strategy:
     * The algorithm uses an adaptive approach to minimize evaluations:
     * - Initially evaluates point x2 = (a + 2*x3)/3
     * - If f(x2) <= f(x3): explore left (1 extra evaluation for x1)
     * - If f(x2) > f(x3): explore right (2 extra evaluations for x4 and x5)
     *
     * \par Invariants:
     * - At start: x3 is the point with the minimum known value
     * - At end: x3 is still the point with the minimum known value
     * - The interval [a, b] reduces by a factor of ~1/3 per iteration
     *
     * \par Reduction Cases:
     * The algorithm considers several cases to reduce the interval:
     *
     * **Case 1** (f(x1) <= f(x2) <= f(x3)):
     * - Minimum in [a, x2]
     * - New center: x1
     * - New interval: [a, x2]
     *
     * **Case 2** (f(x2) < f(x1), f(x3)):
     * - Minimum in [x1, x3]
     * - New center: x2
     * - New interval: [x1, x3]
     *
     * **Case 3** (f(x5) <= f(x4) <= f(x3)):
     * - Minimum in [x4, b]
     * - New center: x5
     * - New interval: [x4, b]
     *
     * **Case 4** (f(x4) < f(x5) and f(x4) < f(x3)):
     * - Minimum in [x3, x5]
     * - New center: x4
     * - New interval: [x3, x5]
     *
     * **Case 5** (f(x3) < f(x2) and f(x3) < f(x4)):
     * - Minimum in [x2, x4]
     * - Center remains: x3
     * - New interval: [x2, x4]
     *
     * \par Point Scheme:
     * \code
     * Interval: [a ---- x2 -- x3 -- x4 ---- b]
     *                |          |          |
     *              x1 (if left)   x5 (if right)
     * \endcode
     *
     * \internal
     */
    bool bracketing()
    {
      // Calcola x2 = (a + 2*x3)/3, divide l'intervallo sinistro
      m_x2 = ( m_a + Real( 2 ) * m_x3 ) / Real( 3 );
      m_f2 = evaluate( m_x2 );

      if ( m_f2 <= m_f3 )
      {
        // Il minimo potrebbe essere nel sottointervallo sinistro [a, x3]
        // Esplora ulteriormente con x1 = (a + x2)/2
        m_x1 = ( m_a + m_x2 ) / Real( 2 );
        m_f1 = evaluate( m_x1 );

        if ( m_f1 <= m_f2 )
        {
          // Caso 1: f(x1) <= f(x2) <= f(x3)
          // Il minimo è in [a, x2], usa x1 come nuovo centro
          m_x3 = m_x1;
          m_f3 = m_f1;
          m_b  = m_x2;
          m_fb = m_f2;
        }
        else
        {
          // Caso 2: f(x2) < f(x1), quindi f(x2) < f(x3)
          // Il minimo è in [x1, x3], usa x2 come nuovo centro
          m_a  = m_x1;
          m_fa = m_f1;
          m_b  = m_x3;
          m_fb = m_f3;
          m_x3 = m_x2;
          m_f3 = m_f2;
        }
      }
      else
      {
        // f(x2) > f(x3): il minimo potrebbe essere nel sottointervallo destro
        // Calcola x4 = (2*x3 + b)/3
        m_x4 = ( m_b + Real( 2 ) * m_x3 ) / Real( 3 );
        m_f4 = evaluate( m_x4 );

        if ( m_f4 <= m_f3 )
        {
          // Il minimo potrebbe essere più a destra, esplora con x5
          m_x5 = ( Real( 2 ) * m_b + m_x3 ) / Real( 3 );
          m_f5 = evaluate( m_x5 );

          if ( m_f5 <= m_f4 )
          {
            // Caso 3: f(x5) <= f(x4) <= f(x3)
            // Il minimo è in [x4, b], usa x5 come nuovo centro
            m_a  = m_x4;
            m_fa = m_f4;
            m_x3 = m_x5;
            m_f3 = m_f5;
          }
          else
          {
            // Caso 4: f(x4) < f(x5) e f(x4) < f(x3)
            // Il minimo è in [x3, x5], usa x4 come nuovo centro
            m_a  = m_x3;
            m_fa = m_f3;
            m_x3 = m_x4;
            m_f3 = m_f4;
            m_b  = m_x5;
            m_fb = m_f5;
          }
        }
        else
        {
          // Caso 5: f(x3) < f(x2) e f(x3) < f(x4)
          // Il minimo è in [x2, x4], x3 rimane il centro
          m_a  = m_x2;
          m_fa = m_f2;
          m_b  = m_x4;
          m_fb = m_f4;
        }
      }

      // Verifica convergenza: intervallo sufficientemente piccolo?
      return ( m_b - m_a ) < m_tolerance;
    }

    /**
     * \brief Main minimization loop
     *
     * Performs successive bracketing iterations until one of the
     * following conditions occurs:
     * 1. Convergence reached (|b-a| < tolerance)
     * 2. Maximum number of iterations reached
     * 3. Maximum number of function evaluations reached
     *
     * \return Estimated minimum location (m_x3)
     *
     * \post m_converged is true if convergence was reached
     * \post m_num_iter_done contains the number of iterations performed
     *
     * \note The loop terminates immediately when convergence is reached,
     *       so m_num_iter_done may be less than m_max_iteration.
     *
     * \internal
     */
    Real minimize()
    {
      m_num_iter_done = 0;
      m_converged     = false;

      while ( m_num_iter_done++ < m_max_iteration )
      {
        m_converged = bracketing();
        if ( m_converged ) break;
        if ( m_num_fun_eval >= m_max_fun_evaluation ) break;
      }

      return m_x3;
    }

    /**
     * \brief Internal implementation of search with interval expansion
     *
     * This method implements a two-phase search:
     *
     * **Phase 1 - Expansion (Minimum Bracketing):**
     * - Start with interval [x-delta, x+delta]
     * - Identify the direction of the minimum by comparing f(x-delta), f(x), f(x+delta)
     * - Expand the interval in the direction of the minimum with geometric steps
     * - Continue until the minimum is bracketed (contained in the interval)
     *
     * **Phase 2 - Refinement:**
     * - Add an intermediate point to improve bracketing
     * - Apply the standard trichotomy algorithm (minimize())
     *
     * \param x Center of the initial interval
     * \param delta Initial radius (must be > 0)
     * \return Estimated minimum location
     *
     * \par Expansion Strategy:
     * - Expansion factor: 2× at each step
     * - Direction: toward decreasing function values
     * - Stop criterion: trend reversal (minimum bracketed)
     *
     * \par Left Expansion Logic:
     * If f(a) < f(b), the minimum appears to be on the left:
     * 1. While f(a) < f(x3), continue expanding left
     * 2. Move b ← x3, x3 ← a
     * 3. Calculate new a = x3 - 2*(b-x3)
     * 4. When f(a) >= f(x3), the minimum is bracketed
     *
     * \par Right Expansion Logic:
     * If f(a) >= f(b), the minimum appears to be on the right:
     * 1. While f(x3) > f(b), continue expanding right
     * 2. Move a ← x3, x3 ← b
     * 3. Calculate new b = x3 + 2*(x3-a)
     * 4. When f(x3) <= f(b), the minimum is bracketed
     *
     * \par Post-Expansion Refinement:
     * After expansion, an intermediate point is added to improve
     * the initial bracketing quality before calling minimize():
     * - On left: x2 = (a+x3)/2
     * - On right: x4 = (x3+b)/2
     *
     * \par Example Sequence:
     * \code
     * Initial: [0.9, 1.0, 1.1]  f(0.9)=5, f(1.0)=3, f(1.1)=4
     * → f(a) > f(b), expand right
     * Step 1: [1.0, 1.1, 1.3]  f(1.3)=2
     * Step 2: [1.1, 1.3, 1.7]  f(1.7)=3
     * → f(x3) < f(b), stop! Minimum bracketed in [1.1, 1.7]
     * Refine: add x4=(1.3+1.7)/2=1.5
     * → Final interval: [1.1, 1.5, 1.7] ready for minimize()
     * \endcode
     *
     * \warning Expansion can require many evaluations for functions
     *          with minima far from the initial point
     *
     * \note IMPORTANT: This version uses m_function->eval() directly
     *       for initial evaluations (like eval_impl) instead of evaluate()
     *       to maintain correct counting.
     *
     * \internal
     */
    Real search_impl( Real x, Real delta )
    {
      m_num_iter_done = 0;
      m_num_fun_eval  = 3;

      // Inizializza con tre punti: x-delta, x, x+delta
      m_a  = x - delta;
      m_fa = m_function->eval( m_a );
      m_x3 = x;
      m_f3 = m_function->eval( m_x3 );
      m_b  = x + delta;
      m_fb = m_function->eval( m_b );

      // Inizializza tolleranza se non impostata
      if ( m_tolerance == Real( 0 ) ) { m_tolerance = pow( machine_eps<Real>(), Real( 2.0 / 3.0 ) ); }

      // Fase di espansione: bracketing del minimo
      if ( m_fa < m_fb )
      {
        // Il minimo sembra essere a sinistra di x
        // Espandi verso sinistra finché f(a) >= f(x3)
        while ( m_fa < m_f3 )
        {
          // Sposta l'intervallo a sinistra con espansione 2×
          m_b  = m_x3;
          m_fb = m_f3;
          m_x3 = m_a;
          m_f3 = m_fa;
          m_a  = m_x3 - Real( 2 ) * ( m_b - m_x3 );
          m_fa = evaluate( m_a );
          ++m_num_iter_done;
        }

        // Raffina il lato sinistro con un punto intermedio
        m_x2 = ( m_a + m_x3 ) / Real( 2 );
        m_f2 = evaluate( m_x2 );

        if ( m_f2 <= m_f3 )
        {
          // x2 è migliore di x3, aggiorna l'intervallo
          m_b  = m_x3;
          m_fb = m_f3;
          m_x3 = m_x2;
          m_f3 = m_f2;
        }
        else
        {
          // x3 rimane il migliore, aggiorna solo il lato sinistro
          m_a  = m_x2;
          m_fa = m_f2;
        }
      }
      else
      {
        // Il minimo sembra essere a destra di x
        // Espandi verso destra finché f(b) <= f(x3)
        while ( m_f3 > m_fb )
        {
          // Sposta l'intervallo a destra con espansione 2×
          m_a  = m_x3;
          m_fa = m_f3;
          m_x3 = m_b;
          m_f3 = m_fb;
          m_b  = m_x3 + Real( 2 ) * ( m_x3 - m_a );
          m_fb = evaluate( m_b );
          ++m_num_iter_done;
        }

        // Raffina il lato destro con un punto intermedio
        m_x4 = ( m_x3 + m_b ) / Real( 2 );
        m_f4 = evaluate( m_x4 );

        if ( m_f4 <= m_f3 )
        {
          // x4 è migliore di x3, aggiorna l'intervallo
          m_a  = m_x3;
          m_fa = m_f3;
          m_x3 = m_x4;
          m_f3 = m_f4;
        }
        else
        {
          // x3 rimane il migliore, aggiorna solo il lato destro
          m_b  = m_x4;
          m_fb = m_f4;
        }
      }

      // Fase di raffinamento: applica l'algoritmo standard di tricotomia
      return minimize();
    }
  };

  /*! @} */  // end of Minimization group

}  // namespace Utils

#endif  // UTILS_TRICHOTOMY_dot_HH

//
// eof: Utils_Trichotomy.hh
//
