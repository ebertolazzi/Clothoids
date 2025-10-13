/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The IPsolver project is distributed under the MIT License.                                    *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_IPSOLVER_HH
#define INCLUDE_IPSOLVER_HH

// Standard libraries
#include <iostream>

// Print IPsolver errors
#ifndef IPSOLVER_ERROR
#define IPSOLVER_ERROR(MSG)             \
  {                                     \
    std::ostringstream os;              \
    os << MSG;                          \
    throw std::runtime_error(os.str()); \
  }
#endif

// Assert for IPsolver
#ifndef IPSOLVER_ASSERT
#define IPSOLVER_ASSERT(COND, MSG) \
  if (!(COND))                     \
  {                                \
    IPSOLVER_ERROR(MSG);           \
  }
#endif

// Warning for IPsolver
#ifndef IPSOLVER_WARNING
#define IPSOLVER_WARNING(MSG)       \
  {                                 \
    std::cout << MSG << std::endl;  \
  }
#endif

// Warning assert for IPsolver
#ifndef IPSOLVER_ASSERT_WARNING
#define IPSOLVER_ASSERT_WARNING(COND, MSG) \
  if (!(COND))                             \
  {                                        \
    IPSOLVER_WARNING(MSG);                 \
  }
#endif


// Define the basic constants for IPsolver
#ifndef IPSOLVER_BASIC_CONSTANTS
#define IPSOLVER_BASIC_CONSTANTS(Real) \
  static constexpr Real EPSILON{std::numeric_limits<Real>::epsilon()};     /**< Machine epsilon epsilon static constant value. */ \
  static constexpr Real EPSILON_HIGH{1.0e-12};                             /**< High precision epsilon static constant value. */ \
  static constexpr Real EPSILON_MEDIUM{1.0e-10};                           /**< Medium precision epsilon static constant value. */ \
  static constexpr Real EPSILON_LOW{1.0e-08};                              /**< Low precision epsilon static constant value. */ \
  static constexpr Real INFTY{std::numeric_limits<Real>::infinity()};      /**< Infinity static constant value. */ \
  static constexpr Real QUIET_NAN{std::numeric_limits<Real>::quiet_NaN()}; /**< Not-a-number static constant value. */
#endif

#ifndef IPSOLVER_DEFAULT_INTEGER_TYPE
#define IPSOLVER_DEFAULT_INTEGER_TYPE int
#endif

/**
* \brief Namespace for the IPsolver library
*
* IPsolver is a simple yet reasonably robust implementation of a primal-dual interior-point solver
* for convex programs with convex inequality constraints (it does not handle equality constraints).
* Precisely speaking, it will compute the solution to the following optimization problem:
* \f[
*  \begin{array}{l}
*    \text{minimize} ~ f(\mathbf{x}) \\
*    \text{subject to} ~ \mathbf{c}(\mathbf{x}) < \mathbf{0}
*  \end{array}
* \f]
* where \f$\mathbf{x} \in \mathbb{R}^n\f$ is the vector of optimization variables, \f$f: \mathbb{R}^n
* \to \mathbb{R}\f$ is the convex function, and \f$\mathbf{c}: \mathbb{R}^n \to \mathbb{R}^m\f$ is a
* convex function, and \f$\mathbf{c}(\mathbf{x})\f$ is a vector-valued function with outputs that are
* convex in \f$\mathbf{x}\f$.
*
* This code is mostly based on the descriptions provided in this reference:
*
* - Paul Armand, Jean C. Gilbert, Sophie Jan-Jegou. A Feasible BFGS Interior Point Algorithm for
*   Solving Convex Minimization Problems. SIAM Journal on Optimization, Vol. 11, No. 1, pp. 199-222.
*
* \note The descent direction can be set to either \c NEWTON for the Newton search direction, \c BFGS
* for the quasi-Newton search direction with the Broyden-Fletcher-Goldfarb-Shanno (BFGS) update, or
* \c  STEEPEST for the steepest descent direction. The steepest descent direction is often quite bad,
* and the solver may fail to converge to the solution if you take this option. Also note that we form
* a quasi-Newton approximation to the objective, not to the Lagrangian (as is usually done). This means
* that you will always have to provide second-order information about the inequality constraint functions.
*
* \note If your optimization problem is large (i.e. it involves a lot of optimization variables or
* inequality constraints) it might speed up the solver to output sparse matrices.
*
* \warning The interior-point solver may not work very well if your problem is very poorly scaled
* (i.e. the Hessian of the objective or the Hessian of one of the constraint functions is poorly
* conditioned). It is up to you to make sure you look at the conditioning of your problem.
*/
namespace IPsolver
{

  /**
  * \brief The Integer type as used for the API.
  *
  * The Integer type, \c \#define the preprocessor symbol \c IPSOLVER_DEFAULT_INTEGER_TYPE. The default
  * value is \c int.
  */
  using Integer = IPSOLVER_DEFAULT_INTEGER_TYPE;

} // namespace IPsolver

#include "IPsolver/Solver.hh"

#endif // INCLUDE_IPSOLVER_HH
