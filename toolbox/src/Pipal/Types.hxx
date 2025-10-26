/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Pipal project is distributed under the MIT License.                                       *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_PIPAL_DEFINES_HXX
#define INCLUDE_PIPAL_DEFINES_HXX

// Standard libraries
#include <string>
#include <iostream>
#include <iomanip>
#include <ostream>

// STL
#include <cmath>
#include <limits>
#include <algorithm>
#include <functional>
#include <type_traits>
#include <numeric>
#include <memory>

// Eigen library
#ifndef PIPAL_EIGEN_EXTERNAL
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#endif

// Print Pipal errors
#ifndef PIPAL_ERROR
#define PIPAL_ERROR(MSG)                \
  {                                     \
    std::ostringstream os;              \
    os << MSG;                          \
    throw std::runtime_error(os.str()); \
  }
#endif

// Assert for Pipal
#ifndef PIPAL_ASSERT
#define PIPAL_ASSERT(COND, MSG) \
  if (!(COND))                  \
  {                             \
    PIPAL_ERROR(MSG);           \
  }
#endif

// Warning for Pipal
#ifndef PIPAL_WARNING
#define PIPAL_WARNING(MSG)         \
  {                                \
    std::cout << MSG << std::endl; \
  }
#endif

// Warning assert for Pipal
#ifndef PIPAL_ASSERT_WARNING
#define PIPAL_ASSERT_WARNING(COND, MSG) \
  if (!(COND))                          \
  {                                     \
    PIPAL_WARNING(MSG);                 \
  }
#endif

#ifndef PIPAL_DEFAULT_INTEGER_TYPE
#define PIPAL_DEFAULT_INTEGER_TYPE int
#endif

/**
 * \brief Namespace for the Pipal library.
 *
 * Penalty-interior-point algorithm (Pipal) is a library for nonlinear constrained optimization with
 * inequality constraints (it does not explicitly handle equality constraints). Precisely speaking,
 * it will compute the solution to the optimization problem
 * \f[
 *  \begin{array}{l}
 *    \text{minimize} ~ f(\mathbf{x}) \\
 *    \text{subject to} ~ \mathbf{c}(\mathbf{x}) \leq \mathbf{0}
 *  \end{array} \text{,}
 * \f]
 * where \f$\mathbf{x} \in \mathbb{R}^n\f$ is the vector of optimization variables, \f$f: \mathbb{R}^n
 * \to \mathbb{R}\f$ is the objective function, and \f$\mathbf{c}: \mathbb{R}^n \to \mathbb{R}^m\f$
 * are the constraints.
 *
 * \note To create an equality constraint \f$ h(\mathbf{x}) = 0 \f$, one can define two inequality
 * constraints \f$ h(\mathbf{x}) \leq 0 \f$ and \f$ -h(\mathbf{x}) \leq 0 \f$.
 *
 * This code is mostly based on the descriptions provided in this reference:
 *
 * - Frank E. Curtis. A penalty-interior-point algorithm for nonlinear constrained optimization.
 *   Mathematical Programming Computation (2012) 4:181-209. DOI: 10.1007/s12532-012-0041-4.
 */
namespace Pipal
{
  /**
   * \brief The Integer type as used for the API.
   *
   * The Integer type, \c \#define the preprocessor symbol \c PIPAL_DEFAULT_INTEGER_TYPE. The default
   * value is \c int.
   */
  using Integer = PIPAL_DEFAULT_INTEGER_TYPE;

  template<typename Real> using Vector       = Eigen::Vector<Real, Eigen::Dynamic>;
  template<typename Real> using Matrix       = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
  template<typename Real> using SparseVector = Eigen::SparseVector<Real>;
  template<typename Real> using SparseMatrix = Eigen::SparseMatrix<Real>;
  template<typename Real> using Array        = Eigen::Array<Real, Eigen::Dynamic, 1>;

  using Indices = Eigen::Array<Integer, Eigen::Dynamic, 1>;
  using Mask    = Eigen::Array<bool, Eigen::Dynamic, 1>;

  /**
   * \brief Enumeration for the algorithm choice.
   *
   * The Algorithm enumeration defines the possible algorithm choices for the solver. The options
   * are \c CONSERVATIVE and \c ADAPTIVE.
   */
  using Algorithm = enum class Algorithm : Integer {CONSERVATIVE = 0, ADAPTIVE = 1};

  /**
   * \brief Select elements from a vector based on a boolean mask.
   * \param[in] mask The boolean mask.
   * \return The selected elements from the input vector.
   */
  static Indices find(Mask const & mask)
  {
    Indices out(mask.count());
    for (Integer i{0}, j{0}; i < mask.size(); ++i) {
      if (mask[i]) {out[j++] = i;}
    }
    return out;
  }

  /**
   * \brief Internal parameters for the solver algorithm.
   * \tparam Real The real number type.
   */
  template<typename Real>
  struct Parameter
  {
    static constexpr Real    rhs_bnd{1.0e+18};               /*!< Maximum absolute value allowed for constraint right-hand side. */
    static constexpr Real    grad_max{1.0e+02};              /*!< Gradient norm limit for scaling. */
    static constexpr Real    infeas_max{1.0e+02};            /*!< Infeasibility limit for penalty parameter update. */
    static constexpr Real    nnz_max{2.0e+04};               /*!< Maximum non-zeros in (upper triangle of) Newton matrix. */
    static constexpr Integer opt_err_mem{6};                 /*!< Optimality error history length. */
    static constexpr Real    ls_factor{5.0e-01};             /*!< Line search reduction factor. */
    static constexpr Real    ls_thresh{1.0e-08};             /*!< Line search threshold value. */
    static constexpr Real    ls_frac{1.0e-02};               /*!< Line search fraction-to-boundary constant. */
    static constexpr Real    slack_min{1.0e-20};             /*!< Slack variable bound. */
    static constexpr Real    shift_min{1.0e-12};             /*!< Hessian shift (non-zero) minimum value. */
    static constexpr Real    shift_factor1{5.0e-01};         /*!< Hessian shift update value (for decreases). */
    static constexpr Real    shift_factor2{6.0e-01};         /*!< Hessian shift update value (for increases). */
    static constexpr Real    shift_max{1.0e+08};             /*!< Hessian shift maximum value. */
    static constexpr Real    rho_init{1.0e-01};              /*!< Penalty parameter initial value. */
    static constexpr Real    rho_min{1.0e-12};               /*!< Penalty parameter minimum value. */
    static constexpr Real    rho_factor{5.0e-01};            /*!< Penalty parameter reduction factor. */
    static constexpr Integer rho_trials{8};                  /*!< Penalty parameter number of trial values per iteration. */
    static constexpr Real    mu_init{1.0e-01};               /*!< Interior-point parameter initial value. */
    static constexpr Real    mu_min{1.0e-12};                /*!< Interior-point parameter minimum value. */
    static constexpr Real    mu_factor{1.0e-01};             /*!< Interior-point parameter reduction factor. */
    static constexpr Real    mu_factor_exp{1.5};             /*!< Interior-point parameter reduction exponent. */
    static constexpr Integer mu_trials{4};                   /*!< Interior-point parameter number of trial values per iteration. */
    static constexpr Real    mu_max{1.0e-01};                /*!< Interior-point parameter maximum value. */
    static constexpr Real    mu_max_exp0{0.0};               /*!< Interior-point parameter maximum exponent in increases (default). */
    static constexpr Real    update_con_1{1.0e-02};          /*!< Steering rule constant 1. */
    static constexpr Real    update_con_2{1.0e-02};          /*!< Steering rule constant 2. */
    static constexpr Real    update_con_3{1.01};             /*!< Adaptive interior-point rule constant. */

    Real      opt_err_tol{1.0e-10};           /*!< Default optimality tolerance. */
    Integer   iter_max{1000};                 /*!< Default iteration limit. */
    Algorithm algorithm{Algorithm::ADAPTIVE}; /*!< Algorithm choice. */
    Real      mu_max_exp{0.0};                /*!< Interior-point parameter maximum exponent in increases. */
  }; // struct Parameter

  /**
   * \brief Internal counters for solver statistics.
   */
  using Counter = struct Counter
  {
    Integer f{0}; // Function evaluation counter
    Integer g{0}; // Gradient evaluation counter
    Integer H{0}; // Hessian evaluation counter
    Integer k{0}; // Iteration counter
    Integer M{0}; // Matrix factorization counter
  }; // struct Counter

  /**
   * \brief Input structure holding all the data defining the optimization problem.
   * \tparam Real The real number type.
   */
  template<typename Real>
  struct Input
  {
    std::string name; /*!< Problem identity. */
    Integer       n0; /*!< Number of original formulation variables. */
    Indices       I1; /*!< Indices of free variables. */
    Indices       I2; /*!< Indices of fixed variables. */
    Indices       I3; /*!< Indices of lower bounded variables. */
    Indices       I4; /*!< Indices of upper bounded variables. */
    Indices       I5; /*!< Indices of lower and upper bounded variables. */
    Indices       I6; /*!< Indices of equality constraints. */
    Indices       I7; /*!< Indices of lower bounded constraints. */
    Indices       I8; /*!< Indices of upper bounded constraints. */
    Indices       I9; /*!< Indices of lower and upper bounded constraints. */
    Vector<Real>  x0; /*!< Initial guess for the primal variables. */
    Vector<Real>  b2; /*!< Right-hand side of fixed variables. */
    Vector<Real>  l3; /*!< Right-hand side of lower bounded variables. */
    Vector<Real>  u4; /*!< Right-hand side of upper bounded variables. */
    Vector<Real>  l5; /*!< Right-hand side of lower half of lower and upper bounded variables. */
    Vector<Real>  u5; /*!< Right-hand side of upper half of lower and upper bounded variables. */
    Vector<Real>  b6; /*!< Right-hand side of equality constraints. */
    Vector<Real>  l7; /*!< Right-hand side of lower bounded constraints. */
    Vector<Real>  u8; /*!< Right-hand side of upper bounded constraints. */
    Vector<Real>  l9; /*!< Right-hand side of lower half of lower and upper bounded constraints. */
    Vector<Real>  u9; /*!< Right-hand side of upper half of lower and upper bounded constraints. */
    Integer       n1; /*!< Number of free variables. */
    Integer       n2; /*!< Number of fixed variables. */
    Integer       n3; /*!< Number of lower bounded variables. */
    Integer       n4; /*!< Number of upper bounded variables. */
    Integer       n5; /*!< Number of lower and upper bounded variables. */
    Integer       n6; /*!< Number of equality constraints. */
    Integer       n7; /*!< Number of lower bounded constraints. */
    Integer       n8; /*!< Number of upper bounded constraints. */
    Integer       n9; /*!< Number of lower and upper bounded constraints. */
    Integer       nV; /*!< Number of variables. */
    Integer       nI; /*!< Number of inequality constraints. */
    Integer       nE; /*!< Number of equality constraints. */
    Integer       nA; /*!< Size of primal-dual matrix. */
    Integer       vi; /*!< Counter for invalid bounds. */
  }; // struct Input

  /**
   * \brief Class for managing the current iterate of the solver.
   * \tparam Real The real number type.
   */
  template<typename Real>
  struct Iterate
  {
    using LDLT = Eigen::SimplicialLDLT<SparseMatrix<Real>, Eigen::Lower>;

    Vector<Real>       x;     /*!< Primal point. */
    Real               rho;   /*!< Penalty parameter value. */
    Real               rho_;  /*!< Penalty parameter last value. */
    Real               mu;    /*!< Interior-point parameter value. */
    Real               f;     /*!< Objective function value (scaled). */
    Real               fu;    /*!< Objective function value (unscaled). */
    Vector<Real>       g;     /*!< Objective gradient value. */
    Array<Real>        r1;    /*!< Equality constraint slack value. */
    Array<Real>        r2;    /*!< Equality constraint slack value. */
    Array<Real>        cE;    /*!< Equality constraint value (scaled). */
    SparseMatrix<Real> JE;    /*!< Equality constraint Jacobian value. */
    Integer            JEnnz; /*!< Equality constraint Jacobian nonzeros. */
    Array<Real>        lE;    /*!< Equality constraint multipliers. */
    Array<Real>        s1;    /*!< Inequality constraint slack value. */
    Array<Real>        s2;    /*!< Inequality constraint slack value. */
    Array<Real>        cI;    /*!< Inequality constraint value (scaled). */
    SparseMatrix<Real> JI;    /*!< Inequality constraint Jacobian value. */
    Integer            JInnz; /*!< Inequality constraint Jacobian nonzeros. */
    Array<Real>        lI;    /*!< Inequality constraint multipliers. */
    SparseMatrix<Real> H;     /*!< Hessian of Lagrangian. */
    Integer            Hnnz;  /*!< Hessian of Lagrangian nonzeros. */
    Real               v;     /*!< Feasibility violation measure value (scaled). */
    Real               vu;    /*!< Feasibility violation measure value (unscaled). */
    Real               v0;    /*!< Feasibility violation measure initial value. */
    Real               phi;   /*!< Merit function value. */
    LDLT               ldlt;  /*!< LDLT factorization of Newton matrix. */
    Integer            Annz;  /*!< Newton matrix (upper triangle) nonzeros. */
    Real               shift; /*!< Hessian shift value. */
    SparseVector<Real> b;     /*!< Newton right-hand side. */
    Vector<Real>       kkt;   /*!< KKT errors. */
    Vector<Real>       kkt_;  /*!< KKT errors last value. */
    Integer            err;   /*!< Function evaluation error flag. */

    Real               fs;      /*!< Objective scaling factor. */
    Array<Real>        cEs;     /*!< Equality constraint scaling factors. */
    Array<Real>        cEu;     /*!< Equality constraint value (unscaled). */
    Array<Real>        cIs;     /*!< Inequality constraint scaling factors. */
    Array<Real>        cIu;     /*!< Inequality constraint value (unscaled). */
    SparseMatrix<Real> A;       /*!< Newton matrix. */
    Real               shift22; /*!< Newton matrix (2,2)-block shift value. */
    Real               v_;      /*!< Feasibility violation measure last value. */
    bool               cut_;    /*!< Boolean value for last backtracking line search. */
  }; // struct Iterate

  /**
   * \brief Class for managing the current search direction of the solver.
   * \tparam Real The real number type.
   */
  template<typename Real>
  struct Direction
  {
    Vector<Real> x;       /*!< Primal direction. */
    Real         x_norm;  /*!< Primal direction norm value. */
    Real         x_norm_; /*!< Primal direction norm last value. */
    Array<Real>  r1;      /*!< Equality constraint slack direction. */
    Array<Real>  r2;      /*!< Equality constraint slack direction. */
    Array<Real>  lE;      /*!< Equality constraint multiplier direction. */
    Array<Real>  s1;      /*!< Inequality constraint slack direction. */
    Array<Real>  s2;      /*!< Inequality constraint slack direction. */
    Array<Real>  lI;      /*!< Inequality constraint multiplier direction. */
    Real         l_norm;  /*!< Constraint multiplier direction norm. */
    Real         lred0;   /*!< Penalty-interior-point linear model value for zero penalty parameter. */
    Real         ltred0;  /*!< Penalty-interior-point linear model reduction value for zero penalty parameter. */
    Real         ltred;   /*!< Penalty-interior-point linear model reduction value. */
    Real         qtred;   /*!< Penalty-interior-point quadratic model reduction value. */
    Real         m;       /*!< Quality function value. */
  }; // struct Direction

  /**
    * \brief Class for managing the acceptance criteria of the solver.
    * \tparam Real The real number type.
    */
  template<typename Real>
  struct Acceptance
  {
    Real p0{0.0};  /*!< Fraction-to-the-boundary steplength. */
    Real p{0.0};   /*!< Primal steplength. */
    Real d{0.0};   /*!< Dual steplength. */
    bool s{false}; /*!< Bool for second-order correction. */
  }; // struct Acceptance

} // namespace Pipal

#endif // INCLUDE_PIPAL_DEFINES_HH
