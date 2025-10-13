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

#ifndef INCLUDE_IPSOLVER_SOLVER_HH
#define INCLUDE_IPSOLVER_SOLVER_HH

// Standard libraries
#include <iostream>
#include <iomanip>
#include <algorithm>

// Eigen library
#include <Eigen/Dense>

// IPsolver includes
#include "IPsolver/Problem.hh"

namespace IPsolver {

  /**
  * \brief Solver class for the IPsolver library
  *
  * The Solver class provides an interface for solving convex optimization problems using
  * interior-point methods.
  *
  * The solver can handle problems of the form:
  * \f[
  *  \begin{array}{l}
  *    \text{minimize} ~ f(\mathbf{x}) \\
  *    \text{subject to} ~ \mathbf{c}(\mathbf{x}) < \mathbf{0}
  *  \end{array}
  * \f]
  * where \f$\mathbf{x} \in \mathbb{R}^n\f$ is the vector of optimization variables, \f$f: \mathbb{R}^n
  * \to \mathbb{R}\f$ is the convex objective function, and \f$\mathbf{c}: \mathbb{R}^n \to \mathbb{R}^m\f$
  * is the vector of convex inequality constraints.
  *
  * The solver supports different descent methods, including Newton's method, BFGS and steepest
  * descent. The user can specify the desired method through the `descent` method.
  * \tparam Real The floating-point type.
  * \tparam N The size of the primal variable vector (default is dynamic).
  * \tparam M The size of the dual variable vector (default is dynamic).
  */
  template<typename Real, Integer N = Eigen::Dynamic, Integer M = Eigen::Dynamic>
  class Solver
  {
  public:
    using Descent = enum class Descent : Integer {
      NEWTON   = 0, /**< Use Newton's method for descent direction */
      BFGS     = 1, /**< Use BFGS method for descent direction */
      STEEPEST = 2  /**< Use steepest descent method for descent direction */
    }; /**< Descent direction enumeration */

    using UniquePtr = std::unique_ptr<Problem<Real, N, M>>;

    using VectorN = typename Problem<Real, N, M>::VectorN;
    using VectorM = typename Problem<Real, N, M>::VectorM;
    using MatrixJ = typename Problem<Real, N, M>::MatrixJ;
    using MatrixH = typename Problem<Real, N, M>::MatrixH;

    using ObjectiveFunc           = typename ProblemWrapper<Real, N, M>::ObjectiveFunc;
    using ObjectiveGradientFunc   = typename ProblemWrapper<Real, N, M>::ObjectiveGradientFunc;
    using ObjectiveHessianFunc    = typename ProblemWrapper<Real, N, M>::ObjectiveHessianFunc;
    using ConstraintsFunc         = typename ProblemWrapper<Real, N, M>::ConstraintsFunc;
    using ConstraintsJacobianFunc = typename ProblemWrapper<Real, N, M>::ConstraintsJacobianFunc;
    using LagrangianHessianFunc   = typename ProblemWrapper<Real, N, M>::LagrangianHessianFunc;

  private:
    std::unique_ptr<Problem<Real, N, M>> m_problem; /**< Problem object */

    Descent m_descent{Descent::NEWTON}; /**< Descent direction method */
    Real    m_tolerance{1.0e-6};        /**< Tolerance for convergence */
    Integer m_max_iterations{100};      /**< Maximum number of iterations */

    // Some algorithm parameters
    bool m_verbose{false};    /**< Verbosity flag */
    Real m_epsilon{1.0e-8};   /**< Small constant to avoid numerical issues */
    Real m_sigma_max{0.5};    /**< Maximum value for the centering parameter */
    Real m_eta_max{0.25};     /**< Maximum value for the step size */
    Real m_mu_min{1.0e-9};    /**< Minimum value for the barrier parameter */
    Real m_alpha_max{0.995};  /**< Maximum value for the line search parameter */
    Real m_alpha_min{1.0e-6}; /**< Minimum value for the line search parameter */
    Real m_beta{0.75};        /**< Value for the backtracking line search */
    Real m_tau{0.01};         /**< Parameter for the sufficient decrease condition */

  public:
    /**
    * \brief Default constructor for the IPSolver class.
    *
    * Initializes the solver with default values for the objective, gradient, constraints, and Jacobian functions.
    */
    Solver() {};

    /**
    * \brief Constructor for the IPSolver class.
    *
    * Initializes the solver with the provided objective, gradient, constraints, and Jacobian functions.
    * \param[in] objective Objective function handle.
    * \param[in] objective_gradient Gradient of the objective function handle.
    * \param[in] constraints Constraints function handle.
    * \param[in] constraints_jacobian Jacobian of the constraints function handle.
    * \param[in] lagrangian_hessian Hessian of the Lagrangian function handle.
    * \warning The default descent direction is set to BFGS (approximation of the Hessian).
    */
    Solver(ObjectiveFunc const & objective, ObjectiveGradientFunc const & objective_gradient,
      ConstraintsFunc const & constraints, ConstraintsJacobianFunc const & constraints_jacobian,
      LagrangianHessianFunc const & lagrangian_hessian)
      : m_problem(std::make_unique<ProblemWrapper<Real, N, M>>(objective, objective_gradient,
          constraints, constraints_jacobian, lagrangian_hessian)), m_descent(Descent::BFGS) {}

    /**
    * \brief Constructor for the IPSolver class (with Hessian).
    *
    * Initializes the solver with the provided objective, gradient, constraints, Jacobian, and Hessian functions.
    * \param[in] objective Objective function handle.
    * \param[in] objective_gradient Gradient of the objective function handle.
    * \param[in] objective_hessian Hessian of the objective function handle.
    * \param[in] constraints Constraints function handle.
    * \param[in] constraints_jacobian Jacobian of the constraints function handle.
    * \param[in] lagrangian_hessian Hessian of the Lagrangian function handle.
    * \warning The default descent direction is set to Newton (exact Hessian).
    */
    Solver(ObjectiveFunc const & objective, ObjectiveGradientFunc const & objective_gradient,
      ObjectiveHessianFunc const & objective_hessian, ConstraintsFunc const & constraints,
      ConstraintsJacobianFunc const & constraints_jacobian, LagrangianHessianFunc const & lagrangian_hessian)
      : m_problem(std::make_unique<ProblemWrapper<Real, N, M>>(objective, objective_gradient, objective_hessian,
          constraints, constraints_jacobian, lagrangian_hessian)), m_descent(Descent::NEWTON) {}

    /**
    * \brief Constructor for the IPSolver class (with a unique pointer to a Problem object).
    *
    * Initializes the solver with the provided unique pointer to a Problem object.
    * \param[in] problem The unique pointer to the Problem object to use.
    * \param[in] descent The descent direction method to use (default is Newton).
    * \warning The pointer is moved into the solver, so it should not be used after this call.
    */
    Solver(std::unique_ptr<Problem<Real, N, M>> && problem, Descent descent = Descent::NEWTON)
      : m_problem(std::move(problem)), m_descent(descent) {}

    /**
    * \brief Deleted copy constructor.
    *
    * This class is not copyable.
    */
    Solver(Solver const &) = delete;

    /**
    * \brief Deleted assignment operator.
    *
    * This class is not assignable.
    */
    Solver& operator=(Solver const &) = delete;

    /**
    * \brief Deleted move constructor.
    *
    * This class is not movable.
    */
    Solver(Solver &&) = delete;

    /**
    * \brief Deleted move assignment operator.
    *
    * This class is not movable.
    */
    Solver& operator=(Solver &&) = delete;

    /**
    * \brief Destructor for the IPSolver class.
    *
    * Cleans up resources used by the solver.
    */
    ~Solver() = default;

    /**
    * \brief Sets the problem to be solved.
    *
    * This method allows the user to specify the problem to be solved.
    * \param[in] problem The problem to set.
    */
    void problem(Problem<Real, N, M> const & problem) {
      this->m_problem = std::make_unique<ProblemWrapper<Real, N, M>>(problem);
    }

    /**
    * \brief Sets the problem to be solved using a unique pointer.
    *
    * This method allows the user to specify the problem to be solved using a unique pointer.
    * \param[in] problem The unique pointer to the problem to set.
    */
    void problem(std::unique_ptr<Problem<Real, N, M>> && problem) {this->m_problem = std::move(problem);}

    /**
    * \brief Gets the current problem being solved.
    * \return A reference to the current problem.
    */
    Problem<Real, N, M> const& problem() const {return *this->m_problem;}

    /**
    * \brief Sets the descent direction method for the solver.
    *
    * This method allows the user to specify the descent direction method to be used by the solver.
    * The available methods are:
    * - Descent::NEWTON: Use Newton's method for descent direction.
    * - Descent::BFGS: Use BFGS method for descent direction.
    * - Descent::STEEPEST: Use steepest descent method for descent direction.
    *
    * \param[in] descent The descent direction method to set.
    */
    void descent(Descent descent) {this->m_descent = descent;}

    /**
    * \brief Gets the current descent direction method.
    * \return The current descent direction method.
    */
    Descent descent() const {return this->m_descent;}

    /**
    * Get the verbose mode.
    * \return The verbose mode.
    */
    bool verbose_mode() {return this->m_verbose;}

    /**
    * Set the verbose mode.
    * \param[in] t_verbose The verbose mode.
    */
    void verbose_mode(bool t_verbose) {this->m_verbose = t_verbose;}

    /**
    * Enable the verbose mode.
    */
    void enable_verbose_mode() {this->verbose_mode(true);}

    /**
    * Disable the verbose mode.
    */
    void disable_verbose_mode() {this->verbose_mode(false);}

    /**
    * \brief Sets the convergence tolerance for the solver.
    *
    * This method allows the user to specify the tolerance for convergence.
    * The solver will stop when the residuals are below this tolerance.
    *
    * \param[in] tolerance The convergence tolerance.
    */
    void tolerance(Real tolerance) {
      IPSOLVER_ASSERT(tolerance > 0.0,
        "IPsolver::Solver::tolerance(...): input value must be positive");
      this->m_tolerance = tolerance;
    }

    /**
    * \brief Gets the current tolerance for convergence.
    * \return The current tolerance for convergence.
    */
    Real tolerance() const {return this->m_tolerance;}

    /**
    * \brief Sets the maximum number of iterations for the solver.
    *
    * This method allows the user to specify the maximum number of iterations
    * the solver will perform before stopping.
    *
    * \param[in] max_iterations The maximum number of iterations.
    */
    void max_iterations(Integer max_iterations) {
      IPSOLVER_ASSERT(max_iterations > 0,
        "IPsolver::Solver::max_iterations(...): input value must be positive");
      this->m_max_iterations = max_iterations;
    }

    /**
    * \brief Gets the current maximum number of iterations.
    * \return The current maximum number of iterations.
    */
    Integer max_iterations() const {return this->m_max_iterations;}

    /**
    * \brief Sets the small constant epsilon to avoid numerical issues.
    * \param[in] epsilon The small positive constant.
    */
    void epsilon(Real epsilon) {
      IPSOLVER_ASSERT(epsilon > 0.0,
        "IPsolver::Solver::epsilon(...): input value must be positive");
      this->m_epsilon = epsilon;
    }

    /**
    * \brief Gets the current epsilon value.
    * \return The current epsilon value.
    */
    Real epsilon() const {return this->m_epsilon;}

    /**
    * \brief Sets the maximum value for the centering parameter sigma.
    * \param[in] sigma_max The maximum value for sigma (must be positive).
    */
    void sigma_max(Real sigma_max) {
      IPSOLVER_ASSERT(sigma_max > 0.0,
        "IPsolver::Solver::sigma_max(...): input value must be positive");
      this->m_sigma_max = sigma_max;
    }

    /**
    * \brief Gets the current maximum value for sigma.
    * \return The current sigma_max value.
    */
    Real sigma_max() const {return this->m_sigma_max;}

    /**
    * \brief Sets the maximum value for the step size eta.
    * \param[in] eta_max The maximum value for eta (must be positive).
    */
    void eta_max(Real eta_max) {
      IPSOLVER_ASSERT(eta_max > 0.0,
        "IPsolver::Solver::eta_max(...): input value must be positive");
      this->m_eta_max = eta_max;
    }

    /**
    * \brief Gets the current maximum value for eta.
    * \return The current eta_max value.
    */
    Real eta_max() const {return this->m_eta_max;}

    /**
    * \brief Sets the minimum value for the barrier parameter mu.
    * \param[in] mu_min The minimum value for mu (must be positive).
    */
    void mu_min(Real mu_min) {
      IPSOLVER_ASSERT(mu_min > 0.0,
        "IPsolver::Solver::mu_min(...): input value must be positive");
      this->m_mu_min = mu_min;
    }

    /**
    * \brief Gets the current minimum value for mu.
    * \return The current mu_min value.
    */
    Real mu_min() const {return this->m_mu_min;}

    /**
    * \brief Sets the maximum value for the line search parameter alpha.
    * \param[in] alpha_max The maximum value for alpha (must be positive).
    */
    void alpha_max(Real alpha_max) {
      IPSOLVER_ASSERT(alpha_max > 0.0,
        "IPsolver::Solver::alpha_max(...): input value must be positive");
      this->m_alpha_max = alpha_max;
    }

    /**
    * \brief Gets the current maximum value for alpha.
    * \return The current alpha_max value.
    */
    Real alpha_max() const {return this->m_alpha_max;}

    /**
    * \brief Sets the minimum value for the line search parameter alpha.
    * \param[in] alpha_min The minimum value for alpha (must be positive).
    */
    void alpha_min(Real alpha_min) {
      IPSOLVER_ASSERT(alpha_min > 0.0,
        "IPsolver::Solver::alpha_min(...): input value must be positive");
      this->m_alpha_min = alpha_min;
    }

    /**
    * \brief Gets the current minimum value for alpha.
    * \return The current alpha_min value.
    */
    Real alpha_min() const {return this->m_alpha_min;}

    /**
    * \brief Sets the value for the backtracking line search parameter beta.
    * \param[in] beta The value for beta (must be positive).
    */
    void beta(Real beta) {
      IPSOLVER_ASSERT(beta > 0.0,
        "IPsolver::Solver::beta(...): input value must be positive");
      this->m_beta = beta;
    }

    /**
    * \brief Gets the current value for beta.
    * \return The current beta value.
    */
    Real beta() const {return this->m_beta;}

    /**
    * \brief Sets the parameter for the sufficient decrease condition tau.
    * \param[in] tau The value for tau (must be positive).
    */
    void tau(Real tau) {
      IPSOLVER_ASSERT(tau > 0.0,
        "IPsolver::Solver::tau(...): input value must be positive");
      this->m_tau = tau;
    }

    /**
    * \brief Gets the current value for tau.
    * \return The current tau value.
    */
    Real tau() const {return this->m_tau;}

    /**
    * \brief Solves the optimization problem using the interior-point method.
    *
    * This method implements the interior-point algorithm to solve the optimization problem defined
    * by the objective function, constraints, and their respective gradients and Jacobians.
    * \param[in] x_guess Initial guess for the optimization variables.
    * \param[out] x_sol Solution vector.
    * \return True if the optimization was successful, false otherwise.
    */
    bool solve(const VectorN & x_guess, VectorN & x_sol)
    {
      #define CMD "IPsolver::Solver::solve(...): "

      using MatrixM = Eigen::Matrix<Real, M, M>;
      using ArrayM = Eigen::Array<Real, M, 1>;
      using MaskM  = Eigen::Array<bool, M, 1>;

      // INITIALIZATION
      // Get the number of primal variables (n), the number of constraints (m), the total number of
      // primal-dual optimization variables (nv), and initialize the Lagrange multipliers and the
      // second-order information
      VectorN x(x_guess);
      VectorM c(this->m_problem->constraints(x));
      Integer n{static_cast<Integer>(x.size())};
      Integer m{static_cast<Integer>(c.size())};
      Integer nv{n + m};
      VectorM z(VectorM::Ones(m));
      MatrixH B(MatrixH::Identity(n, n));

      // Repeat while the convergence criterion has not been satisfied, and we haven't reached the
      // maximum number of iterations
      Real alpha{0.0}, f, norm_r0, eta, sigma, duality_gap, mu, psi, dpsi, psi_new;
      VectorN r_x, g_b, g_old, p_x, g, x_new;
      VectorM r_c, c_epsilon, p_z, z_new;
      MatrixJ J;
      MatrixH H, W;
      MatrixM S;
      MaskM mask;
      ArrayM ratio;
      bool converged{false};
      for (Integer iter{0}; iter < this->m_max_iterations; ++iter)
      {

        // COMPUTE OBJECTIVE, GRADIENT, CONSTRAINTS, ETC
        // Compute the response of the objective function, the gradient of the objective, the
        // response of the inequality constraints, the Jacobian of the inequality constraints, the
        // Hessian of the Lagrangian (minus the Hessian of the objective) and, optionally, the
        // Hessian of the objective.
        f = this->m_problem->objective(x);
        c = this->m_problem->constraints(x);
        g = this->m_problem->objective_gradient(x);
        J = this->m_problem->constraints_jacobian(x, z);
        W = this->m_problem->lagrangian_hessian(x, z);
        if (this->m_descent == Descent::NEWTON) {
          B = this->m_problem->objective_hessian(x);
        }

        // Compute the responses of the unperturbed Karush-Kuhn-Tucker optimality conditions.
        r_x = g + J.transpose() * z; // Dual residual
        r_c = c.array() * z.array(); // Complementarity

        // Set some parameters that affect convergence of the primal-dual interior-point method
        norm_r0     = std::sqrt<Real>(r_x.squaredNorm() + r_c.squaredNorm());
        eta         = std::min<Real>(this->m_eta_max, norm_r0 / nv);
        sigma       = std::min<Real>(this->m_sigma_max, std::sqrt(norm_r0 / nv));
        duality_gap = static_cast<Real>(-c.dot(z));
        mu          = std::max<Real>(this->m_mu_min, sigma * duality_gap / m);

        if (this->m_verbose) {
          std::cout << "Iteration " << iter << std::setprecision(6) << std::scientific
            << ": f = " << f << ", |r| = " << norm_r0 / nv << std::endl;
        }

        // CONVERGENCE CHECK
        // If the norm of the responses is less than the specified tolerance, we are done
        if (norm_r0 / nv < this->m_tolerance) {
          converged = true;
          break;
        }

        // Update the BFGS approximation to the Hessian of the objective
        if (this->m_descent == Descent::BFGS && iter > 0) {
          B = this->bfgs_update(B, alpha * p_x, g - g_old);
        }

        // SOLUTION TO PERTURBED KKT SYSTEM
        // Compute the search direction of x and z
        c_epsilon = c.array() - this->m_epsilon;
        S   = (z.array() / c_epsilon.array()).matrix().asDiagonal();
        g_b = g - mu * J.transpose() * (1.0 / c_epsilon.array()).matrix();
        H   = B + W - J.transpose() * S * J;
        p_x = H.ldlt().solve(-g_b);
        p_z = -(z + mu * (1.0 / c_epsilon.array()).matrix() + S * J * p_x);

        // BACKTRACKING LINE SEARCH
        // To ensure global convergence, execute backtracking line search to determine the step length
        // First, we have to find the largest step size which ensures that z remains feasible
        // Next, we perform backtracking line search
        alpha = this->m_alpha_max;
        mask = (z + p_z).array() < 0.0;
        if (mask.any()) {
          ratio = z.array() / (-p_z.array());
          alpha = this->m_alpha_max * std::min<Real>(1.0, mask.select(ratio, 1.0).minCoeff());
        }

        // Compute the response of the merit function and the directional gradient at the current
        // point and search direction
        psi = this->merit(x, z, f, c, mu);
        dpsi = this->grad_merit(x, z, p_x, p_z, g, c, J, mu);
        while (true) {

          // Compute the candidate point, the constraints, and the response of the objective function
          // and merit function at the candidate point
          x_new = x + alpha * p_x;
          z_new = z + alpha * p_z;
          f = this->m_problem->objective(x_new);
          c = this->m_problem->constraints(x_new);
          psi_new = this->merit(x_new, z_new, f, c, mu);

          // Stop backtracking search if we've found a candidate point that sufficiently decreases
          // the merit function and satisfies all the constraints
          if ((c.array() <= 0.0).all() && psi_new < psi + this->m_tau * eta * alpha * dpsi) {
            x = x_new;
            z = z_new;
            g_old = g;
            break;
          }

          // The candidate point does not meet our criteria, so decrease the step size for 0 < β < 1.
          alpha *= this->m_beta;
          if (alpha < this->m_alpha_min) {
            IPSOLVER_ERROR(CMD "line search step size too small");
            return false;
          }
        }
      }
      x_sol = x;
      return converged;

      #undef CMD
    }

private:
    /**
    * \brief Computes the merit function.
    *
    * This function computes the merit function
    * \f[
    *   \psi(\mathbf{x}, \mathbf{z}) = f(\mathbf{x}) - c(\mathbf{z})^\top \mathbf{z}
    *    -\mu \sum\log(\mathbf{c}^2 \mathbf{z} + \epsilon)
    * \f]
    * It is used to evaluate the quality of the current solution in terms of both the objective function and the constraints.
    * \param[in] x Primal variable vector.
    * \param[in] z Dual variable vector.
    * \param[in] f Objective function value at x.
    * \param[in] c Constraints vector at x.
    * \param[in] mu Barrier parameter.
    * \return The computed merit value.
    */
    Real merit([[maybe_unused]] VectorN const & x, VectorM const & z, Real f, VectorM const & c,
      Real const mu)
    {
      return f - c.dot(z) - mu * ((c.array().square() * z.array() + this->m_epsilon).log().sum());
    }

    /**
    * \brief Computes the directional derivative of the merit function.
    *
    * This function computes the directional derivative of the merit function with respect to the
    * primal and dual variables.
    * \param[in] x Primal variable vector.
    * \param[in] z Dual variable vector.
    * \param[in] p_x Directional derivative of the primal variable.
    * \param[in] p_z Directional derivative of the dual variable.
    * \param[in] g Gradient of the objective function at x.
    * \param[in] c Constraints vector at x.
    * \param[in] J Jacobian matrix of the constraints.
    * \param[in] mu Barrier parameter.
    * \return The computed directional derivative of the merit function.
    */
    Real grad_merit([[maybe_unused]] VectorN const & x, VectorM const & z, VectorN const & p_x,
      VectorM const & p_z, VectorN const & g, VectorM const & c, MatrixJ const & J, Real const mu)
    {
      return p_x.dot(
          g - J.transpose() * z - 2.0 * mu * J.transpose() * (1.0 / (c.array() - this->m_epsilon)).matrix()
        ) - p_z.dot(
          c + mu * (1.0 / (z.array() + this->m_epsilon)).matrix()
        );
    }

    /**
    * \brief Browder-Broyden-Fletcher-Goldfarb-Shanno (BFGS) update for the Hessian approximation.
    *
    * This method updates the Hessian approximation using the Browder-Broyden-Fletcher-Goldfarb-Shanno
    * (BFGS) formula*
    * \f[
    *   \mathbf{B}_{k+1} = \mathbf{B}_{k} - \displaystyle\frac{(\mathbf{B}_{k}\mathbf{s}_{k})(
    *     \mathbf{B}_{k}\mathbf{s}_{k}^\top)}{\mathbf{s}_{k}^\top \mathbf{B}_{k}\mathbf{s}_{k}} +
    *     \displaystyle\frac{\mathbf{y}\mathbf{y}^\top)}{\mathbf{y}^\top\mathbf{s}_{k}}
    * \f]
    * where \f$B\f$ is the current Hessian approximation, \f$s\f$ is the step taken,
    * and \f$y\f$ is the gradient difference.
    * \note The condition \f$y^\top s > 0\f$ must be satisfied for the update to be valid.
    * \param[in] B Current Hessian approximation.
    * \param[in] s Step taken (s_{k} = x_{k+1} - x_{k}).
    * \param[in] y Gradient difference (g_new - g).
    * \return Updated Hessian approximation.
    */
    MatrixH bfgs_update(MatrixH const & B, VectorN const & s, VectorN const & y)
    {
      IPSOLVER_ASSERT(y.dot(s) > 0.0,
        "IPsolver::Solver::bfgs_update(...): update condition yᵀs > 0 not satisfied");
      VectorN x(B * s);
      return B - (x * x.transpose()) / (s.dot(x)) + (y * y.transpose()) / (y.dot(s));
    }

  }; // class IPSolver

} // namespace IPsolver

#endif /* INCLUDE_IPSOLVER_SOLVER_HH */
