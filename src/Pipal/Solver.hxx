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

#ifndef INCLUDE_PIPAL_SOLVER_HXX
#define INCLUDE_PIPAL_SOLVER_HXX

namespace Pipal {

  /**
   * \brief Solver class for the Pipal library.
   *
   * The Solver class provides an interface for solving the optimization problems using Frank E.
   * Curtis Pipal algorithm. It utilizes the Problem class to define the optimization problem and
   * implements various methods for solving it.
   * \tparam Real The floating-point type used for computations (e.g., float, double).
   */
  template <typename Real>
  class Solver
  {
    static_assert(std::is_floating_point_v<Real>,
      "Pipal::Solver<Real>: Real must be a floating-point type.");

  public:
    using ProblemPtr              = typename Problem<Real>::UniquePtr;
    using ObjectiveFunc           = typename ProblemWrapper<Real>::ObjectiveFunc;
    using ObjectiveGradientFunc   = typename ProblemWrapper<Real>::ObjectiveGradientFunc;
    using ConstraintsFunc         = typename ProblemWrapper<Real>::ConstraintsFunc;
    using ConstraintsJacobianFunc = typename ProblemWrapper<Real>::ConstraintsJacobianFunc;
    using LagrangianHessianFunc   = typename ProblemWrapper<Real>::LagrangianHessianFunc;
    using BoundsFunc              = typename ProblemWrapper<Real>::BoundsFunc;

  private:
    Counter          m_counter;    /*!< Internal counters for solver statistics. */
    Acceptance<Real> m_acceptance; /*!< Acceptance criteria for trial points. */
    Input<Real>      m_input;      /*!< Input structure for the solver. */
    Direction<Real>  m_direction;  /*!< Current search direction of the solver. */
    Iterate<Real>    m_iterate;    /*!< Current iterate of the solver. */
    Output<Real>     m_output;     /*!< Output class for managing solver output. */
    Parameter<Real>  m_parameter;  /*!< Internal parameters for the solver algorithm. */
    ProblemPtr       m_problem;    /*!< Problem object pointer. */

    // Some options for the solver
    bool m_verbose{false}; /*!< Verbosity flag. */
    bool m_bfgs{false};    /*!< BFGS update flag. */

    void buildIterate();
    void evalStep();
    void updateParameters();
    void lineSearch();
    void bfgsUpdate(Vector<Real> const & s, Vector<Real> const & y);
    void updateIterate();
    void updatePoint();
    void backtracking();
    void evalFunctions();
    void evalGradients();
    void evalHessian();
    void evalNewtonMatrix();
    void evalScalings();
    void evalModels();
    void fractionToBoundary();
    void evalNewtonRhs();
    void evalSlacks();
    void evalMerit();
    void initNewtonMatrix();
    void resetDirection(Direction<Real> & d) const;
    void evalLambdaOriginal(Vector<Real> & l) const;
    Real evalViolation(Array<Real> const & cE, Array<Real> const & cI) const;
    void evalNewtonStep();
    void evalTrialStep(Direction<Real> & v) const;
    Integer checkTermination() const;
    Integer fullStepCheck();
    Integer secondOrderCorrection();
    void evalXOriginal(Vector<Real> & x);
    void setDirection(Vector<Real> const & dx, Vector<Real> const & dr1, Vector<Real> const & dr2,
      Vector<Real> const & ds1, Vector<Real> const & ds2, Vector<Real> const & dlE,
      Vector<Real> const & dlI, Real const dx_norm, Real const dl_norm);
    void evalTrialSteps(Direction<Real> & d1, Direction<Real> & d2, Direction<Real> & d3);
    void evalTrialStepCut();
    void evalLinearCombination(Direction<Real> const & d1, Direction<Real> const & d2,
      Direction<Real> const & d3, Real const a1, Real const a2, Real const a3);
    Real evalKKTError(Real const rho, Real const mu);
    void evalKKTErrors();
    void setPrimals(Vector<Real> const & x, Array<Real> const & r1, Array<Real> const & r2,
      Array<Real> const & s1, Array<Real> const & s2, Array<Real> const & lE, Array<Real> const & lI,
      Real const f, Array<Real> const & cE, Array<Real> const & cI, Real const phi);
    void buildInput(std::string const & name, Vector<Real> const & x0, Vector<Real> const & bl,
      Vector<Real> const & bu, Vector<Real> const & cl, Vector<Real> const & cu);

    /**
     * \brief Compute scaled and unscaled feasibility violations.
     * \tparam Real Floating-point type used by the algorithm.
     * \param[in] z Current iterate.
     */
    void evalInfeasibility(Iterate<Real> & z) const
    {
      // Evaluate scaled and unscaled feasibility violations
      z.v  = this->evalViolation(z.cE, z.cI)/std::max(1.0, z.v0);
      z.vu = this->evalViolation(z.cEu, z.cIu);
    }

    /**
     * \brief Reset maximum exponent used for mu increases to its default.
     * \tparam Real Floating-point type used by the algorithm.
     */
    void resetMuMaxExp() {this->m_parameter.mu_max_exp = this->m_parameter.mu_max_exp0;}

    /**
     * \brief Initialize algorithm parameters.
     * \param[in] t_algorithm Algorithm selection enumerator.
     */
    void buildParameter(Algorithm t_algorithm) {this->m_parameter.algorithm = t_algorithm;}

    /**
     * \brief Force mu exponent increases to use zero as maximum exponent.
     */
    void setMuMaxExpZero() {this->m_parameter.mu_max_exp = 0.0;}

    /**
     * \brief Set penalty parameter rho.
     * \tparam Real Floating-point type used by the algorithm.
     * \param[in] rho New penalty parameter value.
     */
    void setRho(Real const rho) {this->m_iterate.rho = rho;}

    /**
     * \brief Set last (previous) penalty parameter value.
     * \tparam Real Floating-point type used by the algorithm.
     * \param[in] rho New last penalty parameter value.
     */
    void setRhoLast(Real const rho) {this->m_iterate.rho_ = rho;}

    /**
     * \brief Set interior-point parameter \p mu.
     * \tparam Real Floating-point type used by the algorithm.
     * \param[in] mu New interior-point parameter value.
     */
    void setMu(Real const mu) {this->m_iterate.mu = mu;}

    /**
     * \brief Evaluate quantities that depend on penalty/interior parameters.
     *
     * Runs slack, merit and KKT-error evaluations that are functions of the current penalty and
     * interior-point parameters.
     */
    void evalDependent()
    {
      // Evaluate quantities dependent on penalty and interior-point parameters
      this->evalSlacks();
      this->evalMerit();
      this->evalKKTErrors();
    }

    /**
     * \brief Reset all internal counters to zero.
     */
    void resetCounter()
    {
      this->m_counter.f = this->m_counter.g = this->m_counter.H = this->m_counter.k = this->m_counter.M = 0;
    }

    /**
     * \brief Increment the matrix factorization counter.
     */
    void incrementFactorizationCount() {++this->m_counter.M;}

    /**
     * \brief Increment the function evaluation counter.
     */
    void incrementFunctionCount() {++this->m_counter.f;}

    /**
     * \brief Increment the gradient evaluation counter.
     */
    void incrementGradientCount() {++this->m_counter.g;}

    /**
     * \brief Increment the Hessian evaluation counter.
     */
    void incrementHessianCount() {++this->m_counter.H;}

    /**
     * \brief Increment the iteration counter.
     */
    void incrementIterationCount() {++this->m_counter.k;}

  public:
    /**
     * \brief Default constructor for the Pipal class.
     *
     * Initializes the solver with default values for the objective, gradient, constraints, and Jacobian
     * functions.
     */
    Solver() = default;

    /**
     * \brief Constructor for the Pipal class.
     *
     * Initializes the solver with the provided objective, gradient, constraints, and Jacobian functions.
     * \param[in] name Name of the optimization problem.
     * \param[in] objective Objective function handle.
     * \param[in] objective_gradient Gradient of the objective function handle.
     * \param[in] constraints Constraints function handle.
     * \param[in] constraints_jacobian Jacobian of the constraints function handle.
     * \param[in] lagrangian_hessian Hessian of the Lagrangian function handle.
     * \param[in] primal_lower_bounds Lower bounds on the primal variables handle.
     * \param[in] primal_upper_bounds Upper bounds on the primal variables handle.
     * \param[in] constraints_lower_bounds Lower bounds on the constraints handle.
     * \param[in] constraints_upper_bounds Upper bounds on the constraints handle.
     */
    Solver(
      std::string             const & name,
      ObjectiveFunc           const & objective,
      ObjectiveGradientFunc   const & objective_gradient,
      ConstraintsFunc         const & constraints,
      ConstraintsJacobianFunc const & constraints_jacobian,
      LagrangianHessianFunc   const & lagrangian_hessian,
      BoundsFunc              const & primal_lower_bounds,
      BoundsFunc              const & primal_upper_bounds,
      BoundsFunc              const & constraints_lower_bounds,
      BoundsFunc              const & constraints_upper_bounds
    ) : m_problem(std::make_unique<ProblemWrapper<Real>>(
      name,
      objective,
      objective_gradient,
      constraints,
      constraints_jacobian,
      lagrangian_hessian,
      primal_lower_bounds,
      primal_upper_bounds,
      constraints_lower_bounds,
      constraints_upper_bounds
    )) {}

    /**
     * \brief Constructor for the Pipal class (with a unique pointer to a Problem object).
     *
     * Initializes the solver with the provided unique pointer to a Problem object.
     * \param[in] problem The unique pointer to the Problem object to use.
     * \warning The pointer is moved into the solver, so it should not be used after this call.
     */
    Solver(ProblemPtr && problem) : m_problem(std::move(problem)) {}

    /**
     * \brief Deleted copy constructor.
     * \note This class is not copyable.
     */
    Solver(Solver const &) = delete;

    /**
     * \brief Deleted assignment operator.
     * \note This class is not assignable.
     */
    Solver & operator=(Solver const &) = delete;

    /**
     * \brief Deleted move constructor.
     * \note This class is not movable.
     */
    Solver(Solver &&) = delete;

    /**
     * \brief Deleted move assignment operator.
     * \note This class is not movable.
     */
    Solver & operator=(Solver &&) = delete;

    /**
     * \brief Destructor for the Pipal class.
     *
     * Cleans up resources used by the solver.
     */
    ~Solver() = default;

    /**
     * \brief Set the problem to be solved using a unique pointer.
     *
     * This method allows the user to specify the problem to be solved using a unique pointer.
     * \param[in] problem The unique pointer to the problem to set.
     */
    void problem(ProblemPtr && problem) {this->m_problem = std::move(problem);}

    /**
     * \brief Get the problem being solved.
     * \return A reference to the problem.
     */
    Problem<Real> const & problem() const {return this->m_problem.get();}

    /**
     * \brief Get the verbose mode.
     * \return The verbose mode.
     */
    bool verbose_mode() const {return this->m_verbose;}

    /**
     * \brief Set the verbose mode.
     * \param[in] t_verbose The verbose mode.
     */
    void verbose_mode(bool const t_verbose) {this->m_verbose = t_verbose;}

    /**
     * \brief Get the BFGS mode.
     * \return The BFGS mode.
     */
    bool bfgs() const {return this->m_bfgs;}

    /**
     * \brief Set the BFGS mode.
     * \param[in] t_bfgs The BFGS mode.
     */
    void bfgs(bool const t_bfgs) {this->m_bfgs = t_bfgs;}

    /**
     * \brief Get the algorithm mode.
     * \return The algorithm mode.
     */
    Algorithm algorithm() const {return this->m_parameter.algorithm;}

    /**
     * \brief Set the algorithm mode.
     *
     * This method allows the user to specify the algorithm mode, which can be either \c CONSERVATIVE
     * or \c ADAPTIVE.
     * \param[in] t_algorithm The algorithm mode.
     */
    void algorithm(Algorithm const t_algorithm) {this->m_parameter.algorithm = t_algorithm;}

    /**
     * \brief Set the convergence tolerance for the solver.
     *
     * This method allows the user to specify the tolerance for convergence.
     * The solver will stop when the residuals are below this tolerance.
     *
     * \param[in] t_tolerance The convergence tolerance.
     */
    void tolerance(Real const t_tolerance)
    {
      PIPAL_ASSERT(t_tolerance > 0.0,
        "Pipal::Solver::tolerance(...): input value must be positive");
      this->m_parameter.opt_err_tol = t_tolerance;
    }

    /**
     * \brief Get the tolerance for convergence.
     * \return The tolerance for convergence.
     */
    Real tolerance() const {return this->m_parameter.opt_err_tol;}

    /**
     * \brief Set the maximum number of iterations for the solver.
     *
     * This method allows the user to specify the maximum number of iterations
     * the solver will perform before stopping.
     *
     * \param[in] t_max_iterations The maximum number of iterations.
     */
    void max_iterations(Integer const t_max_iterations) {
      PIPAL_ASSERT(t_max_iterations > 0,
        "Pipal::Solver::max_iterations(...): input value must be positive");
      this->m_parameter.iter_max = t_max_iterations;
    }

    /**
     * \brief Get the maximum number of iterations.
     * \return The maximum number of iterations.
     */
    Integer max_iterations() const {return this->m_parameter.iter_max;}

    /**
     * \brief Solves the optimization problem using the interior-point method.
     *
     * This method implements the interior-point algorithm to solve the optimization problem defined
     * by the objective function, constraints, and their respective gradients and Jacobians.
     * \param[in] x_guess Initial guess for the optimization variables.
     * \param[out] x_sol Solution vector.
     * \return True if the optimization was successful, false otherwise.
     */
    bool optimize(Vector<Real> const & x_guess, Vector<Real> & x_sol)
    {

      #define CMD "Pipal::Solver::optimize(...): "

      // Create alias for easier access
      Counter          & c{this->m_counter};
      Input<Real>      & i{this->m_input};
      Iterate<Real>    & z{this->m_iterate};
      Direction<Real>  & d{this->m_direction};
      Acceptance<Real> & a{this->m_acceptance};

      // Check that the problem is set
      PIPAL_ASSERT(this->m_problem.get() != nullptr,
        CMD "problem not set, use 'problem(...)' method to set it");

      // Get variable bounds
      Vector<Real> bl, bu;
      PIPAL_ASSERT(this->m_problem->primal_lower_bounds(bl),
        CMD "error in evaluating lower bounds on primal variables");
      PIPAL_ASSERT(this->m_problem->primal_upper_bounds(bu),
        CMD "error in evaluating upper bounds on primal variables");

      // Get constraint bounds
      Vector<Real> cl, cu;
      PIPAL_ASSERT(this->m_problem->constraints_lower_bounds(cl),
        CMD "error in evaluating lower bounds on constraints");
      PIPAL_ASSERT(this->m_problem->constraints_upper_bounds(cu),
        CMD "error in evaluating upper bounds on constraints");

      // Reset counters
      resetCounter();

      // Fill input structure
      buildInput(this->m_problem->name(), x_guess, bl, bu, cl, cu);
      buildIterate();
      resetDirection(d);

      // Print header and break line
      if (this->m_verbose) {this->m_output.printHeader(i, z); this->m_output.printBreak(c);}

      // Iterations loop
      while (!this->checkTermination()) {

        // Print iterate
        if (this->m_verbose) {this->m_output.printIterate(c, z);}

        // Evaluate the step
        this->evalStep();

        // Print direction
        if (this->m_verbose) {this->m_output.printDirection(z, d);}

        this->lineSearch();

        // Print accepted
        if (this->m_verbose) {this->m_output.printAcceptance(a);}

        this->updateIterate();

        // Increment iteration counter
        this->incrementIterationCount();

        // Print break
        if (this->m_verbose) {this->m_output.printBreak(c);}
      }
      // Print footer and terminate
      if (this->m_verbose) {this->m_output.printFooter(c, z, this->checkTermination());}

      // Get solution in original variables
      this->evalXOriginal(x_sol);

      // Return success if is finite
      return x_sol.allFinite();

      #undef CMD
    }

    /**
     * \brief Extract a primal-dual solution in original ordering.
     * \tparam Real Floating-point type used by the algorithm.
     * \param[out] x Vector to store primal variables in original space.
     * \param[out] l Vector to store multipliers in original space.
     */
    void getSolution(Vector<Real> & x, Vector<Real> & l)
    {
      this->evalXOriginal(x);
      this->evalLambdaOriginal(l);
    }

  }; // class Solver

} // namespace Solver

#endif // INCLUDE_PIPAL_SOLVER_HXX
