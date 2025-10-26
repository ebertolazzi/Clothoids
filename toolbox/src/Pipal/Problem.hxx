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

#ifndef INCLUDE_PIPAL_PROBLEM_HXX
#define INCLUDE_PIPAL_PROBLEM_HXX

namespace Pipal {
  /**
   * \brief Problem class for the Pipal library.
   *
   * The Problem class serves as a base class for defining optimization problems in the Pipal library.
   * It provides a structure for encapsulating the problem's objective function, constraints, and
   * other necessary components.
    * \tparam Real The real number type.
   */
  template <typename Real>
  class Problem  {
    static_assert(
      std::is_floating_point_v<Real>,
      "Pipal::Problem<Real>: Real must be a floating-point type."
    );

    std::string m_name{"(Unnamed Pipal Problem)"}; /*!< Name of the optimization problem. */

  public:
    using UniquePtr = std::unique_ptr<Problem>;

    /**
     * \brief Default constructor.
     *
     * Initializes the problem with default values for the objective, gradient, hessian, constraints
     * and Jacobian functions.
     */
    Problem() = default;

    /**
     * \brief Problem constructor.
     * \param[in] t_name Name of the optimization problem.
     */
    Problem(std::string t_name) : m_name(std::move(t_name)) {};

    /**
     * \brief Deleted copy constructor.
     * \note This class is not copyable.
     */
    Problem(Problem const &) = delete;

    /**
     * \brief Deleted assignment operator.
     * \note This class is not assignable.
     */
    Problem & operator=(Problem const &) = delete;

    /**
     * \brief Deleted move constructor.
     * \note This class is not movable.
     */
    Problem(Problem &&) = delete;

    /**
     * \brief Deleted move assignment operator.
     * \note This class is not movable.
     */
    Problem & operator=(Problem &&) = delete;

    /**
     * \brief Default destructor.
     */
    virtual ~Problem() = default;

    /**
     * \brief Get the name of the optimization problem.
     * \return The name of the optimization problem.
     */
    [[nodiscard]] std::string const & name() const {return this->m_name;}

    /**
     * \brief Set the name of the optimization problem.
     * \param[in] t_name The name of the optimization problem.
     */
    void name(std::string const & t_name) {this->m_name = t_name;}

    /**
     * \brief Evaluate the objective function.
     * \param[in] x Primal variables.
     * \param[out] out The objective function.
     * \return True if the evaluation was successful, false otherwise.
     */
    virtual bool objective(Vector<Real> const & x, Real & out) const = 0;

    /**
     * \brief Evaluate the gradient of the objective function.
     * \param[in] x Primal variables.
     * \param[out] out The gradient of the objective function.
     * \return True if the evaluation was successful, false otherwise.
     */
    virtual bool objective_gradient(Vector<Real> const & x, Vector<Real> & out) const = 0;

    /**
     * \brief Evaluate the constraints function.
     * \param[in] x Primal variables.
     * \param[out] out The value of the constraints function.
     * \return True if the evaluation was successful, false otherwise.
     */
    virtual bool constraints(Vector<Real> const & x, Vector<Real> & out) const = 0;

    /**
     * \brief Evaluate the Jacobian of the constraints function with respect to the primal variables.
     * \param[in] x Primal variables.
     * \param[out] out The Jacobian matrix of the constraints function.
     * \return True if the evaluation was successful, false otherwise.
     */
    virtual bool constraints_jacobian(Vector<Real> const & x, SparseMatrix<Real> & out) const = 0;

    /**
     * \brief Evaluate the Hessian of the Lagrangian function with respect to the primal variables.
     * \param[in] x Primal variables.
     * \param[in] z Dual variables.
     * \param[out] out The Hessian matrix of the Lagrangian function.
     * \return True if the evaluation was successful, false otherwise.
     */
    virtual bool lagrangian_hessian(Vector<Real> const & x, Vector<Real> const & z, SparseMatrix<Real> & out) const = 0;

    /**
     * \brief Lower bounds on the primal variables.
     * \param[out] out The lower bounds on the primal variables.
     * \return True if the evaluation was successful, false otherwise.
     */
    virtual bool primal_lower_bounds(Vector<Real> & out) const = 0;

    /**
     * \brief Upper bounds on the primal variables.
     * \param[out] out The upper bounds on the primal variables.
     * \return True if the evaluation was successful, false otherwise.
     */
    virtual bool primal_upper_bounds(Vector<Real> & out) const = 0;

    /**
     * \brief Lower bounds on the constraints.
     * \param[out] out The lower bounds on the constraints.
     * \return True if the evaluation was successful, false otherwise.
     */
    virtual bool constraints_lower_bounds(Vector<Real> & out) const = 0;

    /**
     * \brief Upper bounds on the constraints.
     * \param[out] out The upper bounds on the constraints.
     * \return True if the evaluation was successful, false otherwise.
     */
    virtual bool constraints_upper_bounds(Vector<Real> & out) const = 0;

  }; // class Problem

  /**
  * \brief Wrapper class for the Problem class.
  *
  * The ProblemWrapper class is a simple wrapper around the Problem class. It inherits from the
  * Problem class and can be used to define optimization problems through function handles.
  * \tparam Real The real number type.
  */
  template <typename Real>
  class ProblemWrapper : public Problem<Real>
  {
  public:
    using ObjectiveFunc           = std::function<bool(Vector<Real> const &, Real &)>;
    using ObjectiveGradientFunc   = std::function<bool(Vector<Real> const &, Vector<Real> &)>;
    using ConstraintsFunc         = std::function<bool(Vector<Real> const &, Vector<Real> &)>;
    using ConstraintsJacobianFunc = std::function<bool(Vector<Real> const &, SparseMatrix<Real> &)>;
    using LagrangianHessianFunc   = std::function<bool(Vector<Real> const &, Vector<Real> const &, SparseMatrix<Real> &)>;
    using BoundsFunc              = std::function<bool(Vector<Real> &)>;

  private:
    // Problem functions
    ObjectiveFunc           m_objective{nullptr};            /*!< Objective function. */
    ObjectiveGradientFunc   m_objective_gradient{nullptr};   /*!< Gradient of the objective function. */
    ConstraintsFunc         m_constraints{nullptr};          /*!< Constraints function. */
    ConstraintsJacobianFunc m_constraints_jacobian{nullptr}; /*!< Jacobian of the constraints. */
    LagrangianHessianFunc   m_lagrangian_hessian{nullptr};   /*!< Hessian of the Lagrangian. */

    // Bounds functions
    BoundsFunc m_primal_lower_bounds{nullptr};      /*!< Lower bounds on the primal variables. */
    BoundsFunc m_primal_upper_bounds{nullptr};      /*!< Upper bounds on the primal variables. */
    BoundsFunc m_constraints_lower_bounds{nullptr}; /*!< Lower bounds on the constraints. */
    BoundsFunc m_constraints_upper_bounds{nullptr}; /*!< Upper bounds on the constraints. */

  public:
    /**
     * \brief Constructor for the ProblemWrapper class (without the Hessian of the Lagrangian).
     *
     * Initializes the problem with the provided objective, gradient, constraints, and Jacobian
     * functions.
     * \param[in] t_name Name of the optimization problem.
     * \param[in] t_objective Objective function handle.
     * \param[in] t_objective_gradient Gradient of the objective function handle.
     * \param[in] t_constraints Constraints function handle.
     * \param[in] t_constraints_jacobian Jacobian of the constraints function handle.
     * \param[in] t_lagrangian_hessian Hessian of the Lagrangian function handle.
     * \param[in] t_primal_lower_bounds Lower bounds on the primal variables function handle.
     * \param[in] t_primal_upper_bounds Upper bounds on the primal variables function handle.
     * \param[in] t_constraints_lower_bounds Lower bounds on the constraints function handle.
     * \param[in] t_constraints_upper_bounds Upper bounds on the constraints function handle.
     */
    ProblemWrapper(
      std::string             const & t_name,
      ObjectiveFunc           const & t_objective,
      ObjectiveGradientFunc   const & t_objective_gradient,
      ConstraintsFunc         const & t_constraints,
      ConstraintsJacobianFunc const & t_constraints_jacobian,
      LagrangianHessianFunc   const & t_lagrangian_hessian,
      BoundsFunc              const & t_primal_lower_bounds,
      BoundsFunc              const & t_primal_upper_bounds,
      BoundsFunc              const & t_constraints_lower_bounds,
      BoundsFunc              const & t_constraints_upper_bounds
    )
    : Problem<Real>(t_name)
    , m_objective(t_objective)
    , m_objective_gradient(t_objective_gradient)
    , m_constraints(t_constraints)
    , m_constraints_jacobian(t_constraints_jacobian)
    , m_lagrangian_hessian(t_lagrangian_hessian)
    , m_primal_lower_bounds(t_primal_lower_bounds)
    , m_primal_upper_bounds(t_primal_upper_bounds)
    , m_constraints_lower_bounds(t_constraints_lower_bounds)
    , m_constraints_upper_bounds(t_constraints_upper_bounds)
    {}

    /**
     * \brief Default destructor for the ProblemWrapper class.
     */
    ~ProblemWrapper() override {};

    /**
     * \brief Get the objective function.
     * \return The objective function.
     */
    ObjectiveFunc &
    objective() { return this->m_objective; }

    /**
     * \brief Set the objective function.
     * \param[in] objective The objective function to set.
     */
    void
    objective( ObjectiveFunc const & objective )
    { this->m_objective = objective; }

    /**
     * \brief Get the gradient of the objective function.
     * \return The gradient of the objective function.
     */
    ObjectiveGradientFunc &
    objective_gradient() { return this->m_objective_gradient; }

    /**
     * \brief Set the gradient of the objective function.
     * \param[in] objective_gradient The gradient of the objective function to set.
     */
    void
    objective_gradient( ObjectiveGradientFunc const & objective_gradient )
    { this->m_objective_gradient = objective_gradient; }

    /**
     * \brief Get the constraints function.
     * \return The constraints function.
     */
    ConstraintsFunc &
    constraints() { return this->m_constraints; }

    /**
     * \brief Set the constraints function.
     * \param[in] constraints The constraints function to set.
     */
    void
    constraints( ConstraintsFunc const & constraints )
    { this->m_constraints = constraints; }

    /**
     * \brief Get the Jacobian of the constraints function.
     * \return The Jacobian of the constraints function.
     */
    ConstraintsJacobianFunc &
    constraints_jacobian()
    { return this->m_constraints_jacobian; }

    /**
     * \brief Set the Jacobian of the constraints function.
     * \param[in] constraints_jacobian The Jacobian of the constraints function to set.
     */
    void
    constraints_jacobian( ConstraintsJacobianFunc const & constraints_jacobian )
    { this->m_constraints_jacobian = constraints_jacobian; }

    /**
     * \brief Get the lower bounds on the primal variables function.
     * \return The lower bounds on the primal variables function.
     */
    BoundsFunc &
    primal_lower_bounds() { return this->m_primal_lower_bounds; }

    /**
     * \brief Set the lower bounds on the primal variables function.
     * \param[in] primal_lower_bounds The lower bounds on the primal variables function to set.
     */
    void
    primal_lower_bounds( BoundsFunc const & primal_lower_bounds )
    { this->m_primal_lower_bounds = primal_lower_bounds; }

    /**
     * \brief Get the upper bounds on the primal variables function.
     * \return The upper bounds on the primal variables function.
     */
    BoundsFunc &
    primal_upper_bounds()
    { return this->m_primal_upper_bounds; }

    /**
     * \brief Set the upper bounds on the primal variables function.
     * \param[in] primal_upper_bounds The upper bounds on the primal variables function to set.
     */
    void
    primal_upper_bounds( BoundsFunc const & primal_upper_bounds )
    { this->m_primal_upper_bounds = primal_upper_bounds; }

    /**
     * \brief Get the lower bounds on the constraints function.
     * \return The lower bounds on the constraints function.
     */
    BoundsFunc &
    constraints_lower_bounds()
    { return this->m_constraints_lower_bounds; }

    /**
     * \brief Set the lower bounds on the constraints function.
     * \param[in] constraints_lower_bounds The lower bounds on the constraints function to set.
     */
    void
    constraints_lower_bounds( BoundsFunc const & constraints_lower_bounds )
    { this->m_constraints_lower_bounds = constraints_lower_bounds; }

    /**
     * \brief Get the upper bounds on the constraints function.
     * \return The upper bounds on the constraints function.
     */
    BoundsFunc &
    constraints_upper_bounds()
    { return this->m_constraints_upper_bounds; }

    /**
     * \brief Set the upper bounds on the constraints function.
     * \param[in] constraints_upper_bounds The upper bounds on the constraints function to set.
     */
    void
    constraints_upper_bounds( BoundsFunc const & constraints_upper_bounds )
    { this->m_constraints_upper_bounds = constraints_upper_bounds; }

    /**
     * \brief Get the Hessian of the Lagrangian function.
     * \return The Hessian of the Lagrangian function.
     */
    LagrangianHessianFunc &
    lagrangian_hessian() { return this->m_lagrangian_hessian; }

    /**
     * \brief Set the Hessian of the Lagrangian function.
     * \param[in] lagrangian_hessian The Hessian of the Lagrangian function to set.
     */
    void
    lagrangian_hessian( LagrangianHessianFunc const & lagrangian_hessian )
    { this->m_lagrangian_hessian = lagrangian_hessian; }

    /**
     * \brief Evaluate the objective function.
     * \param[in] x Primal variables.
     * \param[out] out The objective function.
     * \return True if the evaluation was successful, false otherwise.
     */
    bool
    objective( Vector<Real> const & x, Real & out ) const override {
      return this->m_objective(x, out);
    }

    /**
     * \brief Evaluate the gradient of the objective function.
     * \param[in] x Primal variables.
     * \param[out] out The gradient of the objective function.
     * \return True if the evaluation was successful, false otherwise.
     */
    bool
    objective_gradient( Vector<Real> const & x, Vector<Real> & out ) const override {
      return this->m_objective_gradient(x, out);
    }

    /**
     * \brief Evaluate the constraints function.
     * \param[in] x Primal variables.
     * \param[out] out The value of the constraints function.
     * \return True if the evaluation was successful, false otherwise.
     */
    bool
    constraints( Vector<Real> const & x, Vector<Real> & out ) const override {
      return this->m_constraints(x, out);
    }

    /**
     * \brief Evaluate the Jacobian of the constraints function.
     * \param[in] x Primal variables.
     * \param[out] out The Jacobian matrix of the constraints function.
     * \return True if the evaluation was successful, false otherwise.
     */
    bool
    constraints_jacobian(
      Vector<Real> const & x,
      SparseMatrix<Real> & out
    ) const override {
      return this->m_constraints_jacobian(x, out);
    }

    /**
     * \brief Evaluate the Hessian of the Lagrangian function.
     * \param[in] x Primal variables.
     * \param[in] z Dual variables.
     * \param[out] out The Hessian matrix of the Lagrangian function.
     * \return True if the evaluation was successful, false otherwise.
     */
    bool
    lagrangian_hessian(
      Vector<Real> const & x,
      Vector<Real> const & z,
      SparseMatrix<Real> & out
    ) const override {
      return this->m_lagrangian_hessian(x, z, out);
    }

    /**
     * \brief Lower bounds on the primal variables.
     * \param[out] out The lower bounds on the primal variables.
     * \return True if the evaluation was successful, false otherwise.
     */
    bool primal_lower_bounds(Vector<Real> & out) const override
    {
      return this->m_primal_lower_bounds(out);
    }

    /**
     * \brief Upper bounds on the primal variables.
     * \param[out] out The upper bounds on the primal variables.
     * \return True if the evaluation was successful, false otherwise.
     */
    bool primal_upper_bounds(Vector<Real> & out) const override
    {
      return this->m_primal_upper_bounds(out);
    }

    /**
     * \brief Lower bounds on the constraints.
     * \param[out] out The lower bounds on the constraints.
     * \return True if the evaluation was successful, false otherwise.
     */
    bool constraints_lower_bounds(Vector<Real> & out) const override
    {
      return this->m_constraints_lower_bounds(out);
    }

    /**
     * \brief Upper bounds on the constraints.
     * \param[out] out The upper bounds on the constraints.
     * \return True if the evaluation was successful, false otherwise.
     */
    bool constraints_upper_bounds(Vector<Real> & out) const override
    {
      return this->m_constraints_upper_bounds(out);
    }

  }; // class ProblemWrapper

} // namespace Pipal

#endif /* INCLUDE_PIPAL_PROBLEM_HH */
