/**
 * PYTHON Wrapper for Clothoids
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#ifdef IPOPT_CLOTHOID_SPLINE

#include "python-ClothoidSpline-IpoptSolver.hh"

#include <IpIpoptApplication.hpp>
#include <IpSmartPtr.hpp>
#include <algorithm>
#include <stdexcept>

namespace G2lib {
  namespace Interpolation {

    using Ipopt::ApplicationReturnStatus;
    using Ipopt::Index;
    using Ipopt::IpoptApplication;
    using Ipopt::Number;
    using Ipopt::SmartPtr;
    using Ipopt::TNLP;

    IpoptSolver::ClothoidSplineProblem::ClothoidSplineProblem(IpoptSolver & solver)
        : m_solver(solver), m_solver_return(SolverReturn::INVALID_OPTION) {}

    bool IpoptSolver::ClothoidSplineProblem::get_nlp_info(
        Index &                theta_size,
        Index &                constraints_size,
        Index &                jacobian_pattern_size,
        Index &                hessian_pattern_size,
        TNLP::IndexStyleEnum & pattern_index_style) {
      theta_size            = m_solver.theta_size();
      constraints_size      = m_solver.constraints_size();
      jacobian_pattern_size = m_solver.jacobian_pattern_size();
      hessian_pattern_size  = m_solver.lagrangian_hessian_size();
      pattern_index_style   = TNLP::C_STYLE;
      return true;
    }

    bool IpoptSolver::ClothoidSplineProblem::get_bounds_info(
        Index    theta_size,
        Number * theta_min,
        Number * theta_max,
        Index    constraints_size,
        Number * constraints_min,
        Number * constraints_max) {
      if (theta_size != m_solver.theta_size())
        return false;

      if (constraints_size != m_solver.constraints_size())
        return false;

      if ((theta_min == nullptr) || (theta_max == nullptr) || (constraints_min == nullptr) ||
          (constraints_max == nullptr))
        return false;

      std::fill_n(constraints_min, constraints_size, 0.0);
      std::fill_n(constraints_max, constraints_size, 0.0);
      std::copy(m_solver.theta_min().begin(), m_solver.theta_min().end(), theta_min);
      std::copy(m_solver.theta_max().begin(), m_solver.theta_max().end(), theta_max);
      return true;
    }

    bool IpoptSolver::ClothoidSplineProblem::get_starting_point(
        Index    theta_size,
        bool     init_theta,
        Number * theta_ics,
        bool     init_z,
        Number * z_L,
        Number * z_U,
        Index    constraints_size,
        bool     init_lambda,
        Number * lambda) {
      if (theta_size != m_solver.theta_size())
        return false;

      if (constraints_size != m_solver.constraints_size())
        return false;

      if (init_z)
        return false;

      if (init_lambda)
        return false;

      if (init_theta) {
        std::copy(m_solver.theta_solution().begin(), m_solver.theta_solution().end(), theta_ics);
      }
      return true;
    }

    bool IpoptSolver::ClothoidSplineProblem::eval_f(
        Index theta_size, const Number * theta, bool new_theta, Number & obj_value) {
      if (theta_size != m_solver.theta_size())
        return false;
      return m_solver.spline().objective(theta, obj_value);
    }

    bool IpoptSolver::ClothoidSplineProblem::eval_grad_f(
        Index theta_size, const Number * theta, bool new_theta, Number * grad_f) {
      if (theta_size != m_solver.theta_size())
        return false;
      return m_solver.spline().gradient(theta, grad_f);
    }

    bool IpoptSolver::ClothoidSplineProblem::eval_g(
        Index theta_size, const Number * theta, bool new_theta, Index constraints_size, Number * g) {
      if (theta_size != m_solver.theta_size())
        return false;
      if (constraints_size != m_solver.constraints_size())
        return false;
      return m_solver.spline().constraints(theta, g);
    }

    bool IpoptSolver::ClothoidSplineProblem::eval_jac_g(
        Index          theta_size,
        const Number * theta,
        bool           new_theta,
        Index          constraints_size,
        Index          jacobian_pattern_size,
        Index *        jacobian_rows,
        Index *        jacobian_cols,
        Number *       jacobian_values) {
      if (theta_size != m_solver.theta_size())
        return false;
      if (constraints_size != m_solver.constraints_size())
        return false;
      if (jacobian_pattern_size != m_solver.jacobian_pattern_size())
        return false;

      bool index_ok = true;
      if ((jacobian_rows != NULL) && (jacobian_cols != NULL)) {
        index_ok = m_solver.spline().jacobian_pattern(jacobian_rows, jacobian_cols);
      }

      bool jac_eval_ok = true;
      if (jacobian_values != NULL) {
        jac_eval_ok = m_solver.spline().jacobian(theta, jacobian_values);
      }

      return index_ok & jac_eval_ok;
    }

    void IpoptSolver::ClothoidSplineProblem::finalize_solution(
        SolverReturn                       status,
        Index                              theta_size,
        const Number *                     theta,
        const Number *                     z_L,
        const Number *                     z_U,
        Index                              constraints_size,
        const Number *                     g,
        const Number *                     lambda,
        Number                             obj_value,
        const Ipopt::IpoptData *           ip_data,
        Ipopt::IpoptCalculatedQuantities * ip_cq) {
      m_solver_return = status;
      std::copy(theta, theta + theta_size, m_solver.theta_solution().begin());
    }

    bool IpoptSolver::solve() {
      SmartPtr<IpoptSolver::ClothoidSplineProblem> spline_problem = new IpoptSolver::ClothoidSplineProblem(*this);
      SmartPtr<IpoptApplication>                   app            = IpoptApplicationFactory();

      app->Options()->SetStringValue("jac_d_constant", "no");
      app->Options()->SetStringValue("hessian_constant", "no");
      app->Options()->SetStringValue("mu_strategy", "adaptive");
      app->Options()->SetStringValue("derivative_test", "none");
      app->Options()->SetStringValue("hessian_approximation", "limited-memory");
      app->Options()->SetStringValue("limited_memory_update_type", "bfgs");

      app->Options()->SetIntegerValue("max_iter", 400);

      app->Options()->SetNumericValue("tol", 1e-10);
      app->Options()->SetNumericValue("derivative_test_tol", 1e-5);

      // Disable all output
      app->Options()->SetIntegerValue("print_level", 0);
      app->Options()->SetStringValue("sb", "yes");

      ApplicationReturnStatus status = app->Initialize();
      if (status != ApplicationReturnStatus::Solve_Succeeded)
        return false;

      status = app->OptimizeTNLP(spline_problem);
      switch (status) {
        case (ApplicationReturnStatus::Solve_Succeeded):
          return true;
          break;
        case (ApplicationReturnStatus::Solved_To_Acceptable_Level):
          return true;
          break;
        case (ApplicationReturnStatus::Infeasible_Problem_Detected):
          std::runtime_error("IPOPT Error :: Infeasible_Problem_Detected");
          break;
        case (ApplicationReturnStatus::Search_Direction_Becomes_Too_Small):
          std::runtime_error("IPOPT Error :: Search_Direction_Becomes_Too_Small");
          break;
        case (ApplicationReturnStatus::Diverging_Iterates):
          std::runtime_error("IPOPT Error :: Diverging_Iterates");
          break;
        case (ApplicationReturnStatus::User_Requested_Stop):
          std::runtime_error("IPOPT Error :: User_Requested_Stop");
          break;
        case (ApplicationReturnStatus::Feasible_Point_Found):
          std::runtime_error("IPOPT Error :: Feasible_Point_Found");
          break;
        case (ApplicationReturnStatus::Maximum_Iterations_Exceeded):
          std::runtime_error("IPOPT Error :: Maximum_Iterations_Exceeded");
          break;
        case (ApplicationReturnStatus::Restoration_Failed):
          std::runtime_error("IPOPT Error :: Restoration_Failed");
          break;
        case (ApplicationReturnStatus::Error_In_Step_Computation):
          std::runtime_error("IPOPT Error :: Error_In_Step_Computation");
          break;
        case (ApplicationReturnStatus::Maximum_CpuTime_Exceeded):
          std::runtime_error("IPOPT Error :: Maximum_CpuTime_Exceeded");
          break;
        case (ApplicationReturnStatus::Not_Enough_Degrees_Of_Freedom):
          std::runtime_error("IPOPT Error :: Not_Enough_Degrees_Of_Freedom");
          break;
        case (ApplicationReturnStatus::Invalid_Problem_Definition):
          std::runtime_error("IPOPT Error :: Invalid_Problem_Definition");
          break;
        case (ApplicationReturnStatus::Invalid_Option):
          std::runtime_error("IPOPT Error :: Invalid_Option");
          break;
        case (ApplicationReturnStatus::Invalid_Number_Detected):
          std::runtime_error("IPOPT Error :: Invalid_Number_Detected");
          break;
        case (ApplicationReturnStatus::Unrecoverable_Exception):
          std::runtime_error("IPOPT Error :: Unrecoverable_Exception");
          break;
        case (ApplicationReturnStatus::NonIpopt_Exception_Thrown):
          std::runtime_error("IPOPT Error :: NonIpopt_Exception_Thrown");
          break;
        case (ApplicationReturnStatus::Insufficient_Memory):
          std::runtime_error("IPOPT Error :: Insufficient_Memory");
          break;
        case (ApplicationReturnStatus::Internal_Error):
          std::runtime_error("IPOPT Error :: Internal_Error");
          break;
      }
      return false;
    }

  } /* namespace Interpolation */
} /* namespace G2lib */

#endif