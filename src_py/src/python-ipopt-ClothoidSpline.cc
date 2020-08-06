/**
 * PYTHON Wrapper for Clothoids
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#ifdef IPOPT_CLOTHOID_SPLINE
#include "python-ipopt-ClothoidSpline.hh"

#include <IpIpoptApplication.hpp>
#include <IpSmartPtr.hpp>
#include <algorithm>
#include <stdexcept>
#include <vector>

namespace G2lib {

  namespace ipopt {

    using Ipopt::TNLP;
    using Ipopt::ApplicationReturnStatus;
    using Ipopt::Index;
    using Ipopt::IpoptApplication;
    using Ipopt::Number;
    using Ipopt::SmartPtr;

    using G2lib::ClothoidList;
    using G2lib::ClothoidSplineG2;
    using G2lib::int_type;
    using G2lib::real_type;

    ClothoidSplineProblem::ClothoidSplineProblem(const ClothoidSplineG2 &spline) : 
      TNLP(), m_spline(spline), m_theta_sol(std::vector< real_type >(m_spline.numTheta(), 0.0)),
      m_theta_guess(std::vector<real_type>(m_spline.numTheta())), m_theta_min(std::vector<real_type>(m_spline.numTheta())),
      m_theta_max(std::vector<real_type>(m_spline.numTheta())) {
      m_spline.guess(&m_theta_guess.front(), &m_theta_min.front(), &m_theta_max.front());
    }

    bool ClothoidSplineProblem::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag,
                                             TNLP::IndexStyleEnum &index_style) {
      n = m_spline.numTheta();
      m = m_spline.numConstraints();
      nnz_jac_g = m_spline.jacobian_nnz();
      nnz_h_lag = (n + m) * (n + m);  // Worst case scenario; we will se in future if we need other dimensions
      index_style = TNLP::C_STYLE;
      return true;
    }

    bool ClothoidSplineProblem::get_bounds_info(Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u) {
      if (n != m_spline.numTheta())
        return false;

      if (m != m_spline.numConstraints())
        return false;

      if ((x_l == nullptr) || (x_u == nullptr) || (g_l == nullptr) || (g_u == nullptr))
        return false;

      std::fill_n(g_l, m, 0.0);
      std::fill_n(g_u, m, 0.0);
      std::copy(m_theta_min.begin(), m_theta_min.end(), x_l);
      std::copy(m_theta_max.begin(), m_theta_max.end(), x_u);
      return true;
    }

    bool ClothoidSplineProblem::get_starting_point(Index n, bool init_x, Number *x, bool init_z, Number *z_L,
                                                   Number *z_U, Index m, bool init_lambda, Number *lambda) {
      if (n != m_spline.numTheta())
        return false;

      if (m != m_spline.numConstraints())
        return false;

      if (init_z)
        return false;

      if (init_lambda)
        return false;

      if (init_x) {
        std::copy(m_theta_guess.begin(), m_theta_guess.end(), x);
      }
      return true;
    }

    bool ClothoidSplineProblem::eval_f(Index n, const Number *x, bool new_x, Number &obj_value) {
      if (n != m_spline.numTheta())
        return false;
      return m_spline.objective(x, obj_value);
    }

    bool ClothoidSplineProblem::eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f) {
      if (n != m_spline.numTheta())
        return false;
      return m_spline.gradient(x, grad_f);
    }

    bool ClothoidSplineProblem::eval_g(Index n, const Number *x, bool new_x, Index m, Number *g) {
      if (n != m_spline.numTheta())
        return false;
      if (m != m_spline.numConstraints())
        return false;
      return m_spline.constraints(x, g);
    }

    bool ClothoidSplineProblem::eval_jac_g(Index n, const Number *x, bool new_x, Index m, Index nele_jac, Index *iRow,
                                           Index *jCol, Number *values) {
      if (n != m_spline.numTheta())
        return false;
      if (m != m_spline.numConstraints())
        return false;
      if (nele_jac != m_spline.jacobian_nnz())
        return false;

      bool index_ok = true;
      if ((iRow != NULL) && (jCol != NULL)) {
        index_ok = m_spline.jacobian_pattern(iRow, jCol);
      }

      bool jac_eval_ok = true;
      if (values != NULL) {
        jac_eval_ok = m_spline.jacobian(x, values);
      }

      return index_ok & jac_eval_ok;
    }

    void ClothoidSplineProblem::finalize_solution(SolverReturn status, Index n, const Number *x, const Number *z_L, 
                                                  const Number *z_U, Index m, const Number *g, const Number *lambda, 
                                                  Number obj_value, const Ipopt::IpoptData *ip_data,
                                                  Ipopt::IpoptCalculatedQuantities *ip_cq) {
      std::copy(x, x + n, m_theta_sol.begin());
    }

    std::vector<real_type> interpolate_clothoid_list(const ClothoidSplineG2 &spline) {
      SmartPtr< ClothoidSplineProblem > spline_problem = new ClothoidSplineProblem(spline);
      SmartPtr< IpoptApplication > app = IpoptApplicationFactory();

      app->Options()->SetStringValue("jac_d_constant", "no");
      app->Options()->SetStringValue("hessian_constant", "no");
      app->Options()->SetStringValue("mu_strategy", "adaptive");
      app->Options()->SetStringValue("derivative_test", "none");
      app->Options()->SetStringValue("hessian_approximation", "limited-memory");
      app->Options()->SetStringValue("limited_memory_update_type", "bfgs");

      app->Options()->SetIntegerValue("max_iter", 400);

      app->Options()->SetNumericValue("tol", 1e-10);
      app->Options()->SetNumericValue("derivative_test_tol", 1e-5);

      if (app->Initialize() != ApplicationReturnStatus::Solve_Succeeded)
        throw std::runtime_error("Cannot initialize solver");

      app->OptimizeTNLP(spline_problem);

      return spline_problem->get_solution();

    }

  };  // namespace ipopt

};  // namespace G2lib

#endif
