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

    using ::Ipopt::ApplicationReturnStatus;
    using ::Ipopt::Index;
    using ::Ipopt::IndexStyleEnum;
    using ::Ipopt::IpoptApplication;
    using ::Ipopt::IpoptApplicationFactory;
    using ::Ipopt::Number;
    using ::Ipopt::SmartPtr;
    using ::Ipopt::TLNP;

    using ::G2lib::ClothoidList;
    using ::G2lib::ClothoidSplineG2;
    using ::G2lib::int_type;
    using ::G2lib::real_type;

    ClothoidSplineProblem::ClothoidSplineProblem(const ClothoidSplineG2 &spline) : TNLP(), m_spline(spline) {
      m_theta_sol = std::vector< rela_type >(m_spline.numTheta(), 0.0);
    }

    bool ClothoidSplineProblem::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag,
                                             IndexStyleEnum &index_style) {
      n = m_spline.numTheta();
      m = m_spline.numConstraints();
      nnz_jac_g = m_spline.jacobian_nnz();
      nnz_h_lag = (n + m) * (n + m);  // Worst case scenario;
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

      std::vector< Number > theta_guess(n);
      m_spline.guess(theta_guess.front(), x_l, x_u);
      std::fill_n(g_l, m, 0.0);
      std::fill_n(g_u, m, 0.0);
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
        std::vector< Number > theta_min(n);
        std::vector< Number > theta_max(n);
        m_spline.guess(x, theta_min.front(), theta_max.front());
      }
      return true
    }

    bool ClothoidSplineProblem::eval_f(Index n, const Number *x, bool new_x, Number &obj_value) {
      if (n != m_spline.numTheta())
        return false;
      return m_spline.objective(x, obj_value);
    }

    bool ClothoidSplineProblem::eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f) {
      if (n != m_spline.numTheta())
        return false;
      return m_spline.gradient(x, grad_f) :
    }

    bool ClothoidSplineProblem::eval_g(Index n, const Number *x, bool new_x, Index m, Number *g) {
      if (n != m_spline.numTheta())
        return false;
      if (m != m_spline.numConstraints())
        return false;
      return m_spline.contraints(x, g);
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

      return index_ok & m_spline.jacobian(x, values);
    }

    void finalize_solution(SolverReturn status, Index n, const Number *x, const Number *z_L, const Number *z_U, Index m,
                           const Number *g, const Number *lambda, Number obj_value, const Ipopt::IpoptData *ip_data,
                           Ipopt::IpoptCalculatedQuantities *ip_cq) {
      std::copy(x, x + n, m_theta_sol.begin());
    }

    std::vector<real_type> interpolate_clothoid_list(const ClothoidSplineG2 &spline) {
      SmartPtr< ClothoidSplineProblem > spline_problem = new ClothoidSplineProblem(spline);
      SmartPtr< IpoptApplication > app = IpoptApplicationFactory();

      app->Options()->setStringValue("jac_d_constant", "no");
      app->Options()->setStringValue("hessian_constant", "no");
      app->Options()->setStringValue("mu_strategy", "adaptive");
      app->Options()->setStringValue("derivative_test", "none");
      app->Options()->setStringValue("hessian_approximation", "limited-memory");
      app->Options()->setStringValue("limited_memory_update_type", "bfgs");

      app->Options()->setIntegerValue("max_iter", 400);

      app->Options()->setNumericValue("tol", 1e-10);
      app->Options()->setNumericValue("derivative_test_tol", 1e-5);

      if (app->Initialize() != ApplicationReturnStatus::Solve_Suceeded)
        throw runtime_exception("Cannot initialize solver");

      app->OptimizeTNLP(spline_problem);

      return spline_problem->get_solution();

    }

  };  // namespace ipopt

};  // namespace G2lib

#endif
