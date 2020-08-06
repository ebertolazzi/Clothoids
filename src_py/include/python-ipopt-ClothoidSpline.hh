/**
 * PYTHON Wrapper for Clothoids
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#pragma once

#ifdef IPOPT_CLOTHOID_SPLINE

#include <Clothoid.hh>
#include <ClothoidList.hh>
#include <G2lib.hh>
#include <IpTNLP.hpp>
#include <vector>

namespace G2lib {

  namespace ipopt {
    using ::G2lib::ClothoidList;
    using ::G2lib::ClothoidSplineG2;

    using ::Ipopt::Index;
    using ::Ipopt::IndexStyleEnum;
    using ::Ipopt::Number;
    using ::Ipopt::TLNP;
    using ::Ipopt::SolverReturn;

    class ClothoidSplineProblem : public TNLP {
      const ClothoidSplineG2& m_spline;
      std::vector< real_type > m_theta_sol;

     public:
      ClothoidSplineProblem(const ClothoidSplineG2& spline);

      virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag,
                                IndexStyleEnum& index_style) override;

      virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) override;

      virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m,
                                      bool init_lambda, Number* lambda) override;

      virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) override;

      virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) override;

      virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) override;

      virtual bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index* jCol,
                              Number* values) override;

      virtual void finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L,
                                     const Number* z_U, Index m, const Number* g, const Number* lambda,
                                     Number obj_value, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq) override;

      const std::vector<real_type> & get_solution() { return m_theta_sol; }
    };

    ClothoidList interpolate_clothoid_list(const ClothoidSplineG2& spline);

  };  // namespace ipopt

};  // namespace G2lib

#endif