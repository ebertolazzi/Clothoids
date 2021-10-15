/**
 * PYTHON Wrapper for Clothoids
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#pragma once

#ifdef IPOPT_CLOTHOID_SPLINE

#include "python-ClothoidSpline-Solver.hh"

#include "python-G2libHeaders.hh"
#include <IpTNLP.hpp>
#include <vector>

namespace G2lib {
  namespace Interpolation {

    using G2lib::ClothoidSplineG2;
    using G2lib::int_type;
    using G2lib::real_type;

    using Ipopt::Index;
    using Ipopt::Number;
    using Ipopt::SolverReturn;
    using Ipopt::TNLP;

    class IpoptSolver : public Solver {
      class ClothoidSplineProblem : public TNLP {
        IpoptSolver & m_solver;
        SolverReturn  m_solver_return;

       public:
        ClothoidSplineProblem(IpoptSolver & solver);

        virtual bool get_nlp_info(
            Index &                theta_size,
            Index &                contraints_size,
            Index &                jacobian_pattern_size,
            Index &                hessian_pattern_size,
            TNLP::IndexStyleEnum & pattern_index_style);

        virtual bool get_bounds_info(
            Index    theta_size,
            Number * theta_min,
            Number * theta_max,
            Index    constraints_size,
            Number * constraints_min,
            Number * constraints_max);

        virtual bool get_starting_point(
            Index    theta_size,
            bool     init_theta,
            Number * theta_ics,
            bool     init_z,
            Number * z_L,
            Number * z_U,
            Index    constraints_size,
            bool     init_lambda,
            Number * lambda);

        virtual bool eval_f(Index theta_size, const Number * theta, bool new_theta, Number & obj_value);

        virtual bool eval_grad_f(Index theta_size, const Number * theta, bool new_theta, Number * grad_f);

        virtual bool eval_g(Index theta_size, const Number * theta, bool new_theta, Index constraints_size, Number * g);

        virtual bool eval_jac_g(
            Index          theta_size,
            const Number * theta,
            bool           new_theta,
            Index          constraints_size,
            Index          jacobian_pattern_size,
            Index *        jacobian_rows,
            Index *        jacobian_cols,
            Number *       jacobian_values);

        virtual void finalize_solution(
            SolverReturn                       status,
            Index                              n,
            const Number *                     x,
            const Number *                     z_L,
            const Number *                     z_U,
            Index                              m,
            const Number *                     g,
            const Number *                     lambda,
            Number                             obj_value,
            const Ipopt::IpoptData *           ip_data,
            Ipopt::IpoptCalculatedQuantities * ip_cq);

        SolverReturn solver_return() { return m_solver_return; }
      };

     public:
      IpoptSolver(const ClothoidSplineG2 & spline) : Solver(spline) {};
      virtual bool solve() override;
    };
  } /* namespace Interpolation */
} /* namespace G2lib */

#endif