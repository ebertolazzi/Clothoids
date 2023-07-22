/**
 * PYTHON Wrapper for Clothoids
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#pragma once

#ifdef LMSOLVE_CLOTHOID_SPLINE

#include "python-ClothoidSpline-Solver.hh"

#include <vector>
#include <Eigen/Eigen>
#include <unsupported/Eigen/LevenbergMarquardt>

namespace G2lib {
  namespace Interpolation {

    using G2lib::ClothoidSplineG2;
    using G2lib::integer;
    using G2lib::real_type;
    using SparseFunctor = Eigen::SparseFunctor<real_type, integer>;

    /* Levemberg Marquardt Solver */
    class LMSolver : public Solver {
      struct ClothoidSplineProblem : SparseFunctor {
        LMSolver &             m_solver;
        std::vector<integer>   m_jacobian_rows;
        std::vector<integer>   m_jacobian_cols;
        std::vector<real_type> m_jacobian_result;

        ClothoidSplineProblem(LMSolver & solver);
        int operator()(const SparseFunctor::InputType & theta, SparseFunctor::ValueType & constraints_value) const;
        int df(const SparseFunctor::InputType & theta, SparseFunctor::JacobianType & jacobian_value);
      };

     public:
      LMSolver(const ClothoidSplineG2 & spline) : Solver(spline) {};
      virtual bool solve() override;
    };
  } /* namespace Interpolation */
} /* namespace G2lib */

#endif
