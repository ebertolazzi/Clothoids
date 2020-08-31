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
    using G2lib::int_type;
    using G2lib::real_type;
    typedef Eigen::SparseFunctor<real_type, int_type> SparseFunctor;

    /* Levemberg Marquardt Solver */
    class LMSolver : public Solver {
      struct ClothoidSplineProblem : SparseFunctor {
        LMSolver &             m_solver;
        std::vector<int_type>  m_jacobian_rows;
        std::vector<int_type>  m_jacobian_cols;
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