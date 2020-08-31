/**
 * PYTHON Wrapper for Clothoids
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-ClothoidSpline-Solver.hh"
#include <iostream>

namespace G2lib {
  namespace Interpolation {

    Solver::Solver(const ClothoidSplineG2 & spline)
        : m_spline(spline), m_theta_size(spline.numTheta()), m_constraints_size(spline.numConstraints()),
          m_jacobian_pattern_size(spline.jacobian_nnz()),
          m_theta_solution(std::vector<real_type>(spline.numTheta(), 0.0)),
          m_theta_min(std::vector<real_type>(spline.numTheta(), 0.0)),
          m_theta_max(std::vector<real_type>(spline.numTheta(), 0.0)) {
      m_jacobian_size           = theta_size() * constraints_size();
      m_lagrangian_hessian_size = (theta_size() + constraints_size()) * (theta_size() + constraints_size());
    }

    void Solver::guess() { m_spline.guess(&m_theta_solution.front(), &m_theta_min.front(), &m_theta_max.front()); }
  } /* namespace Interpolation */
} /* namespace G2lib */