/**
 * PYTHON Wrapper for Clothoids
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-ClothoidSpline-Solver.hh"

namespace G2lib {
namespace Interpolation {

  Solver::Solver(const ClothoidSplineG2 & spline) :
    m_spline(spline), m_theta_size(spline.numTheta()),
    m_constraints_size(spline.numConstraints()),
    m_jacobian_pattern_size(spline.jacobian_nnz()) {
      m_jacobian_size = theta_size() * constraints_size();
      // We don't have a pattern right now for lagrangian
      m_lagrangian_size = jacobian_size() * jacobian_size();
      m_theta_solution.reserve(theta_size());
      m_theta_min.reserve(theta_size());
      m_theta_max.reserve(theta_size());
      spline.guess(&m_theta_solution.front(), &m_theta_min.front(), &m_theta_max.front());
  }
} /* namespace Interpolation */
} /* namespace G2lib */