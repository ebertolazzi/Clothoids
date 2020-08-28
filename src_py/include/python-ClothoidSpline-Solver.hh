/**
 * PYTHON Wrapper for Clothoids
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#pragma once

#include <vector>

#include <G2lib.hh>
#include <Clothoid.hh>
#include <ClothoidList.hh>


namespace G2lib {
namespace Interpolation {

using G2lib::real_type;
using G2lib::int_type;
using G2lib::ClothoidSplineG2;

class Solver {
 private:
  const ClothoidSplineG2 & m_spline;
  int_type m_theta_size;
  int_type m_constraints_size;
  int_type m_jacobian_pattern_size;
  int_type m_jacobian_size;
  int_type m_lagrangian_size;
  std::vector<real_type> m_theta_solution;
  std::vector<real_type> m_theta_min;
  std::vector<real_type> m_theta_max;

 public:
  Solver(const ClothoidSplineG2 & spline);
  
  virtual bool solve() = 0;

  int theta_size() const { return m_theta_size; }
  int constraints_size() const { return m_constraints_size; }
  int jacobian_pattern_size() const { return m_jacobian_pattern_size; }
  int jacobian_size() const { return m_jacobian_size; }
  int lagrangian_size() const { return m_lagrangian_size; }
  
  const ClothoidSplineG2 & spline() { return m_spline; }
  std::vector<real_type> & theta_solution() { return m_theta_solution; }
  const std::vector<real_type> & theta_min() const { return m_theta_min; }
  const std::vector<real_type> & theta_max() const { return m_theta_max; }
};

} /* namespace Interpolation */
} /* namespace G2lib */