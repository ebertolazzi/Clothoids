/**
 * PYTHON Wrapper for Clothoids
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-ClothoidSpline-Interpolation.hh"
#include "python-ClothoidSpline-IpoptSolver.hh"
#include "python-ClothoidSpline-LMSolver.hh"

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace G2lib {
namespace Interpolation {

bool Interpolator::build_clothoid_spline() {
  if (!check_input()) {
    throw std::runtime_error("check_input_failed while building clothoid spline");
    return false;
  }
  m_spline.build(xs().data(), ys().data(), xs().size());
  return true;
}

bool Interpolator::check_input() {
  /* Check equal size for input */
  if (xs().size() != ys().size()) {
    throw std::runtime_error("Input vectors must be of same length");
    return false;
  }

  /* Check enough input */
  if (xs().size() < 2) {
    throw std::runtime_error("Too small input size");
    return false;
  }
  
  const int_type size = xs().size();
  std::vector<real_type> chk;
  for (int_type i = 1; i < size; i++) {
    const real_type x_diff2 = std::pow(xs()[i] - xs()[i-1], 2);
    const real_type y_diff2 = std::pow(ys()[i] - ys()[i-1], 2);
    chk.push_back(x_diff2 + y_diff2);
  }
  const auto min_chk = std::min_element(chk.begin(), chk.end());
  const auto max_chk = std::max_element(chk.begin(), chk.end());
  if (*min_chk == 0) {
    throw std::runtime_error("Minimal distance too short");
    return false;
  }
  if (*min_chk < (1e-10 * *max_chk)) {
    throw std::runtime_error("Problem with too much deviation");
    return false;
  }
  return true;
}

bool Interpolator::build_clothoid_list(const std::vector<real_type> & theta, ClothoidList & result) {
  if (theta.size() < 2) {
    throw std::runtime_error("Result has only two values??");
    return false;
  }
  result.init();
  result.reserve(theta.size() - 1);
  for (int_type i = 0; i < static_cast<int_type>(theta.size()) - 1; i++)
    result.push_back_G1(xs()[i], ys()[i], theta[i], xs()[i + 1], ys()[i + 1], theta[i + 1]);
  return true;
}

#ifdef LMSOLVE_CLOTHOID_SPLINE
bool Interpolator::buildP1(real_type theta_0, real_type theta_1, ClothoidList & result) {
  m_spline.setP1(theta_0, theta_1);
  if (!build_clothoid_spline()) 
    return false;
  return build_using_lm_solver(result);
}

bool Interpolator::buildP2(ClothoidList & result) {
  m_spline.setP2();
  if (!build_clothoid_spline()) 
    return false;
  return build_using_lm_solver(result);
}

bool Interpolator::build_using_lm_solver(ClothoidList & result) {
  LMSolver solver(m_spline);
  if (!solver.solve()) {
    throw std::runtime_error("Optimization failed");
    return false;
  }
  build_clothoid_list(solver.theta_solution(), result);
  return true;
}
#endif

#ifdef IPOPT_CLOTHOID_SPLINE
bool Interpolator::buildP4(ClothoidList & result) {
  m_spline.setP4();
  if (!build_clothoid_spline()) 
    return false;
  return build_using_ipopt_solver(result);
}

bool Interpolator::buildP5(ClothoidList & result) {
  m_spline.setP5();
  if (!build_clothoid_spline()) 
    return false;
  return build_using_ipopt_solver(result);
}

bool Interpolator::buildP6(ClothoidList & result) {
  m_spline.setP6();
  if (!build_clothoid_spline()) 
    return false;
  return build_using_ipopt_solver(result);
}

bool Interpolator::buildP7(ClothoidList & result) {
  m_spline.setP7();
  if (!build_clothoid_spline()) 
    return false;
  return build_using_ipopt_solver(result);
}

bool Interpolator::buildP8(ClothoidList & result) {
  m_spline.setP8();
  if (!build_clothoid_spline()) 
    return false;
  return build_using_ipopt_solver(result);
}

bool Interpolator::buildP9(ClothoidList & result) {
  m_spline.setP9();
  if (!build_clothoid_spline()) 
    return false;
  return build_using_ipopt_solver(result);
}

bool Interpolator::build_using_ipopt_solver(ClothoidList & result) {
  IpoptSolver solver(m_spline);
  if (!solver.solve()) {
    throw std::runtime_error("Optimization failed");
    return false;
  }
  return build_clothoid_list(solver.theta_solution(), result);
}
#endif

} /* namespace Interpolation */
} /* namespace G2lib */
