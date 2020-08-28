/**
 * PYTHON Wrapper for Clothoids
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#pragma once

#include <stdexcept>

#include <G2lib.hh>
#include <Clothoid.hh>
#include <ClothoidList.hh>

namespace G2lib {
namespace Interpolation {
  using G2lib::real_type;
  using G2lib::int_type;

  using G2lib::ClothoidList;
  using G2lib::ClothoidSplineG2;

  class Interpolator {
    const std::vector<real_type> m_xs;
    const std::vector<real_type> m_ys;
    ClothoidSplineG2 m_spline;

   public:
    Interpolator(const std::vector<real_type> & xs, const std::vector<real_type> & ys) : 
      m_xs(xs), m_ys(ys), m_spline(ClothoidSplineG2()) {}

#ifdef LMSOLVE_CLOTHOID_SPLINE
    bool buildP1(real_type theta_0, real_type theta_1, ClothoidList & result);
    bool buildP2(ClothoidList & result);
#else
    bool buildP1(real_type theta_0, real_type theta_1, ClothoidList & result) {
      throw std::runtime_error("Not supported. Recompile with libeigen3-dev library installed!");
      return false;
    }
    bool buildP2(ClothoidList & result) {
      throw std::runtime_error("Not supported. Recompile with libeigen3-dev library installed!");
      return false;
    }
#endif
    // bool buildP3(real_type theta_0, real_type kappa_0, ClothoidList & result);
#ifdef IPOPT_CLOTHOID_SPLINE
    bool buildP4(ClothoidList & result);
    bool buildP5(ClothoidList & result);
    bool buildP6(ClothoidList & result);
    bool buildP7(ClothoidList & result);
    bool buildP8(ClothoidList & result);
    bool buildP9(ClothoidList & result);
#else
    bool buildP4(ClothoidList & result) {
      throw std::runtime_error("Not supported. Recompile with lipipopt-dev library installed!");
      return false;
    }
    bool buildP5(ClothoidList & result) {
      throw std::runtime_error("Not supported. Recompile with lipipopt-dev library installed!");
      return false;
    }
    bool buildP6(ClothoidList & result) {
      throw std::runtime_error("Not supported. Recompile with lipipopt-dev library installed!");
      return false;
    }
    bool buildP7(ClothoidList & result) {
      throw std::runtime_error("Not supported. Recompile with lipipopt-dev library installed!");
      return false;
    }
    bool buildP8(ClothoidList & result) {
      throw std::runtime_error("Not supported. Recompile with lipipopt-dev library installed!");
      return false;
    }
    bool buildP9(ClothoidList & result) {
      throw std::runtime_error("Not supported. Recompile with lipipopt-dev library installed!");
      return false;
    }
#endif

    const std::vector<real_type> & xs() { return m_xs; }
    const std::vector<real_type> & ys() { return m_ys; }

   private:
    bool build_clothoid_spline();
    bool check_input();

#ifdef LMSOLVE_CLOTHOID_SPLINE
    bool build_using_lm_solver(ClothoidList & result);
#else
    bool build_using_lm_solver(ClothoidList & result) {
      throw std::runtime_error("Not supported. Recompile with libeigen3-dev library installed!");
      return false;
    }
#endif

#ifdef IPOPT_CLOTHOID_SPLINE
    bool build_using_ipopt_solver(ClothoidList & result);
#else
    bool build_using_ipopt_solver(ClothoidList & result) {
      throw std::runtime_error("Not supported. Recompile with lipipopt-dev library installed!");
      return false;
    }
#endif

    bool build_clothoid_list(const std::vector<real_type> & theta, ClothoidList & result);
  };

} /* namespace Inteprolation */
} /* namespace G2lib */