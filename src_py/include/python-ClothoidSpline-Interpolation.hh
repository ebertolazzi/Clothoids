/**
 * PYTHON Wrapper for Clothoids
 *
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#pragma once

#include <stdexcept>

#include "python-G2libHeaders.hh"

namespace G2lib {
  namespace Interpolation {
    using G2lib::int_type;
    using G2lib::real_type;

    using G2lib::ClothoidList;
    using G2lib::ClothoidSplineG2;

    class Interpolator {
      const std::vector<real_type> m_xs;
      const std::vector<real_type> m_ys;
      ClothoidSplineG2             m_spline;

     public:
      Interpolator(const std::vector<real_type> & xs, const std::vector<real_type> & ys)
          : m_xs(xs), m_ys(ys), m_spline() {}

#ifdef LMSOLVE_CLOTHOID_SPLINE
      void buildP1(real_type theta_0, real_type theta_1, ClothoidList & result);
      void buildP2(ClothoidList & result);
#else
      void buildP1(real_type theta_0, real_type theta_1, ClothoidList & result) {
        throw std::runtime_error("Not supported. Recompile with libeigen3-dev library installed!");
      }
      void buildP2(ClothoidList & result) {
        throw std::runtime_error("Not supported. Recompile with libeigen3-dev library installed!");
      }
#endif
      // void buildP3(real_type theta_0, real_type kappa_0, ClothoidList & result);
#ifdef IPOPT_CLOTHOID_SPLINE
      void buildP4(ClothoidList & result);
      void buildP5(ClothoidList & result);
      void buildP6(ClothoidList & result);
      void buildP7(ClothoidList & result);
      void buildP8(ClothoidList & result);
      void buildP9(ClothoidList & result);
#else
      void buildP4(ClothoidList & result) {
        throw std::runtime_error("Not supported. Recompile with lipipopt-dev library installed!");
      }
      void buildP5(ClothoidList & result) {
        throw std::runtime_error("Not supported. Recompile with lipipopt-dev library installed!");
      }
      void buildP6(ClothoidList & result) {
        throw std::runtime_error("Not supported. Recompile with lipipopt-dev library installed!");
      }
      void buildP7(ClothoidList & result) {
        throw std::runtime_error("Not supported. Recompile with lipipopt-dev library installed!");
      }
      void buildP8(ClothoidList & result) {
        throw std::runtime_error("Not supported. Recompile with lipipopt-dev library installed!");
      }
      void buildP9(ClothoidList & result) {
        throw std::runtime_error("Not supported. Recompile with lipipopt-dev library installed!");
      }
#endif

      const std::vector<real_type> & xs() { return m_xs; }
      const std::vector<real_type> & ys() { return m_ys; }

     private:
      void build_clothoid_spline();
      void check_input();

#ifdef LMSOLVE_CLOTHOID_SPLINE
      void build_using_lm_solver(ClothoidList & result);
#else
      void build_using_lm_solver(ClothoidList & result) {
        throw std::runtime_error("Not supported. Recompile with libeigen3-dev library installed!");
      }
#endif

#ifdef IPOPT_CLOTHOID_SPLINE
      void build_using_ipopt_solver(ClothoidList & result);
#else
      void build_using_ipopt_solver(ClothoidList & result) {
        throw std::runtime_error("Not supported. Recompile with lipipopt-dev library installed!");
      }
#endif

      void build_clothoid_list(const std::vector<real_type> & theta, ClothoidList & result);
    };

  }  // namespace Interpolation
} /* namespace G2lib */