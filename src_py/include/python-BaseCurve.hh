/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#pragma once

#include <tuple>

#include <pybind11/pybind11.h>
#include "python-G2libHeaders.hh"

namespace py = pybind11;

namespace G2lib {
  namespace python {
    class PythonicBaseCurve : public G2lib::BaseCurve {
      public:

      PythonicBaseCurve(G2lib::CurveType const & __type) : G2lib::BaseCurve::BaseCurve(__type) {} 

      real_type 
      length() const override {
        PYBIND11_OVERLOAD_PURE(real_type, G2lib::BaseCurve, length);
      }

      real_type 
      length_ISO(real_type offs) const override {
        PYBIND11_OVERLOAD_PURE(real_type, G2lib::BaseCurve, length, offs);
      }

      real_type 
      theta(real_type s) const override {
        PYBIND11_OVERLOAD_PURE(real_type, G2lib::BaseCurve, theta, s);
      }

      real_type 
      theta_D(real_type s) const override {
        PYBIND11_OVERLOAD_PURE(real_type, G2lib::BaseCurve, theta_D, s);
      }

      real_type 
      theta_DD(real_type s) const override {
        PYBIND11_OVERLOAD_PURE(real_type, G2lib::BaseCurve, theta_DD, s);
      }

      real_type 
      theta_DDD(real_type s) const override {
        PYBIND11_OVERLOAD_PURE(real_type, G2lib::BaseCurve, theta_DDD, s);
      }

      real_type 
      X(real_type s) const override {
        PYBIND11_OVERLOAD_PURE(real_type, G2lib::BaseCurve, X, s);
      }

      real_type 
      Y(real_type s) const override {
        PYBIND11_OVERLOAD_PURE(real_type, G2lib::BaseCurve, Y, s);
      }

      void 
      bbox(real_type & xmin, real_type & ymin, real_type & xmax, real_type & ymax) const override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, bbox, xmin, ymin, xmax, ymax);
      }

      void 
      bbox_ISO(real_type offs, real_type & xmin, real_type & ymin, real_type & xmax, real_type & ymax) const override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, bbox_ISO, offs, xmin, ymin, xmax, ymax);
      }

      void 
      eval(real_type s, real_type & x, real_type & y) const override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, eval, s, x, y);
      }

      real_type 
      X_D(real_type s) const override {
        PYBIND11_OVERLOAD_PURE(real_type, G2lib::BaseCurve, X_D, s);
      }

      real_type 
      X_DD(real_type s) const override {
        PYBIND11_OVERLOAD_PURE(real_type, G2lib::BaseCurve, X_DD, s);
      }

      real_type 
      X_DDD(real_type s) const override {
        PYBIND11_OVERLOAD_PURE(real_type, G2lib::BaseCurve, X_DDD, s);
      }

      real_type 
      Y_D(real_type s) const override {
        PYBIND11_OVERLOAD_PURE(real_type, G2lib::BaseCurve, Y_D, s);
      }

      real_type 
      Y_DD(real_type s) const override {
        PYBIND11_OVERLOAD_PURE(real_type, G2lib::BaseCurve, Y_DD, s);
      }

      real_type 
      Y_DDD(real_type s) const override {
        PYBIND11_OVERLOAD_PURE(real_type, G2lib::BaseCurve, Y_DDD, s);
      }

      void 
      eval_D(real_type s, real_type & x_D, real_type & y_D) const override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, eval, s, x_D, y_D);
      }

      void 
      eval_DD(real_type s, real_type & x_DD, real_type & y_DD) const override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, eval, s, x_DD, y_DD);
      }

      void 
      eval_DDD(real_type s, real_type & x_DDD, real_type & y_DDD) const override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, eval, s, x_DDD, y_DDD);
      }

      void 
      translate(real_type tx, real_type ty) override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, translate, tx, ty);
      }

      void 
      rotate(real_type angle, real_type cx, real_type cy) override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, rotate, angle, cx, cy);
      }

      void 
      scale(real_type sc) override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, scale, sc);
      }

      void 
      reverse() override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, reverse);
      }

      void 
      changeOrigin(real_type newx0, real_type newy0) override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, changeOrigin, newx0, newy0);
      }

      void 
      trim(real_type s_begin, real_type s_end) override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, trim, s_begin, s_end);
      }

      int_type
      closestPoint_ISO(real_type qx, real_type qy, real_type & x, real_type & y, real_type & s, real_type & t, real_type & dst) const override {
        PYBIND11_OVERLOAD_PURE(int_type, G2lib::BaseCurve, closestPoint_ISO, qx, qy, x, y, s, t, dst);
      }

      int_type
      closestPoint_ISO(real_type qx, real_type qy, real_type offs, real_type & x, real_type & y, real_type & s, real_type & t, real_type & dst) const override {
        PYBIND11_OVERLOAD_PURE(int_type, G2lib::BaseCurve, closestPoint_ISO, qx, qy, offs, x, y, s, t, dst);
      }

      void
      info(ostream_type & stream) const override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, info, stream);
      }

      void
      bbTriangles(vector<Triangle2D> & tvec, real_type max_angle, real_type max_size, int_type icurve) const override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, bbTriangles, tvec, max_angle, max_size, icurve);
      }
      
      void
      bbTriangles_ISO(real_type offs, vector<Triangle2D> & tvec, real_type max_angle, real_type max_size, int_type icurve) const override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, bbTriangles_ISO, offs, tvec, max_angle, max_size, icurve);
      }

      void
      bbTriangles_SAE(real_type offs, vector<Triangle2D> & tvec, real_type max_angle, real_type max_size, int_type icurve) const override {
        PYBIND11_OVERLOAD_PURE(void, G2lib::BaseCurve, bbTriangles_SAE, offs, tvec, max_angle, max_size, icurve);
      }
    };

    void wrap_BaseCurve(py::module & m);
  }
}