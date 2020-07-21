/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-Triangle2D.hh"
#include <pybind11/stl.h>

namespace G2lib {
  namespace python {
    void wrap_Triangle2D(py::module & m) {
      py::class_<Triangle2D>(m, "Triangle2D")
        .def(py::init<>())
        .def(py::init<const Triangle2D &>())
        .def(py::init<real_type, real_type, real_type, real_type, real_type, 
                      real_type, real_type, real_type, int_type>())
        .def(py::init([](std::tuple<real_type, real_type> p0, std::tuple<real_type, real_type> p1, 
                         std::tuple<real_type, real_type> p2, real_type s0, real_type s1, int_type icurve) {
          return Triangle2D(std::get<0>(p0), std::get<1>(p0), std::get<0>(p1), std::get<1>(p1),
                            std::get<0>(p2), std::get<1>(p2), s0, s1, icurve);
        }))
        .def("build", [](Triangle2D * self, std::tuple<real_type, real_type> p0, std::tuple<real_type, real_type> p1, 
                         std::tuple<real_type, real_type> p2, real_type s0, real_type s1, int_type icurve) {
          self->build(std::get<0>(p0), std::get<1>(p0), std::get<0>(p1), std::get<1>(p1),
                      std::get<0>(p2), std::get<1>(p2), s0, s1, icurve);
        })
        .def("build", [](Triangle2D * self, real_type x0, real_type y0, real_type x1, real_type y1,
                         real_type x2, real_type y2, real_type s0, real_type s1, int_type icurve) {
          self->build(x0, y0, x1, y1, x2, y2, s0, s1, icurve);
        })
        .def("Icurve", &Triangle2D::Icurve)
        .def("x1", &Triangle2D::x1)
        .def("y1", &Triangle2D::y1)
        .def("x2", &Triangle2D::x2)
        .def("y2", &Triangle2D::y2)
        .def("x3", &Triangle2D::x3)
        .def("y3", &Triangle2D::y3)
        .def("S0", &Triangle2D::S0)
        .def("S1", &Triangle2D::S1)
        .def("translate", &Triangle2D::translate)
        .def("rotate", &Triangle2D::rotate)
        .def("scale", &Triangle2D::scale)
        .def("bbox", [](Triangle2D * self) {
          real_type x_min, y_min, x_max, y_max;
          self->bbox(x_min, y_min, x_max, y_max);
          return std::make_tuple(
            std::make_tuple(x_min, y_min), 
            std::make_tuple(x_max, y_max));
        })
        .def("baricenterX", &Triangle2D::baricenterX)
        .def("baricenterY", &Triangle2D::baricenterY)
        .def("P1", [](Triangle2D * self) {
          real_type const * p = self->P1();
          return std::make_tuple(p[0], p[1]);
        })
        .def("P2", [](Triangle2D * self) {
          real_type const * p = self->P2();
          return std::make_tuple(p[0], p[1]);
        })
        .def("P3", [](Triangle2D * self) {
          real_type const * p = self->P3();
          return std::make_tuple(p[0], p[1]);
        })
        .def("overlap", &Triangle2D::overlap)
        .def("isCounterClockwise", &Triangle2D::isCounterClockwise)
        .def("isInside", py::overload_cast<real_type, real_type>(&Triangle2D::isInside, py::const_))
        .def("isInside", [](Triangle2D * self, std::tuple<real_type, real_type> p) {
          real_type _p[2] = {std::get<0>(p), std::get<1>(p)};
          return self->isInside(_p);
        })
        .def("distMin", &Triangle2D::distMin)
        .def("distMax", &Triangle2D::distMax)
        .def("__str__", [](Triangle2D * self) {
          std::ostringstream str;
          self->info(str);
          return str.str();
        });
    }
  }
}