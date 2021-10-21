/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "pybind11/pybind11.h"
#include "python-BaseCurve.hh"
#include "python-LineSegment.hh"
#include "python-CircleArc.hh"
#include "python-Triangle2D.hh"
#include "python-AABBtree.hh"
#include "python-Biarc.hh"
#include "python-Clothoid.hh"


namespace py = pybind11;
using namespace G2lib;
using namespace G2lib::python;

PYBIND11_MODULE(G2lib, m) {
  wrap_BaseCurve(m);
  wrap_AABBtree(m);
  wrap_Triangle2D(m);
  wrap_LineSegment(m);
  wrap_CircleArc(m);
  wrap_Biarc(m);
  wrap_ClothoidCurve(m);
  wrap_PolyLine(m);
  wrap_BiarcList(m);
  wrap_ClothoidList(m);
  wrap_ClothoidSplineG2(m);

  m.attr("__version__") = py::str("2.0.9");

  py::class_<Solve2x2>(m, "Solve2x2")
    .def(py::init())
    .def("factorize", [](Solve2x2 * self, std::tuple<std::tuple<real_type, real_type>, std::tuple<real_type, real_type>> A) {
      real_type A_[2][2];
      A_[0][0] = std::get<0>(std::get<0>(A));
      A_[1][0] = std::get<0>(std::get<1>(A));
      A_[0][1] = std::get<1>(std::get<0>(A));
      A_[1][1] = std::get<1>(std::get<1>(A));
      return self->factorize(A_);
    })
    .def("solve", [](Solve2x2 * self, std::tuple<real_type, real_type> b) {
      real_type b_[2], x_[2];
      b_[0] = std::get<0>(b);
      b_[1] = std::get<1>(b);
      self->solve(b_, x_);
      return std::make_tuple(x_[0], x_[1]);
    });
}