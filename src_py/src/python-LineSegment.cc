/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-LineSegment.hh"
#include <pybind11/stl.h>

namespace G2lib {
  namespace python {
    void wrap_LineSegment(py::module & m) {
      py::class_<G2lib::LineSegment, G2lib::BaseCurve>(m, "LineSegment")
        .def(py::init<G2lib::BaseCurve const &>())
        .def(py::init<G2lib::LineSegment const &>())
        .def(py::init<real_type, real_type, real_type, real_type>())
        .def("copy", &LineSegment::copy)
        .def("build", &LineSegment::build)
        .def("build_2P", py::overload_cast<real_type, real_type, real_type, real_type>(&LineSegment::build_2P))
        .def("build_2P", [](LineSegment * self, std::tuple<real_type, real_type> p0, std::tuple<real_type, real_type> p1) {
          real_type _p0[2] {std::get<0>(p0), std::get<1>(p0)};
          real_type _p1[2] {std::get<0>(p1), std::get<1>(p1)};
          self->build_2P(_p0, _p1);
        })
        .def("p1p2", [](LineSegment * self) {
          real_type _p0[2], _p1[2];
          self->build_2P(_p0, _p1);
          return std::make_tuple(std::make_tuple(_p0[0], _p0[1]), std::make_tuple(_p1[0], _p1[1]));
        })
        .def("intersect", [](LineSegment * self, LineSegment const & S) {
          real_type s1, s2;
          self->intersect(S, s1, s2);
          return std::make_tuple(s1, s2);
        })
        .def("intersect", [](LineSegment * self, LineSegment const & S, bool swap_s_vals) {
          IntersectList ilist;
          self->intersect(S, ilist, swap_s_vals);
          return ilist;
        })
        .def("intersect_ISO", [](LineSegment * self, real_type offs, LineSegment const & S, real_type Soffs) {
          real_type s1, s2;
          self->intersect_ISO(offs, S, Soffs, s1, s2);
          return std::make_tuple(s1, s2);
        })
        .def("intersect", [](LineSegment * self, real_type offs, LineSegment const & S, real_type Soffs, bool swap_s_vals) {
          IntersectList ilist;
          self->intersect_ISO(offs, S, Soffs, ilist, swap_s_vals);
          return ilist;
        })
        .def("collision", &LineSegment::collision)
        .def("collision_ISO", &LineSegment::collision_ISO)
        .def("paramNURBS", [](LineSegment * self) {
          int_type n_pnts, n_knots;
          self->paramNURBS(n_knots, n_pnts);
          return std::make_tuple(n_knots, n_pnts);
        })
        .def("toNURBS", [](LineSegment * self){
          using Point = real_type[3];
          int_type n_pnts, n_knots;
          self->paramNURBS(n_knots, n_pnts);
          real_type * knots = new real_type[n_knots];
          Point * Poly = new Point[n_pnts];
          self->toNURBS(knots, Poly);

          std::vector<std::tuple<real_type, real_type, real_type>> Poly_vec;
          for (int_type i = 0; i < n_pnts; i++) {
            Poly_vec.push_back(std::make_tuple(Poly[i][0], Poly[i][1], Poly[i][2]));
          }
          std::vector<real_type> knots_vec;
          for (int_type i = 0; i < n_knots; i++) {
            knots_vec.push_back(knots[i]);
          }
          auto ret = std::make_tuple(knots_vec, Poly_vec);
          delete[] Poly;
          delete[] knots;
          return ret;
        })
        .def("toBS", [](LineSegment * self){
          real_type knots[4], Poly[2][2];
          self->toBS(knots, Poly);
          return std::make_tuple(
            std::make_tuple(knots[0], knots[1], knots[2], knots[3]),
            std::make_tuple(
              std::make_tuple(Poly[0][0], Poly[0][1]),
              std::make_tuple(Poly[1][0], Poly[1][1])));
        });
    }

    void wrap_PolyLine(py::module & m) {
      py::class_<PolyLine, BaseCurve>(m, "PolyLine")
        .def(py::init<>())
        .def(py::init<BaseCurve const &>())
        .def(py::init<PolyLine const &>())
        .def(py::init<LineSegment const &>())
        .def(py::init<CircleArc const &, real_type>())
        .def(py::init<Biarc const &, real_type>())
        .def(py::init<ClothoidCurve const &, real_type>())
        .def(py::init<ClothoidList const &, real_type>())
        .def("getSegment", &PolyLine::getSegment)
        .def("numSegment", &PolyLine::numSegment)
        .def("numPoints", &PolyLine::numPoints)
        .def("polygon", [](PolyLine * self) {
          std::vector<std::tuple<real_type, real_type>> ret;
          int_type numSegment = self->numSegment();
          real_type * x = new real_type[numSegment];
          real_type * y = new real_type[numSegment];
          self->polygon(x, y);
          for (int_type i = 0; i < numSegment; i++) {
            ret.push_back(std::make_tuple(x[i], y[i]));
          }
          delete[] x;
          delete[] y;
          return ret;
        })
        .def("init", &PolyLine::init)
        .def("push_back", py::overload_cast<real_type, real_type>(&PolyLine::push_back))
        .def("push_back", py::overload_cast<LineSegment const &>(&PolyLine::push_back))
        .def("push_back", py::overload_cast<CircleArc const &, real_type>(&PolyLine::push_back))
        .def("push_back", py::overload_cast<Biarc const &, real_type>(&PolyLine::push_back))
        .def("push_back", py::overload_cast<ClothoidCurve const &, real_type>(&PolyLine::push_back))
        .def("push_back", py::overload_cast<ClothoidList const &, real_type>(&PolyLine::push_back))
        .def("build", [](PolyLine * self, std::vector<std::tuple<real_type, real_type>> points) {
          int_type nPoints = static_cast<int_type>(points.size());
          real_type * x = new real_type[nPoints];
          real_type * y = new real_type[nPoints];
          for (int_type i = 0; i < nPoints; i++) {
            x[i] = std::get<0>(points[i]);
            y[i] = std::get<1>(points[i]);
          }
          self->build(x, y, nPoints);
          delete[] x;
          delete[] y;
        })
        .def("build", py::overload_cast<LineSegment const &>(&PolyLine::build))
        .def("build", py::overload_cast<CircleArc const &, real_type>(&PolyLine::build))
        .def("build", py::overload_cast<Biarc const &, real_type>(&PolyLine::build))
        .def("build", py::overload_cast<ClothoidCurve const &, real_type>(&PolyLine::build))
        .def("build", py::overload_cast<ClothoidList const &, real_type>(&PolyLine::build))
        .def("intersect", [](PolyLine * self, PolyLine const & S) {
          std::vector<real_type> s1, s2;
          self->intersect(S, s1, s2);
          return std::make_tuple(s1, s2);
        })
        .def("intersect", [](PolyLine * self, PolyLine const & S, bool swap_s_vals) {
          IntersectList ilist;
          self->intersect(S, ilist, swap_s_vals);
          return ilist;
        })
        .def("build_AABBtree", [](PolyLine * self) {
         std::unique_ptr<AABBtree> tree = std::make_unique<AABBtree>();
         self->build_AABBtree(*tree);
         return tree;
        })
        .def("build_AABBtree", py::overload_cast<>(&PolyLine::build_AABBtree, py::const_));
    }
  }
}