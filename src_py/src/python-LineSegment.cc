/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-LineSegment.hh"
#include "pybind11/stl.h"

namespace G2lib {
  namespace python {
    void wrap_LineSegment(py::module & m) {
      py::class_<G2lib::LineSegment, G2lib::BaseCurve>(m, "LineSegment",
      R"S(
        Class that manages a line segment. There are several possible
        constructors for this class:

         * constructor from a Base Curve
         * constructor from a Line Segment
         * constructor from raw data
          
        for this last constructor:

        :param float x0: starting position **x** coordinate
        :param float y0: starting position **y** coordinate
        :param float theta0: initial angle
        :param float l: length
      )S")
        .def(py::init<>())
        .def(py::init<G2lib::BaseCurve const &>(), py::arg("c"))
        .def(py::init<G2lib::LineSegment const &>(), py::arg("s"))
        .def(py::init<real_type, real_type, real_type, real_type>(),
          py::arg("x0"), py::arg("y0"), py::arg("theta0"), py::arg("l"))

        .def("copy", [](const LineSegment & self) {
          LineSegment other;
          other.copy(self);
          return other;
        },
        R"S(
          Create a copy of the current line segment

          :return: a copy of the current line segment
          :rtype: LineSegment
        )S")

        .def("build", &LineSegment::build,
          py::arg("x0"), py::arg("y0"), py::arg("theta0"), py::arg("l"),
        R"S(
          Builds a line segment with the standard parameters

          :param float x0: starting position **x** coordinate
          :param float y0: starting position **y** coordinate
          :param float theta0: initial angle
          :param float l: length
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("build_2P", py::overload_cast<real_type, real_type, real_type, real_type>(&LineSegment::build_2P),
          py::arg("x0"), py::arg("y0"), py::arg("x1"), py::arg("y1"),
        R"S(
          Builds a line segment with the passed points coordinates

          :param float x0: starting position **x** coordinate
          :param float y0: starting position **y** coordinate
          :param float x1: ending position **x** coordinate
          :param float y1: ending position **y** coordinate
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("build_2P", [](LineSegment & self, const std::tuple<real_type, real_type> p0, const std::tuple<real_type, real_type> p1) {
          real_type _p0[2] {std::get<0>(p0), std::get<1>(p0)};
          real_type _p1[2] {std::get<0>(p1), std::get<1>(p1)};
          self.build_2P(_p0, _p1);
        },
        R"S(
          Builds a line segment with the passed points

          :param Tuple[float, float] p0: tuple containing x and y coordinate of the starting point
          :param Tuple[float, float] p1: tuple containing x and y coordinate of the ending point
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("p1p2", [](LineSegment & self) {
          real_type _p0[2], _p1[2];
          self.build_2P(_p0, _p1);
          return std::make_tuple(std::make_tuple(_p0[0], _p0[1]), std::make_tuple(_p1[0], _p1[1]));
        },
        R"S(
          Returns the first and last points of the line segment

          :return: first and last points of the line segment
          :rtype: Tuple[Tuple[float, float], Tuple[float, float]]
        )S")

        .def("paramNURBS", [](LineSegment & self) {
          int_type n_pnts, n_knots;
          self.paramNURBS(n_knots, n_pnts);
          return std::make_tuple(n_knots, n_pnts);
        },
        R"S(
          Return the number of knots and points for the nurbs of the line segment
          
          :return: knots count and point count
          :rtype: Tuple[int, int]
        )S")

        .def("toNURBS", [](const LineSegment & self){
          using Point = real_type[3];
          using TPoint = std::tuple<float, float, float>;

          int_type n_pnts, n_knots;
          self.paramNURBS(n_knots, n_pnts);
          
          std::vector<real_type> knots(n_knots);
          std::vector<Point> poly(n_pnts);
          std::vector<TPoint> tpoly;
          tpoly.reserve(n_pnts);

          self.toNURBS(knots.data(), poly.data());
          std::for_each(poly.cbegin(), poly.cend(), [&](const Point & p) {
            tpoly.push_back(std::make_tuple(p[0], p[1], p[2]));
          });

          return std::make_tuple(knots, tpoly);
        },
        R"S(
          Returns the nurbs parameters of the line segment, as a tuple with
          knots and point list (as a tuple of 3 values)

          :return: nurbs parameters
          :rtype: Tuple[List[float], List[Tuple[float, float, float]]]
        )S")

        .def("toBS", [](const LineSegment & self){
          real_type knots[4], Poly[2][2];
          self.toBS(knots, Poly);
          return std::make_tuple(
            std::make_tuple(knots[0], knots[1], knots[2], knots[3]),
            std::make_tuple(
              std::make_tuple(Poly[0][0], Poly[0][1]),
              std::make_tuple(Poly[1][0], Poly[1][1])));
        },
        R"S(
          :return: BS parameters
          :rtype: Tuple[Tuple[float, float, float, float], Tuple[Tuple[float, float], Tuple[float, float]]]
        )S");
    }

    void wrap_PolyLine(py::module & m) {
      py::class_<PolyLine, BaseCurve>(m, "PolyLine",
      R"S(
        Class that manages a poly line. There are several possible
        constructors for this class

         * constructor from a Base Curve
         * constructor from a Poly Line
         * constructor from a Line Segment
         * constructor from a Circle Arc
         * constructor from a Biarc
         * constructor from a Clothoid Curve
         * constructor from a Clothoid List   
      )S")
      
        .def(py::init<>())
        .def(py::init<BaseCurve const &>(), py::arg("c"))
        .def(py::init<PolyLine const &>(), py::arg("pl"))
        .def(py::init<LineSegment const &>(), py::arg("ls"))
        .def(py::init<CircleArc const &, real_type>(),
          py::arg("c"), py::arg("tol"))
        .def(py::init<Biarc const &, real_type>(),
          py::arg("b"), py::arg("tol"))
        .def(py::init<ClothoidCurve const &, real_type>(),
          py::arg("b"), py::arg("tol"))
        .def(py::init<ClothoidList const &, real_type>(),
          py::arg("b"), py::arg("tol"))
      
        .def("getSegment", &PolyLine::getSegment,
          py::arg("n"),
        R"S(
          Returns the n-th segment of the poly line 

          :param int n: number of the segment to get
          :return: the requested segment
          :rtype: LineSegment
        )S")

        .def("__getitem__", &PolyLine::getSegment,
          py::arg("n"),
        R"S(
          Returns the n-th segment of the poly line 

          :param int n: number of the segment to get
          :return: the requested segment
          :rtype: LineSegment
        )S")

        .def("numSegments", &PolyLine::numSegments,
        R"S(
          Returns the amount of segments of the poly line

          :return: amount of segments of the poly line
          :rtype: int
        )S")

        .def("__len__", &PolyLine::numSegments,
        R"S(
          Returns the amount of segments of the poly line

          :return: amount of segments of the poly line
          :rtype: int
        )S")

        .def("numPoints", &PolyLine::numPoints,
        R"S(
          Returns the amount of points of the poly line

          :return: amount of points of the poly line
          :rtype: int
        )S")

        .def("polygon", [](const PolyLine & self) {
          int_type n = self.numSegments();
          std::vector<std::tuple<real_type, real_type>> ret;

          std::vector<real_type> x(n);
          std::vector<real_type> y(n);
          ret.reserve(n);
          self.polygon(x.data(), y.data());
          for (int_type i = 0; i < n; i++) {
            ret.push_back(std::make_tuple(x[i], y[i]));
          }
          return ret;
        },
        R"S(
          Returns a list of points representing the polygon described by
          this poly line. Each point is a tuple containg the **x** and
          **y** coordinates of the point. 

          :return: list of points of the polygon
          :rtype: List[Tuple[float, float]]
        )S")

        .def("init", py::overload_cast<>(&PolyLine::init),
        R"S(
          Empties the current list
  
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("push_back", py::overload_cast<real_type, real_type>(&PolyLine::push_back),
          py::arg("x"), py::arg("y"),
        R"S(
          Pushes back another segment in the polyline via ending point

          :param float x: **x** coordinate of the point
          :param float y: **y** coordinate of the point
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("push_back", py::overload_cast<LineSegment const &>(&PolyLine::push_back),
          py::arg("c"),
        R"S(
          Pushes back another segment in the polyline via ending point

          :param LineSegment c: line segment to push back
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("push_back", py::overload_cast<CircleArc const &, real_type>(&PolyLine::push_back),
          py::arg("c"), py::arg("tol"),
        R"S(
          Pushes back another segment in the polyline via ending point

          :param CircleArc c: circle arc to push back
          :param float tol: tolerance
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("push_back", py::overload_cast<Biarc const &, real_type>(&PolyLine::push_back),
          py::arg("c"), py::arg("tol"),
        R"S(
          Pushes back another segment in the polyline via ending point

          :param Biarc c: biarc to push back
          :param float tol: tolerance
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("push_back", py::overload_cast<ClothoidCurve const &, real_type>(&PolyLine::push_back),
          py::arg("c"), py::arg("tol"),
        R"S(
          Pushes back another segment in the polyline via ending point

          :param ClothoidCurve c: clothoid curve to push back
          :param float tol: tolerance
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("push_back", py::overload_cast<ClothoidList const &, real_type>(&PolyLine::push_back),
          py::arg("l"), py::arg("tol"),
        R"S(
          Pushes back another segment in the polyline via ending point

          :param ClothoidList l: clothoid list to push back
          :param float tol: tolerance
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("build", [](PolyLine & self, const std::vector<std::tuple<real_type, real_type>> points) {
          const size_t n = points.size();
          std::vector<real_type> x_, y_;
          x_.reserve(n);
          y_.reserve(n);
          std::for_each(points.cbegin(), points.cend(), [&](const auto & point) {
            x_.push_back(std::get<0>(point));
            y_.push_back(std::get<1>(point));
          });
          return self.build(x_.data(), y_.data(), n);
        },
        R"S(
          Builds a poly line from the list of passed points

          :param List[Tuple[float, float]] points: points from which to build the poly line
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("build", py::overload_cast<LineSegment const &>(&PolyLine::build),
          py::arg("l"),
        R"S(
          Builds a poly line from the passed line segment

          :param LineSegment l: line segment to use to build the poly line
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("build", py::overload_cast<CircleArc const &, real_type>(&PolyLine::build),
          py::arg("c"), py::arg("tol"),
        R"S(
          Builds a poly line from the passed circle arc

          :param CircleArc c: circle arc to use to build the poly line
          :param float tol: tolerance
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("build", py::overload_cast<Biarc const &, real_type>(&PolyLine::build),
          py::arg("b"), py::arg("tol"),
        R"S(
          Builds a poly line from the passed biarc

          :param Biarc b: biarc to use to build the poly line
          :param float tol: tolerance
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("build", py::overload_cast<ClothoidCurve const &, real_type>(&PolyLine::build),
          py::arg("c"), py::arg("tol"),
        R"S(
          Builds a poly line from the passed clothoid curve

          :param ClothoidCurve c: clothiod curve to use to build the poly line
          :param float tol: tolerance
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("build", py::overload_cast<ClothoidList const &, real_type>(&PolyLine::build),
          py::arg("cl"), py::arg("tol"),
        R"S(
          Builds a poly line from the passed clothoid list

          :param ClothoidList cl: clothoid list to use to build the poly line
          :param float tol: tolerance
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("build_AABBtree", [](PolyLine & self) {
         std::unique_ptr<AABBtree> tree = std::make_unique<AABBtree>();
         self.build_AABBtree(*tree);
         return tree;
        },
        R"S(
          Builds and returns an AABB tree on the current poly line

          :return: the built AABB Tree
          :rtype: AABBtree
        )S")

        .def("build_AABBtree", py::overload_cast<>(&PolyLine::build_AABBtree, py::const_),
        R"S(
          Builds an AABB tree on the current poly line

          :return: nothing, works in place
          :rtype: NoneType
        )S"
        );
    }
  }
}