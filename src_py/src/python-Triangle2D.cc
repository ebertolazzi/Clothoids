/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-Triangle2D.hh"
#include "pybind11/stl.h"

namespace G2lib {
  namespace python {
    void wrap_Triangle2D(py::module & m) {
      py::class_<Triangle2D>(m, "Triangle2D",
      R"S(
        Class that manages a 2D Triangle. There are several possible
        constructors for this class:

         * constructor from a Triangle2D
         * constructor from three points, s0, s1 and icurve
         * constructor from raw data

        for this last constructor:

        :param float x1: **x** coordinate of the first point
        :param float y1: **y** coordinate of the first point
        :param float x2: **x** coordinate of the second point
        :param float y2: **y** coordinate of the second point
        :param float x3: **x** coordinate of the third point
        :param float y3: **y** coordinate of the third point
        :param float s0: s0 of this triangle on the curve
        :param float s1: s1 of this triangle on the curve
        :param int icurve: number of the trinagle on the curve
      )S")
        .def(py::init<>())
        .def(py::init<const Triangle2D &>())
        .def(py::init<real_type, real_type, real_type, real_type, real_type, 
                      real_type, real_type, real_type, int_type>())
        .def(py::init([](std::tuple<real_type, real_type> p0, std::tuple<real_type, real_type> p1, 
                         std::tuple<real_type, real_type> p2, real_type s0, real_type s1, int_type icurve) {
          return Triangle2D(std::get<0>(p0), std::get<1>(p0), std::get<0>(p1), std::get<1>(p1),
                            std::get<0>(p2), std::get<1>(p2), s0, s1, icurve);
        }))

        .def("build", [](Triangle2D & self, std::tuple<real_type, real_type> p0, std::tuple<real_type, real_type> p1, 
                         std::tuple<real_type, real_type> p2, real_type s0, real_type s1, int_type icurve) {
          self.build(std::get<0>(p0), std::get<1>(p0), std::get<0>(p1), std::get<1>(p1),
                      std::get<0>(p2), std::get<1>(p2), s0, s1, icurve);
        },
        R"S(
          Builds a triangle based on the passed points

          :param Tuple[float, float] p0: first point of the triangle
          :param Tuple[float, float] p1: second point of the triangle
          :param Tuple[float, float] p2: third point of the triangle
          :param float s0: s0 of the triangle
          :param float s1: s1 of the triangle
          :param int icurve: number of the trinagle on the curve
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("build", [](Triangle2D & self, real_type x0, real_type y0, real_type x1, real_type y1,
                         real_type x2, real_type y2, real_type s0, real_type s1, int_type icurve) {
          self.build(x0, y0, x1, y1, x2, y2, s0, s1, icurve);
        },
        R"S(
          Builds a triangle based on the passed coordinates

          :param float x1: **x** coordinate of the first point
          :param float y1: **y** coordinate of the first point
          :param float x2: **x** coordinate of the second point
          :param float y2: **y** coordinate of the second point
          :param float x3: **x** coordinate of the third point
          :param float y3: **y** coordinate of the third point
          :param float s0: s0 of this triangle on the curve
          :param float s1: s1 of this triangle on the curve
          :param int icurve: number of the trinagle on the curve
        )S")

        .def("Icurve", &Triangle2D::Icurve,
        R"S(
          Returns the number of the curve of the triangle

          :return: number of the curve of the triangle
          :rtype: int
        )S")

        .def("x1", &Triangle2D::x1,
        R"S(
          Returns the **x** coordinate of the first point of the triangle

          :return: **x** coordinate of the first point
          :rtype: float
        )S")

        .def("y1", &Triangle2D::y1,
        R"S(
          Returns the **y** coordinate of the first point of the triangle

          :return: **y** coordinate of the first point
          :rtype: float
        )S")

        .def("x2", &Triangle2D::x2,
        R"S(
          Returns the **x** coordinate of the second point of the triangle

          :return: **x** coordinate of the second point
          :rtype: float
        )S")

        .def("y2", &Triangle2D::y2,
        R"S(
          Returns the **y** coordinate of the second point of the triangle

          :return: **y** coordinate of the second point
          :rtype: float
        )S")

        .def("x3", &Triangle2D::x3,
        R"S(
          Returns the **x** coordinate of the third point of the triangle

          :return: **x** coordinate of the third point
          :rtype: float
        )S")

        .def("y3", &Triangle2D::y3,
        R"S(
          Returns the **y** coordinate of the third point of the triangle

          :return: **y** coordinate of the third point
          :rtype: float
        )S")

        .def("S0", &Triangle2D::S0,
        R"S(
          Returns the s0 of the triangle

          :return: s0 of the triangle
          :rtype: float
        )S")

        .def("S1", &Triangle2D::S1,
        R"S(
          Returns the s1 of the triangle

          :return: s1 of the triangle
          :rtype: float
        )S")

        .def("translate", &Triangle2D::translate,
          py::arg("tx"), py::arg("ty"),
        R"S(
          Translates the triangle by the indicated amount

          :param float tx: translation amount on x
          :param float ty: translation amount on y
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("rotate", &Triangle2D::rotate,
          py::arg("angle"), py::arg("cx"), py::arg("cy"),
        R"S(
          Rotates the triangle by the indicated amount

          :param float angle: angle of which to rotate the triangle
          :param float cx: **x** coordinate of the pivot point
          :param float cy: **y** coordinate of the pivot point
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("scale", &Triangle2D::scale,
          py::arg("sc"),
        R"S(
          Scales the triangle by the indicated amount

          :param float sc: scale factor
          :return: nothing, works in place
          :rtype: NoneType
        )S")

        .def("bbox", [](const Triangle2D & self) {
          real_type x_min, y_min, x_max, y_max;
          self.bbox(x_min, y_min, x_max, y_max);
          return std::make_tuple(
            std::make_tuple(x_min, y_min), 
            std::make_tuple(x_max, y_max));
        },
        R"S(
          Returns the bounding box of the triangle. The returned tuple will
          contain two tuples: the first one will contain the x and y coordinates
          of the first point of the bounding box, while the second one will
          contain the x and y coordinates of the second point of the bounding
          box.

          :return: bounding box of the triangle
          :rtype: Tuple[Tuple[float, float], Tuple[float, float]]
        )S")

        .def("baricenterX", &Triangle2D::baricenterX,
        R"S(
          Returns the x coordinate of the baricenter of the triangle

          :return: **x** coordinate of the baricenter of the triangle
          :rtype: float
        )S")

        .def("baricenterY", &Triangle2D::baricenterY,
        R"S(
          Returns the y coordinate of the baricenter of the triangle

          :return: **y** coordinate of the baricenter of the triangle
          :rtype: float
        )S")

        .def("P1", [](const Triangle2D & self) {
          real_type const * p = self.P1();
          return std::make_tuple(p[0], p[1]);
        },
        R"S(
          Returns the first point of the triangle

          :return: first point of the triangle
          :rtype: Tuple[float, float]
        )S")

        .def("P2", [](const Triangle2D & self) {
          real_type const * p = self.P2();
          return std::make_tuple(p[0], p[1]);
        },
        R"S(
          Returns the second point of the triangle

          :return: second point of the triangle
          :rtype: Tuple[float, float]
        )S")

        .def("P3", [](const Triangle2D & self) {
          real_type const * p = self.P3();
          return std::make_tuple(p[0], p[1]);
        },
        R"S(
          Returns the third point of the triangle

          :return: third point of the triangle
          :rtype: Tuple[float, float]
        )S")

        .def("overlap", &Triangle2D::overlap,
          py::arg("t2"),
        R"S(
          Returns whether the passed triangle overlaps with this triangle

          :param Triangle2D t2: triangle to check
          :return: whether the passed triangle overlaps with this triangle
          :rtype: bool
        )S")

        .def("isCounterClockwise", &Triangle2D::isCounterClockwise,
        R"S(
          Returns whether the triangle has points defined in a counterclockwise
          order. Possible return values are:

           * +1: counterclockwise
           * -1: clockwise
           *  0: deenerate triangle

          :return: whether the triangle is counterclockwise
          :rtype: int
        )S")

        .def("isInside", py::overload_cast<real_type, real_type>(&Triangle2D::isInside, py::const_),
          py::arg("x"), py::arg("y"),
        R"S(
          Returns whether the passed point is inside the triangle.
          Possible return values are:
           * +1: inside
           * -1: outside
           *  0: on the border

          :param float x: **x** coordinate of the point to check
          :param float y: **y** coordinate of the point to check
          :return: whether the passed point is inside the triangle
          :rtype: int
        )S")

        .def("isInside", [](const Triangle2D & self, std::tuple<real_type, real_type> p) {
          real_type _p[2] = {std::get<0>(p), std::get<1>(p)};
          return self.isInside(_p);
        },
        R"S(
          Returns whether the passed point is inside the triangle.
          Possible return values are:
           * +1: inside
           * -1: outside
           *  0: on the border
          
          :param Tuple[float, float]: point to check
          :return: whether the passed point is inside the triangle
          :rtype: int
        )S")

        .def("distMin", &Triangle2D::distMin,
          py::arg("x"), py::arg("y"),
        R"S(
          Returns the minimum distance between the passed point and
          the triangle

          :param float x: **x** coordinate of the point to check
          :param float y: **y** coordinate of the point to check
          :return: minimum distance between the point and the triangle
          :rtype: float
        )S")

        .def("distMax", &Triangle2D::distMax,
          py::arg("x"), py::arg("y"),
        R"S(
          Returns the maximum distance between the passed point and
          the triangle

          :param float x: **x** coordinate of the point to check
          :param float y: **y** coordinate of the point to check
          :return: maximum distance between the point and the triangle
          :rtype: float
        )S")
        
        .def("__str__", [](const Triangle2D & self) {
          std::ostringstream str;
          self.info(str);
          return str.str();
        });
    }
  }
}