/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-Biarc.hh"
#include "pybind11/stl.h"

#ifdef _WIN32
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#endif

namespace G2lib {
  namespace python {
    void wrap_Biarc(py::module & m) {
      py::class_<Biarc, BaseCurve>(m, "Biarc",
      R"S(
        Compute biarc fitting by Hermite data

        There are several possible constructors for this class:

         * Empty constructor
         * Copy contructor
         * Constructor from data
         * Explicit constructor from another base curve (unused)

        For the copy constructor

        :param Biarc curve: curve for copy

        For the data constructor:

        :param float x0: **x** coordinate of the origin
        :param float y0: **y** coordinate of the origin
        :param float theta0: tangent in the origin
        :param float x1: **x** coordinate of the destination
        :param float y1: **y** coordinate of the destination
        :param float theta1: tangent in the destination
      )S")

      .def(py::init())
      .def(py::init<Biarc const &>(), py::arg("curve"))
      .def(py::init<real_type, real_type, real_type, real_type, real_type, real_type>(),
        py::arg("x0"), py::arg("y0"), py::arg("theta0"), 
        py::arg("x1"), py::arg("y1"), py::arg("theta1"))
      .def(py::init<BaseCurve const &>(), py::arg("curve"))
      
      .def("copy", [](const Biarc & self) -> Biarc {
        Biarc other;
        other.copy(self);
        return other;
      },
      R"S(
        Create a copy of the current biarc curve

        :return: a new copy of the current biarc
        :rtype: Biarc
      )S")

      .def("build", &Biarc::build, 
        py::arg("x0"), py::arg("y0"), py::arg("theta0"), 
        py::arg("x1"), py::arg("y1"), py::arg("theta1"),
      R"S(
         Construct a biarc passing from the points
         :math:`(x_0,y_0)` to the point  :math:`(x_1,y_1)`
         with initial angle :math:`\theta_0` and final angle :math:`\theta_1`
        
        :param float x0: **x** coordinate of the origin
        :param float y0: **y** coordinate of the origin
        :param float theta0: tangent in the origin
        :param float x1: **x** coordinate of the destination
        :param float y1: **y** coordinate of the destination
        :param float theta1: tangent in the destination
        :return: false if biarc cannot be computed
        :rtype: bool
      )S")

      .def("c0", &Biarc::C0, 
      R"S(
        Returns the first circle arc of the biarc

        :return: circle arc of the biarc
        :rtype: CircleArc
      )S")

      .def("c1", &Biarc::C1, 
      R"S(
        Returns the second circle arc of the biarc

        :return: circle arc of the biarc
        :rtype: CircleArc
      )S")

      .def("build_3P", &Biarc::build_3P,
        py::arg("x0"), py::arg("y0"), py::arg("x1"), py::arg("y1"), py::arg("x2"), py::arg("y2"),
      R"S(
        Construct a biarc by 3 point at "minimum energy"
        
         * Planar point set fairing and fitting by arc splines
         * Xunnian Yang and Guozhao Wang
         * Computer-Aided Design, vol 33, 2001
        
        :param float x0: **x** coordinate of initial point
        :param float y0: **y** coordinate of initial point
        :param float x1: **x** coordinate of central point
        :param float y1: **y** coordinate of central point
        :param float x2: **x** coordinate of final point
        :param float y2: **y** coordinate of final point
        :return: false if biarc cannot be computed
        :rtype: bool
      )S")
      
      .def("xMiddle", &Biarc::xMiddle, 
      R"S(
        Returns the **x** coordinate of the junction point of the biarc

        :return: **x** coordinates of the junction
        :rtype: float
      )S")
      
      .def("yMiddle", &Biarc::yMiddle,
      R"S(
        Returns the **y** coordinate of the junction point of the biarc

        :return: **y** coordinates of the junction
        :rtype: float
      )S")

      .def("thetaMiddle", &Biarc::thetaMiddle,
      R"S(
        Returns the **x** coordinate of the junction point of the biarc

        :return: **x** coordinates of the junction
        :rtype: float
      )S")
      
      .def("kappa0", &Biarc::kappa0, 
      R"S(
        Curvature of the first circle arc

        :return: curvature of the first circle arc
        :rtype: float
      )S")

      .def("kappa1", &Biarc::kappa1, 
      R"S(
        Curvature of the second circle arc

        :return: curvature of the second circle arc
        :rtype: float
      )S")

      .def("length0", &Biarc::length0, 
      R"S(
        Length of the first circle arc

        :return: length of the first circle arc
        :rtype: float
      )S")

      .def("length1", &Biarc::length1, 
      R"S(
        Length of the second circle arc

        :return: length of the second circle arc
        :rtype: float
      )S")
      
      .def("delta_theta", &Biarc::delta_theta,
      R"S(
        Change of angle of the biarc :math:`\theta_1 - \theta_0`

        :return: the angle difference
        :rtype: float
      )S");
    }

    void wrap_BiarcList(py::module & m) {
      py::class_<BiarcList, BaseCurve>(m, "BiarcList",
      R"S(
        A class to manage a list of biarc curve (not necessarily G1 or G2 connected)

        There exists several constructors for this class:
         
         * empty constructor
         * constructor receiving a curve as parameter (BaseCurve or BiarcList or PolyLine
           or LineSegment or CircleArc or Biarc)
      )S")

      .def(py::init<>())
      .def(py::init<BaseCurve const &>())
      .def(py::init<BiarcList const &>())
      .def(py::init<PolyLine const &>())
      .def(py::init<LineSegment const &>())
      .def(py::init<CircleArc const &>())
      .def(py::init<Biarc const &>())
      
      .def("get", &BiarcList::get, py::arg("idx"),
      R"S(
        Returns the `idx`-th element of the list

        :param int idx: index
        :return: the biarc at that index
        :rtype: Biarc
      )S")

      .def("__getitem__", &BiarcList::get, py::arg("idx"),
      R"S(
        Returns the `idx`-th element of the list
      
        :param int idx: index
        :return: the biarc at that index
        :rtype: Biarc
      )S")

      .def("copy", [](const BiarcList & self) {
        BiarcList other;
        other.copy(self);
        return other;
      }, 
      R"S(
        Returns a copy of the current biarc list

        :return: a copy of the current biarc list
        :rtype: BiarcList
      )S")
      
      .def("init", &BiarcList::init,
      R"S(
        Empties the current biarc list

        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("reserve", &BiarcList::reserve, py::arg("n"),
      R"S(
        Reserve memory for **n** biarc in the list

        :param int n: number of element to reserve
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back", py::overload_cast<LineSegment const &>(&BiarcList::push_back), py::arg("c"),
      R"S(
        Append a line segment to the current biarc list. The curve is transformed
        to a degenerate biarc curve

        :param LineSegment c: the line segment
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back", py::overload_cast<CircleArc const &>(&BiarcList::push_back), py::arg("c"),
      R"S(
        Append a circle arc to the current biarc list. The curve is transformed
        to a degenerate biarc curve

        :param CircleArc c: the circle arc
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back", py::overload_cast<Biarc const &>(&BiarcList::push_back), py::arg("c"),
      R"S(
        Append a biarc to the current biarc list.

        :param Biarc c: the biarc
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back", py::overload_cast<PolyLine const &>(&BiarcList::push_back), py::arg("c"),
      R"S(
        Append a polyline to the current biarc list. The curve is transformed
        to a degenerate biarc curve

        :param PolyLine c: the polyline
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back_G1", py::overload_cast<real_type, real_type, real_type>(&BiarcList::push_back_G1),
      py::arg("x1"), py::arg("y1"), py::arg("theta1"),
      R"S(
        Construct a biarc passing from the points
        :math:`(x_0,y_0)` to the point  :math:`(x_1,y_1)`
        with initial angle :math:`\theta_0` and final angle :math:`\theta_1`
        and append the biarc to the tail of biarc list.
        The initial point and angle is taken from the tail of the biarc list.
        
        :param float x1:     :math:`x_1`
        :param float y1:     :math:`y_1`
        :param float theta1: :math:`\theta_1`
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back_G1", py::overload_cast<real_type, real_type, real_type, real_type, real_type, real_type>(&BiarcList::push_back_G1),
      py::arg("x0"), py::arg("y0"), py::arg("theta0"), py::arg("x1"), py::arg("y1"), py::arg("theta1"),
      R"S(
        Construct a biarc passing from the points
        :math:`(x_0,y_0)` to the point  :math:`(x_1,y_1)`
        with initial angle :math:`\theta_0` and final angle :math:`\theta_1`
        and append the biarc to the tail of biarc list.
        
        :param float x0:     :math:`x_0`
        :param float y0:     :math:`y_0`
        :param float theta0: :math:`\theta_0`
        :param float x1:     :math:`x_1`
        :param float y1:     :math:`y_1`
        :param float theta1: :math:`\theta_1`
        :return: nothing, works in place
        :rtype: NoneType
      )S")
      
      .def("build_G1", [](BiarcList & self, std::vector<real_type> x, std::vector<real_type> y) {
        const size_t n = std::min(x.size(), y.size());
        return self.build_G1(n, x.data(), y.data());
      }, py::arg("x"), py::arg("y"),
      R"S(
        Constructs a G1 biarc list using the points provided. The two vectors 
        are trimmed to the length of the shortest one.

        :param List[float] x: list of x coordinates
        :param List[float] y: list of y coordinates
        :return: true if build succeeds
        :rtype: bool
      )S")

      .def("build_G1", [](BiarcList & self, const std::vector<real_type> & x, 
        const std::vector<real_type> & y, const std::vector<real_type> & theta) {
        const size_t n = std::min({x.size(), y.size(), theta.size()});
        return self.build_G1(n, x.data(), y.data(), theta.data());
      }, py::arg("x"), py::arg("y"), py::arg("theta"),
      R"S(
        Constructs a G1 biarc list using the points provided. The two vectors 
        are trimmed to the length of the shortest one. This version allows to 
        specify the angle at nodes.

        :param List[float] x: list of x coordinates
        :param List[float] y: list of y coordinates
        :param List[float] theta: angles at nodes
        :return: true if build succeeds
        :rtype: bool
      )S")

      .def("getAtS", &BiarcList::getAtS, py::arg("s"),
      R"S(
        Get the biarc at coordinates **s**

        :param float s: **s** coordinates
        :return: biarc at **s**
        :rtype: Biarc
      )S")
      
      .def("numSegments", &BiarcList::numSegments,
      R"S(
        Returns the number of elements in the BiarcList

        :return: the number of elements
        :rtype: int
      )S")

      .def("__len__", &BiarcList::numSegments,
      R"S(
        Returns the number of elements in the BiarcList

        :return: the number of elements
        :rtype: int
      )S")
      
      .def("findAtS", &BiarcList::findAtS, py::arg("s"),
      R"S(
        Get the index of the biarc at coordinate **s**

        :param float s: coordinate **s**
        :return: index of the biarc at that coordinate
        :rtype: int
      )S")
      
      .def("segment_length", &BiarcList::segment_length, py::arg("nseg"),
      R"S(
        Return the length of the segment at index `nseg`

        :param int nseg: index of the segment
        :return: length of the segment at that index
        :rtype: float
      )S")
      
      .def("segment_length_ISO", &BiarcList::segment_length_ISO, py::arg("nseg"), py::arg("offs"),
      R"S(
        Return the length of the segment at index `nseg`. Version with ISO offset.

        :param int nseg: index of the segment
        :param float offs: offset from the curve
        :return: length of the segment at that index
        :rtype: float
      )S")
      
      .def("build_AABBtree_ISO", &BiarcList::build_AABBtree_ISO,
        py::arg("offs"), py::arg("max_angle") = Utils::m_pi/6, py::arg("max_size") = 1e100,
      R"S(
        Build the internal AABB tree of the biarc list with offset (ISO)
        
        :param float offs:      curve offset
        :param float max_angle: maximum angle variation of the arc covered by a triangle
        :param float max_size:  maximum admissible size of the covering tirnagles
        :return: nothing, works in place
        :rtype: NoneType
      )S")
      
      .def("build_AABBtree_SAE", &BiarcList::build_AABBtree_SAE,
        py::arg("offs"), py::arg("max_angle") = Utils::m_pi/6, py::arg("max_size") = 1e100,
      R"S(
        Build the internal AABB tree of the biarc list with offset (SAE)
        
        :param float offs:      curve offset
        :param float max_angle: maximum angle variation of the arc covered by a triangle
        :param float max_size:  maximum admissible size of the covering tirnagles
        :return: nothing, works in place
        :rtype: NoneType
      )S")
      
      .def("getSTK", [](const BiarcList & self) {
        const size_t n = self.numSegments();
        // Avoids segmentation fault for empty list
        if (!n) {
          return make_tuple(std::vector<real_type>(), 
                            std::vector<real_type>(), 
                            std::vector<real_type>());
        }
        std::vector<real_type> s(n), t(n), k(n);
        self.getSTK(s.data(), t.data(), k.data());
        return std::make_tuple(s, t, k);
      },
      R"S(
        Returns the list of nodes angles and curvatures of the biarc

        :return: s, angles and curvatures in nodes as lists
        :rtype: Tuple[List[float], List[float], List[float]]
      )S")
      
      .def("getXY", [](const BiarcList & self) {
        const size_t n = self.numSegments();
        // Avoids segmentation fault for empty list
        if (!n) {
          return make_tuple(std::vector<real_type>(), 
                            std::vector<real_type>());
        }
        std::vector<real_type> x(n+1), y(n+1);
        self.getXY(x.data(), y.data());
        return std::make_tuple(x, y);
      },
      R"S(
        Returns the list of nodes coordinates

        :return: nodes coordinates, a list of x and y values
        :rtype: Tuple[List[float], List[float]]
      )S")
      
      .def("findST1", [](const BiarcList & self, real_type x, real_type y) {
        real_type s, t;
        int_type idx = self.findST1(x, y, s, t);
        return std::make_tuple(idx, s, t);
      }, py::arg("x"), py::arg("y"),
      R"S(
        Find parametric coordinate :math:`(s, t)` given a point. With respect to the 
        classic `findST`, this version returns as first return value the index
        of the segment.

        :param float x: **x** coordinates of the point
        :param float y: **y** coordinates of the point
        :return: a tuple with index, **s** coordinates and **t** coordinates
        :rtype: Tuple[int, float, float]
      )S")

      .def("findST1", [](const BiarcList & self, const std::vector<real_type> & x, const std::vector<real_type> & y) {
        const size_t n = std::min(x.size(), y.size());
        std::vector<real_type> s(n), t(n);
        std::vector<int_type> idx(n);
        for(size_t i = 0; i < n; i++) {
          idx[i] = self.findST1(x[i], y[i], s[i], t[i]);
        }
        return std::make_tuple(idx, s, t);
      }, py::arg("x"), py::arg("y"),
      R"S(
        Find parametric coordinate :math:`(s, t)` given a point. With respect to the 
        classic `findST`, this version returns as first return value the index
        of the segment.
        Vectorial Version.

        :param List[float] x: **x** coordinates of the point
        :param List[float] y: **y** coordinates of the point
        :return: a tuple with index, **s** coordinates and **t** coordinates
        :rtype: Tuple[List[bool], List[float], List[float]]
      )S")
      
      .def("findST1", [](const BiarcList & self, int_type ibegin, int_type iend, real_type x, real_type y) {
        real_type s, t;
        int_type idx = self.findST1(ibegin, iend, x, y, s, t);
        return std::make_tuple(idx, s, t);
      }, py::arg("ibegin"), py::arg("iend"), py::arg("x"), py::arg("y"),
      R"S(
        Find parametric coordinate :math:`(s, t)` given a point. With respect to the 
        classic `findST`, this version returns as first return value the index
        of the segment. The search is done only between elements ``ibegin..iend``.

        :param int ibegin: initial index
        :param int iend: final index
        :param float x: **x** coordinates of the point
        :param float y: **y** coordinates of the point
        :return: a tuple with index, **s** coordinates and **t** coordinates
        :rtype: Tuple[int, float, float]
      )S");
    }
  }
}