/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-BaseCurve.hh"
#include "pybind11/stl.h"


namespace G2lib {
  namespace python {
    void wrap_BaseCurve(py::module & m) {
      py::enum_<CurveType>(m, "CurveType")
        .value("G2LIB_LINE", G2lib::CurveType::G2LIB_LINE)
        .value("G2LIB_POLYLINE", G2lib::CurveType::G2LIB_POLYLINE)
        .value("G2LIB_CIRCLE", G2lib::CurveType::G2LIB_CIRCLE)
        .value("G2LIB_BIARC", G2lib::CurveType::G2LIB_BIARC)
        .value("G2LIB_BIARC_LIST", G2lib::CurveType::G2LIB_BIARC_LIST)
        .value("G2LIB_CLOTHOID", G2lib::CurveType::G2LIB_CLOTHOID)
        .value("G2LIB_CLOTHOID_LIST", G2lib::CurveType::G2LIB_CLOTHOID_LIST)
        .export_values();

      py::class_<BaseCurve, PythonicBaseCurve>(m, "BaseCurve",
      R"S(
        The base curve that rapresents the base for all the curves that are contained
        inside the library. The default constructor select the type from the CurveType
        enum. The class is a pure abstract class and cannot be directly instantiated 
        by the user, but acts as a prototype.

        Several methods are defined in ISO or SAE standard. The difference is in the
        signs for the reference frame attached to the curve.

         * ISO: the positive normal direction is on the left, while the negative is 
           on the right. This standard is used also in the undefined standard methods,
           where the offset is always considered to be 0. This standard follows the 
           right-hand rule.
         * SAE: the positive normal direction is on the right while the negative is
           on the right.

        Select your methods accordinglu

        :param CurveType type: enum representing the type of the curve
      )S")
      
      .def(py::init<CurveType const &>(), py::arg("type"))
      
      .def("type", &BaseCurve::type, 
      R"S(
        Returns the type of the curve

        :return: the type of the curve
        :rtype: CurveType
      )S")
      
      .def("length", &BaseCurve::length,
      R"S(
        Length of the current curve

        :return: length of the curve
        :rtype: float
      )S")
      
      .def("length_ISO", &BaseCurve::length_ISO, py::arg("offs"),
      R"S(
        Length of the current curve, with offset (ISO)
        
        :param float offs: curve offset
        :return: length with offset (ISO)
        :rtype: float
      )S")
      
      .def("length_SAE", &BaseCurve::length_SAE, py::arg("offs"),
      R"S(
        Length of the current curve, with offset (SAE)

        :param float offs: curve offset
        :return: length with offset (SAE)
        :rtype: float
      )S")

      .def("bbox", [](const BaseCurve & self) {
        real_type x_min, y_min, x_max, y_max;
        self.bbox(x_min, y_min, x_max, y_max);
        return std::make_tuple(
          std::make_tuple(x_min, y_min), 
          std::make_tuple(x_max, y_max));
      },
      R"S(
        Compute the bounding box of the curve, as a tuple of 
        tuples with minimum point and maximum point.

        :return: the extrema of the bounding box of the curve
        :rtype: Tuple[Tuple[float, float], Tuple[float, float]]
      )S")
      
      .def("bbox_ISO", [](const BaseCurve & self, real_type offs) {
        real_type x_min, y_min, x_max, y_max;
        self.bbox_ISO(offs, x_min, y_min, x_max, y_max);
        return std::make_tuple(
          std::make_tuple(x_min, y_min), 
          std::make_tuple(x_max, y_max));
      }, py::arg("offs"),
      R"S(
        Compute the bounding box of the curve, as a tuple of 
        tuples with minimum point and maximum point. With ISO offset

        :param float offs: curve offset
        :return: the extrema of the bounding box of the curve, with ISO offset
        :rtype: Tuple[Tuple[float, float], Tuple[float, float]]
      )S")

      .def("bbox_SAE", [](const BaseCurve & self, real_type offs) {
        real_type x_min, y_min, x_max, y_max;
        self.bbox_SAE(offs, x_min, y_min, x_max, y_max);
        return std::make_tuple(
          std::make_tuple(x_min, y_min), 
          std::make_tuple(x_max, y_max));
      }, py::arg("offs"),
      R"S(
        Compute the bounding box of the curve, as a tuple of 
        tuples with minimum point and maximum point. With ISO offset

        :param float offs: curve offset
        :return: the extrema of the bounding box of the curve, with ISO offset
        :rtype: Tuple[Tuple[float, float], Tuple[float, float]]
      )S")

      .def("bbTriangles", [](const BaseCurve & self, real_type max_angle, real_type max_size, int_type icurve) {
        std::vector<Triangle2D> tvec;
        self.bbTriangles(tvec, max_angle, max_size, icurve);
        return tvec;
      }, py::arg("max_angle") = Utils::m_pi/18, py::arg("max_size") = 1e100, py::arg("icurve") = 0,
      R"S(
        Build a cover with triangles of the curve. Returns a list of triangles.
      
        :param float max_angle: maximum angle variation of the curve covered by a triangle
        :param float max_size:  maximum admissible size of the covering tirnagles
        :param int icurve:    index of the covering triangles
        :return: list of covering triangles
        :rtype: List[Triangle2D]
      )S")

      .def("bbTriangles_ISO", [](const BaseCurve & self, real_type offs, real_type max_angle, real_type max_size, int_type icurve) {
        std::vector<Triangle2D> tvec;
        self.bbTriangles_ISO(offs, tvec, max_angle, max_size, icurve);
        return tvec;
      }, py::arg("offs"), py::arg("max_angle") = Utils::m_pi/18, py::arg("max_size") = 1e100, py::arg("icurve") = 0,
      R"S(
        Build a cover with triangles of the curve, using an ISO curve offset. 
        Returns a list of triangles.
      
        :param float offs:      curve offset
        :param float max_angle: maximum angle variation of the curve covered by a triangle
        :param float max_size:  maximum admissible size of the covering tirnagles
        :param int icurve:    index of the covering triangles
        :return: list of covering triangles
        :rtype: List[Triangle2D]
      )S")

      .def("bbTriangles_SAE", [](const BaseCurve & self, real_type offs, real_type max_angle, real_type max_size, int_type icurve) {
        std::vector<Triangle2D> tvec;
        self.bbTriangles_SAE(offs, tvec, max_angle, max_size, icurve);
        return tvec;
      }, py::arg("offs"), py::arg("max_angle") = Utils::m_pi/18, py::arg("max_size") = 1e100, py::arg("icurve") = 0,
      R"S(
        Build a cover with triangles of the curve, using a SAE curve offset. 
        Returns a list of triangles.
      
        :param float offs:      curve offset
        :param float max_angle: maximum angle variation of the curve covered by a triangle
        :param float max_size:  maximum admissible size of the covering tirnagles
        :param int icurve:    index of the covering triangles
        :return: list of covering triangles
        :rtype: List[Triangle2D]
      )S")

      .def("thetaBegin", &BaseCurve::thetaBegin,
      R"S(
        Initial angle of the curve

        :return: the initial angle of the curve
        :rtype: float
      )S")
      
      .def("theta", &BaseCurve::theta, py::arg("s"),
      R"S(
        Angle of the curve at curvilinear value **s**

        :param float s: curvilinear abscissa
        :return: the angle of the curve at **s**
        :rtype: float
      )S")

      .def("theta", [](const BaseCurve & self, const std::vector<real_type> & s) {
        const size_t n = s.size();
        std::vector<real_type> ret(n);
        for (size_t i = 0; i < n; i++) ret[i] = self.theta(s[i]);
        return ret;
      }, py::arg("s"),
      R"S(
        Angle of the curve at curvilinear value **s**. Vectorial Version.

        :param List[float] s: curvilinear abscissa
        :return: the angle of the curve at **s**
        :rtype: List[float]
      )S")

      .def("thetaEnd", &BaseCurve::thetaEnd,
      R"S(
        Final angle of the curve

        :return: the final angle of the curve
        :rtype: float
      )S")

      .def("kappaBegin", &BaseCurve::kappaBegin,
      R"S(
        Initial curvature of the curve

        :return: the initial curvature of the curve
        :rtype: float
      )S")

      .def("kappa", &BaseCurve::kappa, py::arg("s"),
      R"S(
        Curvature of the curve at curvilinear value **s**

        :param float s: curvilinear abscissa
        :return: the curvature angle of the curve
        :rtype: float
      )S")

      .def("kappa", [](const BaseCurve & self, const std::vector<real_type> & s) {
        const size_t n = s.size();
        std::vector<real_type> ret(n);
        for (size_t i = 0; i < n; i++) ret[i] = self.kappa(s[i]);
        return ret;
      }, py::arg("s"),
      R"S(
        Curvature of the curve at curvilinear value **s**. Vectorial Version.

        :param List[float] s: curvilinear abscissa
        :return: the curvature of the curve at **s**
        :rtype: List[float]
      )S")

      .def("kappaEnd", &BaseCurve::kappaEnd, 
      R"S(
        Final curvature of the curve

        :return: the final curvature of the curve
        :rtype: float
      )S")
      
      .def("xBegin", &BaseCurve::xBegin,
      R"S(
        Initial coordinate **x** of the curve

        :return: the initial coordinate **x** of the curve
        :rtype: float
      )S")

      .def("X", &BaseCurve::X, py::arg("s"), 
      R"S(
        Coordinate **x** of the curve at curvilinear value **s**

        :param float s: curvilinear abscissa
        :return: coordinate **x** of the curve
        :rtype: float
      )S")

      .def("X", [](const BaseCurve & self, const std::vector<real_type> & s) {
        const size_t n = s.size();
        std::vector<real_type> ret(n);
        for (size_t i = 0; i < n; i++) ret[i] = self.X(s[i]);
        return ret;
      }, py::arg("s"), 
      R"S(
        Coordinate **x** of the curve at curvilinear value **s**. Vectorial Version.

        :param List[float] s: curvilinear abscissa
        :return: coordinate **x** of the curve
        :rtype: List[float]
      )S")

      .def("xEnd", &BaseCurve::xEnd, 
      R"S(
        Final coordinate **x** of the curve

        :return: the final coordinate **x** of the curve
        :rtype: float
      )S")
      
      .def("yBegin", &BaseCurve::yBegin, 
      R"S(
        Initial coordinate **y** of the curve

        :return: the initial coordinate **y** of the curve
        :rtype: float
      )S")

      .def("Y", &BaseCurve::Y, py::arg("s"),
      R"S(
        Coordinate **y** of the curve at curvilinear value **s**

        :param float s: curvilinear abscissa
        :return: coordinate **y** of the curve
        :rtype: float
      )S")

      .def("Y", [](const BaseCurve & self, const std::vector<real_type> & s) {
        const size_t n = s.size();
        std::vector<real_type> ret(n);
        for (size_t i = 0; i < n; i++) ret[i] = self.Y(s[i]);
        return ret;
      }, py::arg("s"), 
      R"S(
        Coordinate **y** of the curve at curvilinear value **s**. Vectorial Version.

        :param List[float] s: curvilinear abscissa
        :return: coordinate **y** of the curve
        :rtype: List[float]
      )S")

      .def("yEnd", &BaseCurve::yEnd, 
      R"S(
        Final coordinate **y** of the curve

        :return: the final coordinate **y** of the curve
        :rtype: float
      )S")

      .def("xBegin_ISO", &BaseCurve::xBegin_ISO, py::arg("offs"),
      R"S(
        Initial coordinate **x** of the curve, with ISO offset

        :param float offs: curve ISO offset
        :return: the initial coordinate **x** of the curve
        :rtype: float
      )S")

      .def("X_ISO", &BaseCurve::X_ISO, py::arg("s"), py::arg("offs"),
      R"S(
        Coordinate **x** of the curve at curvilinear value **s**

        :param float s: curvilinear abscissa
        :param float offs: curve ISO offset
        :return: coordinate **x** of the curve
        :rtype: float
      )S")

      .def("X_ISO", [](const BaseCurve & self, std::vector<real_type> s, std::vector<real_type> offs) {
        const size_t n = std::min(s.size(), offs.size());
        std::vector<real_type> x(n);
        for (size_t i = 0; i < n; i++) {
          x[i] = self.X_ISO(s[i], offs[i]);
        }
        return x;
      }, py::arg("s"), py::arg("offs"),
      R"S(
        Coordinate **x** of the curve at curvilinear value **s**.
        Vectorial version.

        :param List[float] s: curvilinear abscissa
        :param List[float] offs: curve ISO offset
        :return: coordinate **x** of the curve
        :rtype: List[float]
      )S")

      .def("xEnd_ISO", &BaseCurve::xEnd_ISO, py::arg("offs"),
      R"S(
        Final coordinate **x** of the curve, with ISO offset

        :param float offs: curve ISO offset
        :return: the final coordinate **x** of the curve
        :rtype: float
      )S")
      
      .def("yBegin_ISO", &BaseCurve::yBegin_ISO, py::arg("offs"), 
      R"S(
        Initial coordinate **y** of the curve, with ISO offset

        :param float offs: curve ISO offset
        :return: the initial coordinate **y** of the curve
        :rtype: float
      )S")

      .def("Y_ISO", &BaseCurve::Y_ISO, py::arg("s"), py::arg("offs"),
      R"S(
        Coordinate **y** of the curve at curvilinear value **s**,
        with ISO offset

        :param s: curvilinear abscissa
        :param float offs: curve ISO offset
        :return: coordinate **y** of the curve
        :rtype: float
      )S")

      .def("Y_ISO", [](const BaseCurve & self, std::vector<real_type> s, std::vector<real_type> offs) {
        const size_t n = std::min(s.size(), offs.size());
        std::vector<real_type> y(n);
        for (size_t i = 0; i < n; i++) {
          y[i] = self.Y_ISO(s[i], offs[i]);
        }
        return y;
      }, py::arg("s"), py::arg("offs"),
      R"S(
        Coordinate **y** of the curve at curvilinear value **s**.
        Vectorial version.

        :param List[float] s: curvilinear abscissa
        :param List[float] offs: curve ISO offset
        :return: coordinate **y** of the curve
        :rtype: List[float]
      )S")

      .def("yEnd_ISO", &BaseCurve::yEnd_ISO, py::arg("offs"), 
      R"S(
        Final coordinate **y** of the curve, with ISO offset

        :param float offs: curve ISO offset
        :return: the final coordinate **y** of the curve
        :rtype: float
      )S")
      
      .def("xBegin_SAE", &BaseCurve::xBegin_SAE, py::arg("offs"),
      R"S(
        Initial coordinate **x** of the curve, with SAE offset

        :param float offs: curve SAE offset
        :return: the initial coordinate **x** of the curve
        :rtype: float
      )S")

      .def("X_SAE", &BaseCurve::X_SAE, py::arg("s"), py::arg("offs"),
      R"S(
        Coordinate **x** of the curve at curvilinear value **s**

        :param float s: curvilinear abscissa
        :param float offs: curve SAE offset
        :return: coordinate **x** of the curve
        :rtype: float
      )S")

      .def("xEnd_SAE", &BaseCurve::xEnd_SAE, py::arg("offs"),
      R"S(
        Final coordinate **x** of the curve, with SAE offset

        :param float offs: curve SAE offset
        :return: the final coordinate **x** of the curve
        :rtype: float
      )S")
      
      .def("yBegin_SAE", &BaseCurve::yBegin_SAE, py::arg("offs"), 
      R"S(
        Initial coordinate **y** of the curve, with SAE offset

        :param float offs: curve SAE offset
        :return: the initial coordinate **y** of the curve
        :rtype: float
      )S")

      .def("Y_SAE", &BaseCurve::Y_SAE, py::arg("s"), py::arg("offs"),
      R"S(
        Coordinate **y** of the curve at curvilinear value **s**,
        with SAE offset

        :param s: curvilinear abscissa
        :param float offs: curve SAE offset
        :return: coordinate **y** of the curve
        :rtype: float
      )S")

      .def("yEnd_SAE", &BaseCurve::yEnd_SAE, py::arg("offs"), 
      R"S(
        Final coordinate **y** of the curve, with SAE offset

        :param float offs: curve SAE offset
        :return: the final coordinate **y** of the curve
        :rtype: float
      )S")

      .def("tx_Begin", &BaseCurve::tx_Begin,
      R"S(
        Initial tangent **x** coordinate

        :return: the initial tangent **x** coordinate
        :rtype: float
      )S")

      .def("tx", &BaseCurve::tx, py::arg("s"),
      R"S(
        Tangent **x** coordinate for the curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the tangent **x** coordinate
        :rtype: float
      )S")

      .def("tx", [](const BaseCurve & self, const std::vector<real_type> & s) {
        const size_t n = s.size();
        std::vector<real_type> ret(n);
        for (size_t i = 0; i < n; i++) ret[i] = self.tx(s[i]);
        return ret;
      }, py::arg("s"),
      R"S(
        Tangent **x** coordinate for the curvilinear abscissa **s**.
        Vectorial version.

        :param List[float] s: curvilinear abscissa
        :return: the tangent **x** coordinate
        :rtype: List[float]
      )S")

      .def("tx_End", &BaseCurve::tx_End, 
      R"S(
        Final tangent **x** coordinate

        :return: the final tangent **x** coordinate
        :rtype: float
      )S")

      .def("ty_Begin", &BaseCurve::ty_Begin,
      R"S(
        Initial tangent **y** coordinate

        :return: the initial tangent **y** coordinate
        :rtype: float
      )S")

      .def("ty", &BaseCurve::ty, py::arg("s"),
      R"S(
        Tangent **y** coordinate for the curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the tangent **y** coordinate
        :rtype: float
      )S")

      .def("ty", [](const BaseCurve & self, const std::vector<real_type> & s) {
        const size_t n = s.size();
        std::vector<real_type> ret(n);
        for (size_t i = 0; i < n; i++) ret[i] = self.ty(s[i]);
        return ret;
      }, py::arg("s"),
      R"S(
        Tangent **y** coordinate for the curvilinear abscissa **s**.
        Vectorial version.

        :param List[float] s: curvilinear abscissa
        :return: the tangent **y** coordinate
        :rtype: List[float]
      )S")

      .def("ty_End", &BaseCurve::ty_End,
      R"S(
        Final tangent **y** coordinate

        :return: the final tangent **y** coordinate
        :rtype: float
      )S")

      .def("nx_Begin_ISO", &BaseCurve::nx_Begin_ISO,
      R"S(
        Initial normal **x** coordinate, in ISO standard

        :return: the initial normal **x** coordinate
        :rtype: float
      )S")
      
      .def("ny_Begin_ISO", &BaseCurve::ny_Begin_ISO,
      R"S(
        Initial normal **y** coordinate, in ISO standard

        :return: the initial normal **y** coordinate
        :rtype: float
      )S")
      
      .def("nx_End_ISO", &BaseCurve::nx_End_ISO, 
      R"S(
        Final normal **x** coordinate, in ISO standard

        :return: the final normal **x** coordinate
        :rtype: float
      )S")

      .def("ny_End_ISO", &BaseCurve::ny_End_ISO, 
      R"S(
        Final normal **y** coordinate, in ISO standard

        :return: the final normal **y** coordinate
        :rtype: float
      )S")

      .def("nx_Begin_SAE", &BaseCurve::nx_Begin_SAE,
      R"S(
        Initial normal **x** coordinate, in SAE standard

        :return: the initial normal **x** coordinate
        :rtype: float
      )S")
      
      .def("ny_Begin_SAE", &BaseCurve::ny_Begin_SAE,
      R"S(
        Initial normal **y** coordinate, in SAE standard

        :return: the initial normal **y** coordinate
        :rtype: float
      )S")
      
      .def("nx_End_SAE", &BaseCurve::nx_End_SAE, 
      R"S(
        Final normal **x** coordinate, in SAE standard

        :return: the final normal **x** coordinate
        :rtype: float
      )S")

      .def("ny_End_SAE", &BaseCurve::ny_End_SAE, 
      R"S(
        Final normal **y** coordinate, in SAE standard

        :return: the final normal **y** coordinate
        :rtype: float
      )S")

      .def("theta_D", &BaseCurve::theta_D, py::arg("s"),
      R"S(
        Angle derivative (curvature) at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the angle derivative at **s**
        :rtype: float
      )S")

      .def("theta_D", [](const BaseCurve & self, const std::vector<real_type> & s) {
        const size_t n = s.size();
        std::vector<real_type> ret(n);
        for (size_t i = 0; i < n; i++) ret[i] = self.theta_D(s[i]);
        return ret;
      }, py::arg("s"),
      R"S(
        Angle derivative (curvature) at curvilinear abscissa **s**.
        Vectorial version.

        :param List[float] s: curvilinear abscissa
        :return: the angle derivative at **s**
        :rtype: List[float]
      )S")

      .def("theta_DD", &BaseCurve::theta_DD, py::arg("s"),
      R"S(
        Angle second derivative (derivative of curvature) at 
        curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the derivative of curvature at **s**
        :rtype: float
      )S")

      .def("theta_DD", [](const BaseCurve & self, const std::vector<real_type> & s) {
        const size_t n = s.size();
        std::vector<real_type> ret(n);
        for (size_t i = 0; i < n; i++) ret[i] = self.theta_DD(s[i]);
        return ret;
      }, py::arg("s"),
      R"S(
        Angle second derivative (curvature derivative) at curvilinear abscissa **s**.
        Vectorial version.

        :param List[float] s: curvilinear abscissa
        :return: the curvature derivative at **s**
        :rtype: List[float]
      )S")

      .def("theta_DDD", &BaseCurve::theta_DDD, py::arg("s"),
      R"S(
        Angle third derivative (second derivative of curvature) at 
        curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the second derivative of curvature at **s**
        :rtype: float
      )S")

      .def("theta_DDD", [](const BaseCurve & self, const std::vector<real_type> & s) {
        const size_t n = s.size();
        std::vector<real_type> ret(n);
        for (size_t i = 0; i < n; i++) ret[i] = self.theta_DDD(s[i]);
        return ret;
      }, py::arg("s"),
      R"S(
        Angle third derivative (curvature second derivative) at 
        curvilinear abscissa **s**. Vectorial version.

        :param List[float] s: curvilinear abscissa
        :return: the curvature second derivative at **s**
        :rtype: List[float]
      )S")

      .def("kappa_D", &BaseCurve::kappa_D, py::arg("s"),
      R"S(
        Curvature derivative at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the derivative of curvature at **s**
        :rtype: float
      )S")

      .def("kappa_DD", &BaseCurve::kappa_DD, py::arg("s"),
      R"S(
        Curvature second derivative at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the second derivative of curvature at **s**
        :rtype: float
      )S")

      .def("tx_D", &BaseCurve::tx_D, py::arg("s"),
      R"S(
        Tangent **x** coordinate derivative at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the tangent **x** coordinate derivative at **s**
        :rtype: float
      )S")
      
      .def("tx_DD", &BaseCurve::tx_DD, py::arg("s"),
      R"S(
        Tangent **x** coordinate second derivative at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the tangent **x** coordinate second derivative at **s**
        :rtype: float
      )S")
      
      .def("tx_DDD", &BaseCurve::tx_DDD, py::arg("s"),
      R"S(
        Tangent **x** coordinate third derivative at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the tangent **x** coordinate third derivative at **s**
        :rtype: float
      )S")

      .def("ty_D", &BaseCurve::ty_D, py::arg("s"),
      R"S(
        Tangent **y** coordinate derivative at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the tangent **y** coordinate derivative at **s**
        :rtype: float
      )S")
      
      .def("ty_DD", &BaseCurve::ty_DD, py::arg("s"),
      R"S(
        Tangent **y** coordinate second derivative at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the tangent **y** coordinate second derivative at **s**
        :rtype: float
      )S")
      
      .def("ty_DDD", &BaseCurve::ty_DDD, py::arg("s"),
      R"S(
        Tangent **y** coordinate third derivative at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the tangent **y** coordinate third derivative at **s**
        :rtype: float
      )S")

      .def("nx_ISO", &BaseCurve::nx_ISO, py::arg("s"),
      R"S(
        Normal **x** coordinate at curvilinear abscissa **s**, 
        in ISO standard

        :param float s: curvilinear abscissa
        :return: normal **x** coordinate at **s**
        :rtype: float
      )S")
      
      .def("nx_ISO_D", &BaseCurve::nx_ISO_D, py::arg("s"),
      R"S(
        Normal **x** coordinate derivative at curvilinear abscissa **s**, 
        in ISO standard

        :param float s: curvilinear abscissa
        :return: normal **x** coordinate derivative at **s**
        :rtype: float
      )S")

      .def("nx_ISO_DD", &BaseCurve::nx_ISO_DD, py::arg("s"),
      R"S(
        Normal **x** coordinate second derivative at curvilinear abscissa **s**, 
        in ISO standard

        :param float s: curvilinear abscissa
        :return: normal **x** coordinate second derivative at **s**
        :rtype: float
      )S")

      .def("nx_ISO_DDD", &BaseCurve::nx_ISO_DDD, py::arg("s"),
      R"S(
        Normal **x** coordinate third derivative at curvilinear abscissa **s**, 
        in ISO standard

        :param float s: curvilinear abscissa
        :return: normal **x** coordinate third derivative at **s**
        :rtype: float
      )S")
      
      .def("ny_ISO", &BaseCurve::ny_ISO, py::arg("s"),
      R"S(
        Normal **y** coordinate at curvilinear abscissa **s**, 
        in ISO standard

        :param float s: curvilinear abscissa
        :return: normal **y** coordinate at **s**
        :rtype: float
      )S")

      .def("ny_ISO_D", &BaseCurve::ny_ISO_D, py::arg("s"),
      R"S(
        Normal **y** coordinate derivative at curvilinear abscissa **s**, 
        in ISO standard

        :param float s: curvilinear abscissa
        :return: normal **y** coordinate derivative at **s**
        :rtype: float
      )S")

      .def("ny_ISO_DD", &BaseCurve::ny_ISO_DD, py::arg("s"),
      R"S(
        Normal **y** coordinate second derivative at curvilinear abscissa **s**, 
        in ISO standard

        :param float s: curvilinear abscissa
        :return: normal **x** coordinate second derivative at **s**
        :rtype: float
      )S")

      .def("ny_ISO_DDD", &BaseCurve::ny_ISO_DDD, py::arg("s"),
      R"S(
        Normal **x** coordinate third derivative at curvilinear abscissa **s**, 
        in ISO standard

        :param float s: curvilinear abscissa
        :return: normal **x** coordinate third derivative at **s**
        :rtype: float
      )S")

      .def("nx_SAE", &BaseCurve::nx_SAE, py::arg("s"),
      R"S(
        Normal **x** coordinate at curvilinear abscissa **s**, 
        in SAE standard

        :param float s: curvilinear abscissa
        :return: normal **x** coordinate at **s**
        :rtype: float
      )S")
      
      .def("nx_SAE_D", &BaseCurve::nx_SAE_D, py::arg("s"),
      R"S(
        Normal **x** coordinate derivative at curvilinear abscissa **s**, 
        in SAE standard

        :param float s: curvilinear abscissa
        :return: normal **x** coordinate derivative at **s**
        :rtype: float
      )S")

      .def("nx_SAE_DD", &BaseCurve::nx_SAE_DD, py::arg("s"),
      R"S(
        Normal **x** coordinate second derivative at curvilinear abscissa **s**, 
        in SAE standard

        :param float s: curvilinear abscissa
        :return: normal **x** coordinate second derivative at **s**
        :rtype: float
      )S")

      .def("nx_SAE_DDD", &BaseCurve::nx_SAE_DDD, py::arg("s"),
      R"S(
        Normal **x** coordinate third derivative at curvilinear abscissa **s**, 
        in SAE standard

        :param float s: curvilinear abscissa
        :return: normal **x** coordinate third derivative at **s**
        :rtype: float
      )S"
      )
      
      .def("ny_SAE", &BaseCurve::ny_SAE, py::arg("s"),
      R"S(
        Normal **y** coordinate at curvilinear abscissa **s**, 
        in SAE standard

        :param float s: curvilinear abscissa
        :return: normal **y** coordinate at **s**
        :rtype: float
      )S")

      .def("ny_SAE_D", &BaseCurve::ny_SAE_D, py::arg("s"),
      R"S(
        Normal **y** coordinate derivative at curvilinear abscissa **s**, 
        in SAE standard

        :param float s: curvilinear abscissa
        :return: normal **y** coordinate derivative at **s**
        :rtype: float
      )S")

      .def("ny_SAE_DD", &BaseCurve::ny_SAE_DD, py::arg("s"),
      R"S(
        Normal **y** coordinate second derivative at curvilinear abscissa **s**, 
        in SAE standard

        :param float s: curvilinear abscissa
        :return: normal **x** coordinate second derivative at **s**
        :rtype: float
      )S")

      .def("ny_SAE_DDD", &BaseCurve::ny_SAE_DDD, py::arg("s"),
      R"S(
        Normal **x** coordinate third derivative at curvilinear abscissa **s**, 
        in SAE standard

        :param float s: curvilinear abscissa
        :return: normal **x** coordinate third derivative at **s**
        :rtype: float
      )S")

      .def("tg", [](const BaseCurve & self, real_type s) {
        real_type r_x, r_y;
        self.tg(s, r_x, r_y);
        return std::make_tuple(r_x, r_y);
      }, py::arg("s"),
      R"S(
        Tangent at curvilinear coordinates **s**

        :param float s:
        :return: the tangent at curvilinear coordinate **s**
        :rtype: Tuple[float, float]
      )S")

      .def("tg_D", [](const BaseCurve & self, real_type s) {
        real_type r_x, r_y;
        self.tg_D(s, r_x, r_y);
        return std::make_tuple(r_x, r_y);
      }, py::arg("s"),
      R"S(
        Tangent derivative at curvilinear coordinates **s**

        :param float s:
        :return: the tangent derivative at curvilinear coordinate **s**
        :rtype: Tuple[float, float]
      )S")

      .def("tg_DD", [](const BaseCurve & self, real_type s) {
        real_type r_x, r_y;
        self.tg_DD(s, r_x, r_y);
        return std::make_tuple(r_x, r_y);
      }, py::arg("s"),
      R"S(
        Tangent second derivative at curvilinear coordinates **s**

        :param float s:
        :return: the tangent second derivative at curvilinear coordinate **s**
        :rtype: Tuple[float, float]
      )S")

      .def("tg_DDD", [](const BaseCurve & self, real_type s) {
        real_type r_x, r_y;
        self.tg_DDD(s, r_x, r_y);
        return std::make_tuple(r_x, r_y);
      }, py::arg("s"),
      R"S(
        Tangent third derivative at curvilinear coordinates **s**

        :param float s:
        :return: the tangent third derivative at curvilinear coordinate **s**
        :rtype: Tuple[float, float]
      )S")

      .def("nor_ISO", [](const BaseCurve & self, real_type s) {
        real_type r_x, r_y;
        self.nor_ISO(s, r_x, r_y);
        return std::make_tuple(r_x, r_y);
      }, py::arg("s"),
      R"S(
        Normal at curvilinear coordinates **s**, in ISO standard

        :param float s:
        :return: the normal at curvilinear coordinate **s**
        :rtype: Tuple[float, float]
      )S")

      .def("nor_ISO_D", [](const BaseCurve & self, real_type s) {
        real_type r_x, r_y;
        self.nor_ISO_D(s, r_x, r_y);
        return std::make_tuple(r_x, r_y);
      }, py::arg("s"),
      R"S(
        Normal first derivative at curvilinear coordinates **s**, in ISO standard

        :param float s:
        :return: the normal first derivative at curvilinear coordinate **s**
        :rtype: Tuple[float, float]
      )S")

      .def("nor_ISO_DD", [](const BaseCurve & self, real_type s) {
        real_type r_x, r_y;
        self.nor_ISO_DD(s, r_x, r_y);
        return std::make_tuple(r_x, r_y);
      }, py::arg("s"),
      R"S(
        Normal second derivative at curvilinear coordinates **s**, in ISO standard

        :param float s:
        :return: the normal second derivative at curvilinear coordinate **s**
        :rtype: Tuple[float, float]
      )S")

      .def("nor_ISO_DDD", [](const BaseCurve & self, real_type s) {
        real_type r_x, r_y;
        self.nor_ISO_DDD(s, r_x, r_y);
        return std::make_tuple(r_x, r_y);
      }, py::arg("s"),
      R"S(
        Normal third derivative at curvilinear coordinates **s**, in ISO standard

        :param float s:
        :return: the normal third derivative at curvilinear coordinate **s**
        :rtype: Tuple[float, float]
      )S")

      .def("nor_SAE", [](const BaseCurve & self, real_type s) {
        real_type r_x, r_y;
        self.nor_SAE(s, r_x, r_y);
        return std::make_tuple(r_x, r_y);
      }, py::arg("s"),
      R"S(
        Normal at curvilinear coordinates **s**, in SAE standard

        :param float s:
        :return: the normal at curvilinear coordinate **s**
        :rtype: Tuple[float, float]
      )S")

      .def("nor_SAE_D", [](const BaseCurve & self, real_type s) {
        real_type r_x, r_y;
        self.nor_SAE_D(s, r_x, r_y);
        return std::make_tuple(r_x, r_y);
      }, py::arg("s"),
      R"S(
        Normal first derivative at curvilinear coordinates **s**, in SAE standard

        :param float s:
        :return: the normal first derivative at curvilinear coordinate **s**
        :rtype: Tuple[float, float]
      )S")

      .def("nor_SAE_DD", [](const BaseCurve & self, real_type s) {
        real_type r_x, r_y;
        self.nor_SAE_DD(s, r_x, r_y);
        return std::make_tuple(r_x, r_y);
      }, py::arg("s"),
      R"S(
        Normal second derivative at curvilinear coordinates **s**, in SAE standard

        :param float s:
        :return: the normal second derivative at curvilinear coordinate **s**
        :rtype: Tuple[float, float]
      )S")

      .def("nor_SAE_DDD", [](const BaseCurve & self, real_type s) {
        real_type r_x, r_y;
        self.nor_SAE_DDD(s, r_x, r_y);
        return std::make_tuple(r_x, r_y);
      }, py::arg("s"),
      R"S(
        Normal third derivative at curvilinear coordinates **s**, in SAE standard

        :param float s:
        :return: the normal third derivative at curvilinear coordinate **s**
        :rtype: Tuple[float, float]
      )S")

      .def("evaluate", [](const BaseCurve & self, real_type s) {
        real_type th, k, x, y;
        self.evaluate(s, th, k, x, y);
        return std::make_tuple(th, k, x, y);
      }, py::arg("s"),
      R"S(
        Evaluate curve at curvilinear coordinate `s`.
        
        :param float s:  curvilinear coordinate
        :return: a tuple with angle, curvature, x-coordinate and y-coordinate at **s**
        :rtype: Tuple[float, float, float, float]
      )S")

      .def("evaluate", [](const BaseCurve & self, const std::vector<real_type> & s) {
        const size_t n = s.size();
        std::vector<real_type> th(s), k(s), x(s), y(s);
        for (size_t i = 0; i < n; i++) {
          self.evaluate(s[i], th[i], k[i], x[i], y[i]);
        }
        return std::make_tuple(th, k, x, y);
      }, py::arg("s"),
      R"S(
        Evaluate curve at curvilinear coordinate `s`.
        
        :param float s:  curvilinear coordinate
        :return: a tuple with angle, curvature, x-coordinate and y-coordinate at **s**
        :rtype: Tuple[float, float, float, float]
      )S")

      .def("eval", [](const BaseCurve & self, real_type s) {
        real_type x, y;
        self.eval(s, x, y);
        return std::make_tuple(x, y);
      }, py::arg("s"),
      R"S(
        Evaluate curve at curvilinear coordinate `s`.
        
        :param float s:  curvilinear coordinate
        :return: a tuple with x-coordinate and y-coordinate at **s**
        :rtype: Tuple[float, float]
      )S")

      .def("eval", [](const BaseCurve & self, const std::vector<real_type> & s) {
        const size_t n = s.size();
        std::vector<real_type> x(n), y(n);
        for (size_t i = 0; i < n; i++) {
          self.eval(s[i], x[i], y[i]);
        }
        return std::make_tuple(x, y);
      }, py::arg("s"),
      R"S(
        Evaluate curve at curvilinear coordinate `s`.
        
        :param List[float] s: list of curvilinear coordinate
        :return: a tuple of lists with x-coordinate and y-coordinate at **s**
        :rtype: Tuple[List[float], List[float]]
      )S")

      .def("evaluate_ISO", [](const BaseCurve & self, real_type s, real_type offs) {
        real_type th, k, x, y;
        self.evaluate_ISO(s, offs, th, k, x, y);
        return std::make_tuple(th, k, x, y);
      }, py::arg("s"), py::arg("offs"),
      R"S(
        Evaluate curve at curvilinear coordinate `s`, with an additional offset
        that follows the ISO standard.
        
        :param float s:  curvilinear coordinate
        :param float offs:
        :return: a tuple with angle, curvature, x-coordinate and y-coordinate at **s**
        :rtype: Tuple[float, float, float, float]
      )S")

      .def("evaluate_ISO", [](const BaseCurve & self, const std::vector<real_type> & s, const std::vector<real_type> & offs) {
        const size_t n = std::min(s.size(), offs.size());
        std::vector<real_type> th(n), k(n), x(n), y(n);
        for (size_t i = 0; i < n; i++) {
          self.evaluate_ISO(s[i], offs[i], th[i], k[i], x[i], y[i]);
        }
        return std::make_tuple(th, k, x, y);
      }, py::arg("s"), py::arg("offs"),
      R"S(
        Evaluate curve at curvilinear coordinate `s`, with an additional offset
        that follows the ISO standard.
        
        :param List[float] s:  curvilinear coordinate
        :param List[float] offs:
        :return: a tuple with angle, curvature, x-coordinate and y-coordinate at **s**
        :rtype: Tuple[List[float], List[float], List[float], List[float]]
      )S")

      .def("eval_ISO", [](const BaseCurve & self, real_type s, real_type offs) {
        real_type x, y;
        self.eval_ISO(s, offs, x, y);
        return std::make_tuple(x, y);
      }, py::arg("s"), py::arg("offs"),
      R"S(
        Evaluate curve at curvilinear coordinate `s`, with an additional offset
        that follows the ISO standard.
        
        :param float s:  curvilinear coordinate
        :param float offs:
        :return: a tuple with x-coordinate and y-coordinate at **s**
        :rtype: Tuple[float, float]
      )S")

      .def("eval_ISO", [](const BaseCurve & self, const std::vector<real_type> & s, const std::vector<real_type> & offs) {
        const size_t n = std::min(s.size(), offs.size());
        std::vector<real_type> x(n), y(n);
        for (size_t i = 0; i < n; i++) {
          self.eval_ISO(s[i], offs[i], x[i], y[i]);
        }
        return std::make_tuple(x, y);
      }, py::arg("s"), py::arg("offs"),
      R"S(
        Evaluate curve at curvilinear coordinate `s`, with an additional offset
        that follows the ISO standard.
        Vector version.
        
        :param List[float] s:  curvilinear coordinate
        :param List[float] offs:
        :return: a tuple with x-coordinate and y-coordinate at **s**
        :rtype: Tuple[List[float], List[float]]
      )S")

      .def("evaluate_SAE", [](const BaseCurve & self, real_type s, real_type offs) {
        real_type th, k, x, y;
        self.evaluate_SAE(s, offs, th, k, x, y);
        return std::make_tuple(th, k, x, y);
      }, py::arg("s"), py::arg("offs"),
      R"S(
        Evaluate curve at curvilinear coordinate `s`, with an additional offset
        that follows the SAE standard.
        
        :param float s:  curvilinear coordinate
        :param float offs:
        :return: a tuple with angle, curvature, x-coordinate and y-coordinate at **s**
        :rtype: Tuple[float, float, float, float]
      )S")

      .def("eval_SAE", [](const BaseCurve & self, real_type s, real_type offs) {
        real_type x, y;
        self.eval_SAE(s, offs, x, y);
        return std::make_tuple(x, y);
      }, py::arg("s"), py::arg("offs"),
      R"S(
        Evaluate curve at curvilinear coordinate `s`, with an additional offset
        that follows the SAE standard.
        
        :param float s:  curvilinear coordinate
        :param float offs:
        :return: a tuple with x-coordinate and y-coordinate at **s**
        :rtype: Tuple[float, float]
      )S")
      
      .def("X_D", &BaseCurve::X_D, py::arg("s"),
      R"S(
        **x** coordinate derivative at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: **x** coordinate derivative
        :rtype: float
      )S")
      
      .def("X_DD", &BaseCurve::X_DD, py::arg("s"),
      R"S(
        **x** coordinate second derivative at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: **x** coordinate second derivative
        :rtype: float
      )S")

      .def("X_DDD", &BaseCurve::X_DDD, py::arg("s"),
      R"S(
        **x** coordinate third derivative at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: **x** coordinate third derivative
        :rtype: float
      )S")

      .def("Y_D", &BaseCurve::Y_D, py::arg("s"),
      R"S(
        **y** coordinate derivative at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: **y** coordinate derivative
        :rtype: float
      )S")

      .def("Y_DD", &BaseCurve::Y_DD, py::arg("s"),
      R"S(
        **y** coordinate second derivative at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: **y** coordinate second derivative
        :rtype: float
      )S")

      .def("Y_DDD", &BaseCurve::Y_DDD, py::arg("s"),
      R"S(
        **y** coordinate third derivative at curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: **y** coordinate third derivative
        :rtype: float
      )S")

      .def("eval_D", [](const BaseCurve & self, real_type s) {
        real_type x, y;
        self.eval_D(s, x, y);
        return std::make_tuple(x, y);
      }, py::arg("s"),
      R"S(
        Evaluate curve first derivative at curvilinear coordinate `s`.
        
        :param float s:  curvilinear coordinate
        :return: a tuple with derivative x-coordinate and y-coordinate at **s**
        :rtype: Tuple[float, float]
      )S")

      .def("eval_D", [](const BaseCurve & self, const std::vector<real_type> & s) {
        const size_t n = s.size();
        std::vector<real_type> x(n), y(n);
        for (size_t i = 0; i < n; i++) {
          self.eval_D(s[i], x[i], y[i]);
        }
        return std::make_tuple(x, y);
      }, py::arg("s"),
      R"S(
        Evaluate curve first derivative at curvilinear coordinate `s`.
        Vector version.
        
        :param List[float] s:  curvilinear coordinate
        :return: a tuple with derivative x-coordinate and y-coordinate at **s**
        :rtype: Tuple[List[float], List[float]]
      )S")

      .def("eval_DD", [](const BaseCurve & self, real_type s) {
        real_type x, y;
        self.eval_DD(s, x, y);
        return std::make_tuple(x, y);
      }, py::arg("s"),
      R"S(
        Evaluate curve second derivative at curvilinear coordinate `s`.
        
        :param float s:  curvilinear coordinate
        :return: a tuple with second derivative x-coordinate and y-coordinate at **s**
        :rtype: Tuple[float, float]
      )S")

      .def("eval_DD", [](const BaseCurve & self, const std::vector<real_type> & s) {
        const size_t n = s.size();
        std::vector<real_type> x(n), y(n);
        for (size_t i = 0; i < n; i++) {
          self.eval_DD(s[i], x[i], y[i]);
        }
        return std::make_tuple(x, y);
      }, py::arg("s"),
      R"S(
        Evaluate curve second derivative at curvilinear coordinate `s`.
        Vector version
        
        :param List[float] s:  curvilinear coordinate
        :return: a tuple with second derivative x-coordinate and y-coordinate at **s**
        :rtype: Tuple[List[float], List[float]]
      )S")
      
      .def("eval_DDD", [](const BaseCurve & self, real_type s) {
        real_type x, y;
        self.eval_DDD(s, x, y);
        return std::make_tuple(x, y);
      }, py::arg("s"),
      R"S(
        Evaluate curve third derivative at curvilinear coordinate `s`.
        
        :param float s:  curvilinear coordinate
        :return: a tuple with third derivative x-coordinate and y-coordinate at **s**
        :rtype: Tuple[float, float]
      )S")

      .def("eval_DDD", [](const BaseCurve & self, const std::vector<real_type> & s) {
        const size_t n = s.size();
        std::vector<real_type> x(n), y(n);
        for (size_t i = 0; i < n; i++) {
          self.eval_DDD(s[i], x[i], y[i]);
        }
        return std::make_tuple(x, y);
      }, py::arg("s"),
      R"S(
        Evaluate curve third derivative at curvilinear coordinate `s`.
        Vector version.
        
        :param List[float] s:  curvilinear coordinate
        :return: a tuple with third derivative x-coordinate and y-coordinate at **s**
        :rtype: Tuple[List[float], List[float]]
      )S")

      .def("X_ISO_D", &BaseCurve::X_ISO_D, py::arg("s"), py::arg("offs"),
      R"S(
        **x** coordinate derivative at curvilinear abscissa **s**, in ISO
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: **x** coordinate derivative
        :rtype: float
      )S")

      .def("X_ISO_DD", &BaseCurve::X_ISO_DD, py::arg("s"), py::arg("offs"),
      R"S(
        **x** coordinate second derivative at curvilinear abscissa **s**, in ISO
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: **x** coordinate second derivative
        :rtype: float
      )S")
      
      .def("X_ISO_DDD", &BaseCurve::X_ISO_DDD, py::arg("s"), py::arg("offs"),
      R"S(
        **x** coordinate third derivative at curvilinear abscissa **s**, in ISO
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: **x** coordinate third derivative
        :rtype: float
      )S")

      .def("Y_ISO_D", &BaseCurve::Y_ISO_D, py::arg("s"), py::arg("offs"),
      R"S(
        **y** coordinate derivative at curvilinear abscissa **s**, in ISO
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: **y** coordinate derivative
        :rtype: float
      )S")

      .def("Y_ISO_DD", &BaseCurve::Y_ISO_DD, py::arg("s"), py::arg("offs"),
      R"S(
        **y** coordinate second derivative at curvilinear abscissa **s**, in ISO
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: **y** coordinate second derivative
        :rtype: float
      )S")

      .def("Y_ISO_DDD", &BaseCurve::Y_ISO_DDD, py::arg("s"), py::arg("offs"),
      R"S(
        **y** coordinate third derivative at curvilinear abscissa **s**, in ISO
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: **y** coordinate third derivative
        :rtype: float
      )S")

      // SAE

      .def("X_SAE_D", &BaseCurve::X_SAE_D, py::arg("s"), py::arg("offs"),
      R"S(
        **x** coordinate derivative at curvilinear abscissa **s**, in SAE
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: **x** coordinate derivative
        :rtype: float
      )S")

      .def("X_SAE_DD", &BaseCurve::X_SAE_DD, py::arg("s"), py::arg("offs"),
      R"S(
        **x** coordinate second derivative at curvilinear abscissa **s**, in SAE
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: **x** coordinate second derivative
        :rtype: float
      )S")
      
      .def("X_SAE_DDD", &BaseCurve::X_SAE_DDD, py::arg("s"), py::arg("offs"),
      R"S(
        **x** coordinate third derivative at curvilinear abscissa **s**, in SAE
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: **x** coordinate third derivative
        :rtype: float
      )S")

      .def("Y_SAE_D", &BaseCurve::Y_SAE_D, py::arg("s"), py::arg("offs"),
      R"S(
        **y** coordinate derivative at curvilinear abscissa **s**, in SAE
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: **y** coordinate derivative
        :rtype: float
      )S")

      .def("Y_SAE_DD", &BaseCurve::Y_SAE_DD, py::arg("s"), py::arg("offs"),
      R"S(
        **y** coordinate second derivative at curvilinear abscissa **s**, in SAE
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: **y** coordinate second derivative
        :rtype: float
      )S")

      .def("Y_SAE_DDD", &BaseCurve::Y_SAE_DDD, py::arg("s"), py::arg("offs"),
      R"S(
        **y** coordinate third derivative at curvilinear abscissa **s**, in SAE
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: **y** coordinate third derivative
        :rtype: float
      )S")
  
      .def("eval_ISO_D", [](const BaseCurve & self, real_type s, real_type offs) {
        real_type x, y;
        self.eval_ISO_D(s, offs, x, y);
        return std::make_tuple(x, y);
      }, py::arg("s"), py::arg("offs"),
      R"S(
        Coordinate derivatives at curvilinear abscissa **s**, in ISO
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: coordinate derivatives
        :rtype: Tuple[float, float]
      )S")
      
      .def("eval_ISO_DD", [](const BaseCurve & self, real_type s, real_type offs) {
        real_type x, y;
        self.eval_ISO_DD(s, offs, x, y);
        return std::make_tuple(x, y);
      }, py::arg("s"), py::arg("offs"),
      R"S(
        Coordinate second derivatives at curvilinear abscissa **s**, in ISO
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: coordinate second derivatives
        :rtype: Tuple[float, float]
      )S")

      .def("eval_ISO_DDD", [](const BaseCurve & self, real_type s, real_type offs) {
        real_type x, y;
        self.eval_ISO_DDD(s, offs, x, y);
        return std::make_tuple(x, y);
      }, py::arg("s"), py::arg("offs"),
      R"S(
        Coordinate third derivatives at curvilinear abscissa **s**, in ISO
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: coordinate third derivatives
        :rtype: Tuple[float, float]
      )S")

      .def("eval_SAE_D", [](const BaseCurve & self, real_type s, real_type offs) {
        real_type x, y;
        self.eval_SAE_D(s, offs, x, y);
        return std::make_tuple(x, y);
      }, py::arg("s"), py::arg("offs"),
      R"S(
        Coordinate derivatives at curvilinear abscissa **s**, in SAE
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: coordinate derivatives
        :rtype: Tuple[float, float]
      )S")
      
      .def("eval_SAE_DD", [](const BaseCurve & self, real_type s, real_type offs) {
        real_type x, y;
        self.eval_SAE_DD(s, offs, x, y);
        return std::make_tuple(x, y);
      }, py::arg("s"), py::arg("offs"),
      R"S(
        Coordinate second derivatives at curvilinear abscissa **s**, in SAE
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: coordinate second derivatives
        :rtype: Tuple[float, float]
      )S")

      .def("eval_SAE_DDD", [](const BaseCurve & self, real_type s, real_type offs) {
        real_type x, y;
        self.eval_SAE_DDD(s, offs, x, y);
        return std::make_tuple(x, y);
      }, py::arg("s"), py::arg("offs"),
      R"S(
        Coordinate third derivatives at curvilinear abscissa **s**, in SAE
        reference frame

        :param float s: curvilinear abscissa
        :param float offs: offset from the curve
        :return: coordinate third derivatives
        :rtype: Tuple[float, float]
      )S")
    
      .def("translate", &BaseCurve::translate, py::arg("tx"), py::arg("ty"),
      R"S(
        Translate the curve by :math:`(t_x, t_y)`. This method works in place,
        and forces the recalculation of the AABBtree of the curve

        :param float tx: **x** coordinates offset
        :param float ty: **y** coordinates offset
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("rotate", &BaseCurve::rotate, py::arg("angle"), py::arg("cx"), py::arg("cy"),
      R"S(
        Rotate curve by angle :math:`\theta` centered at point  :math:`(c_x, c_y)`.
        This method works in place and forces the recalculation of the AABBtree 
        of the curve
        
        :param float angle: angle :math:`\theta`
        :param float cx: center origin coordinate :math:`c_x`
        :param float cy: center origin coordinate :math:`c_y`
        :return: nothing, works in place
        :rtype: NoneType
      )S")
      
      .def("scale", &BaseCurve::scale, py::arg("scale"),
      R"S(
        Scale a curve by the factor provided as input. This method works in place
        and forces the recalculation of the AABBtree of the curve.

        :param float sc: the scaling factor
        :return: nothing, works in place
        :rtype: NoneType
      )S")
      
      .def("reverse", &BaseCurve::reverse,
      R"S(
        Reverses the curve parametrization. This method works in place and
        forces the recalculation of the AABBtree of the curve

        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("changeOrigin", &BaseCurve::changeOrigin, py::arg("newx0"), py::arg("newy0"),
      R"S(
        Translate the curve so that the origin will be :math:`(x_{0,new}, y_{0,new})`.
        This method works in place and forces the recalculation of the AABBtree of the
        curve.

        :param float newx0: new origin **x** coordinate
        :param float newy0: new origin **y** coordinate
        :return: nothing, works in place
        :rtype: NoneType
      )S")
      
      .def("trim", &BaseCurve::trim, py::arg("s_begin"), py::arg("s_end"),
      R"S(
        Cuts the curve between parametric curvilinear abscissa range 
        :math:`(s_{begin}, s_{end})`. This method works in place and forces
        the recalculation of the AABBtree of the curve

        :param float s_begin: starting curvilinear abscissa
        :param float e_sen: ending curvilinear abscissa
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("collision", &BaseCurve::collision, py::arg("curve"),
      R"S(
        Checks collisions with another curve

        :param BaseCurve curve: the other curve to check
        :return: true if curves collide
        :rtype: bool
      )S")

      .def("collision_ISO", &BaseCurve::collision_ISO, py::arg("offs"), py::arg("curve"), py::arg("curve_offs"),
      R"S(
        Checks collisions with another curve, considering offset on both 
        curves. Uses ISO standard references

        :param float offs: offset on the current curve
        :param BaseCurve curve: the other curve to check
        :param float curve_offs: offset on the other curve
        :return: true if curves collide
        :rtype: bool
      )S")

      .def("collision_SAE", &BaseCurve::collision_SAE, py::arg("offs"), py::arg("curve"), py::arg("curve_offs"),
      R"S(
        Checks collisions with another curve, considering offset on both 
        curves. Uses SAE standard references

        :param float offs: offset on the current curve
        :param BaseCurve curve: the other curve to check
        :param float curve_offs: offset on the other curve
        :return: true if curves collide
        :rtype: bool
      )S")

      .def("intersect", [](const BaseCurve & self, const BaseCurve & curve, bool swap_s_vals) {
        IntersectList ilist;
        self.intersect(curve, ilist, swap_s_vals);
        return ilist;
      }, py::arg("curve"), py::arg("swap_s_vals") = false,
      R"S(
        Intersect the curve with another curve. The result is a list of pairs, where 
        each element of the pair represents the ascissa coordinate at which intersection
        occurs
        
        :param BaseCurve curve: second curve intersect
        :param bool swap_s_vals: if true store `(s2,s1)` instead of `(s1,s2)` for each
                                 intersection. Default false.
        :return: a list of pair with intersection coordinates on both curves
        :rtype: List[Tuple[float, float]]
      )S")

      .def("intersect_ISO", [](const BaseCurve & self, real_type offs, const BaseCurve & curve, 
        real_type curve_offs, bool swap_s_vals) {
        IntersectList ilist;
        self.intersect_ISO(offs, curve, curve_offs, ilist, swap_s_vals);
        return ilist;
      }, py::arg("offs"), py::arg("curve"), py::arg("curve_offs"), py::arg("swap_s_vals") = false,
      R"S(
        Intersect the curve with another curve. The result is a list of pairs, where 
        each element of the pair represents the ascissa coordinate at which intersection
        occurs. Version with offset in ISO standard reference
        
        :param float offs: offset on the curve
        :param BaseCurve curve: second curve intersect
        :param float curve_offs: offset on the second curve
        :param bool swap_s_vals: if true store `(s2,s1)` instead of `(s1,s2)` for each
                                 intersection. Default false.
        :return: a list of pair with intersection coordinates on both curves
        :rtype: List[Tuple[float, float]]
      )S")
      
      .def("intersect_SAE", [](const BaseCurve & self, real_type offs, const BaseCurve & curve, 
        real_type curve_offs, bool swap_s_vals) {
        IntersectList ilist;
        self.intersect_SAE(offs, curve, curve_offs, ilist, swap_s_vals);
        return ilist;
      }, py::arg("offs"), py::arg("curve"), py::arg("curve_offs"), py::arg("swap_s_vals") = false,
      R"S(
        Intersect the curve with another curve. The result is a list of pairs, where 
        each element of the pair represents the ascissa coordinate at which intersection
        occurs. Version with offset in ISO standard reference
        
        :param float offs: offset on the curve
        :param BaseCurve curve: second curve intersect
        :param float curve_offs: offset on the second curve
        :param bool swap_s_vals: if true store `(s2,s1)` instead of `(s1,s2)` for each
                                 intersection. Default false.
        :return: a list of pair with intersection coordinates on both curves
        :rtype: List[Tuple[float, float]]
      )S")

      .def("closestPoint", [](const BaseCurve & self, real_type qx, real_type qy) {
        int_type ret;
        real_type x, y, s, t, dst;
        ret = self.closestPoint_ISO(qx, qy, x, y, s, t, dst);
        return std::make_tuple(ret, x, y, s, t, dst);
      }, py::arg("qx"), py::arg("qy"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the closest point on the curve.

        There are 6 values in the response tuple:

         1. result of the projection. If the value is 1 the projection is unique and 
            orthogonal, if the value is 0 there is more than one orthogonal projection
            and only the first is returned, if the value is -1 the minimum distance point 
            has a projection which is not orthogonal
         2. **x** coordinate of the point on the clothoid
         3. **y** coordinate of the point on the clothoid
         4. **s** curvilinear abscissa of the point
         5. **t** normal distance of the point on curvilinear abscissa
         6. **dst** distance of the point from the curve

        :param float qx: x coordinates of the point
        :param float qy: y coordinates of the point
        :return: a tuple of results as described
        :rtype: Tuple[int, float, float, float, float, int]
      )S")

      .def("closestPoint", [](const BaseCurve & self, const std::vector<real_type> & qx, const std::vector<real_type> & qy) {
        const size_t n = std::min(qx.size(), qy.size());
        std::vector<int_type> ret(n);
        std::vector<real_type> x(n), y(n), s(n), t(n), dst(n);
        for (size_t i = 0; i < n; i++) {
          ret[i] = self.closestPoint_ISO(qx[i], qy[i], x[i], y[i], s[i], t[i], dst[i]);
        }
        return std::make_tuple(ret, x, y, s, t, dst);
      }, py::arg("qx"), py::arg("qy"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the closest point on the curve.
        Vectorial version.

        There are 6 values in the response tuple:

         1. result of the projection. If the value is 1 the projection is unique and 
            orthogonal, if the value is 0 there is more than one orthogonal projection
            and only the first is returned, if the value is -1 the minimum distance point 
            has a projection which is not orthogonal
         2. **x** coordinate of the point on the clothoid
         3. **y** coordinate of the point on the clothoid
         4. **s** curvilinear abscissa of the point
         5. **t** normal distance of the point on curvilinear abscissa
         6. **dst** distance of the point from the curve

        :param List[float] qx: x coordinates of the point
        :param List[float] qy: y coordinates of the point
        :return: a tuple of results as described
        :rtype: Tuple[List[int], List[float], List[float], List[float], List[float], List[int]]
      )S")

      .def("closestPoint_ISO", [](const BaseCurve & self, real_type qx, real_type qy) {
        int_type ret;
        real_type x, y, s, t, dst;
        ret = self.closestPoint_ISO(qx, qy, x, y, s, t, dst);
        return std::make_tuple(ret, x, y, s, t, dst);
      }, py::arg("qx"), py::arg("qy"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the closest point on the curve. This
        function uses the ISO standard.

        There are 6 values in the response tuple:

         1. result of the projection. If the value is 1 the projection is unique and 
            orthogonal, if the value is 0 there is more than one orthogonal projection
            and only the first is returned, if the value is -1 the minimum distance point 
            has a projection which is not orthogonal
         2. **x** coordinate of the point on the clothoid
         3. **y** coordinate of the point on the clothoid
         4. **s** curvilinear abscissa of the point
         5. **t** normal distance of the point on curvilinear abscissa
         6. **dst** distance of the point from the curve

        :param float qx: x coordinates of the point
        :param float qy: y coordinates of the point
        :return: a tuple of results as described
        :rtype: Tuple[int, float, float, float, float, int]
      )S")

      .def("closestPoint_SAE", [](const BaseCurve & self, real_type qx, real_type qy) {
        int_type ret;
        real_type x, y, s, t, dst;
        ret = self.closestPoint_SAE(qx, qy, x, y, s, t, dst);
        return std::make_tuple(ret, x, y, s, t, dst);
      }, py::arg("qx"), py::arg("qy"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the closest point on the curve. This
        function uses SAE standard.

        There are 6 values in the response tuple:

         1. result of the projection. If the value is 1 the projection is unique and 
            orthogonal, if the value is 0 there is more than one orthogonal projection
            and only the first is returned, if the value is -1 the minimum distance point 
            has a projection which is not orthogonal
         2. **x** coordinate of the point on the clothoid
         3. **y** coordinate of the point on the clothoid
         4. **s** curvilinear abscissa of the point
         5. **t** normal distance of the point on curvilinear abscissa
         6. **dst** distance of the point from the curve

        :param float qx: x coordinates of the point
        :param float qy: y coordinates of the point
        :return: a tuple of results as described
        :rtype: Tuple[int, float, float, float, float, int]
      )S")

      .def("distance", &BaseCurve::distance, py::arg("qx"), py::arg("qy"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the distance between the point and the
        curve.

        :param float qx: x coordinates of the point
        :param float qy: y coordinates of the point
        :return: the distance of the point
        :rtype: float
      )S")

      .def("distance_ISO", &BaseCurve::distance_ISO, py::arg("qx"), py::arg("qy"), py::arg("offs"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the distance between the point and the
        curve, even with an offset. It uses ISO standard reference frame.

        :param float qx: x coordinates of the point
        :param float qy: y coordinates of the point
        :param float offs: offset from the curve
        :return: the distance of the point
        :rtype: float
      )S")

      .def("distance_SAE", &BaseCurve::distance_SAE, py::arg("qx"), py::arg("qy"), py::arg("offs"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the distance between the point and the
        curve, even with an offset. It uses SAE standard reference frame.

        :param float qx: x coordinates of the point
        :param float qy: y coordinates of the point
        :param float offs: offset from the curve
        :return: the distance of the point
        :rtype: float
      )S")
        
      .def("findST", [](const BaseCurve & self, real_type x, real_type y) {
        real_type s, t;
        bool ret = self.findST_ISO(x, y, s, t);
        return std::make_tuple(ret, s, t);
      }, py::arg("x"), py::arg("y"), 
      R"S(
        Find the curvilinear coordinate of point :math:`P = (x, y)`
        with respect to the curve, such that: :math:`P = C(s) + t\,N(s)`
        where :math:`C(s)` is the curve position respect to the curvilinear 
        coordinates and :math:`t` is the normal :math:`N(s)` at the point.
       
        :param float x: **x** component
        :param float y: **y** component
        :return: a tuple with a boolean value (projection found or not) and
                 the **s** and **t** coordinates on the curve.
        :rtype: Tuple[bool, float, float]
      )S")

      .def("findST", [](const BaseCurve & self, const std::vector<real_type> & x, const std::vector<real_type> & y) {
        const size_t n = std::min(x.size(), y.size());
        std::vector<real_type> s(n), t(n);
        std::vector<bool> ret(n);
        for (size_t i = 0; i < n; i++) {
          ret[i] = self.findST_ISO(x[i], y[i], s[i], t[i]);
        }
        return std::make_tuple(ret, s, t);
      }, py::arg("x"), py::arg("y"), 
      R"S(
        Find the curvilinear coordinate of point :math:`P = (x, y)`
        with respect to the curve, such that: :math:`P = C(s) + t\,N(s)`
        where :math:`C(s)` is the curve position respect to the curvilinear 
        coordinates and :math:`t` is the normal :math:`N(s)` at the point.
        Vectorial version.
       
        :param List[float] x: **x** component
        :param List[float] y: **y** component
        :return: a tuple with a boolean value (projection found or not) and
                 the **s** and **t** coordinates on the curve.
        :rtype: Tuple[List[bool], List[float], List[float]]
      )S")

      .def("findST_ISO", [](const BaseCurve & self, real_type x, real_type y) {
        real_type s, t;
        bool ret = self.findST_ISO(x, y, s, t);
        return std::make_tuple(ret, s, t);
      }, py::arg("x"), py::arg("y"), 
      R"S(
        Find the curvilinear coordinate of point :math:`P = (x, y)`
        with respect to the curve, such that: :math:`P = C(s) + t\,N(s)`
        where :math:`C(s)` is the curve position respect to the curvilinear 
        coordinates and :math:`t` is the normal :math:`N(s)` at the point.
        It uses the ISO reference frame.
       
        :param float x: **x** component
        :param float y: **y** component
        :return: a tuple with a boolean value (projection found or not) and
                 the **s** and **t** coordinates on the curve.
        :rtype: Tuple[bool, float, float]
      )S")

      .def("findST_SAE", [](const BaseCurve & self, real_type x, real_type y) {
        real_type s, t;
        bool ret = self.findST_SAE(x, y, s, t);
        return std::make_tuple(ret, s, t);
      }, py::arg("x"), py::arg("y"), 
      R"S(
        Find the curvilinear coordinate of point :math:`P = (x, y)`
        with respect to the curve, such that: :math:`P = C(s) + t\,N(s)`
        where :math:`C(s)` is the curve position respect to the curvilinear 
        coordinates and :math:`t` is the normal :math:`N(s)` at the point.
        It uses the SAE reference frame.
       
        :param float x: **x** component
        :param float y: **y** component
        :return: a tuple with a boolean value (projection found or not) and
                 the **s** and **t** coordinates on the curve.
        :rtype: Tuple[bool, float, float]
      )S")
      
      .def("__str__", [](const BaseCurve & self) {
        std::ostringstream str;
        self.info(str);
        return str.str();
      });
    }
  }
}