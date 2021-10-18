/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-BaseCurve.hh"


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
        :return: the curvature at **s**
        :rtype: float
      )S")

      .def("theta_DD", &BaseCurve::theta_DD, py::arg("s"),
      R"S(
        Angle second derivative (derivative of curvature) at 
        curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the derivative of curvature at **s**
        :rtype: float
      )S")

      .def("theta_DDD", &BaseCurve::theta_DDD, py::arg("s"),
      R"S(
        Angle third derivative (second derivative of curvature) at 
        curvilinear abscissa **s**

        :param float s: curvilinear abscissa
        :return: the second derivative of curvature at **s**
        :rtype: float
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

      .def("eval", [](const BaseCurve & self, real_type s) {
        real_type x, y;
        self.eval(s, x, y);
        return std::make_tuple(x, y);
      }, py::arg("s"),
      R"S(
        Evaluate curve at curvilinear coordinate `s`.
        
        :param float s:  curvilinear coordinate
        :return: a tuple with angle, curvature, x-coordinate and y-coordinate at **s**
        :rtype: Tuple[float, float, float, float]
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
        :rtype: Tuple[float, float, float, float]
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
        :rtype: Tuple[float, float, float, float]
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
        :rtype: Tuple[float, float, float, float]
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
        Translate the curve by :math:`(t_x, t_y)`. This method works in place

        :param float tx: **x** coordinates offset
        :param float ty: **y** coordinates offset
        :return: Nothing
        :rtype: NoneType
      )S")

      .def("rotate", &BaseCurve::rotate, py::arg("angle"), py::arg("cx"), py::arg("cy"),
      R"S(

      )S")
      
      .def("scale", &BaseCurve::scale)
      .def("reverse", &BaseCurve::reverse)
      .def("changeOrigin", &BaseCurve::changeOrigin)
      .def("trim", &BaseCurve::trim)
      .def("collision", &BaseCurve::collision)
      .def("collision_ISO", &BaseCurve::collision_ISO)
      .def("collision_SAE", &BaseCurve::collision_SAE)

      // intersect
      // intersect_ISO
      // itersect_SAE

      .def("closestPoint_ISO", [](BaseCurve * self, real_type qx, real_type qy) {
        real_type x, y, s, t, dst;
        self->closestPoint_ISO(qx, qy, x, y, s, t, dst);
        return std::make_tuple(x, y, s, t, dst);
      })  
      .def("closestPoint_ISO", [](BaseCurve * self, real_type qx, real_type qy, real_type offs) {
        real_type x, y, s, t, dst;
        self->closestPoint_ISO(qx, qy, offs, x, y, s, t, dst);
        return std::make_tuple(x, y, s, t, dst);
      })   
      .def("closestPoint_SAE", [](BaseCurve * self, real_type qx, real_type qy) {
        real_type x, y, s, t, dst;
        self->closestPoint_SAE(qx, qy, x, y, s, t, dst);
        return std::make_tuple(x, y, s, t, dst);
      })  
      .def("closestPoint_SAE", [](BaseCurve * self, real_type qx, real_type qy, real_type offs) {
        real_type x, y, s, t, dst;
        self->closestPoint_SAE(qx, qy, offs, x, y, s, t, dst);
        return std::make_tuple(x, y, s, t, dst);
      })

      .def("distance", &BaseCurve::distance)
      .def("distance_ISO", &BaseCurve::distance_ISO)
      .def("distance_SAE", &BaseCurve::distance_SAE)
      .def("info", &BaseCurve::info)
        
      

      
        
        
      .def("findST_ISO", [](BaseCurve * self, real_type x, real_type y) {
        real_type s, t;
        self->findST_ISO(x, y, s, t);
        return std::make_tuple(s, t);
      })
      .def("findST_SAE", [](BaseCurve * self, real_type x, real_type y) {
        real_type s, t;
        self->findST_SAE(x, y, s, t);
        return std::make_tuple(s, t);
      })
      
      .def("__str__", [](BaseCurve * self) {
        std::ostringstream str;
        self->info(str);
        return str.str();
      });
    }
  }
}